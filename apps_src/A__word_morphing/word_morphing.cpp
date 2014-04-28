#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dimage.h"
#include "dtimer.h"
#include "dmorphink.h"
#include "dthresholder.h"
#include "dfeaturevector.h"
#include "dwordfeatures.h"
#include <math.h>

#ifndef D_NOTHREADS
#include "dthreads.h"
#endif


#define DO_FAST_PASS_FIRST 1

/*comparison function for qsort()ing doubles in nondecreasing order */
int compareDoubles(const void *p1, const void *p2){
  if ( (*(double*)(p1)) < (*(double*)(p2)) )
    return -1;
  else if ( (*(double*)(p1)) > (*(double*)(p2)) )
    return 1;
  return 0;
}



typedef struct {
  int wordIdx;
  double morphCost;
} MORPHCOST_T;

int compareMORPHCOST_T(const void *p1, const void *p2){
  MORPHCOST_T *m1,*m2;
  m1 = (MORPHCOST_T*)p1;
  m2 = (MORPHCOST_T*)p2;

  if((m1->morphCost) < (m2->morphCost))
    return -1;
  else if((m1->morphCost) > (m2->morphCost))
    return 1;
  return 0;
}


typedef struct{
  int numThreads;//how many threads are doing comparisons
  int threadNum;//which thread number this is (0..numThreads-1)
  int testWordIdx; //0-based (even if first test word >0) index into arrays
  int *rgNumLexReductionWordsSkipped;//each thread uses 1 slot of shared array
  DImage *pimgTest;//pointer shared by all threads
  DImage *rgTrainingImages;//pointer shared by all threads
  int numTrain;
  double *rgCostsMorph;//shared by all threads
  double *rgCostsDP; //shared by all threads
  double lengthPenalty;
  int meshSpacingStatic;//-1 if auto-calculate like originally done
  int numRefinesStatic;//-1 if auto-calculate like originally done
  double meshDiv;//what to divide word image height by to get mesh size [4.0]
  int bandWidthDP;//Sakoe-Chiba band DP constraint (rath/manmatha used 15)
#if DO_FAST_PASS_FIRST
  bool fFastPass; // true if this is a fast pass, false otherwise
  int fastPassN; // what is N for the top-N from fast pass?
  int *rgFastPassTopNidxs;//passed into thread func for slow pass
  double *rgSlowPassTopNMorphCosts;//filled in by thread func during slow pass
#endif
} WORDWARP_THREAD_PARMS;


void* word_morphing_thread_func(void *params){
 //Note:This looks like where the real comparison is done
  WORDWARP_THREAD_PARMS *pparms;
  int numTrain;
  DMorphInk mobj;
  int testWordIdx;

  pparms = (WORDWARP_THREAD_PARMS*)params;
  numTrain = pparms->numTrain;
  testWordIdx = pparms->testWordIdx;

  //now compare to all training values
  int numLexReductionWordsSkipped;
  numLexReductionWordsSkipped = 0;
  for(int tr=(pparms->threadNum); tr < numTrain; tr+=(pparms->numThreads)){
  	double morphCost=0.;
  	double DPcost=0.;
	
#if DO_FAST_PASS_FIRST
    	if(pparms->fFastPass)
      	morphCost = mobj.getWordMorphCostFast(*(pparms->pimgTest),
					  pparms->rgTrainingImages[tr],
					  pparms->bandWidthDP,/*15 */
					  0./*nonDiagonalCostDP*/,
					  pparms->meshSpacingStatic,
					  pparms->numRefinesStatic,
					  pparms->meshDiv,
					  pparms->lengthPenalty);
    	else
#else
      	morphCost = mobj.getWordMorphCost(*(pparms->pimgTest),
					pparms->rgTrainingImages[tr],
					pparms->bandWidthDP,/*15 bandWidthDP*/
					0./*nonDiagonalCostDP*/,
					pparms->meshSpacingStatic,
					pparms->numRefinesStatic,
					pparms->meshDiv,
					pparms->lengthPenalty);
#endif
	    DPcost = mobj.warpCostDP;
	    pparms->rgCostsMorph[testWordIdx*(long)numTrain+tr] = morphCost;
	    pparms->rgCostsDP[testWordIdx*(long)numTrain+tr] = DPcost;
	    pparms->rgNumLexReductionWordsSkipped[pparms->threadNum] =
		 numLexReductionWordsSkipped;
  	}
  	return NULL;
}


#if DO_FAST_PASS_FIRST
void* word_morphing_thread_func_slow_pass(void *params){
  WORDWARP_THREAD_PARMS *pparms;
  int numTrain;
  DMorphInk mobj;
  int testWordIdx;

  pparms = (WORDWARP_THREAD_PARMS*)params;
  numTrain = pparms->numTrain;
  testWordIdx = pparms->testWordIdx;

  //now do full morphing comparison to topN training values from fast pass
  for(int tr=(pparms->threadNum); tr < (pparms->fastPassN);
      tr+=(pparms->numThreads)){
    double morphCost=0.;
    int trIdx;
    trIdx = pparms->rgFastPassTopNidxs[tr];
    morphCost =
      mobj.getWordMorphCost(*(pparms->pimgTest),
			    pparms->rgTrainingImages[tr],
			    pparms->bandWidthDP,/*15 bandWidthDP*/
			    0./*nonDiagonalCostDP*/,
			    pparms->meshSpacingStatic,
			    pparms->numRefinesStatic,
			    pparms->meshDiv,
			    pparms->lengthPenalty);
    //pparms->rgCostsMorph[testWordIdx*(long)numTrain+tr] = morphCost;
    pparms->rgSlowPassTopNMorphCosts[tr] = morphCost;
  }
  return NULL;
}
#endif

// ---------------This function is not currently being used
//Return clipped version of a word image so that any whitespace around
//it is removed.  Assumes that the image is black and white (not
//grayscale).
DImage clipWordImageToInk(DImage &img){
  int x0, x1, y0, y1;
  int w, h;
  D_uint8 *p8;
  
  if(DImage::DImage_u8 != img.getImageType()){
    fprintf(stderr,"clipWordImageToInk() only accepts 8-bit black-and-white images\n");
    exit(1);
  }
  w = img.width();
  h = img.height();
  if((w<1) || (h<1)){
    fprintf(stderr,"clipWordImageToInk() width(%d) and height(%d) must be > 0!\n", w,h);
    exit(1);
  }
  x0 = w;
  x1 = 0;
  y0 = h;
  y1 = 0;
  p8 = img.dataPointer_u8();
  for(int y=0, idx=0; y < h; ++y){
    for(int x = 0; x < w; ++x, ++idx){
      if(0x00 == p8[idx]){//ink
	if(x < x0)
	  x0 = x;
	if(x > x1)
	  x1 = x;
	if(y < y0)
	  y0 = y;
	if(y > y1)
	  y1 = y;
      }
      else if(0xff != p8[idx]){// not 25word_morphing ../datasets/smith_jhl_vol1_lasso/smith_jhl_vol1_lasso.prj_intermediates/ 0 999 1000 1999 0.1 -1 -1 -1 4.0 14 10 /tmp/smith_results.dat5 or 0!
	fprintf(stderr,"clipWordImageToInk() expects all pixel values to be ink(0) or background(255).  saw %d at x,y=%d,%d\n", p8[idx], x, y);
	exit(1);
      }
    }
  }
  if((x1 < x0)||(y1 < y0)){
    fprintf(stderr,"clipWordImageToInk() found no ink pixels in the image!\n");
    DImage imgEmpty;
    return imgEmpty;
  }

  
  return (img.copy(x0,y0, x1-x0+1,y1-y0+1)).padEdges(1,1,1,1,
						     DImage::DImagePadValue, 255.);
}


int main(int argc, char **argv);
//Uhh, ya. This is not readable. Like, at all...
int main(int argc, char **argv){
  char stPathIn[1025];
  char stOutfile[1025];
  char stTmp[1025];
  int trainFirst, trainLast;
  int testFirst, testLast;
  int numTrain, numTest;
  DImage *rgTrainingImages;
  DImage testImage;
  double *rgCostsMorph;
  double *rgCostsDP;
  DTimer t1;
  DTimer t3;
  std::string *rgLabelsTrain;
  std::string *rgLabelsTest;
  int *rgPagesTrain;
  int *rgPagesTest;
  int *rgAuthorIdsTrain;
  int *rgAuthorIdsTest;
  int *rgThresholdsTrain;
  WORDWARP_THREAD_PARMS *rgThreadParms;
  int *rgNumLexReductionWordsSkipped;
  int numThreads = 1;
  double lengthPenalty = 0.;
  int meshSpacingStatic;
  int numRefinesStatic;
  double meshDiv;
  int bandWidthDP = 15;
  int slowPassN;//number of best matches from fast pass to do slow pass for
#ifndef D_NOTHREADS
  pthread_t *rgThreadID;
#else
  numThreads = 1;
#endif


  if(argc != 14){
    fprintf(stderr, "usage: %s <dataset_path> <first_training_num> <last_training_num> <first_test_num> <last_test_num> <lengthPenalty=0.> <numThreads=-1> <meshSpacingStatic=-1> <numRefinesStatic=-1> <meshDiv=4> <bandWidthDP=15> <slowPassN=10> <output_file>\n",
	    argv[0]);
    exit(1);
  }
  t3.start();

  sprintf(stPathIn, "%s", argv[1]);
  trainFirst = atoi(argv[2]);
  trainLast = atoi(argv[3]);
  testFirst = atoi(argv[4]);
  testLast = atoi(argv[5]);
  lengthPenalty =  atof(argv[6]);
  numThreads = atoi(argv[7]);
  meshSpacingStatic = atoi(argv[8]);
  numRefinesStatic = atoi(argv[9]);
  meshDiv = atof(argv[10]);
  bandWidthDP = atoi(argv[11]);
  slowPassN = atoi(argv[12]);
  sprintf(stOutfile, "%s", argv[13]);

  //Validate input values.
  
  if((meshDiv < 1) || (meshDiv > 50)){
    fprintf(stderr,"meshDiv should be 1.0 - 50.0 (default=4). was %lf\n",
	    meshDiv);
    exit(1);
  }

  if((slowPassN<-1)||(slowPassN>1000)){
    fprintf(stderr,"check slowPassN. -1 means don't do fast pass. 0 means only do fast pass. (N>0 means do fast pass for everything, then do best-N in slow pass\n");
    exit(1);
  }
  if((meshSpacingStatic==0) ||
     ((meshSpacingStatic!=-1)&&(meshSpacingStatic>200))){
    fprintf(stderr,"meshSpacingStatic expected to be 1-200 or -1! (was %d)\n",
	    meshSpacingStatic);
    exit(1);
  }
  if((numRefinesStatic < -1) || (numRefinesStatic > 4)){
    fprintf(stderr,"numRefinesStatic expected to be 0-4 or -1! (was %d)\n",
	    numRefinesStatic);
    exit(1);
  }

  if((numThreads!=-1) && ((numThreads<1)||(numThreads>64))){
    fprintf(stderr,"numThreads(%d) should be 1-64, or -1 to use all processors\n",numThreads);
    exit(1);
  }
  if((bandWidthDP < 4) || (bandWidthDP > 250)){
    fprintf(stderr,"bandWidthDP should be in range 10 to 250 (was %d)\n",
	    bandWidthDP);
    exit(1);
  }
  if(numThreads < 1){
#ifndef D_NOTHREADS
    numThreads = getNumCPUs();
#else
    numThreads = 1;
#endif
  }

#ifndef D_NOTHREADS
  rgThreadID = new pthread_t[numThreads];
  D_CHECKPTR(rgThreadID);
#else
  numThreads = 1;
#endif

  rgThreadParms = (WORDWARP_THREAD_PARMS *)malloc(sizeof(WORDWARP_THREAD_PARMS)*
						  numThreads);
  D_CHECKPTR(rgThreadParms);
  rgNumLexReductionWordsSkipped = new int[numThreads];
  D_CHECKPTR(rgNumLexReductionWordsSkipped);
  
  numTrain = trainLast - trainFirst + 1;
  numTest = testLast - testFirst + 1;

  if(numTrain < 1){
    fprintf(stderr, "numTrain must be greater than 0\n");
    exit(1);
  }
  if(numTest < 1){
    fprintf(stderr, "numTest must be greater than 0\n");
    exit(1);
  }
#if DO_FAST_PASS_FIRST
  if(slowPassN >= numTrain){
    fprintf(stderr,"slowPassN must be greater than numTrain\n");
    exit(1);
  }
#endif

  rgTrainingImages = new DImage[numTrain];
  D_CHECKPTR(rgTrainingImages);
  rgLabelsTrain = new std::string[numTrain];
  D_CHECKPTR(rgLabelsTrain);
  rgLabelsTest = new std::string[numTest];
  D_CHECKPTR(rgLabelsTest);
  rgThresholdsTrain = (int*)malloc(sizeof(int)*numTrain);
  D_CHECKPTR(rgThresholdsTrain);
  rgPagesTrain = new int[numTrain];
  D_CHECKPTR(rgPagesTrain);
  rgPagesTest = new int[numTest];
  D_CHECKPTR(rgPagesTest);
  rgAuthorIdsTrain = new int[numTrain];
  D_CHECKPTR(rgAuthorIdsTrain);
  rgAuthorIdsTest = new int[numTest];
  D_CHECKPTR(rgAuthorIdsTest);
  rgCostsMorph = (double*)malloc(sizeof(double)*(long)numTrain*(long)numTest);
  D_CHECKPTR(rgCostsMorph);
  rgCostsDP = (double*)malloc(sizeof(double)*(long)numTrain*(long)numTest);
  D_CHECKPTR(rgCostsDP);


  
  t1.start();
  int maxTrainWidth = 0;
  int maxTrainHeight = 0;
  /////////////////////////////////////
  printf("loading training data (%d words)...\n", numTrain);
  fflush(stdout);

  int heightHist[1000];
  for(int hh=0; hh < 1000; ++hh)
    	heightHist[hh] = 0;
  for(int tt=trainFirst, i=0; tt <= trainLast; ++tt,++i){
    	sprintf(stTmp,"%s/w_%08d.pgm",stPathIn,tt);
    	if(!rgTrainingImages[i].load(stTmp)){
      	fprintf(stderr,"couldn't load training image '%s'\n",stTmp);
      	exit(1);
    	}

    	if(rgTrainingImages[i].getNumProperties() >= 4){
      	std::string strTmp;
      	strTmp = rgTrainingImages[i].getPropertyVal(std::string("threshold"));
      	if(strTmp.size()<1){
			fprintf(stderr,"couldn't find threshold property in training img %d\n",
				i);
			exit(1);
      	}
      	rgThresholdsTrain[i] = atoi(strTmp.c_str());

      	strTmp = rgTrainingImages[i].getPropertyVal(std::string("label"));
      	if(strTmp.size()<1){
			fprintf(stderr,"couldn't find label property in training img %d\n",i);
			exit(1);
     	}
      	rgLabelsTrain[i] = strTmp;
      
      	strTmp = rgTrainingImages[i].getPropertyVal(std::string("page"));
      	if(strTmp.size()<1){
			fprintf(stderr,"couldn't find page property in training img %d\n",i);
			exit(1);
      	}
      	rgPagesTrain[i] = atoi(strTmp.c_str());

      	strTmp = rgTrainingImages[i].getPropertyVal(std::string("authorId"));
      	if(strTmp.size()<1){
			fprintf(stderr,"couldn't find authorId property in training img %d\n",i);
			exit(1);
      	}
      	rgAuthorIdsTrain[i] = atoi(strTmp.c_str());
    	}
    	else{
     	if(3 != rgTrainingImages[i].getNumComments()){
			fprintf(stderr, "comment in image '%s' should be: #threshval\\n#label\\n#pageNum\\n",
				stTmp);
		exit(1);
      	}
      	rgThresholdsTrain[i] =
			atoi(rgTrainingIimgTestmages[i].getCommentByIndex(0).c_str());
      	rgLabelsTrain[i] = rgTrainingImages[i].getCommentByIndex(1);
    	}
    	DThresholder::threshImage_(rgTrainingImages[i],rgTrainingImages[i],
     			       rgThresholdsTrain[i]);
    	if(rgTrainingImages[i].width() > maxTrainWidth)
     	maxTrainWidth = rgTrainingImages[i].width();
    	if(rgTrainingImages[i].height() > maxTrainHeight)
      	maxTrainHeight = rgTrainingImages[i].height();

    	DImage imgTmp;
    	// char stTmp2[1025];
    	// sprintf(stTmp2,"/tmp/clip/noclip%04d.pgm",tt);
    	// rgTrainingImages[i].save(stTmp2);
    	imgTmp = rgTrainingImages[i];
    	imgTmp.copyProperties(rgTrainingImages[i]);
    	imgTmp.copyComments(rgTrainingImages[i]);
    	// rgTrainingImages[i] = clipWordImageToInk(imgTmp);
    	// rgTrainingImages[i].copyProperties(imgTmp);
    	// rgTrainingImages[i].copyComments(imgTmp);
    	int hh;
    	hh = rgTrainingImages[i].height();
    	if(hh >= 1000)
      	hh = 1000;
    	++(heightHist[hh]);
    	// sprintf(stTmp2,"/tmp/clip/clip%04d.pgm",tt);
    	// rgTrainingImages[i].save(stTmp2);
  }
  
  /////////////////////////////////////////////////
  t1.stop();
  printf("took %.02f seconds\n", t1.getAccumulated());



#if DO_FAST_PASS_FIRST
  int *rgFastPassTopNidxs;
  double *rgSlowPassTopNMorphCosts;
  MORPHCOST_T *rgSlowPassMORPHCOST_T = NULL;
  if(slowPassN > 0){
    rgFastPassTopNidxs = new int[slowPassN];
    D_CHECKPTR(rgFastPassTopNidxs);
    rgSlowPassTopNMorphCosts = new double[slowPassN];
    D_CHECKPTR(rgSlowPassTopNMorphCosts);
    rgSlowPassMORPHCOST_T = new MORPHCOST_T[numTrain];
    D_CHECKPTR(rgSlowPassMORPHCOST_T);
  }
#endif


  t1.start();
  printf("comparing test images to training images...\n");fflush(stdout);
  //////////////////////////////////////////////////////////////////////
  for(int tt=testFirst, i=0; tt <= testLast; ++tt,++i){
	    DImage imgTest;
	    int tval;
	    DTimer t2;
	    
	    t2.start();
	    if(0==(i%1))
			printf(" image %d of %d (%d%%)\n",i,numTest, 100*(i+1)/numTest);
	    		sprintf(stTmp,"%s/w_%08d.pgm",stPathIn,tt);
	    if(!imgTest.load(stTmp)){
		 	fprintf(stderr,"couldn't load test image '%s'\n",stTmp);
		 	exit(1);
	    }

	    if(imgTest.getNumProperties() >= 4){
		 	std::string strTmp;
		 	strTmp = imgTest.getPropertyVal(std::string("threshold"));
		 	if(strTmp.size()<1){
				fprintf(stderr,"couldn't find threshold property in training img %d\n",i);
				exit(1);
		 	}
		 	tval = atoi(strTmp.c_str());
		 	strTmp = imgTest.getPropertyVal(std::string("label"));
		 	if(strTmp.size()<1){
				fprintf(stderr,"couldn't find label property in training img %d\n",i);
				exit(1);
		 	}
		 	rgLabelsTest[i] = strTmp;

		 	strTmp = imgTest.getPropertyVal(std::string("page"));
		 	if(strTmp.size()<1){
				fprintf(stderr,"couldn't find page property in training img %d\n",i);
				exit(1);
		 	}
		 	rgPagesTest[i] = atoi(strTmp.c_str());

		 	strTmp = imgTest.getPropertyVal(std::string("authorId"));
		 	if(strTmp.size()<1){
				fprintf(stderr,"couldn't find authorId property in training img %d\n",i);
				exit(1);
		 	}
		 	rgAuthorIdsTest[i] = atoi(strTmp.c_str());


	    }
	    else{
		 	if(3 != imgTest.getNumComments()){
				fprintf(stderr, "comment in image '%s' should be: #threshval\\n#label\\n#pageNum\\n",stTmp);
				exit(1);
		 	}
		 	tval = atoi(imgTest.getCommentByIndex(0).c_str());
		 	rgLabelsTest[i] = imgTest.getCommentByIndex(1);
	    }
	    //    printf(" tval=%d label=%s\n",tval, rgLabelsTest[i].c_str());
	    DThresholder::threshImage_(imgTest,imgTest, tval);
	    //DThresholder::otsuThreshImage_(imgTest,imgTest);





	    //now compare to all training values
	    for(int tnum=numThreads-1; tnum >=0; --tnum){//launch threads reverse order
		 rgThreadParms[tnum].numThreads = numThreads;
		 rgThreadParms[tnum].threadNum = tnum;
		 rgThreadParms[tnum].testWordIdx = i;
		 rgThreadParms[tnum].rgNumLexReductionWordsSkipped =
		rgNumLexReductionWordsSkipped;
		 rgThreadParms[tnum].pimgTest = &imgTest;
		 rgThreadParms[tnum].rgTrainingImages = rgTrainingImages;
		 rgThreadParms[tnum].numTrain = numTrain;
		 rgThreadParms[tnum].rgCostsMorph = rgCostsMorph;
		 rgThreadParms[tnum].rgCostsDP = rgCostsDP;
		 rgThreadParms[tnum].lengthPenalty = lengthPenalty;
		 rgThreadParms[tnum].meshSpacingStatic = meshSpacingStatic;
		 rgThreadParms[tnum].numRefinesStatic = numRefinesStatic;
		 rgThreadParms[tnum].meshDiv = meshDiv;
		 rgThreadParms[tnum].bandWidthDP = bandWidthDP;
#if DO_FAST_PASS_FIRST
		 rgThreadParms[tnum].fFastPass = true;
		 if(slowPassN==-1)
		rgThreadParms[tnum].fFastPass = false;
#endif
#ifdef D_NOTHREADS
		 word_morphing_thread_func(rgThreadParms);
#else
		 if(0 == tnum){//don't spawn thread zero. use the current thread.
		word_morphing_thread_func(rgThreadParms);
		 }
		 else{//spawn all other threads besides zero
		if(0 !=pthread_create(&rgThreadID[tnum], NULL,
					 word_morphing_thread_func, &rgThreadParms[tnum])){
		  fprintf(stderr,"failed to spawn thread #%d. Exiting.\n",tnum);
		  exit(1);
		}
		 }
#endif
	    }
#ifndef D_NOTHREADS
	    // wait for all threads to finish
	    for(int tnum = 1; tnum < numThreads; ++tnum){
		 if(pthread_join(rgThreadID[tnum],NULL)){
		fprintf(stderr, "Thread #%d failed to join. Exiting.\n", tnum);
		exit(1);
		 }
	    }
#endif
	    int numLexReductionWordsSkipped;
	    numLexReductionWordsSkipped = 0;
	    for(int tnum = 0; tnum < numThreads; ++tnum){
		 numLexReductionWordsSkipped += rgNumLexReductionWordsSkipped[tnum];
	    }

#if DO_FAST_PASS_FIRST
	    //figure out which are the top N from the fast pass
	    DMorphInk mobj;
	    if(slowPassN > 0){
			for(int fpi=0; fpi < numTrain; ++fpi){
				rgSlowPassMORPHCOST_T[fpi].wordIdx=fpi;
				rgSlowPassMORPHCOST_T[fpi].morphCost=rgCostsMorph[i*(long)numTrain+fpi];
		 	}
		 	qsort((void*)rgSlowPassMORPHCOST_T, numTrain, sizeof(MORPHCOST_T),
		    		compareMORPHCOST_T);
		 	//right now we are just doing the slow pass sequentially
		 	if(slowPassN > 100){//if it's more than a few, thread this
				fprintf(stderr,"You should probably thread this part of the code\n");
				exit(1);
			 }
			 else{
				for(int fpi=0; fpi < slowPassN; ++fpi){
		  			int trIdx;
		  			double morphCostNew;
		 			 trIdx = rgSlowPassMORPHCOST_T[fpi].wordIdx;
		 			 morphCostNew =
		 				mobj.getWordMorphCost(imgTest,
						rgTrainingImages[trIdx],
					  	bandWidthDP,
					  	0./*nonDiagonalCostDP*/,
					  	meshSpacingStatic,
					  	numRefinesStatic,
					  	meshDiv,
					  	lengthPenalty);
		 			 if(morphCostNew < rgCostsMorph[i*(long)numTrain+trIdx])
		    				rgCostsMorph[i*(long)numTrain+trIdx] = morphCostNew;
				}//end for(fpi...
		 	}//end else
	    
	    }//end if (slowPassN > 0)


#endif
	    t2.stop();
	    printf("took %.02f seconds  numLexReductionWordsSkipped=%d\n", t2.getAccumulated(), numLexReductionWordsSkipped);fflush(stdout);
  }
  /////////////////////////////////////////////////////
  t1.stop();
  printf("took %.02f seconds\n", t1.getAccumulated());

  delete [] rgNumLexReductionWordsSkipped;


  //now output the costs for analysis
  FILE *fout;
  printf("saving results to '%s'\n",stOutfile);
  fout = fopen(stOutfile, "wb");
  if(!fout){
    fprintf(stderr, "couldn't open output file '%s\n",stOutfile);
    exit(1);
  }
  fprintf(fout,"#%s trained on w_%08d.pgm through w_%08d.pgm\n",stPathIn,
	  trainFirst, trainLast);
  t3.stop();
  fprintf(fout,"walltime: %.2lf\n",t3.getAccumulated());
  fprintf(fout,"trainFirst: %08d trainLast: %08d\n",trainFirst,trainLast);
  fprintf(fout,"testFirst: %08d testLast: %08d\n",testFirst,testLast);
  for(int tt=testFirst, i=0; tt <= testLast; ++tt,++i){
    fprintf(fout,"#%d ------\n#%s\n",tt, rgLabelsTest[i].c_str());
    double bestMorphCost, bestDPCost;
    int bestMorphIdx, bestDPIdx;
    //bestMorphCost = rgCostsMorph[0];// This is wrong! (but just for comment)
    bestMorphCost = rgCostsMorph[i*(long)numTrain+0];//corrected
    bestDPCost = rgCostsDP[i*(long)numTrain+0];//corrected
    bestMorphIdx = 0;
    bestDPIdx = 0;
    for(int tr=0; tr < numTrain; ++tr){
      if(rgCostsMorph[i*(long)numTrain+tr] < bestMorphCost){
	bestMorphCost = rgCostsMorph[i*(long)numTrain+tr];
	bestMorphIdx = tr;
      }
      if(rgCostsDP[i*(long)numTrain+tr] < bestDPCost){
	bestDPCost = rgCostsDP[i*(long)numTrain+tr];
	bestDPIdx = tr;
      }
    }
    fprintf(fout,"# bestMorph=w_%08d.pgm[%s]%f   bestDP=w_%08d.pgm[%s]%f\n",
	   bestMorphIdx+trainFirst, rgLabelsTrain[bestMorphIdx].c_str(),
	   bestMorphCost,
	   bestDPIdx+trainFirst, rgLabelsTrain[bestDPIdx].c_str(),bestDPCost);
  }
  for(int tt=testFirst, i=0; tt <= testLast; ++tt,++i){
    fprintf(fout,"%d :\n",tt);
    for(int tr=0; tr < numTrain; ++tr){
      fprintf(fout,"   %f %f\n", rgCostsMorph[i*(long)numTrain+tr],
	      rgCostsDP[i*(long)numTrain+tr]);
    }
    fprintf(fout,"\n");
  }
  fprintf(fout,"#training labels:\n");
  for(int tr=0; tr < numTrain; ++tr){
    fprintf(fout,"#%d\n%s\n",tr+trainFirst, rgLabelsTrain[tr].c_str());
  }
  fclose(fout);

#ifndef D_NOTHREADS
  delete [] rgThreadID;
#endif
  free(rgThreadParms);
  free(rgCostsMorph);
  free(rgCostsDP);
  free(rgThresholdsTrain);
  delete [] rgTrainingImages;
  delete [] rgLabelsTrain;
  delete [] rgLabelsTest;
  delete [] rgPagesTrain;
  delete [] rgPagesTest;
  delete [] rgAuthorIdsTrain;
  delete [] rgAuthorIdsTest;

  return 0;
}
