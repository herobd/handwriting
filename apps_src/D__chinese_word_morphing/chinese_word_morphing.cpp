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
#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>

#ifndef D_NOTHREADS
#include "dthreads.h"
#endif

#define SCALE_IMAGES 0

#define ONLY_COARSE_ALIGNMENT 0

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
  WORDWARP_THREAD_PARMS *pparms;
  int numTrain;
  DMorphInk mobj;
  int testWordIdx;
#if ONLY_COARSE_ALIGNMENT 
  mobj.fOnlyDoCoarseAlignment = true;
#endif
  pparms = (WORDWARP_THREAD_PARMS*)params;
  numTrain = pparms->numTrain;
  testWordIdx = pparms->testWordIdx;

  //now compare to all training values
  for(int tr=(pparms->threadNum); tr < numTrain; tr+=(pparms->numThreads)){
    double morphCost, DPcost;

#if DO_FAST_PASS_FIRST
    if(pparms->fFastPass)
      morphCost =
	mobj.getWordMorphCostFast(*(pparms->pimgTest),
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
    pparms->rgCostsMorph[tr] = morphCost;
    pparms->rgCostsDP[tr] = DPcost;
  }
  return NULL;
}


#if DO_FAST_PASS_FIRST
void* word_morphing_thread_func_slow_pass(void *params){
  WORDWARP_THREAD_PARMS *pparms;
  int numTrain;
  DMorphInk mobj;
  int testWordIdx;

#if ONLY_COARSE_ALIGNMENT 
  mobj.fOnlyDoCoarseAlignment = true;
#endif

  pparms = (WORDWARP_THREAD_PARMS*)params;
  numTrain = pparms->numTrain;
  testWordIdx = pparms->testWordIdx;

  //now do full morphing comparison to topN training values from fast pass
  for(int tr=(pparms->threadNum); tr < (pparms->fastPassN);
      tr+=(pparms->numThreads)){
    double morphCost, DPcost;
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
      else if(0xff != p8[idx]){// not 255 or 0!
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

//startAt is the position in the array to start filling in
int loadImagesFromGnt(const char *stGntName, DImage *rgImgs, int startAt,
		      std::string *rgLabels){
  FILE *fin;
  fin = fopen(stGntName, "rb");
  if(!fin){
    fprintf(stderr,"file '%s' could not be opened for reading\n",stGntName);
    exit(1);
  }
  size_t imgNum;
  imgNum = 0;

  while(!feof(fin)){
    uint32_t blockSize32;
    uint16_t tagCode16;
    uint16_t width16;
    uint16_t height16;
    unsigned char *p8;

    if(1!=fread(&blockSize32, sizeof(uint32_t), 1, fin)){
      //we must have reached end of file
      break;
    }
    if(1!=fread(&tagCode16, sizeof(uint16_t), 1, fin)){
      fprintf(stderr,"an error occurred reading the tagcode for image #%ld\n",
	      imgNum);
      exit(1);
    }
    if(1!=fread(&width16, sizeof(uint16_t), 1, fin)){
      fprintf(stderr,"an error occurred reading the width for image #%ld\n",
	      imgNum);
      exit(1);
    }
    if(1!=fread(&height16, sizeof(uint16_t), 1, fin)){
      fprintf(stderr,"an error occurred reading the height for image #%ld\n",
	      imgNum);
      exit(1);
    }
    //printf("image %ld w=%d h=%d\n",imgNum, (int)width16, (int)height16);
    size_t datalen;
    datalen = (unsigned int)width16 * (unsigned int)height16;
    DImage *pimg;

#if SCALE_IMAGES
    DImage imgUnscaled;
    pimg = &imgUnscaled;
#else
    pimg = &rgImgs[startAt + imgNum];
#endif

    pimg->create((int)width16,(int)height16, DImage::DImage_u8, 1);
    p8 = pimg->dataPointer_u8();
    if(1!=fread(p8, datalen, 1, fin)){
      fprintf(stderr,"error reading image data\n");
      exit(1);
    }
    
#if SCALE_IMAGES
    pimg->scaled_(rgImgs[startAt + imgNum],0.75,0.75,DImage::DImageTransSmooth);
    pimg = &(rgImgs[startAt + imgNum]);
#endif

    //find threshold using otsu
    DThresholder::otsuThreshImage_(*pimg, *pimg);
    char stLabelTmp[20];
    sprintf(stLabelTmp,"%03u%03u",
	    (unsigned char)tagCode16, (unsigned char)(tagCode16>>8));
    rgLabels[startAt + imgNum] = std::string(stLabelTmp);
    
    ++imgNum;
  }
  fclose(fin);
  return imgNum;
}

std::vector<int> readFileList(const char *stFile){
  std::vector<int> vect;
  FILE *fin;
  fin = fopen(stFile,"rb");
  if(!fin){
    fprintf(stderr,"couldn't open file list '%s'\n",stFile);
    exit(1);
  }
  int gntNum;
  while(!feof(fin)){
    if(1!=fscanf(fin,"%d", &gntNum))
      break;
    vect.push_back(gntNum);
  }
  fclose(fin);
  return vect;
}


int main(int argc, char **argv);

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
  double lexiconReductionThresh = 0.;
  WORDWARP_THREAD_PARMS *rgThreadParms;
  int numThreads = 1;
  double lengthPenalty = 0.;
  int meshSpacingStatic;
  int numRefinesStatic;
  double meshDiv;
  int bandWidthDP = 15;
  int slowPassN;//numer of best matches from fast pass to do slow pass for
  int topN = 10;//number of best matches to output (for later use with Ngrams)
#ifndef D_NOTHREADS
  pthread_t *rgThreadID;
#else
  numThreads = 1;
#endif
  char stTrainFileList[2048];
  char stTestFileList[2048];
  std::vector<int> vectTrainGntFiles;
  std::vector<int> vectTestGntFiles;
  long totalCorrect = 0;
  long totalInTopN = 0;
  FILE *fout;
  
  size_t first_gnum = 0;
  int first_tti = 0;
  size_t last_gnum = 2147483647;
  int last_tti = 2147483647;

  fout = stdout;

  if((argc != 13) && (argc != 17)){
    fprintf(stderr, "usage: %s <dataset_path> <training_list>  <test_list> <lengthPenalty=0.> <numThreads=-1> <meshSpacingStatic=-1> <numRefinesStatic=-1> <meshDiv=4> <bandWidthDP=15> <slowPassN=100> <topN=10> <output_file> [<first_gnum> <first_tti> <last_gnum> <last_tti>]\n",   
	    argv[0]);
    fprintf(stderr," first_gnum: skip all .gnt files before the <first_gnum+1>th file in test list\n");
    fprintf(stderr," first_tti: skip all test words before <first_tti> in first_gnum (0-based)\n");
    fprintf(stderr," last_gnum: skip all .gnt files after the <last_gnum+1>th\n");
    fprintf(stderr," last_tti: skip all test words after <last_tti> in last_gnum (0-based)\n");
    fprintf(stderr," ***defaults: 2147483647 for last_gnum/last_tti, 0 for first_gnum/first_tti\n");
    exit(1);
  }
  t3.start();

  sprintf(stPathIn, "%s", argv[1]);
  strncpy(stTrainFileList, argv[2],2048);
  strncpy(stTestFileList, argv[3],2048);
  // trainFirst = atoi(argv[2]);
  // trainLast = atoi(argv[3]);
  // testFirst = atoi(argv[4]);
  // testLast = atoi(argv[5]);

  lengthPenalty =  atof(argv[4]);
  numThreads = atoi(argv[5]);
  meshSpacingStatic = atoi(argv[6]);
  numRefinesStatic = atoi(argv[7]);
  meshDiv = atof(argv[8]);
  bandWidthDP = atoi(argv[9]);
  slowPassN = atoi(argv[10]);
  topN = atoi(argv[11]);
  sprintf(stOutfile, "%s", argv[12]);
  if(17==argc){
    first_gnum = atoi(argv[13]);
    first_tti = atoi(argv[14]);
    last_gnum = atoi(argv[15]);
    last_tti = atoi(argv[16]);
  }

  if(topN > slowPassN){
    fprintf(stderr,"topN(%d) must be <= slowPassN(%d)\n",topN,slowPassN);
    exit(1);
  }

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
  char stOutputLog[2048];
  sprintf(stOutputLog,"%s_%d.log",stOutfile,getpid());
  fout = fopen(stOutputLog,"wb");
  if(!fout){
    fprintf(stderr,"couldn't open '%s' for log output\n",stOutputLog);
    exit(1);
  }
  fprintf(stderr,"using log file '%s'\n",stOutputLog);
  printf("using log file '%s'\n",stOutputLog);
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

  vectTrainGntFiles = readFileList(stTrainFileList);
  vectTestGntFiles = readFileList(stTestFileList);

  printf("vectTrainGntFiles.size()=%ld\n",vectTrainGntFiles.size());
  printf("vectTestGntFiles.size()=%ld\n",vectTestGntFiles.size());

  for(size_t i=0; i < vectTestGntFiles.size(); ++i)
    printf(" [%ld]=%d\n",i, vectTestGntFiles[i]);

  rgThreadParms = (WORDWARP_THREAD_PARMS *)malloc(sizeof(WORDWARP_THREAD_PARMS)*
						  numThreads);
  D_CHECKPTR(rgThreadParms);
  



#define MAX_TRAIN 2144749/*number train in hwdb1.0 + number train in hwdb1.1*/
#define MAX_TEST 10000/*max number that can be read from a single .gnt file*/
  rgTrainingImages = new DImage[MAX_TRAIN];
  D_CHECKPTR(rgTrainingImages);
  rgLabelsTrain = new std::string[MAX_TRAIN];
  D_CHECKPTR(rgLabelsTrain);

  DImage *rgTestImages = NULL;
  rgTestImages = new DImage[MAX_TEST];
  D_CHECKPTR(rgTestImages);
  rgLabelsTest = new std::string[MAX_TEST];
  D_CHECKPTR(rgLabelsTest);


  printf("loading training data (%ld gnt files)...\n",vectTrainGntFiles.size());
  fflush(stdout);
  fprintf(fout,"loading training data (%ld gnt files)...\n",
	  vectTrainGntFiles.size());
  fflush(fout);

  char stGntName[2048];
  numTrain = 0;
  for(size_t gnum=0; gnum <vectTrainGntFiles.size(); ++gnum){
    int numRead;
    if(vectTrainGntFiles[gnum] < 1000)
      sprintf(stGntName,"%s/%03d-f.gnt",stPathIn,vectTrainGntFiles[gnum]);
    else
      sprintf(stGntName,"%s/%4d-f.gnt",stPathIn,vectTrainGntFiles[gnum]);

    // printf("%d: loading '%s'\n",gnum,stGntName);
    numRead = loadImagesFromGnt(stGntName, rgTrainingImages, numTrain,
				rgLabelsTrain);
    printf("read %d training images from '%03u-f.gnt' (%ld of %ld)\n",
	   numRead,(unsigned int)vectTrainGntFiles[gnum], gnum,
	   vectTrainGntFiles.size());
    numTrain += numRead;
  }
  printf("read a total of %d training images\n",numTrain);

  rgCostsMorph = (double*)malloc(sizeof(double)*(long)numTrain);
  D_CHECKPTR(rgCostsMorph);
  rgCostsDP = (double*)malloc(sizeof(double)*(long)numTrain);
  D_CHECKPTR(rgCostsDP);

#if DO_FAST_PASS_FIRST
  int *rgFastPassTopNidxs;
  double *rgSlowPassTopNMorphCosts;
  if(slowPassN > 0){
    rgFastPassTopNidxs = new int[slowPassN];
    D_CHECKPTR(rgFastPassTopNidxs);
    rgSlowPassTopNMorphCosts = new double[slowPassN];
    D_CHECKPTR(rgSlowPassTopNMorphCosts);
  }
#endif
  //we will use this to sort results to get topN for output even if not
  //doing the fast pass first
  MORPHCOST_T *rgSlowPassMORPHCOST_T = NULL;
  rgSlowPassMORPHCOST_T = new MORPHCOST_T[numTrain];
  D_CHECKPTR(rgSlowPassMORPHCOST_T);

  printf("loading and recognizing test data (%ld gnt files)...\n", vectTestGntFiles.size());
  fflush(stdout);
 


  numTest = 0;
  for(size_t gnum=first_gnum; (gnum <=last_gnum) &&
	(gnum < vectTestGntFiles.size()); ++gnum){
    int numRead;
    sprintf(stGntName,"%s/%03d-f.gnt",stPathIn,vectTestGntFiles[gnum]);
    numRead = loadImagesFromGnt(stGntName, rgTestImages, 0, rgLabelsTest);
    printf("read %d test images from '%03d-f.gnt'\n",
	   numRead,vectTestGntFiles[gnum]);
    //numTest += numRead;

    // recognize the words from this .gnt file and report the results
    for(int tti=first_tti; (tti<=last_tti)&&(tti < numRead); ++tti){
      printf(" image %d of %d in this gnt file\n",tti, numRead);
      DTimer t2;
      t2.start();
      //now compare to all training values
      for(int tnum=numThreads-1; tnum >=0;--tnum){//launch threads reverse order
	rgThreadParms[tnum].numThreads = numThreads;
	rgThreadParms[tnum].threadNum = tnum;
	rgThreadParms[tnum].testWordIdx = tti;
	rgThreadParms[tnum].pimgTest = &rgTestImages[tti];
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
				word_morphing_thread_func,
				&rgThreadParms[tnum])){
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


#if DO_FAST_PASS_FIRST
      //figure out which are the top N from the fast pass
      DMorphInk mobj;
#if ONLY_COARSE_ALIGNMENT 
      mobj.fOnlyDoCoarseAlignment = true;
#endif
      if(slowPassN > 0){
	for(int fpi=0; fpi < numTrain; ++fpi){
	  rgSlowPassMORPHCOST_T[fpi].wordIdx=fpi;
	  rgSlowPassMORPHCOST_T[fpi].morphCost = rgCostsMorph[fpi];
	}
	qsort((void*)rgSlowPassMORPHCOST_T, numTrain, sizeof(MORPHCOST_T),
	      compareMORPHCOST_T);
	//is it faster to spawn threads or just do the top N sequentially?
	if(slowPassN > 100){//if it's more than a few, thread this
	  fprintf(stderr,"May want to thread this code\n");
	  exit(1);
	}
	else{
	  for(int fpi=0; fpi < slowPassN; ++fpi){
	    int trIdx;
	    double morphCostNew;
	    trIdx = rgSlowPassMORPHCOST_T[fpi].wordIdx;
	    morphCostNew =
	      mobj.getWordMorphCost(rgTestImages[tti],
				    rgTrainingImages[trIdx],
				    bandWidthDP,
				    0./*nonDiagonalCostDP*/,
				    meshSpacingStatic,
				    numRefinesStatic,
				    meshDiv,
				    lengthPenalty);
	    if(morphCostNew < rgCostsMorph[trIdx]){
	      rgCostsMorph[trIdx] = morphCostNew;//we don't use in this program
	      rgSlowPassMORPHCOST_T[fpi].morphCost = morphCostNew;//for re-sort
	    }
	  }//end for(fpi...
	}//end else
	//now re-sort just the 
	qsort((void*)rgSlowPassMORPHCOST_T, slowPassN, sizeof(MORPHCOST_T),
	      compareMORPHCOST_T);
      }//end if (slowPassN > 0)
#endif
      t2.stop();
      if(slowPassN<1){//we didn't already sort the results, so do it now
	for(int fpi=0; fpi < numTrain; ++fpi){
	  rgSlowPassMORPHCOST_T[fpi].wordIdx=fpi;
	  rgSlowPassMORPHCOST_T[fpi].morphCost = rgCostsMorph[fpi];
	}
	qsort((void*)rgSlowPassMORPHCOST_T, numTrain, sizeof(MORPHCOST_T),
	      compareMORPHCOST_T);
      }
      // output topN recognition results (rgSlowPassMORPHCOST is now in sorted
      // order for the topN indexes)
      bool fRecognizedCorrectly;
      bool fInTopN;
      fRecognizedCorrectly = false;
      fInTopN = false;
      if(0==strcmp(rgLabelsTest[tti].c_str(),
		   rgLabelsTrain[rgSlowPassMORPHCOST_T[0].wordIdx].c_str()))
	fRecognizedCorrectly = true;
      fprintf(fout,"#%d (gnum=%ld %d) groundtruth='%s' %s\n",numTest,gnum,tti,
	      rgLabelsTest[tti].c_str(), fRecognizedCorrectly ? "CORRECT" :
	      "WRONG!");
      fprintf(stderr,"#%d (gnum=%ld %d) groundtruth='%s' %s\n",numTest,gnum,tti,
	      rgLabelsTest[tti].c_str(), fRecognizedCorrectly ? "CORRECT" :
	      "WRONG!");
      for(int n=0; n < topN; ++n){
	fprintf(fout,"  %d trainingImg%d  cost=%lf  label='%s'\n", n,
		rgSlowPassMORPHCOST_T[n].wordIdx,
		rgSlowPassMORPHCOST_T[n].morphCost,
		rgLabelsTrain[rgSlowPassMORPHCOST_T[n].wordIdx].c_str());
	if(0==strcmp(rgLabelsTest[tti].c_str(),
		     rgLabelsTrain[rgSlowPassMORPHCOST_T[n].wordIdx].c_str()))
	  fInTopN = true;
      }
      if(fRecognizedCorrectly)
	++totalCorrect;
      if(fInTopN)
	++totalInTopN;
      ++numTest;
      fprintf(fout,"  took %.02f seconds\n",t2.getAccumulated());
      fprintf(fout,"  so far %ld/%d correct (%.2lf%%)  in top%d=%.2lf%%\n",
	      totalCorrect, numTest, 100.*totalCorrect/(double)numTest,topN,
	      100.*totalInTopN/(double)numTest);
      fprintf(stderr,"  took %.02f seconds\n",t2.getAccumulated());
      fprintf(stderr,"  so far %ld/%d correct (%.2lf%%)  in top%d=%.2lf%%\n",
	      totalCorrect, numTest, 100.*totalCorrect/(double)numTest,topN,
	      100.*totalInTopN/(double)numTest);
      fflush(fout);
      if(0==(tti%20))
	fprintf(stderr,"***note:details are in log file '%s'\n",stOutputLog);
    }//for(tti=0; tti < numRead; ++tti)
    t1.stop();
    fprintf(fout,"took %.02f seconds for gnum=%ld\n", t1.getAccumulated(),gnum);
    fflush(fout);
    fprintf(stderr,"took %.02f seconds for gnum=%ld\n",
	    t1.getAccumulated(),gnum);
  }//for(size_t gnum=0; gnum <vectTestGntFiles.size(); ++gnum)

  fclose(fout);//the output log file

  return 0;


#if 0//------------------------------------------------------------------
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
#endif//------------------------------------------------------------------

#ifndef D_NOTHREADS
  delete [] rgThreadID;
#endif
  free(rgThreadParms);
  free(rgCostsMorph);
  free(rgCostsDP);
  delete [] rgTrainingImages;
  delete [] rgLabelsTrain;
  delete [] rgLabelsTest;
  delete [] rgPagesTrain;
  delete [] rgPagesTest;
  delete [] rgAuthorIdsTrain;
  delete [] rgAuthorIdsTest;

  delete [] rgTestImages;

  delete [] rgSlowPassMORPHCOST_T;
  if(slowPassN > 0){
    delete [] rgFastPassTopNidxs;
    delete [] rgSlowPassTopNMorphCosts;
  }
  return 0;

}
