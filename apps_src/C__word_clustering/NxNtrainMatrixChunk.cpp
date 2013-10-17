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

//#define D_NOTHREADS

#ifndef D_NOTHREADS
#include "dthreads.h"
#endif

#define USE_FAST_PASS_FOR_NxN_TRAINING 1

//n is numTrain, maxCost is th value mapped to 255 (anything > is clipped)
void saveMatrixImage(char *stFileName, double *rgMatrix, int n,
		     double maxCost){
  FILE *fout;
  fout = fopen(stFileName,"wb");
  if(!fout){
    fprintf(stderr,"saveMatrixImage() couldn't open '%s' for writing\n",
	    stFileName);
    exit(1);
  }
  fprintf(fout,"P5\n%d %d\n255\n",n,n);
  double scale;
  scale = 255./maxCost;
  for(int r=0,idx=0; r < n; ++r){
    for(int c=0; c < n; ++c,++idx){
      int val;
      val = (int)(scale*rgMatrix[r*n+c]);
      if(val > 255)
	val = 255;
      fprintf(fout,"%c",(char)val);
    }
  }
  fclose(fout);
}

typedef struct{
  int numThreads;
  int threadNum;
  DImage *rgTrainingImages;//pointer shared by all threads
  int numTrain;
  int chunkFirst;
  int chunkLast;
  double *rgTrainCostMatrix;
  double weightMovement; // 0. to 1. (how much to weight the movement in cost)
  double lengthPenalty;
  int meshSpacingStatic;//-1 if auto-calculate like originally done
  int numRefinesStatic;//-1 if auto-calculate like originally done
  double meshDiv;//what to divide word image height by to get mesh size [4.0]
  int bandWidthDP;//sakeo-chiba bandwidth for DP
} TRAIN_NXN_THREAD_PARMS;



void* word_morphing_NxN_training_thread_func(void *params){
  TRAIN_NXN_THREAD_PARMS *pparms;
  int numTrain;
  DMorphInk mobj;
  double *rgCosts;
  DTimer t1;
  t1.start();
  int numRows;

  pparms = (TRAIN_NXN_THREAD_PARMS*)params;
  numRows = 1 + (pparms->chunkLast) - (pparms->chunkFirst);
  numTrain = pparms->numTrain;
  rgCosts = pparms->rgTrainCostMatrix;

  for(long int idx = (pparms->threadNum)+(long)numTrain*(long)(pparms->chunkFirst),idx2=pparms->threadNum;
      idx < (numTrain)*(long)(pparms->chunkLast+1);
      idx+=(pparms->numThreads), idx2+=(pparms->numThreads)){
    double morphCost;
    int r, c;
    r = idx / numTrain;
    c = idx % numTrain;
    if(r<c){
      // morphCost = mobj.getWordMorphCost(pparms->rgTrainingImages[r],
      // 					pparms->rgTrainingImages[c],
      // 					15,/*bandWidthDP*/
      // 					0./*nonDiagonalCostDP*/);
#if USE_FAST_PASS_FOR_NxN_TRAINING
      morphCost = 
	mobj.getWordMorphCostFast(pparms->rgTrainingImages[r],
				  pparms->rgTrainingImages[c],
				  pparms->bandWidthDP,
				  0./*nonDiagonalCostDP*/,
				  pparms->meshSpacingStatic,
				  pparms->numRefinesStatic,
				  pparms->meshDiv,
				  pparms->lengthPenalty);
#else
      morphCost = 
	mobj.getWordMorphCost(pparms->rgTrainingImages[r],
			      pparms->rgTrainingImages[c],
			      pparms->bandWidthDP,/*bandWidthDP*/
			      0./*nonDiagonalCostDP*/,
			      pparms->meshSpacingStatic,
			      pparms->numRefinesStatic,
			      pparms->meshDiv,
			      pparms->lengthPenalty);
#endif
      rgCosts[idx2] = morphCost;
    }
    else if(r==c){
      rgCosts[idx2] = 0.;
    }
    else if(r>c){
      rgCosts[idx2] = 999998.;
    }
    if(0==c){
      if(0==(r%1)){
	double pctRowsComplete;
	pctRowsComplete = 100.*(r-(pparms->chunkFirst)) / (double)numRows;
	printf(" NxN %.2lf%% complete (%.2lf seconds have elapsed total)\n",
	       pctRowsComplete, t1.getAccumulated());fflush(stdout);
      }
    }
  }
  return NULL;
}

int main(int argc, char **argv);

int main(int argc, char **argv){
  char stPathIn[1025];
  char stOutfile[1025];
  char stTrainCostMatrix[1025];//file to load or write
  char stTmp[1025];
  int trainFirst, trainLast;
  int chunkFirst, chunkLast;
  int numTrain, numChunk;
  DImage *rgTrainingImages;
  double *rgTrainCostMatrix;
  DTimer t1;
  std::string *rgLabelsTrain;
  TRAIN_NXN_THREAD_PARMS *rgTrainThreadParms;
  int numThreads = 1;
  double weightMovement = 0.;
  double lengthPenalty = 0.;
  int meshSpacingStatic;
  int numRefinesStatic;
  double meshDiv;
  DMorphInk mobj;
  int bandWidth = 15;
#ifndef D_NOTHREADS
  pthread_t *rgThreadID;
#else
  numThreads = 1;
#endif


  if(argc < 15){
    fprintf(stderr, "usage: %s <dataset_path> <trainCostMatrix> <first_training_num> <last_training_num> <first_chunk_num> <last_chunk_num> <weightMovement=0.> <lengthPenalty=0.> <numThreads=-1> <meshSpacingStatic=-1> <numRefinesStatic=-1> <meshDiv=4> <SakoeChibaBandwidth=15> <output_file>\n",
	    argv[0]);
    exit(1);
  }

  sprintf(stPathIn, "%s", argv[1]);
  sprintf(stTrainCostMatrix, "%s", argv[2]);
  trainFirst = atoi(argv[3]);
  trainLast = atoi(argv[4]);
  chunkFirst = atoi(argv[5]);
  chunkLast = atoi(argv[6]);
  weightMovement = atof(argv[7]);
  lengthPenalty =  atof(argv[8]);
  numThreads = atoi(argv[9]);
  meshSpacingStatic = atoi(argv[10]);
  numRefinesStatic = atoi(argv[11]);
  meshDiv = atof(argv[12]);
  bandWidth = atoi(argv[13]);
  sprintf(stOutfile, "%s", argv[14]);

  if((meshDiv < 1) || (meshDiv > 50)){
    fprintf(stderr,"meshDiv should be 1.0 - 50.0 (default=4). was %lf\n",
	    meshDiv);
    exit(1);
  }

  if((bandWidth < 2) || (bandWidth > 50)){
    fprintf(stderr,"bandWidth should be 2-50 (default=15). was %d\n",
	    bandWidth);
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
  if((weightMovement < 0.) || (weightMovement > 1.)){
    fprintf(stderr,"weightMovement should be 0. to 1. (was %lf)\n",
	    weightMovement);
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

  rgTrainThreadParms =
    (TRAIN_NXN_THREAD_PARMS*)malloc(sizeof(TRAIN_NXN_THREAD_PARMS)*numThreads);
  D_CHECKPTR(rgTrainThreadParms);
  
  numTrain = trainLast - trainFirst + 1;
  numChunk = chunkLast - chunkFirst + 1;

  if(numTrain < 1){
    fprintf(stderr, "numTrain must be greater than 0\n");
    exit(1);
  }
  if(numChunk < 1){
    fprintf(stderr, "numChunk must be greater than 0\n");
    exit(1);
  }

  rgTrainingImages = new DImage[numTrain];
  D_CHECKPTR(rgTrainingImages);
  rgLabelsTrain = new std::string[numTrain];
  D_CHECKPTR(rgLabelsTrain);

  
  t1.start();
  int maxTrainWidth = 0;
  int maxTrainHeight = 0;
  printf("loading training data....\n");

  for(int tt=trainFirst, i=0; tt <= trainLast; ++tt,++i){
    // sprintf(stTmp,"%s/w_%08d.pgm",stPathIn,tt);
    sprintf(stTmp,"%s/thresh_w_%08d.pgm",stPathIn,tt);
    if(!rgTrainingImages[i].load(stTmp)){
      fprintf(stderr,"couldn't load training image '%s'\n",stTmp);
      exit(1);
    }
    std::string strTmp;
    strTmp = rgTrainingImages[i].getPropertyVal(std::string("label"));
    if(strTmp.size()<1){//couldn't find label property, try comments
      if(3 != rgTrainingImages[i].getNumComments()){
	fprintf(stderr, "image '%s' needs 'label' property or else comments should be: #threshval\\n#label\\n#pageNum\\n",stTmp);
	exit(1);
      }
      rgLabelsTrain[i] = rgTrainingImages[i].getCommentByIndex(1);
    }
    else
      rgLabelsTrain[i] = strTmp;
    strTmp =
      rgTrainingImages[i].getPropertyVal(std::string("textlineExtractionOK"));
    if(strTmp.size()>0){
      if(0 == strcmp("0",strTmp.c_str())){
	fprintf(stderr,"textlineExtractionOK=0 for training image '%s'\n",
		stTmp);
	exit(1);
      }
    }
    if(rgTrainingImages[i].width() > maxTrainWidth)
      maxTrainWidth = rgTrainingImages[i].width();
    if(rgTrainingImages[i].height() > maxTrainHeight)
      maxTrainHeight = rgTrainingImages[i].height();
  }
  

  t1.stop();
  printf("took %.02f seconds\n", t1.getAccumulated());

  printf("numTrain=%d\n",numTrain);

  rgTrainCostMatrix = new double[(long)numTrain*(long)numChunk];
  D_CHECKPTR(rgTrainCostMatrix);
  for(long int i=0; i < numTrain*(long)numChunk; ++i)
    rgTrainCostMatrix[i] = -42.0;

  // for(int r=0,idx=0; r < numChunk; ++r){
  //   for(int c=0; c < numTrain; ++c, ++idx){
  //     printf(" %9.2lf",rgTrainCostMatrix[idx]);
  //   }
  //   printf("\n");
  // }


  t1.start();
  printf("doing NxN comparison of training words\n");
  for(int tnum=numThreads-1; tnum >=0; --tnum){//launch threads reverse order
    rgTrainThreadParms[tnum].numThreads = numThreads;
    rgTrainThreadParms[tnum].threadNum = tnum;
    rgTrainThreadParms[tnum].rgTrainingImages = rgTrainingImages;
    rgTrainThreadParms[tnum].numTrain = numTrain;
    rgTrainThreadParms[tnum].rgTrainCostMatrix = rgTrainCostMatrix;
    rgTrainThreadParms[tnum].weightMovement = weightMovement;
    rgTrainThreadParms[tnum].lengthPenalty = lengthPenalty;
    rgTrainThreadParms[tnum].meshSpacingStatic = meshSpacingStatic;
    rgTrainThreadParms[tnum].numRefinesStatic = numRefinesStatic;
    rgTrainThreadParms[tnum].meshDiv = meshDiv;
    rgTrainThreadParms[tnum].bandWidthDP = bandWidth;
    rgTrainThreadParms[tnum].chunkFirst = chunkFirst;
    rgTrainThreadParms[tnum].chunkLast = chunkLast;
#ifdef D_NOTHREADS
    word_morphing_NxN_training_thread_func(rgTrainThreadParms);
#else
    if(0 == tnum){//don't spawn thread zero. use the current thread.
      word_morphing_NxN_training_thread_func(rgTrainThreadParms);
    }
    else{//spawn all other threads besides zero
      if(0 !=pthread_create(&rgThreadID[tnum], NULL,
			    word_morphing_NxN_training_thread_func,
			    &rgTrainThreadParms[tnum])){
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
  

  double elapsedTime;
  t1.stop();
  elapsedTime = t1.getAccumulated();
  printf("NxM comparison of training data took %.02f seconds\n", t1.getAccumulated());
  //  int numCompares = (numTrain*numChunk-numChunk-1)/2;
    // printf("(n*n-n)/2=%d comparisons = %.6lf sec per compare\n", numCompares,
    // 	   t1.getAccumulated()/numCompares);

    // now save the cost matrix
  t1.start();
  printf("SAVING matrix chunk file '%s'\n",stTrainCostMatrix);
  FILE *ftrin;
  ftrin = fopen(stTrainCostMatrix,"wb");
  if(!ftrin){
    fprintf(stderr,"couldn't open '%s' to save training cost matrix\n",
	    stTrainCostMatrix);
    exit(1);
  }
  double rgHeader[128];//firstTrain,lastTrain, then whatever text
  char *stHeaderText;
  char stText[4096];
  rgHeader[0] = trainFirst;
  rgHeader[1] = trainLast;
  rgHeader[2] = chunkFirst;
  rgHeader[3] = chunkLast;
  stHeaderText = (char*)&(rgHeader[4]);
  sprintf(stText,"\n%d compares in %.2lf seconds\n"
	  "trainFirst=%d trainLast=%d\n"
	  "weightMovement=%.2lf lengthPenalty=%.2lf meshSpacingStatic=%d numRefinesStatic=%d meshDiv=%.2lf bandWidth=%d\ndir=%s\nchunkFirst=%d chunkLast=%d\n",
	  -1, elapsedTime, trainFirst, trainLast, weightMovement,
	  lengthPenalty, meshSpacingStatic, numRefinesStatic,
	  meshDiv, bandWidth, stPathIn, chunkFirst, chunkLast);
  for(int i=4; i < 128; ++i)
    rgHeader[i] = 0.;
  strncpy(stHeaderText,stText,124*sizeof(double));
  stHeaderText[124*sizeof(double)-1]='\0';
  if(128 != fwrite(rgHeader,sizeof(double),128,ftrin)){
    fprintf(stderr,"couldn't fwrite header to '%s'\n",stTrainCostMatrix);
    exit(1);
  }
  if(numTrain*(long)numChunk != (long)
     fwrite(rgTrainCostMatrix, sizeof(double),numTrain*(long)numChunk,
		 ftrin)){
    fprintf(stderr,"couldn't fwrite the cost matrix block to '%s'\n",
	    stTrainCostMatrix);
    exit(1);
  }
  fclose(ftrin);
  t1.stop();
  printf("finished writing the cost matrix (took %.02f seconds)\n",
	 t1.getAccumulated());

  // for(int r=0,idx=0; r < numChunk; ++r){
  //   for(int c=0; c < numTrain; ++c, ++idx){
  //     printf(" %9.2lf",rgTrainCostMatrix[idx]);
  //   }
  //   printf("\n");
  // }

  // //debug: print the cost matrix
#if 0
  for(int r=0,idx=0;r < numTrain; ++r){
    for(int c=0; c < numTrain; ++c,++idx){
      printf(" %9.2lf",rgTrainCostMatrix[idx]);
    }
    printf("\n");
  }
#endif  

  
  // saveMatrixImage("/tmp/costs.pgm", rgTrainCostMatrix, numTrain, 25);


#ifndef D_NOTHREADS
  delete [] rgThreadID;
#endif
  free(rgTrainThreadParms);
  delete [] rgTrainingImages;
  delete [] rgLabelsTrain;
  delete [] rgTrainCostMatrix;
  return 0;
}
