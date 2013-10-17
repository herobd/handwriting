#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

#define TEN 10
#define SAVE_IMAGES 0
#define SHOW_DP_STATS 0

typedef struct {
  std::string *label;
  double costMorph;
  double costDP;
  int trainNum;
} MATCH_S;

int compareMATCH_S_morph(const void *p1, const void *p2){
  MATCH_S *pm1, *pm2;
  pm1 = (MATCH_S*)p1;
  pm2 = (MATCH_S*)p2;
  if((pm1->costMorph) < (pm2->costMorph))
    return -1;
  else if ((pm1->costMorph) > (pm2->costMorph))
    return 1;
  return 0;
}

int compareMATCH_S_dp(const void *p1, const void *p2){
  MATCH_S *pm1, *pm2;
  pm1 = (MATCH_S*)p1;
  pm2 = (MATCH_S*)p2;
  if((pm1->costDP) < (pm2->costDP))
    return -1;
  else if ((pm1->costDP) > (pm2->costDP))
    return 1;
  return 0;
}


#define D_CHECKPTR(x) {if(!(x)){fprintf(stderr,"mem err!\n");abort();}}

int main(int argc, char **argv){
  int trainFirst, trainLast;
  int testFirst, testLast;
  int numTrain, numTest;
  std::string *rgLabelsTrain;
  std::string *rgLabelsTest;
  double *rgCostsMorph;
  double *rgCostsDP;
  bool *rgFoundData;
  FILE *fin;
  char stComment[1025];
  char stTmp[1025];
  char stOutputPath[1025];
	
  if(argc < 6){
    fprintf(stderr, "usage: %s trainFirst trainLast testFirst testLast "
	    "datafile0.dat [...]\n", argv[0]);
    exit(1);
  }
  
  trainFirst = atoi(argv[1]);
  trainLast = atoi(argv[2]);
  testFirst = atoi(argv[3]);
  testLast = atoi(argv[4]);
  
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
  
  
  rgLabelsTrain = new std::string[numTrain];
  D_CHECKPTR(rgLabelsTrain);
  rgLabelsTest = new std::string[numTest];
  D_CHECKPTR(rgLabelsTest);
  rgCostsMorph = (double*)malloc(sizeof(double)*numTrain*numTest);
  D_CHECKPTR(rgCostsMorph);
  rgCostsDP = (double*)malloc(sizeof(double)*numTrain*numTest);
  D_CHECKPTR(rgCostsDP);

  rgFoundData = (bool*)malloc(sizeof(bool)*numTest);
  D_CHECKPTR(rgFoundData);
  for(int i=0; i < numTest; ++i)
    rgFoundData[i] = false;
  double *rgWallTimes;
  double wallTimeTot = 0.;
  rgWallTimes = new double[argc-5];
  D_CHECKPTR(rgWallTimes);
  for(int iii=5; iii<argc; ++iii){//read each datafile
    printf("reading '%s'\n", argv[iii]);
    fin = fopen(argv[iii],"rb");
    if(!fin){
      fprintf(stderr, "couldn't read datafile '%s'\n",argv[iii]);
      exit(1);
    }
    char *stSkip;
    stSkip = fgets(stTmp,1024,fin);
    if(iii == 5){
      sprintf(stComment,"%s",stTmp);
      sscanf(&(stComment[1]),"%s",stOutputPath);
      if(stOutputPath[strlen(stOutputPath)-1] != '/'){
	stOutputPath[strlen(stOutputPath)] = '/';
      }
      stOutputPath[strlen(stOutputPath)]='\0';
    }
    if(0 != strcmp(stComment, stTmp)){
      fprintf(stderr, "the comment line (first line of file '%s') doesn't "
	      "match!\n  first line:'%s'\n  expected  :'%s'\n", argv[iii],
	      stTmp,stComment);
      //      exit(1);
    }
    int a, b, c, d;
    stSkip = fgets(stTmp,1024,fin);
    
    if(0==strncmp("walltime:",stTmp,strlen("walltime:"))){
      sscanf(stTmp,"%*s%lf",&(rgWallTimes[iii-5]));
      wallTimeTot += rgWallTimes[iii-5];
      stSkip = fgets(stTmp,1024,fin);
    }
    sscanf(stTmp, "%*s%d%*s%d",&a, &b);
    stSkip = fgets(stTmp,1024,fin);
    sscanf(stTmp, "%*s%d%*s%d",&c, &d);
    if((a!=trainFirst) || (b!=trainLast) || (c<testFirst) || (c>testLast) ||
       (d<testFirst) || (d > testLast)){
      fprintf(stderr, "the training or test ranges in '%s' seem incorrect\n",
	      argv[iii]);
      fprintf(stderr, "a=%d b=%d c=%d d=%d\ntrainFirst=%d trainLast=%d testFirst=%d testLast=%d\n",a,b,c,d,trainFirst,trainLast,testFirst, testLast);
      exit(1);
    }
    //read the test labels from the human-readable comments
    printf("test data (tt=%d to %d)\n",c,d);
    for(int tt=c; tt <=d; ++tt){
		 int tmptt;
		 stSkip = fgets(stTmp,1024,fin);
		 sscanf(stTmp,"%*c%d",&tmptt);
		 if(tmptt != tt){
			fprintf(stderr, "tmptt=%d tt=%d\n",tmptt, tt);
			exit(1);
		 }
		 stSkip = fgets(stTmp,1024,fin);
		 if(stTmp[strlen(stTmp)-1]=='\n')
			stTmp[strlen(stTmp)-1]='\0';
		 char *ptmp;
		 ptmp = &(stTmp[1]);
		 rgLabelsTest[tt-testFirst]=std::string(ptmp);

		 stSkip = fgets(stTmp,1024,fin);
    }
    //read the actual matching costs
    //    printf("looping to read the matching costs (tt=%d to %d)\n",c,d);
    for(int tt=c; tt <=d; ++tt){
		stSkip = fgets(stTmp,1024,fin);
     	int ttTmp;
      	sscanf(stTmp,"%d",&ttTmp);
      	if(tt!=ttTmp){
			fprintf(stderr, "*****tt!=ttTmp\n");
			exit(1);
      	}
      	for(int tr=0; tr < numTrain; ++tr){
			stSkip = fgets(stTmp,1024,fin);
			sscanf(stTmp,"%lf%lf",&(rgCostsMorph[(tt-testFirst)*numTrain+tr]),
	       		&(rgCostsDP[(tt-testFirst)*numTrain+tr]));
      	}
      	stSkip = fgets(stTmp,1024,fin);//skip empty line
      	rgFoundData[tt-testFirst] = true;
    }
    //read the training labels
    //    printf("looping to read training labels!!!!!!!!\n");
    stSkip = fgets(stTmp, 1024, fin);
    for(int tr=a; tr <= b; ++tr){
      int matchtr;
      stSkip = fgets(stTmp,1024,fin);
      sscanf(stTmp,"%*c%d",&matchtr);
      if(tr != matchtr){
	fprintf(stderr, "matchtr(%d)!=tr(%d)\n",matchtr,tr);
	exit(1);
      }
      stSkip = fgets(stTmp, 1024, fin);
      if(stTmp[strlen(stTmp)-1] == '\n')
	stTmp[strlen(stTmp)-1] = '\0';//get rid of newline at end
      if(iii==5){
	rgLabelsTrain[tr-trainFirst] = std::string(stTmp);
      }
      else{
	if(0 != strcmp(stTmp,rgLabelsTrain[tr-trainFirst].c_str())){
	  fprintf(stderr, "the training labels don't seem to match! (%s!=%s)\n",
		  stTmp,rgLabelsTrain[tr-trainFirst].c_str());
	  exit(1);
	}
      }
    }
    fclose(fin);
  }
    
  //now figure out the stats, etc.
  printf("done reading all the files\n");
  bool fNotFound = false;
  for(int i=0; i < numTest; ++i){
    if(!rgFoundData[i]){
      fprintf(stderr, "some of the test data is not found in the data files (%d)\n",i+testFirst);
      fNotFound = true;
      //exit(1);
    }
  }
  if(fNotFound)
    exit(1);
  int numOoV; // number out of vocabulary (there isn't an example in train data)
  int numCorrectDP1;//# test words matched correctly at best 1 match
  int numCorrectDP10;//# test words matched correctly at best 10 matches
  int numCorrectMorph1;//# test words matched correctly at best 1 match
  int numCorrectMorph10;//# test words matched correctly at best 10 matches
  double avgCostDifferenceToActualCorrect;
  double avgCostDifferenceToIncorrect;
  int numDifferenceToActualCorrect;
  int numDifferenceToIncorrect;

  MATCH_S *rgMatches;
  rgMatches = (MATCH_S*)malloc(sizeof(MATCH_S)*numTrain);
  D_CHECKPTR(rgMatches);


  numOoV = 0;
  numCorrectDP1 = 0;
  numCorrectDP10 = 0;
  numCorrectMorph1 = 0;
  numCorrectMorph10 = 0;
  avgCostDifferenceToActualCorrect = 0.;
  avgCostDifferenceToIncorrect = 0.;
  numDifferenceToActualCorrect = 0;
  numDifferenceToIncorrect = 0;

#if SAVE_IMAGES
  {
    system("echo \"P2\n10 1 255\n255 255 255 255 255 255 255 255 255 255\" > /tmp/wrongwords.pgm");
    system("echo \"P2\n10 1 255\n255 255 255 255 255 255 255 255 255 255\" > /tmp/space.pgm");
  }
#endif

  for(int w=0; w < numTest; ++w){
    // int rgMorphBest10idx[10];
    // double rgMorphBest10Cost[10];
    // int rgDPBest10idx[10];
    // double rgDPBest10Cost[10];
    
    //copy the matches for this test word into the array for sorting
    for(int tr=0; tr < numTrain; ++tr){
      rgMatches[tr].label = &(rgLabelsTrain[tr]);
      rgMatches[tr].costMorph = rgCostsMorph[w*numTrain+tr];
      rgMatches[tr].costDP = rgCostsDP[w*numTrain+tr];
      rgMatches[tr].trainNum = tr+trainFirst;
    }
    qsort((void*)rgMatches, numTrain, sizeof(MATCH_S), compareMATCH_S_morph);
    printf("%d[%s]",w+testFirst,rgLabelsTest[w].c_str());
    //check for OoV
    bool fOOV = true;
    for(int i=0; i < numTrain; ++i){
      if(0 == strcmp(rgMatches[i].label->c_str(), rgLabelsTest[w].c_str())){
	fOOV = false;
	break;
      }
    }
    if(fOOV){
      ++numOoV;
      printf("(OoV)");
    }
    //best match correct? (Morph)
    bool fCorrect =false;
    if(0 == strcmp(rgMatches[0].label->c_str(), rgLabelsTest[w].c_str())){
      ++numCorrectMorph1;
      fCorrect = true;
    }
    else if(!fOOV)
      printf("(wrong)");
    //show best 10 matches (Morph)
    printf("\n  best %d morph:\n", TEN);
    bool fBest10;
    fBest10 = false;
    for(int b=0; b < TEN; ++b){
      if(0 == strcmp(rgMatches[b].label->c_str(), rgLabelsTest[w].c_str()))
	fBest10 = true;
      printf("    %s cost=%lf trainNum=%d\n",
	     rgMatches[b].label->c_str(), rgMatches[b].costMorph,
	     rgMatches[b].trainNum);
    }
    if(fBest10)
      ++numCorrectMorph10;
    
    //if correct, get the difference in cost between match and best non-match
    //otherwise difference in cost between best match and best non-match
    if(fCorrect){
      for(int b=1; b < numTest; ++b){
	if(0 != strcmp(rgMatches[b].label->c_str(), rgLabelsTest[w].c_str())){
	  avgCostDifferenceToIncorrect += 
	    rgMatches[b].costMorph - rgMatches[0].costMorph;
	  ++numDifferenceToIncorrect;
	  break;
	}
      }
    }
    else if(!fOOV){
#if SAVE_IMAGES
      {
	char stImg[1025];
	char stCmd[2048];
	sprintf(stCmd,"pnmcat -lr -white %sthresh_w_%08d.pgm",
		stOutputPath,testFirst+w);
	for(int b=0; b < 4; ++b){
	  char stCmdTmp[1025];
	  sprintf(stCmdTmp," /tmp/space.pgm %sthresh_w_%08d.pgm",
		  stOutputPath,rgMatches[b].trainNum);
	  strcat(stCmd,stCmdTmp);
	}
	strcat(stCmd," > /tmp/tmpword.pgm");
	system(stCmd);
	sprintf(stCmd,"pnmcat -tb -jleft -white /tmp/wrongwords.pgm /tmp/tmpword.pgm > /tmp/wrongwords2.pgm");
	system(stCmd);
	system("mv /tmp/wrongwords2.pgm /tmp/wrongwords.pgm");
	system(stCmd);
      }
#endif


      for(int b=1; b < numTest; ++b){
	if(0 == strcmp(rgMatches[b].label->c_str(), rgLabelsTest[w].c_str())){
	  avgCostDifferenceToActualCorrect += 
	    rgMatches[b].costMorph - rgMatches[0].costMorph;
	  ++numDifferenceToActualCorrect;
	  break;
	}
      }
    }

#if SHOW_DP_STATS
    qsort((void*)rgMatches, numTrain, sizeof(MATCH_S), compareMATCH_S_dp);
    //best match correct? (DP)
    if(0 == strcmp(rgMatches[0].label->c_str(), rgLabelsTest[w].c_str()))
      ++numCorrectDP1;
    //show best 10 matches (DP)
    printf(" best %d DP:\n", TEN);
    fBest10 = false;
    for(int b=0; b < TEN; ++b){
      if(0 == strcmp(rgMatches[b].label->c_str(), rgLabelsTest[w].c_str()))
	fBest10 = true;
      printf("    %s cost=%lf trainNum=%d\n",
	     rgMatches[b].label->c_str(), rgMatches[b].costDP,
	     rgMatches[b].trainNum);
    }
    if(fBest10)
      ++numCorrectDP10;
#endif /*SHOW_DP_STATS*/
  }

  for(int i=0; i < (argc-5); ++i){
    printf("walltime[%d]=%.2lf\n",i,rgWallTimes[i]);
  }
  printf("totalWallTime=%.2lf\n",wallTimeTot);
  
  printf("---------------------------------------\n");
  printf("Training words: %d\n", numTrain);
  printf("Test words: %d\n", numTest);
  printf("OoV test words: %d (%.2f%%)\n", numOoV, 100.*numOoV/(double)numTest);
  printf("In-Vocab test words: %d (%.2f%%)\n", numTest-numOoV, 100.*(numTest-numOoV)/(double)numTest);
  printf("numCorrectMorph1: %d (%.2f%% of total, %.2f%% of non-OoV)\n",
	 numCorrectMorph1, 100.*numCorrectMorph1/(double)numTest,
	 100.*numCorrectMorph1/(double)(numTest-numOoV));
  printf("numCorrectMorph10: %d (%.2f%% of total, %.2f%% of non-OoV)\n",
	 numCorrectMorph10, 100.*numCorrectMorph10/(double)numTest,
	 100.*numCorrectMorph10/(double)(numTest-numOoV));

#if SHOW_DP_STATS
  printf("numCorrectDP1: %d (%.2f%% of total, %.2f%% of non-OoV)\n",
	 numCorrectDP1, 100.*numCorrectDP1/(double)numTest,
	 100.*numCorrectDP1/(double)(numTest-numOoV));
  printf("numCorrectDP10: %d (%.2f%% of total, %.2f%% of non-OoV)\n",
	 numCorrectDP10, 100.*numCorrectDP10/(double)numTest,
	 100.*numCorrectDP10/(double)(numTest-numOoV));
#endif /*SHOW_DP_STATS*/
  if(numDifferenceToActualCorrect > 0)
    avgCostDifferenceToActualCorrect /= numDifferenceToActualCorrect;
  if(numDifferenceToIncorrect > 0)
    avgCostDifferenceToIncorrect /= numDifferenceToIncorrect;
  printf("avgCostDifferenceToActualCorrect=%f (num=%d)\n",
	 avgCostDifferenceToActualCorrect, numDifferenceToActualCorrect);
  printf("avgCostDifferenceToIncorrect=%f (num=%d)\n",
	 avgCostDifferenceToIncorrect, numDifferenceToIncorrect);

  free(rgFoundData);
  free(rgCostsDP);
  free(rgCostsMorph);
  delete [] rgLabelsTrain;
  delete [] rgLabelsTest;
#if SAVE_IMAGES
    printf(">>>stOutputPath='%s'\n",stOutputPath);
#endif
  return 0;
}
