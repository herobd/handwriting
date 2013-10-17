#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv){
  char stOutFile[2048];
  int trainFirst = -1;
  int trainLast = -1;
  long int numTrain = -1;
  bool *rgDataRead = NULL;

  double rgHeader[128];
  char *stHeaderText;
  bool fAllOk;
  double *rgMatrix = NULL;

  if(argc < 3){
    fprintf(stderr,"usage: %s <outputFile> <chunkFile1> [<chunkFileN>...]\n",
	    argv[0]);
    exit(1);
  }

  strncpy(stOutFile,argv[1],2047);
  stOutFile[2047] = '\0';



  fAllOk = true;


  for(int ii=2; ii <argc;++ii){
    FILE *fin;

    fin=fopen(argv[ii],"rb");
    if(!fin){
      fprintf(stderr,"couldn't open chunk file '%s' for inpute\n",argv[ii]);
      fAllOk = false;
      continue;
    }

    if(128!=fread(rgHeader,sizeof(double),128,fin)){
      fprintf(stderr,"fread did not return 128 header doubles for chunk file '%s'\n", argv[ii]);
      fAllOk = false;
      fclose(fin);
      continue;
    }
    int firstTrainTmp, lastTrainTmp, chunkFirstTmp, chunkLastTmp;
    long int numTrainTmp, chunkSizeTmp;
    firstTrainTmp = (int)(rgHeader[0]);
    lastTrainTmp = (int)(rgHeader[1]);
    chunkFirstTmp = (int)(rgHeader[2]);
    chunkLastTmp = (int)(rgHeader[3]);
    numTrainTmp = lastTrainTmp-firstTrainTmp+1;
    chunkSizeTmp = chunkLastTmp-chunkFirstTmp+1;

    if(-1 == trainFirst){
      trainFirst = firstTrainTmp;
      trainLast = lastTrainTmp;
      numTrain = numTrainTmp;
      printf("allocating %ld bytes\n",numTrain*sizeof(bool));
      rgDataRead = new bool[numTrain];
      if(!rgDataRead){
	fprintf(stderr,"allocation error of rgDataRead numTrain=%ld\n",numTrain);
	exit(1);
      }
      for(long int i=0; i < numTrain;++i)
	rgDataRead[i] = false;
      printf("allocating %ld bytes\n",numTrain*numTrain*sizeof(double));
      rgMatrix = new double[numTrain*numTrain];
      if(!rgMatrix){
	fprintf(stderr,"allocation error.  couldn't allocate rgMatrix (%ld bytes)\n",
		numTrain*numTrain*(long int)sizeof(double));
	exit(1);
      }
      else{
	printf("successfully allocated trainMatrix %ldx%ld (%.2lf Gig)\n",
	       numTrain,numTrain,
	       sizeof(double)*(numTrain*numTrain)/1000000000.);
      }
      for(long int i=0; i < numTrain*numTrain; ++i)
	rgMatrix[i] = 999999.;
    }

    if((firstTrainTmp != trainFirst) || (lastTrainTmp != trainLast)){
      fprintf(stderr,"file '%s' has trainFirst=%d trainLast=%d, %d and %d were specified in first chunk file\n",argv[ii],firstTrainTmp,lastTrainTmp,
	      trainFirst,trainLast);
      fAllOk = false;
      fclose(fin);
      continue;
    }
    if(chunkFirstTmp > chunkLastTmp){
      fprintf(stderr,"chunkFirstTmp(%d) > chunkLastTmp(%d)\n",
	      chunkFirstTmp, chunkLastTmp);
      fAllOk = false;
      fclose(fin);
      continue;
    }
    if(chunkFirstTmp < trainFirst){
      fprintf(stderr,"chunkFirst=%d but trainFirst=%d -- file '%s'\n",
	      chunkFirstTmp, trainFirst, argv[ii]);
      fAllOk = false;
      fclose(fin);
      continue;
    }
    if(chunkLastTmp > trainLast){
      fprintf(stderr,"chunkLast=%d but trainLast=%d -- file '%s'\n",
	      chunkLastTmp, trainLast, argv[ii]);
      fAllOk = false;
      fclose(fin);
      continue;
    }
    bool fOverlap;
    fOverlap = false;
    for(long int i=(chunkFirstTmp-trainFirst); i<=(chunkLastTmp-trainFirst); ++i){
      if(rgDataRead[i]){
	fprintf(stderr,"overlap at %ld\n",i);
	fOverlap = true;
      }
    }
    if(fOverlap){
      fprintf(stderr,"chunk '%s' overlaps with a previously read chunk!\n",
	      argv[ii]);
      fAllOk = false;
      fclose(fin);
      continue;
    }

    if(chunkSizeTmp*numTrainTmp !=
       (int)fread(&(rgMatrix[(chunkFirstTmp-trainFirst)*numTrainTmp]),
		  sizeof(double),chunkSizeTmp*numTrainTmp,fin)){
      fprintf(stderr,"error trying to read chunk size of %ld rows from "
	      "chunk file '%s'\n",chunkSizeTmp, argv[ii]);
      fAllOk = false;
      fclose(fin);
      continue;
    }
    else{
      printf("successfully read chunk %d-%d (of %ld total) from '%s'\n",
	     chunkFirstTmp,chunkLastTmp,numTrain,argv[ii]);
      for(long int i=(chunkFirstTmp-trainFirst); i <= (chunkLastTmp-trainFirst); ++i)
	rgDataRead[i] = true;
    }
    fclose(fin);
  }
  
  // for(int i=0; i < numTrain; ++i){
  //   printf("rgDataRead[%d]=%d\n",i,(int)rgDataRead[i]);
  // }
  for(long int i=0; i < numTrain; ++i){
    if(!rgDataRead[i]){
      fAllOk = false;
      if((i==0)||((i>0)&&rgDataRead[i-1]))
	fprintf(stderr,"missing chunk %ld-",i);
      if((i==(numTrain-1))||((i<(numTrain-1))&&(rgDataRead[i+1])))
	fprintf(stderr,"-%ld\n",i);
    }
  }

  // printf("debug printing matrix:\n");
  // for(int r=0; r < numTrain; ++r){
  //   for(int c=0; c < numTrain; ++c){
  //     printf(" %10.2lf", rgMatrix[r*numTrain+c]);
  //   }
  //   printf("\n");
  // }
  // printf("\n");
  if(fAllOk){
    for(long int r=0; r < numTrain; ++r){
      for(long int c=0; c < r; ++c){
	rgMatrix[r*numTrain+c] = rgMatrix[c*numTrain+r];
      }
    }
  }


  // printf("debug printing matrix after symmetry copy:\n");
  // for(int r=0; r < numTrain; ++r){
  //   for(int c=0; c < numTrain; ++c){
  //     printf(" %10.2lf", rgMatrix[r*numTrain+c]);
  //   }
  //   printf("\n");
  // }
  // printf("\n");


  if(fAllOk){
    printf("saving matrix file '%s'...\n",stOutFile);
    FILE *fout;
    fout = fopen(stOutFile,"wb");
    if(!fout){
      fprintf(stderr,"couldn't open '%s' for output\n",stOutFile);
      exit(1);
    }
    rgHeader[0]=trainFirst;
    rgHeader[1]=trainLast;
    for(int i=2; i < 126; ++i)
      rgHeader[i] = rgHeader[i+2];
    rgHeader[126] = 0.;
    rgHeader[127] = 0.;
    stHeaderText = (char*)&(rgHeader[2]);
    //erase the 'chunkFirst' and 'chunkLast' from the string
    for(int i=0, len=(int)strlen(stHeaderText); i < len; ++i){
      if(0 == strncmp("chunkFirst=",&(stHeaderText[i]),11)){
	for(int j=i; j < len; ++j)
	  stHeaderText[j]='\0';
	break;
      }
    }
    if(128 != fwrite(rgHeader,sizeof(double),128,fout)){
      fprintf(stderr,"couldn't fwrite header to '%s'\n",stOutFile);
      exit(1);
    }
    if(numTrain*numTrain !=
       (long int)fwrite(rgMatrix, sizeof(double),numTrain*numTrain,fout)){
      fprintf(stderr,"couldn't fwrite the cost matrix block to '%s'\n",
	      stOutFile);
      exit(1);
    }
    fclose(fout);
    printf("successfully combined chunk files into '%s'\n",stOutFile);
  }
  else{
    fprintf(stderr,"ERROR! did not save output matrix file because of previous errors\n");
  }

  if(NULL != rgMatrix)
    delete [] rgMatrix;
  if(NULL != rgDataRead)
    delete [] rgDataRead;

  return 0;
}
