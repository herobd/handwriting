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
#include <vector>
#include <queue>
#include <stack>
#include <map>
#include <sys/stat.h>

#define USE_MEAN_VECTORS 1
#define DONT_USE_SMART_KEY_IMAGE_SELECTION 1
#define USE_MOST_FREQUENT_LABELS_AS_SAMPLES 0

#define USE_FAST_PASS_FIRST 1
#define USE_FAST_PASS_FOR_NxN_TRAINING 1
#define TRACK_THE_BEST_N 1
#define BEST_N_UNIQUE_LABELS 1

#define DP_ONLY 0


//#define D_NOTHREADS

#ifndef D_NOTHREADS
#include "dthreads.h"
#endif


/*comparison function for qsort()ing doubles in nondecreasing order */
int compareDoubles(const void *p1, const void *p2){
  if ( (*(double*)(p1)) < (*(double*)(p2)) )
    return -1;
  else if ( (*(double*)(p1)) > (*(double*)(p2)) )
    return 1;
  return 0;
}


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

//n is numTrain, maxCost is th value mapped to 255 (anything > is clipped)
void saveMatrixImageRGB(char *stFileName, double *rgMatrix, int n,
			double maxCost){
  FILE *fout;
  fout = fopen(stFileName,"wb");
  if(!fout){
    fprintf(stderr,"saveMatrixImage() couldn't open '%s' for writing\n",
	    stFileName);
    exit(1);
  }
  fprintf(fout,"P6\n%d %d\n255\n",n,n);
  double scale;
  scale = 255./maxCost;
  for(int r=0,idx=0; r < n; ++r){
    for(int c=0; c < n; ++c,++idx){
      int val;
      val = (int)(scale*rgMatrix[r*n+c]);
      if(val > 255)
	val = 255;
      if(val < -255.)
	val = -255.;
      if(val < 0){
	fprintf(fout,"%c%c%c",(char)val,0,0);
      }
      else{
	fprintf(fout,"%c%c%c",0,(char)val,0);
      }
    }
  }
  fclose(fout);
}



typedef struct{
  int numThreads;//how many threads are doing comparisons
  int threadNum;//which thread number this is (0..numThreads-1)
  int testWordIdx; //0-based (even if first test word >0) index into arrays
  DImage *pimgTest;//pointer shared by all threads
  DImage *rgTrainingImages;//pointer shared by all threads
  int numTrain;
  double *rgCostsMorph;//shared by all threads
  double weightMovement; // 0. to 1. (how much to weight the movement in cost)
  double lengthPenalty;
  int meshSpacingStatic;//-1 if auto-calculate like originally done
  int numRefinesStatic;//-1 if auto-calculate like originally done
  double meshDiv;//what to divide word image height by to get mesh size [4.0]
  int bandWidthDP;//sakeo-chiba bandwidth for DP
  int r,c;//row and column of the NxN training cost matrix
} WORDWARP_THREAD_PARMS;


typedef struct{
  int numThreads;
  int threadNum;
  DImage *rgTrainingImages;//pointer shared by all threads
  int numTrain;
  double *rgTrainCostMatrix;
  double weightMovement; // 0. to 1. (how much to weight the movement in cost)
  double lengthPenalty;
  int meshSpacingStatic;//-1 if auto-calculate like originally done
  int numRefinesStatic;//-1 if auto-calculate like originally done
  double meshDiv;//what to divide word image height by to get mesh size [4.0]
  int bandWidthDP;//sakeo-chiba bandwidth for DP
} TRAIN_NXN_THREAD_PARMS;


class HAC_TREE_NODE;//forward declaration so we can use pointer to own type
class HAC_TREE_NODE{
public:
  HAC_TREE_NODE();
  ~HAC_TREE_NODE();
  //tree information
  int clustID;
  HAC_TREE_NODE *pParent;
  HAC_TREE_NODE **rgpChildren;
  int numChildren;
  int numDescendantTreeNodes;//doesn't include this node itself

  //word information
  int numWordsInClustAtThisLevel;
  int *rgWordsInClustAtThisLevel;//indexes in rgTrain
  int numWordsIncludingDescendants;
  int *rgWordsIncludingDescendants;//indexes in rgTrain

  //comparison information
  int centerIdx;//index in rgTrain of word considered to be center of this node
  double maxDistFromCenter;//max dist of any (including children) to center
  double sumDistFromCenter;//used to compute avg.  includes descendants too

  DFeatureVector fvCenter;
  DFeatureVector fvMaxDistFromCenter;
  double maxDistFromCenterFv;

  DFeatureVector fvSumDescendants;//used to calculate fvMean using all descendts
  DFeatureVector fvMean;//mean of all descendant feature vectors (from key imgs)
  DFeatureVector fvMaxDistFromDescendantsMean;
  double maxDistFromDescendantsMeanFv;

  int *rgTopNTrainingMatches;// these are indexes (0-based) into training array

  double centerAspectRatio_w_div_h;
  static int numInstances;//this is not thread safe!! only build one at a time
};
int HAC_TREE_NODE::numInstances = 0;

HAC_TREE_NODE::HAC_TREE_NODE(){
  //  clustID=-1;
  pParent = NULL;
  rgpChildren = NULL;
  numChildren = 0;
  numDescendantTreeNodes = 0;
  numWordsInClustAtThisLevel = 0;
  rgWordsInClustAtThisLevel = NULL;
  numWordsIncludingDescendants = 0;
  rgWordsIncludingDescendants = NULL;
  centerIdx = -1;
  maxDistFromCenter = 0.;
  sumDistFromCenter = 0.;
  maxDistFromCenterFv = 0.;
  maxDistFromDescendantsMeanFv = 0.;
  rgTopNTrainingMatches = NULL;
  centerAspectRatio_w_div_h = 0.;
  clustID=HAC_TREE_NODE::numInstances;
  ++HAC_TREE_NODE::numInstances;
}
HAC_TREE_NODE::~HAC_TREE_NODE(){
  clustID =-1;
  pParent = NULL;
  if(NULL != rgpChildren){
    for(int i=0; i < numChildren; ++i){
      if(rgpChildren[i]!= NULL)
	delete rgpChildren[i];
    }
    delete [] rgpChildren;
    rgpChildren = NULL;
  }
  numChildren = 0;
  numWordsInClustAtThisLevel = 0;
  if(NULL != rgWordsInClustAtThisLevel)
    delete [] rgWordsInClustAtThisLevel;
  rgWordsInClustAtThisLevel = NULL;
  numWordsIncludingDescendants = 0;
  if(NULL != rgWordsIncludingDescendants)
    delete [] rgWordsIncludingDescendants;
  if(NULL != rgTopNTrainingMatches)
    delete [] rgTopNTrainingMatches;
  rgWordsIncludingDescendants = NULL;
  centerIdx = -1;
  maxDistFromCenter = 0.;
  sumDistFromCenter = 0.;
  --HAC_TREE_NODE::numInstances;
}


//node type for use by std::priority_queue
typedef struct{
  HAC_TREE_NODE *pHTnode;
  double priority;
  double morphCostFromTestToCenter;
  double distFromTestToCenterFV;
  double distFromTestToMeanFV;
} PQ_NODE_T;

typedef struct{
  int trIdx;//idx of training example
  double cost; //cost from current word to training word trIdx
} FASTPASS_SORT_NODE_T;

int compare_fastpass_sort(const void *pv1, const void *pv2){
  FASTPASS_SORT_NODE_T *p1, *p2;
  p1 = (FASTPASS_SORT_NODE_T*)pv1;
  p2 = (FASTPASS_SORT_NODE_T*)pv2;
  if((p1->cost) < (p2->cost))
    return -1;
  else if((p1->cost) > (p2->cost))
    return 1;
  return 0;
}

typedef struct{
  int trIdx;//idx of training example
  double cost; //cost from current word to training word trIdx
} TOPN_MATCHES_T;

int compare_TOPN_matches(const void *pv1, const void *pv2){
  TOPN_MATCHES_T *p1, *p2;
  p1 = (TOPN_MATCHES_T*)pv1;
  p2 = (TOPN_MATCHES_T*)pv2;
  if((p1->cost) < (p2->cost))
    return -1;
  else if((p1->cost) > (p2->cost))
    return 1;
  return 0;
}



//comparison function for use by std::priority_queue
bool compare_for_min_priority_queue(const PQ_NODE_T *a, const PQ_NODE_T *b){
  return ((a->priority) > (b->priority));
}
class Compare_for_min_pri_queue{
public:
 bool operator()(const PQ_NODE_T *p1, const PQ_NODE_T *p2){
   return ((p1->priority) > (p2->priority));
 }
  // bool comp(const PQ_NODE_T *p1, const PQ_NODE_T *p2){
  //   return ((p1->priority) > (p2->priority));
  // }
};
void calculateTreeKeyVectors(HAC_TREE_NODE *pTreeRoot,
			     int *rgKeyWordIdxs, int numKeyWords,
			     double *rgTrainCostMatrix, int numTrain,
			     HAC_TREE_NODE **rgOrigWordLeafNodes){
  std::stack<HAC_TREE_NODE*> searchStack;
  HAC_TREE_NODE *pCur;
  double *rgData;
  double *rgDataZero;
  rgData = new double[numKeyWords];
  D_CHECKPTR(rgData);
  rgDataZero = new double[numKeyWords];
  D_CHECKPTR(rgDataZero);
  for(int kk=0; kk < numKeyWords; ++kk)
    rgDataZero[kk] = 0.;
  searchStack.push(pTreeRoot);
  while(!searchStack.empty()){
    pCur = searchStack.top();
    if(0 == pCur->fvCenter.dimensions){//this is on the way down
      for(int kk=0; kk < numKeyWords; ++kk){
	rgData[kk] =
	  rgTrainCostMatrix[(pCur->centerIdx)*(long)numTrain+rgKeyWordIdxs[kk]];
      }
      pCur->fvCenter.setData_dbl(rgData, numKeyWords,1,true,true,true);
      pCur->maxDistFromCenterFv = 0.;
      pCur->fvSumDescendants.setData_dbl(rgDataZero,numKeyWords,1,
					 true,true,true);
      if(0 == pCur->numChildren){//this is a leaf node
	//set the sum vector to this word's vector
	pCur->fvSumDescendants.setData_dbl(rgData,numKeyWords,1,true,true,true);
      }
      for(int i=0; i < pCur->numChildren; ++i){
	searchStack.push(pCur->rgpChildren[i]);
      }
    }
    else{//this is on the way back up after children have been processed
      searchStack.pop();
      // pCur->fvSumDescendants.setData_dbl(rgDataZero,numKeyWords,1,
      // 					 true,true,true);
      for(int i=0; i < pCur->numChildren; ++i){
	pCur->fvSumDescendants.add(pCur->rgpChildren[i]->fvSumDescendants);
      }
      pCur->fvMean = pCur->fvSumDescendants;
      if(pCur->numChildren > 0)
	pCur->fvMean.divideByScalar((double)(pCur->numWordsIncludingDescendants));
      //now figure out maxDistFromCenterFv and maxDistFromDescendantsMeanFv
      pCur->maxDistFromCenterFv = 0.;
      pCur->maxDistFromDescendantsMeanFv = 0.;
      for(int i=0; i < pCur->numWordsIncludingDescendants; ++i){
	double distTmpSq;
	HAC_TREE_NODE *pDesc;//pointer to descendant word leaf node
	pDesc = rgOrigWordLeafNodes[pCur->rgWordsIncludingDescendants[i]];
	distTmpSq = 
	  DFeatureVector::getEuclideanDistSquared(pCur->fvCenter,
						  pDesc->fvCenter);
	if(distTmpSq > pCur->maxDistFromCenterFv)
	  pCur->maxDistFromCenterFv = distTmpSq;

	distTmpSq = 
	  DFeatureVector::getEuclideanDistSquared(pCur->fvMean,
						  pDesc->fvCenter);
	if(distTmpSq > pCur->maxDistFromDescendantsMeanFv)
	  pCur->maxDistFromDescendantsMeanFv = distTmpSq;
      }
      pCur->maxDistFromCenterFv = sqrt(pCur->maxDistFromCenterFv);
      pCur->maxDistFromDescendantsMeanFv =
	sqrt(pCur->maxDistFromDescendantsMeanFv);
    }
  }
  delete [] rgData;
  delete [] rgDataZero;
}

//assumes that calculateTreeKeyVectors() has already been called
void save3dWordPlot(char *stFileName, HAC_TREE_NODE *pTreeRoot,
		    std::string *rgLabelsTrain, int numTrain){
  FILE *fout;
  double *rgXs, *rgYs, *rgZs;
  int *rgIdxs;
  std::stack<HAC_TREE_NODE*> searchStack;
  HAC_TREE_NODE *pCur;

  fout=fopen(stFileName,"wb");
  if(!fout){
    fprintf(stderr,"couldn't open '%s' for output of matlab commands\n",
	    stFileName);
    exit(1);
  }
  rgXs = new double[numTrain];
  D_CHECKPTR(rgXs);
  rgYs = new double[numTrain];
  D_CHECKPTR(rgYs);
  rgZs = new double[numTrain];
  D_CHECKPTR(rgZs);
  rgIdxs = new int[numTrain];
  D_CHECKPTR(rgIdxs);

  int num;
  num = 0;
  searchStack.push(pTreeRoot);
  while(!searchStack.empty()){
    pCur = searchStack.top();
    searchStack.pop();
    if(0 == pCur->numChildren){//leaf node
      rgXs[num] = pCur->fvCenter.pDbl[0];
      rgYs[num] = pCur->fvCenter.pDbl[1];
      rgZs[num] = pCur->fvCenter.pDbl[2];
      rgIdxs[num] = pCur->centerIdx;
      ++num;
    }
    for(int i=0; i < pCur->numChildren; ++i){
      searchStack.push(pCur->rgpChildren[i]);
    }
  }//end while
  if(num != numTrain){
    fprintf(stderr,"WEIRD! num=%d but numTrain=%d\n",num,numTrain);
    exit(1);
  }
  fprintf(fout,"x=[");
  for(int i=0; i < numTrain; ++i){
    if(i>0)fprintf(fout,",");
    fprintf(fout,"%lf",rgXs[i]);
  }
  fprintf(fout,"];\n");

  fprintf(fout,"y=[");
  for(int i=0; i < numTrain; ++i){
    if(i>0)fprintf(fout,",");
    fprintf(fout,"%lf",rgYs[i]);
  }
  fprintf(fout,"];\n");

  fprintf(fout,"z=[");
  for(int i=0; i < numTrain; ++i){
    if(i>0)fprintf(fout,",");
    fprintf(fout,"%lf",rgZs[i]);
  }
  fprintf(fout,"];\n");

  fprintf(fout,"labels={");
  for(int i=0; i < numTrain; ++i){
    if(i>0)fprintf(fout,";");
    fprintf(fout,"'%s'",rgLabelsTrain[rgIdxs[i]].c_str());
  }
  fprintf(fout,"};\n");
  fprintf(fout,"\nfigure\nplot3(x,y,z,'r*')\ntext(x,y,z,labels);\n");
  fclose(fout);
}

void saveTreeForGraphviz(char *stFileName,
			 HAC_TREE_NODE *pTreeRoot, std::string *rgLabelsTrain,
			 double *rgTrainCostMatrix, int numTrain,
			 int *rgKeyWordIdxs, int numKeyWords){
  std::stack<HAC_TREE_NODE*> searchStack;
  HAC_TREE_NODE *root;
  FILE *fout;
  
  fout = fopen(stFileName,"wb");
  if(!fout){
    fprintf(stderr,"ERROR trying to open graphviz output file '%s'\n",
	    stFileName);
    return;
  }
  fprintf(fout,"digraph G {\n");
  //nodes
  searchStack.push(pTreeRoot);
  while(!searchStack.empty()){
    root = searchStack.top();
    char stFV[128],stFVtmp[128];//strings to print feature vector
    //sprintf(stFV,"");
    stFV[0]='\0';
    for(int kk=0; kk < numKeyWords; ++kk){
      sprintf(stFVtmp," %.2lf", root->fvMean.pDbl[kk]);
      strcat(stFV,stFVtmp);
    }
    if(0 == root->numChildren){
      fprintf(fout,"I%d [label=\"#%d(%s)%.2lf\",color=blue];\n",root->clustID,
	      root->centerIdx, rgLabelsTrain[root->centerIdx].c_str(),
	      root->maxDistFromCenter);
      // fprintf(fout,"I%d [label=\"#%d(%s)%.2lf\\n%s\",color=blue];\n",root->clustID,
      // 	      root->centerIdx, rgLabelsTrain[root->centerIdx].c_str(),
      // 	      root->maxDistFromDescendantsMeanFv/*maxDistFromCenterFv*/,stFV);
    }
    else{
      fprintf(fout,"I%d [label=\"#%d(%s)%.2lf\"];\n",root->clustID,
	      root->centerIdx, rgLabelsTrain[root->centerIdx].c_str(),
	      root->maxDistFromCenter);
    }
      // fprintf(fout,"I%d [label=\"#%d(%s)%.2lf\\n%s\"];\n",root->clustID,
      // 	      root->centerIdx, rgLabelsTrain[root->centerIdx].c_str(),
      // 	      root->maxDistFromDescendantsMeanFv/*maxDistFromCenterFv*/,stFV);
    searchStack.pop();
    for(int i=0; i < root->numChildren; ++i){
      searchStack.push(root->rgpChildren[i]);
    }
  }//end while
  //edges
  searchStack.push(pTreeRoot);
  while(!searchStack.empty()){
    root = searchStack.top();
    // fprintf(fout,"I%d [label=\"%d#%d(%s)\"];\n",root->clustID,root->clustID,
    // 	    root->centerIdx, rgLabelsTrain[root->centerIdx].c_str());
    searchStack.pop();
    for(int i=0; i < root->numChildren; ++i){
      //    for(int i=root->numChildren-1; i >= 0; --i){
      if(0!=i)
	fprintf(fout," I%d -> I%d [fontsize=10,label=\"%.2lf\",color=blue];\n",
		root->clustID,
		root->rgpChildren[i]->clustID,
		rgTrainCostMatrix[(root->centerIdx)*(long)numTrain+root->rgpChildren[i]->centerIdx]);
      else
	fprintf(fout," I%d -> I%d [fontsize=10,label=\"%.2lf\"];\n",
		root->clustID,
		root->rgpChildren[i]->clustID,
		rgTrainCostMatrix[(root->centerIdx)*(long)numTrain+root->rgpChildren[i]->centerIdx]);
      searchStack.push(root->rgpChildren[i]);
    }
  }//end while
  fprintf(fout,"}\n");
  fclose(fout);
}

//if numCompares is not NULL, the number of morphCompares will be put in it
HAC_TREE_NODE* findBestMatchNodeInTree(DImage &imgTest,
				       HAC_TREE_NODE *pTreeRoot,
				       DImage *rgTrainingImages,
				       int numTrain,
				       int bandWidth,
				       int meshSpacingStatic,
				       int numRefinesStatic,
				       double meshDiv,
				       double lengthPenalty,
				       double alpha,
				       int *numCompares,
				       double *minC,
				       std::string *rgLabelsTrain,
				       std::string strLabelTest,
				       int slowPassTopN,
				       FASTPASS_SORT_NODE_T *rgFastpassSorted,
				       int topNMatches
				       ){
  // now traverse the tree and find the best match using branch and bound
  double minCostSoFar; // this is m in branch and bound
  HAC_TREE_NODE *pMinNode;
  std::priority_queue<PQ_NODE_T*, std::vector<PQ_NODE_T*>,Compare_for_min_pri_queue> priQueue;
  DMorphInk mobj;
#if DP_ONLY
  mobj.fOnlyDoCoarseAlignment = true;
#endif

  // here we could prime the search by choosing a few frequent words or using
  // a dynamic cache of frequent words

  // initialize to search for this word's best training match
  int numMorphCompares;
  double *rgCostFromTestToTrain;
  rgCostFromTestToTrain = new double[numTrain];
  D_CHECKPTR(rgCostFromTestToTrain);

  for(int jjj=0; jjj < numTrain; ++jjj)
    rgCostFromTestToTrain[jjj] = 999999.;
  numMorphCompares = 0;

  int j;
  j = pTreeRoot->centerIdx;
  if(true/*999999. == rgCostFromTestToTrain[j]*/){
#if USE_FAST_PASS_FIRST
    rgCostFromTestToTrain[j] =
      mobj.getWordMorphCostFast(imgTest, rgTrainingImages[j],
				bandWidth,/*bandWidthDP*/
				0./*nonDiagonalCostDP*/,
				meshSpacingStatic,
				numRefinesStatic,
				meshDiv,
				lengthPenalty);
#else
    rgCostFromTestToTrain[j] =
      mobj.getWordMorphCost(imgTest, rgTrainingImages[j],
			    bandWidth,/*bandWidthDP*/
			    0./*nonDiagonalCostDP*/,
			    meshSpacingStatic,
			    numRefinesStatic,
			    meshDiv,
			    lengthPenalty);
#endif //USE_FAST_PASS_FIRST
    ++numMorphCompares;
  }
  pMinNode = pTreeRoot;
  minCostSoFar = rgCostFromTestToTrain[j];

  {//this block is experimental and should be removed------------------------
    //it finds the best match for each word (as an ORACLE) and starts with it
    //as the "best so far" to see how well the search prunes.
#if 0    
    std::stack<HAC_TREE_NODE*> searchStack;
    HAC_TREE_NODE *root;

    searchStack.push(pTreeRoot);
    while(!searchStack.empty()){
      root = searchStack.top();
      searchStack.pop();
      for(int ii=0; ii < root->numWordsInClustAtThisLevel; ++ii){
	if(0==strcmp(rgLabelsTrain[root->rgWordsInClustAtThisLevel[ii]].c_str(),
		     strLabelTest.c_str())){
	  double tmpMorphCost;
#if USE_FAST_PASS_FIRST
	  tmpMorphCost =
	    mobj.getWordMorphCostFast(imgTest,
				      rgTrainingImages[root->rgWordsInClustAtThisLevel[ii]],
				      bandWidth,/*bandWidthDP*/
				      0./*nonDiagonalCostDP*/,
				      meshSpacingStatic,
				      numRefinesStatic,
				      meshDiv,
				      lengthPenalty);
#else
	  tmpMorphCost =
	    mobj.getWordMorphCost(imgTest,
				  rgTrainingImages[root->rgWordsInClustAtThisLevel[ii]],
				  bandWidth,/*bandWidthDP*/
				  0./*nonDiagonalCostDP*/,
				  meshSpacingStatic,
				  numRefinesStatic,
				  meshDiv,
				  lengthPenalty);
#endif //USE_FAST_PASS_FIRST
	  if(tmpMorphCost < minCostSoFar){
	    minCostSoFar = tmpMorphCost;
	    pMinNode = root;
	  }
	}
      }
      for(int i=0; i < root->numChildren; ++i){
	searchStack.push(root->rgpChildren[i]);
      }
    }//end while
    printf("initializing search to cost %.2lf minNode=#%d(%s)\n",
	   minCostSoFar,pMinNode->centerIdx,rgLabelsTrain[pMinNode->centerIdx].c_str());
#endif
  }//end of experimental block-------------------------------------------------


  PQ_NODE_T *pQN;
  pQN = new PQ_NODE_T;
  pQN->pHTnode = pTreeRoot;
  pQN->priority = rgCostFromTestToTrain[j];
  pQN->morphCostFromTestToCenter = rgCostFromTestToTrain[j];
  priQueue.push(pQN);
  while(!priQueue.empty()){
    double mcostToCenter;
    HAC_TREE_NODE *pHTN;
    pQN = priQueue.top();
    priQueue.pop();
    mcostToCenter = pQN->morphCostFromTestToCenter;
    pHTN = pQN->pHTnode;
    delete pQN; pQN = NULL;
    // here we check the bounding function to see if we need to prune
    if( ((mcostToCenter) - (pHTN->maxDistFromCenter)) >
	(alpha * minCostSoFar)){//prune
      //printf("p(%d)",pHTN->numWordsIncludingDescendants);fflush(stdout);
    }
    else{// don't prune, expand the node by adding children to queue
      for(int nn=0; nn < (pHTN->numChildren); ++nn){
	int jj;
	HAC_TREE_NODE *pHTNnew;
	double mcostToCenterNew;
	
	pHTNnew = pHTN->rgpChildren[nn];
	jj=pHTNnew->centerIdx;
	if(999999. == rgCostFromTestToTrain[jj]){
	  //printf("*");fflush(stdout);
#if USE_FAST_PASS_FIRST
	  rgCostFromTestToTrain[jj] =
	    mobj.getWordMorphCostFast(imgTest,
				      rgTrainingImages[jj],
				      bandWidth,/*bandWidthDP*/
				      0./*nonDiagonalCostDP*/,
				      meshSpacingStatic,
				      numRefinesStatic,
				      meshDiv,
				      lengthPenalty);
#else
	  rgCostFromTestToTrain[jj] =
	    mobj.getWordMorphCost(imgTest,
				  rgTrainingImages[jj],
				  bandWidth,/*bandWidthDP*/
				  0./*nonDiagonalCostDP*/,
				  meshSpacingStatic,
				  numRefinesStatic,
				  meshDiv,
				  lengthPenalty);
#endif //USE_FAST_PASS_FIRST

	  ++numMorphCompares;
	}
	else{
	  //printf(".");fflush(stdout);
	}
	mcostToCenterNew = rgCostFromTestToTrain[jj];
	if( ((mcostToCenterNew) - (pHTNnew->maxDistFromCenter)) >
	    (alpha * minCostSoFar)){//prune (don't add to priorityQueue)
	  //printf("P(%d)",pHTNnew->numWordsIncludingDescendants);fflush(stdout);
	}
	else{
	  PQ_NODE_T *pQNnew;
	  pQNnew = new PQ_NODE_T;
	  pQNnew->pHTnode = pHTNnew;
	  pQNnew->priority = mcostToCenterNew;
	  pQNnew->morphCostFromTestToCenter = mcostToCenterNew;
	  priQueue.push(pQNnew);
	  if(mcostToCenterNew < minCostSoFar){
	    minCostSoFar = mcostToCenterNew;
	    pMinNode = pHTNnew;
	  }
	}
      }
    }
  }//end while(!priQueue.empty())
  
  // printf("bestMatch='%s'\n", rgLabelsTrain[pMinNode->centerIdx].c_str());
  if(NULL != numCompares){
    (*numCompares) = numMorphCompares;
  }
  if(NULL != minC){
    (*minC) = minCostSoFar;
  }
#if TRACK_THE_BEST_N
  if((slowPassTopN > 0)||(topNMatches > 0)){
    for(int i=0; i < numTrain; ++i){
      rgFastpassSorted[i].trIdx = i;
      rgFastpassSorted[i].cost = rgCostFromTestToTrain[i];
    }
    qsort((void*)rgFastpassSorted, numTrain, sizeof(FASTPASS_SORT_NODE_T),
	  compare_fastpass_sort);
    // for(int i=0; i < slowPassTopN; ++i){
    //   rgTopNFoundInTree[i] = rgFastPassSortNodes[i].trIdx;
    //   rgTopNFoundInTreeCosts[i] = rgFastPassSortNodes[i].cost;
    // }
  }
#endif
  delete [] rgCostFromTestToTrain;
  return pMinNode;
} 



void* word_morphing_NxN_training_thread_func(void *params){
  TRAIN_NXN_THREAD_PARMS *pparms;
  int numTrain;
  DMorphInk mobj;
#if DP_ONLY
  mobj.fOnlyDoCoarseAlignment = true;
#endif
  double *rgCosts;

  pparms = (TRAIN_NXN_THREAD_PARMS*)params;
  numTrain = pparms->numTrain;
  rgCosts = pparms->rgTrainCostMatrix;

  for(long int idx = (pparms->threadNum); idx < numTrain*(long)numTrain;
      idx+=(pparms->numThreads)){
    double morphCost;
    long int r, c;
    r = idx / numTrain;
    c = idx % numTrain;
    if(r<c){
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
      rgCosts[idx] = morphCost;
    }
    else if(r==c){
      rgCosts[idx] = 0.;
    }
    else if(r>c){
      rgCosts[idx] = 999999.;
    }
    if(0==c){
      if(0==(r%1))
	printf(" NxN %.2lf%% complete row%ld\n",(r*100/(double)numTrain),r);
    }
  }
  return NULL;
}


typedef struct{
  int numThreads;//how many threads are doing comparisons
  int threadNum;//which thread number this is (0..numThreads-1)
  DImage *rgTrainingImages;//pointer shared by all threads
  int numTrain;
  int testFirst;
  int testLast;
  double *rgCostsMorph;//shared by all threads
  HAC_TREE_NODE **rgMatchNodes;//shared by all threads
  int *rgNumMorphCompares;//shared by all threads
  double weightMovement; // 0. to 1. (how much to weight the movement in cost)
  double lengthPenalty;
  int meshSpacingStatic;//-1 if auto-calculate like originally done
  int numRefinesStatic;//-1 if auto-calculate like originally done
  double meshDiv;//what to divide word image height by to get mesh size [4.0]
  int bandWidthDP;//sakeo-chiba bandwidth for DP
  HAC_TREE_NODE *treeRoot;//shared by all threads
  double alpha;//constant multiplier for pruning the tree
  char *stPathIn;//shared by all threads
  std::string *rgLabelsTest;//shared by all threads
  std::string *rgLabelsTrain;//shared by all threads
  bool *rgfCorrect;//shared by all threads
  int *rgNumPrunedFV;//shared by all threads
  int *rgNumPrunedChildFV;//shared by all threads
  int *rgKeyWordIdxs;//shared by all threads
  int numKeyWords;
  int slowPassTopN;//top N to check after fast pass
  HAC_TREE_NODE **rgOrigWordLeafNodes;//shared by all threads
  int topN;//N for word N-grams (currently must be <= slowPassTopN)
  TOPN_MATCHES_T *rgTopNMatches;//shared by all threads
  FASTPASS_SORT_NODE_T *rgFastpassSorted;
} TREE_SEARCH_THREAD_PARMS;

void* tree_search_thread_func(void *params){
  TREE_SEARCH_THREAD_PARMS *pparms;
  int numTrain, numTest;
  char stTmp[2048];//for loading test images
  DMorphInk mobj;
  FASTPASS_SORT_NODE_T *rgFastpassSorted = NULL;

  pparms = (TREE_SEARCH_THREAD_PARMS*)params;
  numTrain = pparms->numTrain;
  numTest = (pparms->testLast) - (pparms->testFirst) + 1;

#if TRACK_THE_BEST_N
  if((pparms->slowPassTopN) > numTrain){
    fprintf(stderr,"slowPassTopN must be less than numTrain!!!\n");
    exit(1);
  }
  // int *rgTopNFoundInTree = NULL;
  // double *rgTopNFoundInTreeCosts = NULL;
  // if((pparms->slowPassTopN) > 0){
  //   rgTopNFoundInTree = new int[pparms->slowPassTopN];
  //   D_CHECKPTR(rgTopNFoundInTree);
  //   rgTopNFoundInTreeCosts = new double[pparms->slowPassTopN];
  //   D_CHECKPTR(rgTopNFoundInTreeCosts);
  // }
  rgFastpassSorted = new FASTPASS_SORT_NODE_T[numTrain];
  D_CHECKPTR(rgFastpassSorted);
#endif
  //now find best training match in tree for each test word
  for(int tt=(pparms->testFirst)+(pparms->threadNum),i=pparms->threadNum;
      tt<=(pparms->testLast); tt+=(pparms->numThreads),i+=(pparms->numThreads)){
    DImage imgTest;
    DTimer t2;
    
    //    printf("thread %d: tt=%d\n",pparms->threadNum,tt);

    t2.start();
    sprintf(stTmp,"%s/thresh_w_%08d.pgm",pparms->stPathIn,tt);
    if(!imgTest.load(stTmp)){
      fprintf(stderr,"couldn't load test image '%s'\n",stTmp);
      exit(1);
    }
    std::string strTmp;
    strTmp = imgTest.getPropertyVal(std::string("label"));
    if(strTmp.size()<1){//couldn't find label property, try comments
      if(3 != imgTest.getNumComments()){
	fprintf(stderr, "image '%s' needs 'label' property or else comments should be: #threshval\\n#label\\n#pageNum\\n",stTmp);
	exit(1);
      }
      pparms->rgLabelsTest[i] = imgTest.getCommentByIndex(1);
    }
    else
      pparms->rgLabelsTest[i] = strTmp;
    //    if(pparms->threadNum==(pparms->numThreads - 1))



    double morphCost;
    HAC_TREE_NODE *pMinNode;
    int numMorphCompares;
    int numPrunedFV, numPrunedChildFV;
    numPrunedFV = numPrunedChildFV = 0;
    pMinNode = 
      findBestMatchNodeInTree(imgTest,
      			      pparms->treeRoot, 
      			      pparms->rgTrainingImages,
      			      numTrain,
      			      pparms->bandWidthDP,
      			      pparms->meshSpacingStatic,
      			      pparms->numRefinesStatic, 
      			      pparms->meshDiv,
      			      pparms->lengthPenalty,
      			      pparms->alpha,
      			      &numMorphCompares, &morphCost,
      			      //next two parameters are for experimental
			      pparms->rgLabelsTrain,pparms->rgLabelsTest[i],
			      // ,rgTopNFoundInTree,rgTopNFoundInTreeCosts,
			      pparms->slowPassTopN,
			      rgFastpassSorted,
			      pparms->topN
			      );
#if USE_FAST_PASS_FIRST
    //use the proper leaf node instead of an internal tree node as min
    int leafIdx;
    leafIdx = pMinNode->centerIdx;
    pMinNode = pparms->rgOrigWordLeafNodes[pMinNode->centerIdx];
    if((pMinNode->centerIdx) != leafIdx){//debug can be removed
      fprintf(stderr,"Found another error!\n");exit(1);//debug
    }//debug
    // do full morph on the top N to see if any are better
    for(int topn=0; topn < (pparms->slowPassTopN); ++topn){
      double newMorphCost;
      int trIdx;
#if TRACK_THE_BEST_N
      //trIdx = rgTopNFoundInTree[topn];
      trIdx = rgFastpassSorted[topn].trIdx;
#else
      trIdx = pMinNode->rgTopNTrainingMatches[topn];
#endif
      if(trIdx >= 0){
	newMorphCost = 
	  mobj.getWordMorphCost(imgTest,
				pparms->rgTrainingImages[trIdx],
				pparms->bandWidthDP,/*bandWidthDP*/
				0./*nonDiagonalCostDP*/,
				pparms->meshSpacingStatic,
				pparms->numRefinesStatic,
				pparms->meshDiv,
				pparms->lengthPenalty);
	++numMorphCompares;
	if(newMorphCost < morphCost){
	  morphCost = newMorphCost;
	  pMinNode = pparms->rgOrigWordLeafNodes[trIdx];
	  if(pMinNode->centerIdx != trIdx){//this debug block can be removed
	    fprintf(stderr,"Found an error!\n");exit(1);//debug
	  }//debug
	  
	}//end if(newMorphCost < morphCost)
#if TRACK_THE_BEST_N
	if(newMorphCost < rgFastpassSorted[topn].cost){
	  rgFastpassSorted[topn].cost = newMorphCost;//so we can sort for topN
	}
#endif
      }//end if(trIdx >= 0)
    }//end for(int topn...
#endif //USE_FAST_PASS_FIRST

#if TRACK_THE_BEST_N
    if((pparms->topN) > 0){
      qsort((void*)rgFastpassSorted, numTrain, sizeof(FASTPASS_SORT_NODE_T),
	    compare_fastpass_sort);
    }

#if BEST_N_UNIQUE_LABELS
    //initialize in case not enough
    for(int topn=0; topn < (pparms->topN); ++topn){
      pparms->rgTopNMatches[(pparms->topN)*(long)i+topn].trIdx = -1;
      pparms->rgTopNMatches[(pparms->topN)*(long)i+topn].cost = 999999.;
    }
    {
      int numSoFar;
      int trIdxTmp;
      trIdxTmp = 0;
      numSoFar = 0;
      while((numSoFar<(pparms->topN)) && (trIdxTmp<numTrain)){
	bool fFound;
	fFound = false;
	for(int jj=0; jj < numSoFar; ++jj){
	   // fprintf(stderr,"jj=%d numSoFar=%d trIdxTmp=%d rgFastpassSorted[trIdxTmp].trIdx=%d\n",
	   // 	  jj, numSoFar, trIdxTmp,
	   // 	  rgFastpassSorted[trIdxTmp].trIdx);
	  if(0==strcmp((pparms->rgLabelsTrain[pparms->rgTopNMatches[(pparms->topN)*(long)i+jj].trIdx]).c_str(),
		       (pparms->rgLabelsTrain[rgFastpassSorted[trIdxTmp].trIdx]).c_str()))
	    fFound = true;
	}
	if(!fFound){
	  //	  fprintf(stderr,"  '%s' not found.. adding\n",pparms->rgLabelsTrain[rgFastpassSorted[trIdxTmp].trIdx].c_str());
	  pparms->rgTopNMatches[(pparms->topN)*(long)i+numSoFar].trIdx =
	    rgFastpassSorted[trIdxTmp].trIdx;
	  pparms->rgTopNMatches[(pparms->topN)*(long)i+numSoFar].cost =
	    rgFastpassSorted[trIdxTmp].cost;
	  ++numSoFar;
	}
	else{
	  //	  fprintf(stderr,"found '%s'.. skipping\n",pparms->rgLabelsTrain[rgFastpassSorted[trIdxTmp].trIdx].c_str());
	}
	++trIdxTmp;
      }
    }
#else
    for(int topn=0; topn < (pparms->topN); ++topn){
      pparms->rgTopNMatches[(pparms->topN)*(long)i+topn].trIdx =
	rgFastpassSorted[topn].trIdx;
      pparms->rgTopNMatches[(pparms->topN)*(long)i+topn].cost =
	rgFastpassSorted[topn].cost;
    }
#endif //BEST_N_UNIQUE_LABELS
#endif //TRACK_THE_BEST_N

    pparms->rgCostsMorph[i] = morphCost;
    pparms->rgMatchNodes[i] = pMinNode;
    pparms->rgNumMorphCompares[i] = numMorphCompares;
    if(0 == strcmp(pparms->rgLabelsTrain[pMinNode->centerIdx].c_str(),
		   pparms->rgLabelsTest[i].c_str())){
      pparms->rgfCorrect[i] = true;
    }
    else{
      pparms->rgfCorrect[i] = false;
    }
    pparms->rgNumPrunedFV[i] = numPrunedFV;
    pparms->rgNumPrunedChildFV[i] = numPrunedChildFV;

    printf("thread%2d: tt=%d image %d(%s) of %d (%d%%) time word=%.2lf sec  cost%lf numCompares=%d  matched#%d(%s) %s\n",
	   pparms->threadNum,tt,i,pparms->rgLabelsTest[i].c_str(),
	   numTest, 100*(i+1)/numTest, t2.getAccumulated(), morphCost,
	   numMorphCompares,
	   pMinNode->centerIdx,
	   pparms->rgLabelsTrain[pMinNode->centerIdx].c_str(),
	   (pparms->rgfCorrect[i]) ? "correct" : "WRONG!");
    fflush(stdout);

  }
#if TRACK_THE_BEST_N
  // delete [] rgTopNFoundInTree;
  // delete [] rgTopNFoundInTreeCosts;
  delete [] rgFastpassSorted;
#endif
  return NULL;
}



// this version computes distance only from the "center" word of each node
inline double getClusterNodeDist(HAC_TREE_NODE *node1,
			  HAC_TREE_NODE *node2,
			  double *rgCostMatrix, int N){
  return rgCostMatrix[(node1->centerIdx) * (long)N + (node2->centerIdx)];
}



// this will create a node node in the tree that combines the nodes at
// indices iMinDist and jMinDist and collapses the array
//
// this version of the function just keeps the center of i as the center of the
// combined node
void mergeClusterNodes(HAC_TREE_NODE **rgHACNodes, int *numHACNodes,
		       int iMinDist, int jMinDist, double minDist,
		       std::string *rgLabels,
		       double *rgTrainCostMatrix, int numTrain, bool fPrint){
  HAC_TREE_NODE *pnodei;
  HAC_TREE_NODE *pnodej;
  HAC_TREE_NODE *proot;

  pnodei = rgHACNodes[iMinDist];
  pnodej = rgHACNodes[jMinDist];
  proot = new HAC_TREE_NODE;
  D_CHECKPTR(proot);
  pnodei->pParent = proot;
  pnodej->pParent = proot;
  proot->rgpChildren = new HAC_TREE_NODE*[2];
  D_CHECKPTR(proot->rgpChildren);
  proot->rgpChildren[0] = pnodei;
  proot->rgpChildren[1] = pnodej;
  proot->numChildren = 2;

  proot->numWordsInClustAtThisLevel = 0;
  proot->numWordsIncludingDescendants =
    pnodei->numWordsIncludingDescendants + pnodej->numWordsIncludingDescendants;
  proot->centerIdx = pnodei->centerIdx;
  proot->centerAspectRatio_w_div_h = pnodei->centerAspectRatio_w_div_h;
  proot->rgWordsIncludingDescendants = new int[proot->numWordsIncludingDescendants];
  D_CHECKPTR(proot->rgWordsIncludingDescendants);
  if(pnodei->numWordsIncludingDescendants > 0){
    memcpy(proot->rgWordsIncludingDescendants,
	   pnodei->rgWordsIncludingDescendants,
	   sizeof(int)*pnodei->numWordsIncludingDescendants);
  }
  if(pnodej->numWordsIncludingDescendants > 0){
    memcpy(&(proot->rgWordsIncludingDescendants[pnodei->numWordsIncludingDescendants]),
	   pnodej->rgWordsIncludingDescendants,
	   sizeof(int)*pnodej->numWordsIncludingDescendants);
    //we can delete the child's copy now if we want since we wont need it again
  }
  // TODO: update maxDistFromCenter and sumDistFromCenter
  proot->maxDistFromCenter = pnodei->maxDistFromCenter;
  proot->sumDistFromCenter = pnodei->sumDistFromCenter;
  for(int jj=0; jj < pnodej->numWordsIncludingDescendants; ++jj){
    double costToCenter;
    costToCenter =
      rgTrainCostMatrix[pnodej->rgWordsIncludingDescendants[jj]*(long)numTrain+
			proot->centerIdx];
    if(costToCenter > proot->maxDistFromCenter)
      proot->maxDistFromCenter = costToCenter;
    proot->sumDistFromCenter += costToCenter;
  }

  //we can delete the children copies now if we want since we wont need again

  //debug print
  if(fPrint)
    printf("numHACNodes=%d merging i=%d j=%d clusters with centers #%d and #%d (%s & %s) max=%lf avg=%lf\n",
	   (*numHACNodes), iMinDist, jMinDist,
	   pnodei->centerIdx,pnodej->centerIdx,
	   rgLabels[pnodei->centerIdx].c_str(),
	   rgLabels[pnodej->centerIdx].c_str(),
	   proot->maxDistFromCenter,
	   proot->sumDistFromCenter/(proot->numWordsIncludingDescendants));
  rgHACNodes[iMinDist] = proot;
  for(int j=jMinDist; j < ((*numHACNodes)-1); ++j)
    rgHACNodes[j] = rgHACNodes[j+1];
  --(*numHACNodes);
}


void printHACTree(HAC_TREE_NODE *node, int depth, std::string *rgLabels,
		    int num){
  // print this node
  for(int i=0; i < depth; ++i)
    printf(" ");
  printf("clust#%d centerIdx=%d(%s) words=%d descendWords=%d\n",
	 node->clustID,node->centerIdx,
	 rgLabels[node->centerIdx].c_str(), node->numWordsInClustAtThisLevel,
	 node->numWordsIncludingDescendants);
  // print words in this node
  for(int i=0; i <node->numWordsInClustAtThisLevel; ++i){
    for(int j=0; j < depth; ++j)
      printf(" ");
    printf(">%d(%s)\n",node->rgWordsInClustAtThisLevel[i],
	   rgLabels[node->rgWordsInClustAtThisLevel[i]].c_str());
  }
  // traverse all children nodes
  for(int i=0; i < node->numChildren; ++i){
    printHACTree(node->rgpChildren[i],depth+1,rgLabels,num);
  }
}

void getLeafNodePtrForEachWord(HAC_TREE_NODE *node, HAC_TREE_NODE **rgPtrs){
  rgPtrs[node->centerIdx] = node;
  for(int i=0; i < node->numChildren; ++i)
    getLeafNodePtrForEachWord(node->rgpChildren[i], rgPtrs);
  printf(".");
}

int getPathToRoot(HAC_TREE_NODE *node, HAC_TREE_NODE **rgPtrs, int idx){
  rgPtrs[idx] = node;
  if(NULL == node)
    return idx;
  return getPathToRoot(node->pParent, rgPtrs, idx+1);
}

bool isNode1ChildOfNode2(HAC_TREE_NODE *node1, HAC_TREE_NODE *node2){
  if(node1 == node2)
    return true;
  if(node1 == NULL)
    return false;
  return isNode1ChildOfNode2(node1->pParent, node2);
}

void printPath(HAC_TREE_NODE **rgPathNodes, int pathLen, bool fReverse,
	       std::string *rgLabelsTrain){
  if(fReverse){
    for(int i=pathLen-1; i >= 0; --i)
      printf("  %p:%d(%s)",rgPathNodes[i],rgPathNodes[i]->centerIdx,
	     rgLabelsTrain[rgPathNodes[i]->centerIdx].c_str());
  }
  else{
    for(int i=0; i < pathLen; ++i)
      printf("  %p:%d(%s)",rgPathNodes[i],rgPathNodes[i]->centerIdx,
	     rgLabelsTrain[rgPathNodes[i]->centerIdx].c_str());
  }

}

// returns the training index of the leaf node found by traversing tree
int traverseTreeForTrainingWord(int trainIdxToFind, HAC_TREE_NODE *node,
				double *rgTrainCosts, int numTrain,
				std::string *rgLabelsTrain,
				HAC_TREE_NODE **rgPath, int *pathLen,
				int nearestIdx, double nearestDist,
				double *cost){
  rgPath[*pathLen] = node;
  ++(*pathLen);
  int minDistChildIdx = -1;
  double minDistChild = 999999.;
  if(0 == (node->numChildren)){
    (*cost) = nearestDist;
    rgPath[*pathLen] = NULL; //pathLen was already incremented above
    return nearestIdx;
  }

  for(int i=0; i < node->numChildren; ++i){
    double morphCost;
    if(node->rgpChildren[i]->centerIdx == nearestIdx){
      // use the same distance as passed in
      morphCost = nearestDist;
    }
    else{
      // compute it new (actually pull from training cost matrix)
      morphCost = rgTrainCosts[(node->rgpChildren[i]->centerIdx)*(long)numTrain+
			       trainIdxToFind];
    }
    if((0==i) || (morphCost < minDistChild)){
      minDistChild = morphCost;
      minDistChildIdx = i;
    }
  }
  return traverseTreeForTrainingWord(trainIdxToFind,
				     node->rgpChildren[minDistChildIdx],
				     rgTrainCosts, numTrain, rgLabelsTrain,
				     rgPath, pathLen, 
				     node->rgpChildren[minDistChildIdx]->centerIdx,
				     minDistChild, cost);
}

int main(int argc, char **argv);

int main(int argc, char **argv){
  char stPathIn[1025];
  char stOutfile[1025];
  char stTrainCostMatrix[1025];//file to load or write
  char stTmp[1025];
  int trainFirst, trainLast;
  int testFirst, testLast;
  int numTrain, numTest;
  DImage *rgTrainingImages;
  DImage testImage;
  double *rgTrainCostMatrix;
  double *rgCostsMorph;
  DTimer t1;
  std::string *rgLabelsTrain;
  std::string *rgLabelsTest;
  WORDWARP_THREAD_PARMS *rgThreadParms;
  TRAIN_NXN_THREAD_PARMS *rgTrainThreadParms;
  TREE_SEARCH_THREAD_PARMS *rgTreeThreadParms;
  int numThreads = 1;
  double weightMovement = 0.;
  double lengthPenalty = 0.;
  int meshSpacingStatic;
  int numRefinesStatic;
  double meshDiv;
  DMorphInk mobj;
#if DP_ONLY
  mobj.fOnlyDoCoarseAlignment = true;
#endif
  int bandWidth = 15;
  double alpha = 1.0;
  char stHACMergeFile[2048];//list of nodes to merge (to avoid recalculating)
  FILE *fmerge = NULL;
  int slowPassTopN = 10;
  int topNMatches = 10;// top N matches to save for N-gram post-processing
  TOPN_MATCHES_T *rgTopNMatches;//top N matches for all test words
#ifndef D_NOTHREADS
  pthread_t *rgThreadID;
#else
  numThreads = 1;
#endif


  if((argc < 18)||(argc > 19)){
    fprintf(stderr, "usage: %s <dataset_path> <trainCostMatrix> <first_training_num> <last_training_num> <first_test_num> <last_test_num> <weightMovement=0.> <lengthPenalty=0.> <numThreads=-1> <meshSpacingStatic=-1> <numRefinesStatic=-1> <meshDiv=4> <SakoeChibaBandwidth=15> <alpha=1.0> <slowPassTopN=10> <topNmatches=10> <output_file> [mergeFileHAC]\n",
	    argv[0]);
    exit(1);
  }

  sprintf(stPathIn, "%s", argv[1]);
  sprintf(stTrainCostMatrix, "%s", argv[2]);
  trainFirst = atoi(argv[3]);
  trainLast = atoi(argv[4]);
  testFirst = atoi(argv[5]);
  testLast = atoi(argv[6]);
  weightMovement = atof(argv[7]);
  lengthPenalty =  atof(argv[8]);
  numThreads = atoi(argv[9]);
  meshSpacingStatic = atoi(argv[10]);
  numRefinesStatic = atoi(argv[11]);
  meshDiv = atof(argv[12]);
  bandWidth = atoi(argv[13]);
  alpha = atof(argv[14]);
  slowPassTopN = atoi(argv[15]);
  topNMatches = atoi(argv[16]);
  sprintf(stOutfile, "%s", argv[17]);
  stHACMergeFile[0] = '\0';
  if(argc > 18)
    sprintf(stHACMergeFile, "%s", argv[18]);

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

  if((alpha < 0.1) || (alpha > 5.0)){
    fprintf(stderr,"alpha should be 0.1-5.0. was %lf\n",alpha);
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

  rgThreadParms = (WORDWARP_THREAD_PARMS *)malloc(sizeof(WORDWARP_THREAD_PARMS)*
						  numThreads);
  D_CHECKPTR(rgThreadParms);
  rgTrainThreadParms =
    (TRAIN_NXN_THREAD_PARMS*)malloc(sizeof(TRAIN_NXN_THREAD_PARMS)*numThreads);
  D_CHECKPTR(rgTrainThreadParms);
  rgTreeThreadParms =
    (TREE_SEARCH_THREAD_PARMS*)malloc(sizeof(TREE_SEARCH_THREAD_PARMS)*
				      numThreads);
  D_CHECKPTR(rgTreeThreadParms);
  
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

  rgTrainingImages = new DImage[numTrain];
  D_CHECKPTR(rgTrainingImages);
  rgLabelsTrain = new std::string[numTrain];
  D_CHECKPTR(rgLabelsTrain);
  rgLabelsTest = new std::string[numTest];
  D_CHECKPTR(rgLabelsTest);
  rgCostsMorph = (double*)malloc(sizeof(double)*numTest);
  D_CHECKPTR(rgCostsMorph);


  rgTopNMatches = new TOPN_MATCHES_T[(long)numTest*topNMatches];
  D_CHECKPTR(rgTopNMatches);
  
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

  printf("allocating %ld bytes for rgTrainCostMatrix\n",
	 sizeof(double)*numTrain*(long)numTrain);
  rgTrainCostMatrix = new double[numTrain*(long)numTrain];
  D_CHECKPTR(rgTrainCostMatrix);

  // load the training NxN cost matrix
  FILE *ftrin;
  ftrin = fopen(stTrainCostMatrix,"rb");
  if(ftrin){
    char *stHeaderText;
    double rgHeader[128];//firstTrain,lastTrain, then whatever text
    int numTrainTmp;
    int firstTrainTmp, lastTrainTmp;

    printf("loading '%s' into rgTrainCostMatrix\n",stTrainCostMatrix);
    fflush(stdout);
    if(128!=fread(rgHeader,sizeof(double),128,ftrin)){
      fprintf(stderr,"fread did not return 128 header doubles\n");
      exit(1);
    }
    t1.start();
    firstTrainTmp = (int)(rgHeader[0]);
    lastTrainTmp = (int)(rgHeader[1]);
    numTrainTmp = lastTrainTmp-firstTrainTmp+1;
    stHeaderText = (char*)(&rgHeader[2]);
    stHeaderText[128*sizeof(double)-1]='\0';
    if((trainFirst < firstTrainTmp)||(trainLast > lastTrainTmp)){
      fprintf(stderr,"requested #%d-#%d as training, but matrix file '%s' only has costs for #%d-%d\n",trainFirst, trainLast, stTrainCostMatrix,
	      firstTrainTmp, lastTrainTmp);
      exit(1);
    }

    //skip forward in file if necessary
    long int skipForward = 0;
    skipForward = trainFirst - firstTrainTmp;
    if(skipForward > 0){
      printf("skipForward=%ld skipping %ld doubles\n",skipForward,
	     skipForward*(long)numTrainTmp);
      if(-1==fseek(ftrin,sizeof(double)*(skipForward*(long)numTrainTmp+
					 skipForward),
		   SEEK_CUR)){
	fprintf(stderr,"error fseeking in matrix file\n");
	exit(1);
      }
    }
    if(numTrainTmp == numTrain){//row sizes are the same, load entire block
      if(numTrain*(long)numTrain != (long)fread(rgTrainCostMatrix,sizeof(double),
					       numTrain*(long)numTrain,ftrin)){
	fprintf(stderr,"error reading the matrix as a block\n");
	exit(1);
      }
    }
    else{//row sizes don't match, read a row at a time
      for(int tr=0; tr < numTrain; ++tr){
	if((long)numTrain!=(long)fread(&(rgTrainCostMatrix[tr*(long)numTrain]),
				sizeof(double), numTrain,ftrin)){
	  fprintf(stderr,"error reading a line from the matrix\n");
	  exit(1);
	}
	if(-1==fseek(ftrin,sizeof(double)*(long)(numTrainTmp-numTrain),
		     SEEK_CUR)){
	  fprintf(stderr,"error fseeking to next line in matrix file\n");
	  exit(1);
	}
      }
    }
    t1.stop();
    printf("successfully loaded the matrix (took %.2f seconds)\n",
	   t1.getAccumulated());
    
  }
  else{
    printf("couldn't find matrix file '%s'\n",stTrainCostMatrix);

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
    // now mirror the upper half of the matrix
    printf("mirroring the upper half of the matrix...");fflush(stdout);
    for(int r=0;r < numTrain; ++r){
      for(int c=0; c < r; ++c){
    	rgTrainCostMatrix[r*(long)numTrain+c] =
	  rgTrainCostMatrix[c*(long)numTrain+r];
      }
    }

  //   for(int r=0,idx=0;r < numTrain; ++r){
  //     for(int c=0; c < numTrain; ++c,++idx){
  // 	printf(" %9.2lf",rgTrainCostMatrix[idx]);
  //     }
  //     printf("\n");
  //   }
  // }




    double elapsedTime;
    t1.stop();
    elapsedTime = t1.getAccumulated();
    printf("NxN comparison of training data took %.02f seconds\n", t1.getAccumulated());
    long int numCompares = (numTrain*(long)numTrain-numTrain)/2;
    printf("(n*n-n)/2=%ld comparisons = %.6lf sec per compare\n", numCompares,
	   t1.getAccumulated()/numCompares);

    // now save the cost matrix
    printf("SAVING matrix file '%s'\n",stTrainCostMatrix);
    t1.start();
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
    stHeaderText = (char*)&(rgHeader[2]);
    sprintf(stText,"\n%ld compares in %.2lf seconds\n"
	    "trainFirst=%d trainLast=%d\n"
	    "weightMovement=%.2lf lengthPenalty=%.2lf meshSpacingStatic=%d numRefinesStatic=%d meshDiv=%.2lf bandWidth=%d\ndir=%s\n",
	    numCompares, elapsedTime, trainFirst, trainLast, weightMovement,
	    lengthPenalty, meshSpacingStatic, numRefinesStatic,
	    meshDiv, bandWidth, stPathIn);
    for(int i=2; i < 128; ++i)
      rgHeader[i] = 0.;
    strncpy(stHeaderText,stText,126*sizeof(double));
    stHeaderText[126*sizeof(double)-1]='\0';
    if(128 != fwrite(rgHeader,sizeof(double),128,ftrin)){
      fprintf(stderr,"couldn't fwrite header to '%s'\n",stTrainCostMatrix);
      exit(1);
    }
    if(numTrain*(long)numTrain !=
       (int)fwrite(rgTrainCostMatrix, sizeof(double),numTrain*(long)numTrain,
		   ftrin)){
      fprintf(stderr,"couldn't fwrite the cost matrix block to '%s'\n",
	      stTrainCostMatrix);
      exit(1);
    }
    fclose(ftrin);
    t1.stop();
    printf("finished writing the cost matrix (took %.2f seconds)\n",
	   t1.getAccumulated());
	   
  }
  // printf("TRAINING WORDS:");
  // for(int i=0; i < numTrain; ++i){
  //   if(0==(i%10))
  //     printf("\n");
  //   printf(" %d(%s)",i,rgLabelsTrain[i].c_str());
  // }
  // printf("\n");
  // //debug: print the cost matrix
#if 0
  for(long int r=0,idx=0;r < numTrain; ++r){
    for(int c=0; c < numTrain; ++c,++idx){
      printf(" %9.2lf",rgTrainCostMatrix[idx]);
    }
    printf("\n");
  }
#endif  

  




  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //     do hierarchical clustering
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  t1.start();

  // int *rgClustAssignments;
  // rgClustAssignments = new int[numTrain];
  // D_CHECKPTR(rgClustAssignments);


  // int rgLevelSizes[3] = {50,1000,10000};
  // int numLevels=3;
  // int rgLevelSizes[3] = {200,10000,100000};
  // int numLevels=3;
  int rgLevelSizes[1] = {16};//use the words at level 1 as key images
  int numLevels=1;
  while((rgLevelSizes[numLevels-1] >= numTrain)&&(numLevels>1))
    --numLevels;

  printf(">> numLevels=%d\n", numLevels);

  HAC_TREE_NODE ***rgLevelHACNodes;
  rgLevelHACNodes = new HAC_TREE_NODE**[numLevels];
  D_CHECKPTR(rgLevelHACNodes);
  for(int lvl=0; lvl < numLevels; ++lvl){
    rgLevelHACNodes[lvl] = new HAC_TREE_NODE*[rgLevelSizes[lvl]];
    D_CHECKPTR(rgLevelHACNodes[lvl]);
  }


  //start with each word as its own node
  HAC_TREE_NODE **rgHACNodes;
  rgHACNodes = new HAC_TREE_NODE*[numTrain];
  int numHACNodes;
  numHACNodes = numTrain;
  for(int i=0; i < numHACNodes; ++i){
    rgHACNodes[i] = new HAC_TREE_NODE;
  }
  for(int i=0; i < numTrain; ++i){
    //    rgHACNodes[i]->clustID = i;
    rgHACNodes[i]->numWordsInClustAtThisLevel = 1;
    rgHACNodes[i]->numWordsIncludingDescendants = 1;
    rgHACNodes[i]->rgWordsInClustAtThisLevel = new int[1];
    rgHACNodes[i]->rgWordsInClustAtThisLevel[0] = i;
    rgHACNodes[i]->rgWordsIncludingDescendants = new int[1];
    rgHACNodes[i]->rgWordsIncludingDescendants[0] = i;

    rgHACNodes[i]->centerIdx = i;
    if(rgTrainingImages[i].height() != 0.)
      rgHACNodes[i]->centerAspectRatio_w_div_h = 
	rgTrainingImages[i].width() / rgTrainingImages[i].height();
    else
      rgHACNodes[i]->centerAspectRatio_w_div_h = 0.;
  }


  HAC_TREE_NODE **rgOrigWordLeafNodes;
  rgOrigWordLeafNodes = new HAC_TREE_NODE*[numTrain];
  D_CHECKPTR(rgOrigWordLeafNodes);
  for(int tr=0; tr < numTrain; ++tr){
    rgOrigWordLeafNodes[tr] = rgHACNodes[tr];
  }

  int curLevel;
  curLevel = numLevels-1;

#if USE_FAST_PASS_FIRST
#if !TRACK_THE_BEST_N
  FASTPASS_SORT_NODE_T *rgFastPassSortNodes;
  rgFastPassSortNodes = new FASTPASS_SORT_NODE_T[numTrain];
  D_CHECKPTR(rgFastPassSortNodes);
  printf("setting up top-N matches for every leaf in tree...(not threaded yet)\n");fflush(stdout);
  DTimer t8;
  t8.start();
  for(int tr=0; tr < numTrain; ++tr){
    int numSoFar;
    numSoFar = 0;
    rgHACNodes[tr]->rgTopNTrainingMatches = new int[slowPassTopN];
    D_CHECKPTR(rgHACNodes[tr]->rgTopNTrainingMatches);
    for(int jj=0; jj < slowPassTopN; ++jj)
      rgHACNodes[tr]->rgTopNTrainingMatches[jj] = -1;
    for(int trIdx=0; trIdx < numTrain; ++trIdx){
      rgFastPassSortNodes[trIdx].trIdx = trIdx;
      rgFastPassSortNodes[trIdx].cost =
	rgTrainCostMatrix[tr*(long)numTrain+trIdx];
    }
    qsort((void*)rgFastPassSortNodes, numTrain, sizeof(FASTPASS_SORT_NODE_T),
	  compare_fastpass_sort);
    if(0 != rgFastPassSortNodes[0].cost){
      fprintf(stderr,"oops!\n");
      exit(1);
    }
    int trIdxTmp;
    trIdxTmp = 0;
    // if(tr < 20)
    //   printf("tr%d:'%s' topN:\n",tr,rgLabelsTrain[tr].c_str());
    while((numSoFar<slowPassTopN) && (trIdxTmp<numTrain)){
      bool fFound;
      fFound = false;
      for(int jj=0; jj < numSoFar; ++jj){
	if(0==strcmp(rgLabelsTrain[rgHACNodes[tr]->rgTopNTrainingMatches[jj]].c_str(),
		     rgLabelsTrain[rgFastPassSortNodes[trIdxTmp].trIdx].c_str()))
	  fFound = true;
      }
#if BEST_N_UNIQUE_LABELS
#else
      fFound=false;
#endif
      if(!fFound){
	// if(tr < 20)
	//   printf("  [%d] %d %d:'%s'\n",numSoFar,trIdxTmp,rgFastPassSortNodes[trIdxTmp].trIdx,
	// 	 rgLabelsTrain[rgFastPassSortNodes[trIdxTmp].trIdx].c_str());
	rgHACNodes[tr]->rgTopNTrainingMatches[numSoFar] =
	  rgFastPassSortNodes[trIdxTmp].trIdx;
	++numSoFar;
      }
      else{
	// if(tr < 20)
	// printf("      skip %d %d:'%s'\n",trIdxTmp,
	//        rgFastPassSortNodes[trIdxTmp].trIdx,
	//        rgLabelsTrain[rgFastPassSortNodes[trIdxTmp].trIdx].c_str());
      }
      ++trIdxTmp;
    }
  }
  t8.stop();
  printf("time to set up top-N matches for tree: %.2lf seconds\n",
	 t8.getAccumulated());
  delete [] rgFastPassSortNodes;
#endif
#endif

  //lineNo0 iMerge jMerge *\n
  //lineNo1 iMerge jMerge *\n
  // ...

  //see if the merge file exists. if so, read all the merges 
  //file format: each line is exactly 32 bytes long and has two numbers
  //(15 digits each) separated by a space and a newline as the 32nd char
  //the 2 numbers on the first line are numMergesInFile and numTrain
  //numMergesInFile gets updated after each additional merge line is added
  //any partial line after the 
  std::vector<long int> vectMergeIJs;
  long int numHACMergeIdxsInFile = 0;
  if(0 != strlen(stHACMergeFile)){
    char stTmpString[1024];
    long int numTrainInFile = 0;
    long int lineNo = 0;
    fmerge = fopen(stHACMergeFile, "rb");
    if(fmerge){//opened it, so read the merges
      printf("reading merge file '%s'...\n",stHACMergeFile);
      if(2 != fscanf(fmerge,"%s%ld",stTmpString, &numTrainInFile)){
	fprintf(stderr,"ERROR! expect 'numTrain %%ld\\n' in merge file '%s'\n",
		stHACMergeFile);
	exit(1);
      }
      if(0 != strcmp(stTmpString,"numTrain")){
	fprintf(stderr,"ERROR!!! expect 'numTrain %%ld\\n' in file '%s'\n",
		stHACMergeFile);
	exit(1);
      }
      if(numTrainInFile != (long)numTrain){
	fprintf(stderr,"ERROR!!! numTrain=%d but it is %ld in merge file '%s'\n",numTrain, numTrainInFile, stHACMergeFile);
	exit(1);
      }
      while(!feof(fmerge)){
	long int mergeI, mergeJ;
	int numRet;
	numRet=fscanf(fmerge,"%ld%ld%ld%s",&lineNo,&mergeI,&mergeJ,stTmpString);
	if(4!=numRet){
	  if(feof(fmerge))//end of file
	    break;
	  fprintf(stderr,
		  "ERROR!! on line %ld of mergefile '%s'\n",
		  numHACMergeIdxsInFile+2, stHACMergeFile);
	  exit(1);
	}
	if(lineNo != numHACMergeIdxsInFile){
	  fprintf(stderr,"expect lineNo=%ld but was %ld on line %ld of '%s'\n",
		  numHACMergeIdxsInFile, lineNo, numHACMergeIdxsInFile+2,
		  stHACMergeFile);
	  exit(1);
	}
	if(0 != strcmp(stTmpString,"*")){
	  fprintf(stderr,"expect '*' at end of line %ld in merge file '%s'\n",
		  numHACMergeIdxsInFile+2, stHACMergeFile);
	  exit(1);
	}
	vectMergeIJs.push_back(mergeI);
	vectMergeIJs.push_back(mergeJ);
	++numHACMergeIdxsInFile;
      }
      fclose(fmerge);
      fmerge = NULL;
      // now merge them, then if not completely merged, open file for append
      for(long int ll=0; ll < numHACMergeIdxsInFile; ++ll){
	mergeClusterNodes(rgHACNodes, &numHACNodes, vectMergeIJs[ll*2],
			  vectMergeIJs[ll*2+1], 999999./*parm ignored*/,
			  rgLabelsTrain, rgTrainCostMatrix, numTrain,false);
	if((curLevel>=0) && (numHACNodes == rgLevelSizes[curLevel])){
	  printf("copying %d nodes for level %d\n",numHACNodes,curLevel);
	  for(int nd=0; nd < numHACNodes; ++nd){
	    printf(" %d=%d(%s)",nd,rgHACNodes[nd]->centerIdx,
		   rgLabelsTrain[rgHACNodes[nd]->centerIdx].c_str());
	    rgLevelHACNodes[curLevel][nd] = rgHACNodes[nd];
	  }
	  printf("\n");
	  --curLevel;
	}
	if(numHACNodes < 1){
	  fprintf(stderr,
		  "ERROR!!! did merges from file '%s' and ran out of nodes!\n",
		  stHACMergeFile);
	  exit(1);
	}
      }
      vectMergeIJs.clear();
      if(numHACNodes > 1){
	printf("there are still %d HAC nodes... appending to file '%s'\n",
	       numHACNodes, stHACMergeFile);
	fmerge = fopen(stHACMergeFile,"ab");
	if(!fmerge){
	  fprintf(stderr,"WARNING!!! couldn't open '%s' for append\n",
		  stHACMergeFile);
	}
      }
    }
    else{//couldn't open, so create it
      printf("merge file '%s' not found, so creating it...\n",stHACMergeFile);
      fmerge = fopen(stHACMergeFile,"wb");
      if(!fmerge){
	fprintf(stderr,"ERROR! couldn't open OR create merge file '%s'\n",
		stHACMergeFile);
	exit(1);
      }
      fprintf(fmerge,"numTrain %ld\n",(long)numTrain);
      fflush(fmerge);
    }
  }
  

  //cluster until all nodes are combined
  while(numHACNodes > 1){
    int iMinDist, jMinDist;
    double minDist;
    iMinDist = 0; jMinDist = 1; minDist = 999999.;
    for(int i=0; i < numHACNodes; ++i){
      for(int j=0; j < numHACNodes; ++j){
	if(i==j)
	  continue;
	double dist;
	dist = getClusterNodeDist(rgHACNodes[i], rgHACNodes[j],
				  rgTrainCostMatrix, numTrain);
	if(((i==0) && (j==1)) || (dist < minDist)){
	  minDist = dist;
	  iMinDist = i;
	  jMinDist = j;
	}
      }
    }
    //now merge the two clusters that are nearest each other (iMinDist,jMinDist)
    //the function will remove jMinDist and shorten the array length by 1
    mergeClusterNodes(rgHACNodes, &numHACNodes, iMinDist, jMinDist, minDist,
		      rgLabelsTrain, rgTrainCostMatrix, numTrain,true);
    if(fmerge){//save to merge file
      fprintf(fmerge,"%ld %ld %ld *\n",numHACMergeIdxsInFile,(long)iMinDist,
	      (long)jMinDist);
      fflush(fmerge);
      ++numHACMergeIdxsInFile;
    }
    if((curLevel>=0) && (numHACNodes == rgLevelSizes[curLevel])){
      printf("copying %d nodes for level %d\n",numHACNodes,curLevel);
      for(int nd=0; nd < numHACNodes; ++nd){
	printf(" %d=%d(%s)",nd,rgHACNodes[nd]->centerIdx,
	       rgLabelsTrain[rgHACNodes[nd]->centerIdx].c_str());
	rgLevelHACNodes[curLevel][nd] = rgHACNodes[nd];
      }
      printf("\n");
      --curLevel;
    }
    if(0 == numHACNodes%10)
      printf("   (total elapsed time so far for clustering is %.2lf seconds\n",
	     t1.getAccumulated());
  }
  t1.stop();
  printf("clustering took %.2lf seconds\n", t1.getAccumulated());
  if(fmerge){
    fclose(fmerge);
  }

  //count up the number of descendant tree nodes for each node in the tree
  printf("counting number of descendant tree nodes for each node in tree...");
  t1.start();
  for(int i=0; i < numTrain; ++i){
    HAC_TREE_NODE *np;
    np = rgOrigWordLeafNodes[i];
    while(np != NULL){
      np = np->pParent;
      if(np != NULL){
	np->numDescendantTreeNodes += 1;
      }
    }
  }
  t1.stop();
  printf("(%.2lf seconds)\n", t1.getAccumulated());

  //use the first level of HAC nodes as key word images for pseudo-features
  int numKeyWords;
  int *rgKeyWordIdxs;
  numKeyWords = rgLevelSizes[0];
  rgKeyWordIdxs = new int[numKeyWords];
  D_CHECKPTR(rgKeyWordIdxs);
#if DONT_USE_SMART_KEY_IMAGE_SELECTION
  for(int kk=0; kk < numKeyWords; ++kk){
    rgKeyWordIdxs[kk] = rgLevelHACNodes[0][kk]->centerIdx;
    printf("rgKeyWordIdxs[%d]=%d(%s)\n",kk,rgKeyWordIdxs[kk],
	   rgLabelsTrain[rgLevelHACNodes[0][kk]->centerIdx].c_str());
  }
  for(int kk=0; kk < numKeyWords; ++kk){
    printf("I%d [color=red];\n",rgLevelHACNodes[0][kk]->clustID);
  }
#else /*smart selection of key images*/
  {
    int maxNumSamples = 200;
    int numSamples = 0;
    int *rgSampleIdxs;
    rgSampleIdxs = new int[maxNumSamples];
    D_CHECKPTR(rgSampleIdxs);
    memset(rgSampleIdxs,0, numSamples*sizeof(int));

 #if USE_MOST_FREQUENT_LABELS_AS_SAMPLES
    std::map<std::string,int> mapLabelCounts;
    for(int i=0; i < numTrain; ++i){
      mapLabelCounts[rgLabelsTrain[i]] += 1;
    }
    die();
    abort();
 #else
    std::map<std::string, int> mapLabelsSampled;
    std::map<std::string, int>::iterator mapiter;
    std::string strTmp;
    mapLabelsSampled.clear();
    for(int i=0; (i < numTrain) && (numSamples < maxNumSamples); ++i){
      mapiter = mapLabelsSampled.find(rgLabelsTrain[i]);
      if(mapLabelsSampled.end() == mapiter){
	strTmp = rgLabelsTrain[i];
	mapLabelsSampled[strTmp] = 1;
	rgSampleIdxs[numSamples] = i;
	printf(" sample%d=rgLabelsTrain[%d](%s)\n",numSamples,i,
	       rgLabelsTrain[i].c_str());
	++numSamples;
      }
    }
 #endif   
    //rgKeyWordIdxs[0]: select sample with greatest variance in dist to other
    //sample images (for each row k, find mean_k then var_k. rgKeyWordIdxs[0]
    //is the k for which var_k is highest
    double bestVar;
    int bestVarSampleIdx;
    bestVar = -1.;
    bestVarSampleIdx = 0;
    
    double *rgDists;
    rgDists = new double[numSamples];
    D_CHECKPTR(rgDists);
    for(int samp=0; samp<numSamples; ++samp){
      double sampDistVar;//variance of distances from samp to other samples (ss)
      for(int ss=0,distnum=0; ss<numSamples; ++ss){
	if(ss==samp)
	  continue;
	rgDists[distnum] =
	  rgTrainCostMatrix[(long int)rgSampleIdxs[ss]*(long)numTrain+
			    (long int)rgSampleIdxs[samp]];
	++distnum;
      }
      sampDistVar = DMath::variance(rgDists,numSamples-1);
      printf(" variance%d=%.2lf\n",samp,sampDistVar);
      if(sampDistVar > bestVar){
	bestVar = sampDistVar;
	bestVarSampleIdx = samp;
      }
    }
    rgKeyWordIdxs[0] = rgSampleIdxs[bestVarSampleIdx];
    printf("***rgKeyWordIdxs[0]=samp%d trainIdx=%d(%s)  variance=%.2lf\n",
	   bestVarSampleIdx, rgKeyWordIdxs[0],
	   rgLabelsTrain[rgKeyWordIdxs[0]].c_str(), bestVar);

    //rgKeyWordIdxs[1..K]: select sample image, samp, that increases variance
    //for all other sample image vectors when samp is used as key image k
    bool *rgfSampUsed;
    rgfSampUsed = new bool[numSamples];
    D_CHECKPTR(rgfSampUsed);
    double **rgFV_data;
    rgFV_data = new double*[numSamples];
    D_CHECKPTR(rgFV_data);
    for(int ss=0; ss < numSamples; ++ss){
      rgFV_data[ss] = new double[numSamples];
      D_CHECKPTR(rgFV_data[ss]);
      rgFV_data[ss][0] = rgTrainCostMatrix[rgKeyWordIdxs[0]*(long)numTrain+
					  rgSampleIdxs[ss]];
      rgfSampUsed[ss] = false;
    }
    rgfSampUsed[bestVarSampleIdx] = true;
    for(int k=1; k < numKeyWords; ++k){//choose the rest of the keywords
      int bestSampIdx;
      double bestSampFVvar;//variance in distance in feature space
      DFeatureVector fvMeanTmp;
      DFeatureVector fvTmp;
      bestSampIdx = -1;
      bestSampFVvar = 0.;

      fvTmp.setData_dbl(rgFV_data[0],k+1,1,true,false/*no copy*/,true);
      fvMeanTmp.setData_dbl(rgFV_data[0],k+1,1,true,true/*copy*/,true);
      for(int samp=0; samp < numSamples; ++samp){
	if(rgfSampUsed[samp])
	  continue;//only consider samp if it isn't already used as key image
	fvMeanTmp.setValuesToZero();
	for(int ss=0; ss < numSamples; ++ss){
	  rgFV_data[ss][k]=rgTrainCostMatrix[rgSampleIdxs[samp]*(long)numTrain+
					     rgSampleIdxs[ss]];
	  fvTmp.pDbl = rgFV_data[ss];
	  fvMeanTmp.add(fvTmp);
	}
	fvMeanTmp.divideByScalar((double)numSamples);
	double varianceTmp;
	varianceTmp = 0.;
	for(int ss=0; ss < numSamples; ++ss){
	  double dblDistFVTmp;
	  fvTmp.pDbl = rgFV_data[ss];
	  dblDistFVTmp = DFeatureVector::getEuclideanDist(fvTmp,fvMeanTmp);
	  varianceTmp += (dblDistFVTmp * dblDistFVTmp);
	}
	varianceTmp /= numSamples;
	if(varianceTmp > bestSampFVvar){
	  bestSampFVvar = varianceTmp;
	  bestSampIdx = samp;
	}
      }//end for samp
      if(-1==bestSampIdx){
	fprintf(stderr,"this shouldn't happen! bestSampIdx ==-1. are there less samples than key word images?\n");
	exit(1);
      }
      rgKeyWordIdxs[k] = rgSampleIdxs[bestSampIdx];
      rgfSampUsed[bestSampIdx] = true;
      for(int ss=0; ss < numSamples; ++ss){
	rgFV_data[ss][k] = rgTrainCostMatrix[rgKeyWordIdxs[k]*(long)numTrain+
					    rgSampleIdxs[ss]];
      }
      printf("***rgKeyWordIdxs[%d]=samp%d trainIdx=%d(%s)  variance=%.2lf\n",
	     k,bestSampIdx, rgKeyWordIdxs[k],
	     rgLabelsTrain[rgKeyWordIdxs[k]].c_str(), bestSampFVvar);
    }//end for(k=1; k < numKeyWords...
    delete [] rgDists;
    delete [] rgfSampUsed;
    for(int i=0; i < numSamples; ++i)
      delete [] rgFV_data[i];
    delete [] rgFV_data;

    delete [] rgSampleIdxs;
  }
#endif // end #else (do use SMART_KEY_IMAGE_SELECTION)
  // exit(1);

  for(int lvl=0; lvl < numLevels; ++lvl){
    delete [] rgLevelHACNodes[lvl];
    rgLevelHACNodes[lvl] = NULL;
  }
  delete [] rgLevelHACNodes;


  //calculate pseudo distances and max pseudo distances for each tree node
  DTimer t4;
  t4.start();
  printf("calling calculateTreeKeyVectors()...");fflush(stdout);
  calculateTreeKeyVectors(rgHACNodes[0],rgKeyWordIdxs, numKeyWords,
			  rgTrainCostMatrix, numTrain, rgOrigWordLeafNodes);
  t4.stop();
  printf("took %.2f seconds\n",t4.getAccumulated());

  //printf("NOT saving tree graphviz diagram or 3d matlab word plot\n");
 // saveTreeForGraphviz("/tmp/tree_graphviz_file.dot",rgHACNodes[0],rgLabelsTrain,
 // 		      		      rgTrainCostMatrix, numTrain, rgKeyWordIdxs, numKeyWords);
 //   printf("to create the tree image do 'dot -Tps /tmp/tree_graphviz_file.dot -o /tmp/tree_graphviz_file.ps'\n");




  HAC_TREE_NODE **rgMatchNodes;
  rgMatchNodes = new HAC_TREE_NODE*[numTest];
  D_CHECKPTR(rgMatchNodes);
  int *rgNumMorphCompares;
  rgNumMorphCompares = new int[numTest];
  D_CHECKPTR(rgNumMorphCompares);
  bool *rgfCorrect;
  rgfCorrect = new bool[numTest];
  D_CHECKPTR(rgfCorrect);
  int *rgNumPrunedFV;
  int *rgNumPrunedChildFV;
  rgNumPrunedFV = new int[numTest];
  D_CHECKPTR(rgNumPrunedFV);
  rgNumPrunedChildFV = new int[numTest];
  D_CHECKPTR(rgNumPrunedChildFV);

  t1.start();
  printf("Doing threaded tree search to classify each test word\n");
  for(int tnum=numThreads-1; tnum >=0; --tnum){//launch threads reverse order
    rgTreeThreadParms[tnum].numThreads = numThreads;
    rgTreeThreadParms[tnum].threadNum = tnum;
    rgTreeThreadParms[tnum].rgTrainingImages = rgTrainingImages;
    rgTreeThreadParms[tnum].numTrain = numTrain;
    rgTreeThreadParms[tnum].testFirst = testFirst;
    rgTreeThreadParms[tnum].testLast = testLast;
    rgTreeThreadParms[tnum].rgCostsMorph = rgCostsMorph;
    rgTreeThreadParms[tnum].rgMatchNodes = rgMatchNodes;
    rgTreeThreadParms[tnum].rgNumMorphCompares = rgNumMorphCompares;
    rgTreeThreadParms[tnum].weightMovement = weightMovement;
    rgTreeThreadParms[tnum].lengthPenalty = lengthPenalty;
    rgTreeThreadParms[tnum].meshSpacingStatic = meshSpacingStatic;
    rgTreeThreadParms[tnum].numRefinesStatic = numRefinesStatic;
    rgTreeThreadParms[tnum].meshDiv = meshDiv;
    rgTreeThreadParms[tnum].bandWidthDP = bandWidth;
    rgTreeThreadParms[tnum].treeRoot = rgHACNodes[0];
    rgTreeThreadParms[tnum].alpha = alpha;
    rgTreeThreadParms[tnum].stPathIn = stPathIn;//don't alloc/copy, just point
    rgTreeThreadParms[tnum].rgLabelsTest = rgLabelsTest;
    rgTreeThreadParms[tnum].rgLabelsTrain = rgLabelsTrain;
    rgTreeThreadParms[tnum].rgfCorrect = rgfCorrect;
    rgTreeThreadParms[tnum].rgNumPrunedFV = rgNumPrunedFV;
    rgTreeThreadParms[tnum].rgNumPrunedChildFV = rgNumPrunedChildFV;
    rgTreeThreadParms[tnum].rgKeyWordIdxs = rgKeyWordIdxs;
    rgTreeThreadParms[tnum].numKeyWords = numKeyWords;
    rgTreeThreadParms[tnum].rgOrigWordLeafNodes = rgOrigWordLeafNodes;
    rgTreeThreadParms[tnum].slowPassTopN = slowPassTopN;
    rgTreeThreadParms[tnum].rgTopNMatches = rgTopNMatches;
    rgTreeThreadParms[tnum].topN = topNMatches;
#ifdef D_NOTHREADS
      tree_search_thread_func(rgTreeThreadParms);
#else
      if(0 == tnum){//don't spawn thread zero. use the current thread.
	tree_search_thread_func(rgTreeThreadParms);
      }
      else{//spawn all other threads besides zero
	if(0 !=pthread_create(&rgThreadID[tnum], NULL,
			      tree_search_thread_func,
			      &rgTreeThreadParms[tnum])){
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
  printf("Tree search to classify all test words took %.02f seconds\n",
	 t1.getAccumulated());
  
  int numCorrect = 0;
  int sumNumMorphCompares;
  int sumNumPrunedFV;
  int sumNumPrunedChildFV;
  sumNumPrunedFV = 0;
  sumNumPrunedChildFV =0;
  numCorrect = 0;
  sumNumMorphCompares = 0;
  for(int i=0; i < numTest; ++i){
    if(rgfCorrect[i])
      ++numCorrect;
    sumNumMorphCompares += rgNumMorphCompares[i];
    sumNumPrunedFV += rgNumPrunedFV[i];
    sumNumPrunedChildFV += rgNumPrunedChildFV[i];
  }
  printf("numCorrect=%d out of %d (%.2lf%%)\n",numCorrect, numTest, (100.*numCorrect)/(double)numTest);
  printf("avg number of morph compares=%.2lf\n",
	 sumNumMorphCompares/(double)numTest);
  printf("avgNumPrunedFV=%.2lf avgNumPrunedChildFV=%.2lf (these include all of their pruned descendants)\n",
	 sumNumPrunedFV/(double)numTest, sumNumPrunedChildFV/(double)numTest);

// typedef struct{
//   int trIdx;//idx of training example
//   double cost; //cost from current word to training word trIdx
// } TOPN_MATCHES_T;

 //  rgTopNMatches = new TOPN_MATCHES_T[(long)numTest*topNMatches];



  //now output the costs for analysis
  FILE *fout;
  printf("saving results to '%s'\n",stOutfile);
  fout = fopen(stOutfile, "wb");
  if(!fout){
    fprintf(stderr, "couldn't open output file '%s\n",stOutfile);
    exit(1);
  }
  fprintf(fout,"Top %d recognition results for each test word. TRACK_THE_BEST_N=%d BEST_N_UNIQUE_LABELS=%d\n",topNMatches,TRACK_THE_BEST_N,BEST_N_UNIQUE_LABELS);

  fprintf(fout,"#%s trained on w_%08d.pgm through w_%08d.pgm\n",stPathIn,
	  trainFirst, trainLast);
  fprintf(fout,"trainFirst: %08d trainLast: %08d\n",trainFirst,trainLast);
  fprintf(fout,"testFirst: %08d testLast: %08d\n",testFirst,testLast);
  for(int tt=testFirst, i=0; tt <= testLast; ++tt,++i){
    fprintf(fout,"#%d ------\n#%s\n",tt, rgLabelsTest[i].c_str());
    fprintf(fout,"%d :\n",tt);
    for(int topn=0; topn < topNMatches; ++topn){
      int trIdx;
      trIdx = rgTopNMatches[i*(long)topNMatches+topn].trIdx;
      fprintf(fout,"   %d %f #%s\n",
	      trIdx+trainFirst,
	      rgTopNMatches[i*(long)topNMatches+topn].cost,
	      rgLabelsTrain[trIdx].c_str());
    }
  }
  fclose(fout);

#ifndef D_NOTHREADS
  delete [] rgThreadID;
#endif
  delete [] rgMatchNodes;
  delete [] rgNumMorphCompares;
  delete [] rgfCorrect;
  delete [] rgNumPrunedFV;
  delete [] rgNumPrunedChildFV;

  free(rgTreeThreadParms);
  free(rgTrainThreadParms);
  free(rgThreadParms);
  free(rgCostsMorph);
  delete [] rgTrainingImages;
  delete [] rgLabelsTrain;
  delete [] rgLabelsTest;
  delete [] rgTrainCostMatrix;
  delete [] rgKeyWordIdxs;
  delete [] rgTopNMatches;
  delete [] rgOrigWordLeafNodes;
  delete rgHACNodes[0];
  delete [] rgHACNodes;
  return 0;
}
