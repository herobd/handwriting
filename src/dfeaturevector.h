#ifndef DFEATUREVECTOR_H
#define DFEATUREVECTOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ddefs.h"
#include "dinstancecounter.h"

/** This class is meant as a container for float/double features to be
    used in clustering, dynamic programming, and classification.  If the
    dimension is one and the length of two DFeatureVectors is equal, a
    Euclidean distance can be taken between them for purposes such as
    clustering, nearest-neighbor classification, etc.

    The feature vector may also be composed of more than one feature
    per slot by specifying a dimension greater than one.  For example,
    I use it this way for Rath/Manmatha style word-level feature
    vectors that include the word profile, the upper/lower word
    profiles, and ink transition counts at each pixel column of the
    word.  While there is one DFeatureVector for each word (and the
    length of the vector depends on how many pixels wide the word
    image is), each column can actually be thought of as a vector
    within that vector, having the values of each feature
    corresponding to that pixel column of the word image.  Instead of
    a direct Euclidean distance between two DFeatureVectors, I can
    then do a Dynamic Programming (Dynamic Time Warping) of one to
    the other to align them (thus accounting for variability/
    stretching) and use that alignment cost as a matching cost, or
    align them and do further processing such as with DMorphInk to
    compute some other matching cost for use in clustering /
    classification.

    Keep in mind that if dimension is greater than 1, then you must make
    sure that you have put it all into a single data array of size
    dimensions*vectLen and it should be packed correctly (all of dim 0,
    then all of dim 1, and so on if fInputDataGroupedByDimension=true,
    otherwise dim0...dimN for column 0, then dim0...dimN for column 2, and
    so on).  If fCopy is true, then the internal storage order will be
    whichever is specified by fStoreDataGroupedByDimension.
    
*/




class DFeatureVector{
public:
  DFeatureVector();
  DFeatureVector(const DFeatureVector &src);
  ~DFeatureVector();
  const DFeatureVector& operator=(const DFeatureVector &src);
  void add(const DFeatureVector &src);
  void divideByScalar(double scalar);
  void multiplyByScalar(double scalar);
  void setValuesToZero();
  void clearData();
  void setData_flt(float *rgData, int vectorLength, int dimensions = 1,
		   bool fInputDataGroupedByDimension = true,
		   bool fCopy = true,
		   bool fStoreDataGroupedByDimension = true);
  void setData_dbl(double *rgData, int vectorLength, int dimensions = 1,
		   bool fInputDataGroupedByDimension = true,
		   bool fCopy = true,
		   bool fStoreDataGroupedByDimension = true);
  void print(int maxlen=-1);
  static double getEuclideanDist(const DFeatureVector &fv1,
				 const DFeatureVector &fv2);
  static double getEuclideanDist_nochecks(const DFeatureVector &fv1,
					  const DFeatureVector &fv2);
  static double getEuclideanDistSquared(const DFeatureVector &fv1,
				 const DFeatureVector &fv2);
  static double getEuclideanDistSquared_nochecks(const DFeatureVector &fv1,
					  const DFeatureVector &fv2);
  bool fDataIsDouble;
  bool fDataGroupedByDim;
  bool fDataCopied;
  double *pDbl;
  float *pFlt;
  int dimensions;
  int vectLen;
};
inline DFeatureVector::DFeatureVector(){
  DInstanceCounter::addInstance("DFeatureVector");
  pDbl = NULL;
  pFlt = NULL;
  fDataIsDouble = fDataCopied  = fDataGroupedByDim = false;
  dimensions = vectLen = 0;
}
///makes a deep copy
inline DFeatureVector::DFeatureVector(const DFeatureVector &src){
  DInstanceCounter::addInstance("DFeatureVector");
  fDataIsDouble = src.fDataIsDouble;
  fDataGroupedByDim = src.fDataGroupedByDim;
  fDataCopied = true;
  dimensions = src.dimensions;
  vectLen = src.vectLen;
  if(src.pDbl != NULL){
    pDbl = (double*)malloc(sizeof(double)*dimensions*vectLen);
    D_CHECKPTR(pDbl);
    memcpy(pDbl, src.pDbl, sizeof(double)*dimensions*vectLen);
  }
  if(src.pFlt != NULL){
    pFlt = (float*)malloc(sizeof(float)*dimensions*vectLen);
    D_CHECKPTR(pFlt);
    memcpy(pFlt, src.pFlt, sizeof(float)*dimensions*vectLen);
  }
}
inline void DFeatureVector::clearData(){
  if(fDataCopied){
    if(NULL!=pDbl)
      free(pDbl);
    if(NULL!=pFlt)
      free(pFlt);
  }
  pFlt = NULL;
  pDbl = NULL;
}
inline DFeatureVector::~DFeatureVector(){
  DInstanceCounter::removeInstance("DFeatureVector");
  clearData();
  dimensions = vectLen = 0;
  fDataIsDouble = fDataGroupedByDim = fDataCopied = false;
}
///does a deep copy of the data (unless src==this)
inline const DFeatureVector& DFeatureVector::operator=(const DFeatureVector &src){
  if(this != &src){
    //first deallocate old data so we don't leak memory
    clearData();
    //now copy data from src
    fDataIsDouble = src.fDataIsDouble;
    fDataGroupedByDim = src.fDataGroupedByDim;
    fDataCopied = true;
    dimensions = src.dimensions;
    vectLen = src.vectLen;
    if(src.pDbl != NULL){
      pDbl = (double*)malloc(sizeof(double)*dimensions*vectLen);
      D_CHECKPTR(pDbl);
      memcpy(pDbl, src.pDbl, sizeof(double)*dimensions*vectLen);
    }
    if(src.pFlt != NULL){
      pFlt = (float*)malloc(sizeof(float)*dimensions*vectLen);
      D_CHECKPTR(pFlt);
      memcpy(pFlt, src.pFlt, sizeof(float)*dimensions*vectLen);
    }
  }
  return *this;
}

///same as getEuclideanDistSquared() but doesn't do any sanity checks
inline double DFeatureVector::getEuclideanDistSquared_nochecks(const DFeatureVector &fv1,const DFeatureVector &fv2){
  double dblSum = 0.;
  if(fv1.fDataIsDouble){
    for(int i = 0, len=fv1.vectLen; i < len; ++i){
      dblSum += (fv1.pDbl[i] - fv2.pDbl[i])*(fv1.pDbl[i] - fv2.pDbl[i]);
    }
  }
  else{
    for(int i = 0, len=fv1.vectLen; i < len; ++i){
      dblSum += (double)(fv1.pFlt[i] - fv2.pFlt[i])*(fv1.pFlt[i] - fv2.pFlt[i]);
    }
  }
  return dblSum;
}

///same as getEuclideanDist() but doesn't do any sanity checks
inline double DFeatureVector::getEuclideanDist_nochecks(const DFeatureVector
							&fv1,
							const DFeatureVector
							&fv2){
  return sqrt(DFeatureVector::getEuclideanDistSquared_nochecks(fv1,fv2));
}
///Returns the Squared Euclidean distance between the two DFeatureVectors
/**This function requires that the feature vectors be dimension 1 and the same length.*/
inline double DFeatureVector::getEuclideanDistSquared(const DFeatureVector &fv1,
						      const DFeatureVector &fv2){
  if(fv1.vectLen != fv2.vectLen){
    fprintf(stderr, "DFeatureVector::getEuclideanDist() vector lengths don't "
	    "match! (%d!=%d)\n",fv1.vectLen, fv2.vectLen);
    return HUGE_VAL;
  }
  if((fv1.dimensions!=1)||(fv2.dimensions!= 1)){
    fprintf(stderr, "DFeatureVector::getEuclideanDist() vector dimensions must "
	    "both be 1! (dimensions are %d and %d)\n",
	    fv1.dimensions, fv2.dimensions);
    return HUGE_VAL;
  }
  // if(fv1.fDataGroupedByDim != fv2.fDataGroupedByDim){
  //   fprintf(stderr, "DFeatureVector::getEuclideanDist() vector data must "
  // 	    "be grouped/interleaved the same way for each vector\n");
  //   return HUGE_VAL;
  // }
  if(fv1.fDataIsDouble != fv2.fDataIsDouble){
    fprintf(stderr, "DFeatureVector::getEuclideanDist() data types must "
	    "match\n");
    return HUGE_VAL;
  }
  if(fv1.vectLen <=0){
    fprintf(stderr, "DFeatureVector::getEuclideanDist() length must be > 0\n ");
    return HUGE_VAL;
  }
  return DFeatureVector::getEuclideanDistSquared_nochecks(fv1,fv2);
}
///Returns the Euclidean distance between the two DFeatureVectors
/**This function requires that the feature vectors be dimension 1 and the same length.*/
inline double DFeatureVector::getEuclideanDist(const DFeatureVector &fv1,
					       const DFeatureVector &fv2){
  return sqrt(DFeatureVector::getEuclideanDistSquared(fv1,fv2));
}

#endif
