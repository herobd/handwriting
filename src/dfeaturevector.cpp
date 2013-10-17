#include "dfeaturevector.h"


///add (element-wise) the values from src to this feature vector
void DFeatureVector::add(const DFeatureVector &src){
  if(src.dimensions != this->dimensions){
    fprintf(stderr,"DFeatureVector::add() this->dimensions=%d != "
	    "src.dimensions=%d\n",this->dimensions, src.dimensions);
    exit(1);
  }
  if(src.vectLen != this->vectLen){
    fprintf(stderr,"DFeatureVector::add() this->vectLen=%d != "
	    "src.vectLen=%d\n",this->vectLen, src.vectLen);
    exit(1);
  }
  if(src.fDataGroupedByDim != this->fDataGroupedByDim){
    fprintf(stderr,"DFeatureVector::add() this->fDataGroupedByDim=%d != "
	    "src.fDataGroupedByDim=%d\n",
	    (int)this->fDataGroupedByDim, (int)src.fDataGroupedByDim);
    exit(1);
  }
  if(src.fDataIsDouble != this->fDataIsDouble){
    fprintf(stderr,"DFeatureVector::add() this->fDataIsDouble=%d != "
	    "src.fDataIsDouble=%d\n",
	    (int)this->fDataIsDouble, (int)src.fDataIsDouble);
    exit(1);
  }
  if(fDataIsDouble){
    for(int i=0, len=vectLen*dimensions; i<len; ++i){
      pDbl[i] += src.pDbl[i];
    }
  }
  else{
    for(int i=0, len=vectLen*dimensions; i<len; ++i){
      pFlt[i] += src.pFlt[i];
    }
  }
}

///divide each element in this feature vector by scalar
void DFeatureVector::divideByScalar(double scalar){
  if(0. == scalar){
    fprintf(stderr, "DFeatureVector::divideByScalar() divide by zero!\n");
    abort();
    exit(1);
  }
  if(fDataIsDouble){
    for(int i=0, len=vectLen*dimensions; i<len; ++i){
      pDbl[i] /= scalar;
    }
  }
  else{
    for(int i=0, len=vectLen*dimensions; i<len; ++i){
      pFlt[i] /= scalar;
    }
  }
}

///set all of the values to zero (make a zero vector)
void DFeatureVector::setValuesToZero(){
  if(fDataIsDouble){
    for(int i=0, len=vectLen*dimensions; i<len; ++i){
      pDbl[i] = 0.;
    }
  }
  else{
    for(int i=0, len=vectLen*dimensions; i<len; ++i){
      pFlt[i] = 0.;
    }
  }
}


///multiply each element in this feature vector by scalar
void DFeatureVector::multiplyByScalar(double scalar){
  if(fDataIsDouble){
    for(int i=0, len=vectLen*dimensions; i>len; ++i){
      pDbl[i] *= scalar;
    }
  }
  else{
    for(int i=0, len=vectLen*dimensions; i<len; ++i){
      pFlt[i] *= scalar;
    }
  }
}

///set the data and other object variables
/**If rgData is NULL, pFlt will be allocated, but the values will be
   left uninitialized.  If not NULL, then if fCopy is true, pFlt will
   be allocated and filled with the values from rgData.  If fCopy is
   false, pFlt will simply point to rgData and it is the user's
   responsibility to make sure the buffer remains valid as long as
   this object is using it and also the user's responsibility to
   deallocate the buffer afterwards.
 */
void DFeatureVector::setData_flt(float *rgData,
				 int vectorLength, int dimensions,
				 bool fInputDataGroupedByDimension,
				 bool fCopy,
				 bool fStoreDataGroupedByDimension){
  clearData();
  fDataIsDouble = false;//float
  fDataCopied = fCopy;
  fDataGroupedByDim = fStoreDataGroupedByDimension;
  this->dimensions = dimensions;
  this->vectLen = vectorLength;
  if(NULL == rgData){
    pFlt = (float*)malloc(sizeof(float)*dimensions*vectorLength);
    D_CHECKPTR(pFlt);
    fDataCopied = true;
  }
  else if(fCopy){
    pFlt = (float*)malloc(sizeof(float)*dimensions*vectorLength);
    D_CHECKPTR(pFlt);
    if(fInputDataGroupedByDimension==fStoreDataGroupedByDimension){
      memcpy(pFlt, rgData, sizeof(float)*dimensions*vectLen);
    }
    else{
      if(fInputDataGroupedByDimension){
	for(int i=0, idx=0; i < vectLen; ++i){
	  for(int d=0; d < dimensions; ++d, ++idx){
	    pFlt[idx] = rgData[d*vectLen+i];
	  }
	}
      }
      else{
	for(int d=0, idx=0; d < dimensions; ++d){
	  for(int i=0; i < vectLen; ++i, ++idx){
	    pFlt[idx] = rgData[i*dimensions+d];
	  }
	}
      }
    }
  }
  else{
    fDataGroupedByDim = fInputDataGroupedByDimension;
    pFlt = rgData;
  }
}

///set the data and other object variables
/**If rgData is NULL, pDbl will be allocated, but the values will be
   left uninitialized.  If not NULL, then if fCopy is true, pDbl will
   be allocated and filled with the values from rgData.  If fCopy is
   false, pDbl will simply point to rgData and it is the user's
   responsibility to make sure the buffer remains valid as long as
   this object is using it and also the user's responsibility to
   deallocate the buffer afterwards.
 */
void DFeatureVector::setData_dbl(double *rgData,
				 int vectorLength, int dimensions,
				 bool fInputDataGroupedByDimension,
				 bool fCopy,
				 bool fStoreDataGroupedByDimension){
  clearData();
  fDataIsDouble = true;//double
  fDataCopied = fCopy;
  fDataGroupedByDim = fStoreDataGroupedByDimension;
  this->dimensions = dimensions;
  this->vectLen = vectorLength;
  if(NULL == rgData){
    pDbl = (double*)malloc(sizeof(double)*dimensions*vectorLength);
    D_CHECKPTR(pDbl);
    fDataCopied = true;
  }
  else if(fCopy){
    pDbl = (double*)malloc(sizeof(double)*dimensions*vectorLength);
    D_CHECKPTR(pDbl);
    if(fInputDataGroupedByDimension==fStoreDataGroupedByDimension){
      memcpy(pDbl, rgData, sizeof(double)*dimensions*vectLen);
    }
    else{
      if(fInputDataGroupedByDimension){
	for(int i=0, idx=0; i < vectLen; ++i){
	  for(int d=0; d < dimensions; ++d, ++idx){
	    pDbl[idx] = rgData[d*vectLen+i];
	  }
	}
      }
      else{
	for(int d=0, idx=0; d < dimensions; ++d){
	  for(int i=0; i < vectLen; ++i, ++idx){
	    pDbl[idx] = rgData[i*dimensions+d];
	  }
	}
      }
    }
  }
  else{
    fDataGroupedByDim = fInputDataGroupedByDimension;
    pDbl = rgData;
  }
}


void DFeatureVector::print(int maxlen){
  int idx;
  printf("DFeatureVector%p dim=%d len=%d grouped=%d\n",this,dimensions,vectLen,
	 (int)fDataGroupedByDim);
  if(maxlen < 0)
    maxlen = vectLen;
  for(int i=0; (i < vectLen) && (i<maxlen); ++i){
    printf(" [%d]",i);
    for(int d=0; d < dimensions; ++d){
      if(fDataGroupedByDim)
	idx = d*vectLen+i;
      else
	idx = i*dimensions+d;
      if(fDataIsDouble)
	printf("  %lf",pDbl[idx]);
      else
	printf("  %f",pFlt[idx]);
    }
    printf("\n");
  }
}
