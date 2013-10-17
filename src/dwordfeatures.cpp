#include <string.h>
#include "dwordfeatures.h"



///Extract word-level features based on those used by Rath, Manmatha
/**This method assumes that img has already been properly prepared
   including: deskewed, sized, background removed, noise removed,slant
   removed, cropped.  The resulting feature vector will have 4
   dimensions: the (grayscale) profile, upper profile, lower profile,
   and black-white transition counts.  The length of the resulting
   feature vector is the width of img.  The transition counts are in
   the opposite direction of those used by Rath/Manmatha if
   fUseMyTransitions is true.  If fInkIsBlack is true then dark colors
   are foreground and light is BG.  If fRangeIs255 is false, then it
   is assumed that the range is 0-1 instead of 0-255.  If the image is
   grayscale instead of bitonal, pixels less than or equal to tval
   (default is 127) will be treated as if zero and pixels above tval
   will be treated as 255.  Features will be of type double.  If
   floats are desired, call extractWordFeatures_flt() instead.  */
DFeatureVector DWordFeatures::extractWordFeatures(DImage &img,
						  bool fUseGrayscaleProf,
						  bool fUseMyTransitions,
						  bool fInkIsBlack,
						  bool fRangeIs255,
						  D_uint8 tval,
						  double weight_prof,
						  double weight_upper,
						  double weight_lower,
						  double weight_trans){
  int w, h;
  DFeatureVector fv;
  w = img.width();
  h = img.height();
  double *pProf, *pUpper, *pLower, *pTrans;//start of each feature within data
  D_uint8 *p8;
  
  if(DImage::DImage_u8 != img.getImageType()){
    fprintf(stderr, "DWordFeatures::extractWordFeatures() only supports 8-bit "
	    "grayscale data\n");
    exit(1);
  }
  //passing NULL in here means to allocate space, but leave it uninitialized
  fv.setData_dbl(NULL/*allocate*/, w, 4, true, true, true);
  pProf = fv.pDbl;
  pUpper = &(pProf[w]);
  pLower = &(pProf[w*2]);
  pTrans = &(pProf[w*3]);
  memset(pProf,0,sizeof(double)*w*4);

  p8 = img.dataPointer_u8();

  //TODO: instead of making four passes through the image (one for each feature,
  //I could combine them all into a single loop.  I'll leave it separate for
  //now for debug purposes and because this isn't a bottleneck anyway.


  // Feature 0: Profile
  if(fUseGrayscaleProf){
    if(!fRangeIs255){
      fprintf(stderr,"DWordFeatures::extractWordFeatures() if "
	      "fUseGrayscaleProf is true, then fRangeIs255 should be too!\n");
      exit(1);
    }
    if(fInkIsBlack){
      for(int y=0, idx=0; y < h; ++y){
	for(int x=0; x < w; ++x, ++idx){
	  pProf[x] += 255.-(double)(p8[idx]);
	}
      }
    }//if fInkIsBlack
    else{
      for(int y=0, idx=0; y < h; ++y){
	for(int x=0; x < w; ++x, ++idx){
	  pProf[x] += (double)(p8[idx]);
	}
      }
    }
  }//if fUseGrayscaleProf
  else{
    if(fRangeIs255){
      if(fInkIsBlack){
	for(int y=0, idx=0; y < h; ++y){
	  for(int x=0; x < w; ++x, ++idx){
	    if(p8[idx] <= tval)
	      pProf[x] += 1;
	  }
	}
      }
      else{
	for(int y=0, idx=0; y < h; ++y){
	  for(int x=0; x < w; ++x, ++idx){
	    if(p8[idx] > tval)
	      pProf[x] += 1;
	  }
	}
      }
    }//end if fRangeIs255
    else{
      fprintf(stderr, "DWordFeatures::extractWordFeatures() NYI for fRangeIs255=false\n");
      exit(1);
    }
  }
  

  // Feature 1: Upper Profile
  for(int x=0; x < w; ++x)
    pUpper[x] = 999999.;

  if(fRangeIs255){
    if(fInkIsBlack){
      for(int y=0, idx=0; y < h; ++y){
	for(int x=0; x < w; ++x, ++idx){
	  if((p8[idx] <= tval)&&(y<pUpper[x]))
	    pUpper[x] = y;
	}
      }
    }
    else{
      for(int y=0, idx=0; y < h; ++y){
	for(int x=0; x < w; ++x, ++idx){
	  if((p8[idx] > tval)&&(y<pUpper[x]))
	    pUpper[x] = y;
	}
      }
    }
  }//end if fRangeIs255
  else{
    fprintf(stderr, "DWordFeatures::extractWordFeatures() NYI for fRangeIs255=false\n");
    exit(1);
  }
  // interpolate values where there is no ink in the column
  {
    int idxPrev, idxNext;
    for(int x=0; x < w; ++x){
      if(pUpper[x] == 999999.){
	idxPrev = x-1;
	idxNext = x+1;
	double prevVal, nextVal;//values to interpolate
	if(x>0)
	  prevVal = pUpper[x-1];
	while((idxNext < w)&& (pUpper[idxNext] == 999999.))
	  ++idxNext;
	if(idxNext >= w){//we went past the right end
	  if(idxPrev<0){//this shouldn't happen.  There must be no ink!
	    fprintf(stderr, "WARNING: DWordFeatures::extractWordFeatures() "
		    "found no ink for upper profile!\n");
	    break;
	  }
	  nextVal = prevVal;
	}
	else
	  nextVal = pUpper[idxNext];
	if(x<=0)
	  prevVal = nextVal;
	for(int i= idxPrev+1; i < idxNext; ++i){
	  double s;// distance along the segment between idxPrev and idxNext
	  s = (i-idxPrev) / (idxNext-idxPrev);
	  pUpper[i] = (1.0-s)*prevVal + s*nextVal;
	}
	x=idxNext-1;
      }
    }
  }




  // Feature 3: Lower Profile
  for(int x=0; x < w; ++x)
    pLower[x] = 999999.;

  if(fRangeIs255){
    if(fInkIsBlack){
      for(int y=0, idx=0; y < h; ++y){
	for(int x=0; x < w; ++x, ++idx){
	  if((p8[idx] <= tval)&&((h-y)<pLower[x]))
	    pLower[x] = (h-y);
	}
      }
    }
    else{
      for(int y=0, idx=0; y < h; ++y){
	for(int x=0; x < w; ++x, ++idx){
	  if((p8[idx] > tval)&&((h-y)<pLower[x]))
	    pLower[x] = (h-y);
	}
      }
    }
  }//end if fRangeIs255
  else{
    fprintf(stderr, "DWordFeatures::extractWordFeatures() NYI for fRangeIs255=false\n");
    exit(1);
  }
  // interpolate values where there is no ink in the column
  {
    int idxPrev, idxNext;
    for(int x=0; x < w; ++x){
      if(pLower[x] == 999999.){
	idxPrev = x-1;
	idxNext = x+1;
	double prevVal, nextVal;//values to interpolate
	if(x>0)
	  prevVal = pLower[x-1];
	while((idxNext < w)&& (pLower[idxNext] == 999999.))
	  ++idxNext;
	if(idxNext >= w){//we went past the right end
	  if(idxPrev<0){//this shouldn't happen.  There must be no ink!
	    fprintf(stderr, "WARNING: DWordFeatures::extractWordFeatures() "
		    "found no ink for upper profile!\n");
	    break;
	  }
	  nextVal = prevVal;
	}
	else
	  nextVal = pLower[idxNext];
	if(x<=0)
	  prevVal = nextVal;
	for(int i= idxPrev+1; i < idxNext; ++i){
	  double s;// distance along the segment between idxPrev and idxNext
	  s = (i-idxPrev) / (idxNext-idxPrev);
	  pLower[i] = (1.0-s)*prevVal + s*nextVal;
	}
	x=idxNext-1;
      }
    }
  }


  // Feature 3: Background to Ink Transition Counts
  if(!fUseMyTransitions){//use Rath and Manmatha's instead
    if(fRangeIs255){
      if(fInkIsBlack){
	for(int x=0; x < w; ++x){//first row
	  if(p8[x] <= tval)
	    pTrans[x] = 1.0;
	}
	for(int y=1, idx=w; y < h; ++y){
	  for(int x=0; x < w; ++x, ++idx){
	    if((p8[idx] <= tval)&&(p8[idx-w]>tval))
	      pTrans[x] += 1.0;
	  }
	}
      }
      else{
	for(int x=0; x < w; ++x){//first row
	  if(p8[x] > tval)
	    pTrans[x] = 1.0;
	}
	for(int y=1, idx=w; y < h; ++y){
	  for(int x=0; x < w; ++x, ++idx){
	    if((p8[idx] > tval)&&(p8[idx-w]<=tval))
	      pTrans[x] += 1.0;
	  }
	}
      }
    }//end if fRangeIs255
    else{
      fprintf(stderr, "DWordFeatures::extractWordFeatures() NYI for fRangeIs255=false\n");
      exit(1);
    }
  }
  else{
    fprintf(stderr, "DWordFeatures::extractWordFeatures() NYI for fUseMyTransitions=true\n");
    exit(1);
  }

  

  //normalize the features and weight them appropriately----------
  //Rath/Manmatha normalized the grayscale profile, the upper profile,
  //and lower profile from 0 to 1 (based on max values), and divided
  //the transition counts by 6.
  double maxProf, maxUpper, maxLower;
  maxProf = maxUpper = maxLower = 0.;
  for(int x=0; x < w; ++x){
    if(pProf[x] > maxProf)
      maxProf = pProf[x];
  }
  if(maxProf > 0.){
    for(int x=0; x < w; ++x){
      pProf[x] /= maxProf;
    }
  }
  for(int x=0; x < w; ++x){
    if(pUpper[x] > maxUpper)
      maxUpper = pUpper[x];
  }
  if(maxUpper > 0.){
    for(int x=0; x < w; ++x){
      pUpper[x] /= maxUpper;
    }
  }
  for(int x=0; x < w; ++x){
    if(pLower[x] > maxLower)
      maxLower = pLower[x];
  }
  if(maxLower > 0.){
    for(int x=0; x < w; ++x){
      pLower[x] /= maxLower;
    }
  }
  for(int x=0; x < w; ++x){
    pTrans[x] /= 6.;
  }


  return fv;
}


///same as extractWordFeatures() but creates DFeatureVectors with float data
/**The main reason one might wish to use float instead of double data is that
   if dealing with large numbers of feature vectors (comparing each word to a big training set, for example), the amount of memory required is smaller.*/
DFeatureVector DWordFeatures::extractWordFeatures_flt(DImage &img,
						      bool fUseGrayscaleProf,
						      bool fUseMyTransitions,
						      bool fInkIsBlack,
						      bool fRangeIs255,
						      float weight_prof,
						      float weight_upper,
						      float weight_lower,
						      float weight_trans){
  fprintf(stderr,"DWordFeatures::extractWordFeatures_flt() NYI!\n");
  exit(1);
}

/**create a feature vector from the first few coefficients of the four word-level features like Manmatha did so the features can be quickly compared and easily clustered instead of having to do dynamic warping between all pairs, etc.*/
DFeatureVector DWordFeatures::extractWordFourierFeatures_flt(){
  fprintf(stderr,"DWordFeatures::extractWordFourierFeatures_flt() NYI!\n");
  exit(1);
}


