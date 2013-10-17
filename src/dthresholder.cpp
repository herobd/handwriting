#include <stdio.h>
#include <math.h>
#include <string.h>
#include "dthresholder.h"
#include "dprogress.h"
#include "dkernel2d.h"
#include "dconvolver.h"
#include "dconnectedcomplabeler.h"

#include "dtimer.h"//for debug.  this can be removed after debugging

#ifndef D_NOTHREADS
#include "dthreads.h"
#define MAX_CCTHRESHVAL_THREADS 64
//structure for passing parameters to the thread function
/// \cond
typedef struct{
  const DImage *pImgSrc;//pointer to the image (shared by all threads)
  int numThreads;//how many threads are processing
  int threadNum;//which thread number this is (0..numThreads-1)
  int *rgNumCCs;//all threads point to same array.
  int tOtsu;
  int otsuCCArea;
  int otsuNumCCs;
} CCTHRESHVAL_THREAD_PARMS;
/// \endcond
#endif /* D_NOTHREADS not defined */



/// perform a simple threshold at tval
/** pixels less than or equal to tval go to zero, pixels greater than
 * tval go to 255.  Source and destination images may both be the same
 * image since this is a simple threshold.
 */
void DThresholder::threshImage_(DImage &imgDst, const DImage &imgSrc,
				double tval){
  int w, h;
  int len;
  D_uint8 *pDst8;
  
  w = imgSrc.width();
  h = imgSrc.height();

  if( (&imgSrc) != (&imgDst))
    imgDst.create(imgSrc.width(), imgSrc.height(), DImage::DImage_u8,
		  1, imgSrc.getAllocMethod());
  else if(imgDst.getImageType() != DImage::DImage_u8){
    fprintf(stderr, "DThresholder::threshImage_() imgDst can't be the same "
	    "image as imgSrc unless it is of type DImage::DImage_u8\n");
    abort();
  }
  pDst8 = imgDst.dataPointer_u8();

  switch(imgSrc.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint8 *pSrc8;
	pSrc8 = imgSrc.dataPointer_u8();
	len = w * h;
	for(int idx = 0; idx < len; ++idx){
	  if(pSrc8[idx] <= tval)
	    pDst8[idx] = 0x00;
	  else
	    pDst8[idx] = 0xff;
	}
      }
      break;
    case DImage::DImage_u16:
      {
	D_uint16 *pSrc16;
	pSrc16 = imgSrc.dataPointer_u16();
	len = w * h;
	for(int idx = 0; idx < len; ++idx){
	  if(pSrc16[idx] <= tval)
	    pDst8[idx] = 0x00;
	  else
	    pDst8[idx] = 0xff;
	}
      }
      break;
    case DImage::DImage_u32:
      {
	D_uint32 *pSrc32;
	pSrc32 = imgSrc.dataPointer_u32();
	len = w * h;
	for(int idx = 0; idx < len; ++idx){
	  if(pSrc32[idx] <= tval)
	    pDst8[idx] = 0x00;
	  else
	    pDst8[idx] = 0xff;
	}
      }
      break;
    case DImage::DImage_RGB:
    case DImage::DImage_RGB_16:
      fprintf(stderr, "DThresholder::threshImage_() does not support RGB\n");
      exit(1);
      break;
    case DImage::DImage_flt_multi:
      {
	float *pflt;

	if(imgSrc.numChannels() > 1){
	  fprintf(stderr, "DThresholder::threshImage_() does not support "
		  "multiple-channel images\n");
	  exit(1);
	}
	pflt = imgSrc.dataPointer_flt();
	len = w * h;
	for(int idx = 0; idx < len; ++idx){
	  if(pflt[idx] <= tval)
	    pDst8[idx] = 0x00;
	  else
	    pDst8[idx] = 0xff;
	}
      }
      break;
    case DImage::DImage_dbl_multi:
      {
	double *pdbl;

	if(imgSrc.numChannels() > 1){
	  fprintf(stderr, "DThresholder::threshImage_() does not support "
		  "multiple-channel images\n");
	  exit(1);
	}
	pdbl = imgSrc.dataPointer_dbl();
	len = w * h;
	for(int idx = 0; idx < len; ++idx){
	  if(pdbl[idx] <= tval)
	    pDst8[idx] = 0x00;
	  else
	    pDst8[idx] = 0xff;
	}
      }
      break;
    case DImage::DImage_cmplx:
      fprintf(stderr, "DThresholder::threshImage_() does not support "
	      "complex images\n");
      exit(1);
      break;
  } // end switch
}

///set pixel values conditional on relationship to the range [t1..t2]
/**For each pixel "p": If p<t1, p gets set to lowVal (or left alone if
 * lowVal==-1).  If t1<=p<=t2, then p gets set to midVal (or left
 * alone if midVal==-1).  If t2<p, then p gets set to highVal (or left
 * alone if highVal==-1).  Source and destination images may both be
 * the same image.  To use only one threshold, set t1==t2.
 */
void DThresholder::threshImageSpecial_(DImage &imgDst, const DImage &imgSrc,
				       double t1, double t2, double lowVal,
				       double midVal, double highVal){
  int w, h;
  int len;
  
  w = imgSrc.width();
  h = imgSrc.height();

  if( (&imgSrc) != (&imgDst))
    imgDst.create(imgSrc.width(), imgSrc.height(), imgSrc.getImageType(),
		  1, imgSrc.getAllocMethod());

  switch(imgSrc.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint8 *pSrc8;
	D_uint8 *pDst8;
	D_uint8 lv,mv,hv; // lowVal, midVal, highVal (casted to type of imgDst)
	pSrc8 = imgSrc.dataPointer_u8();
	pDst8 = imgDst.dataPointer_u8();
	len = w * h;
	lv = (D_uint8)lowVal;
	mv = (D_uint8)midVal;
	hv = (D_uint8)highVal;
	for(int idx = 0; idx < len; ++idx){
	  D_uint8 pxl;
	  pxl = pSrc8[idx];
	  if(pxl < t1){
	    if(-1. != lowVal)
	      pDst8[idx] = lv;
	    else
	      pDst8[idx] = pxl;
	  }
	  else if(pxl > t2){
	    if(-1. != highVal)
	      pDst8[idx] = hv;
	    else
	      pDst8[idx] = pxl;
	  }
	  else{
	    if(-1. != midVal)
	      pDst8[idx] = mv;
	    else
	      pDst8[idx] = pxl;
	  }
	}
      }
      break;
    case DImage::DImage_u16:
      {
	D_uint16 *pSrc16;
	D_uint16 *pDst16;
	D_uint16 lv,mv,hv; //lowVal,midVal,highVal (casted to type of imgDst)
	pSrc16 = imgSrc.dataPointer_u16();
	pDst16 = imgDst.dataPointer_u16();
	len = w * h;
	lv = (D_uint16)lowVal;
	mv = (D_uint16)midVal;
	hv = (D_uint16)highVal;
	for(int idx = 0; idx < len; ++idx){
	  D_uint16 pxl;
	  pxl = pSrc16[idx];
	  if(pxl < t1){
	    if(-1. != lowVal)
	      pDst16[idx] = lv;
	    else
	      pDst16[idx] = pxl;
	  }
	  else if(pxl > t2){
	    if(-1. != highVal)
	      pDst16[idx] = hv;
	    else
	      pDst16[idx] = pxl;
	  }
	  else{
	    if(-1. != midVal)
	      pDst16[idx] = mv;
	    else
	      pDst16[idx] = pxl;
	  }
	}
      }
      break;
    case DImage::DImage_u32:
      {
	D_uint32 *pSrc32;
	D_uint32 *pDst32;
	D_uint32 lv,mv,hv; //lowVal,midVal,highVal (casted to type of imgDst)
	pSrc32 = imgSrc.dataPointer_u32();
	pDst32 = imgDst.dataPointer_u32();
	len = w * h;
	lv = (D_uint32)lowVal;
	mv = (D_uint32)midVal;
	hv = (D_uint32)highVal;
	for(int idx = 0; idx < len; ++idx){
	  D_uint32 pxl;
	  pxl = pSrc32[idx];
	  if(pxl < t1){
	    if(-1. != lowVal)
	      pDst32[idx] = lv;
	    else
	      pDst32[idx] = pxl;
	  }
	  else if(pxl > t2){
	    if(-1. != highVal)
	      pDst32[idx] = hv;
	    else
	      pDst32[idx] = pxl;
	  }
	  else{
	    if(-1. != midVal)
	      pDst32[idx] = mv;
	    else
	      pDst32[idx] = pxl;
	  }
	}
      }
      break;
    case DImage::DImage_RGB:
    case DImage::DImage_RGB_16:
      fprintf(stderr, "DThresholder::threshImageSpecial_() does not "
	      "support RGB\n");
      exit(1);
      break;
    case DImage::DImage_flt_multi:
      {
	float *pflt;
	float *pfdst;
	float lv,mv,hv;//lowVal,midVal,highVal (casted to type of imgDst)

	if(imgSrc.numChannels() > 1){
	  fprintf(stderr, "DThresholder::threshImageSpecial_() does not "
		  "support multiple-channel images\n");
	  abort();
	}
	pflt = imgSrc.dataPointer_flt();
	pfdst = imgDst.dataPointer_flt();
	len = w * h;
	lv = (float)lowVal;
	mv = (float)midVal;
	hv = (float)highVal;
	for(int idx = 0; idx < len; ++idx){
	  float pxl;
	  pxl = pflt[idx];
	  if(pxl < t1){
	    if(-1. != lowVal)
	      pfdst[idx] = lv;
	    else
	      pfdst[idx] = pxl;
	  }
	  else if(pxl > t2){
	    if(-1. != highVal)
	      pfdst[idx] = hv;
	    else
	      pfdst[idx] = pxl;
	  }
	  else{
	    if(-1. != midVal)
	      pfdst[idx] = mv;
	    else
	      pfdst[idx] = pxl;
	  }
	}
      }
      break;
    case DImage::DImage_dbl_multi:
      {
	double *pdbl;
	float *pddst;

	if(imgSrc.numChannels() > 1){
	  fprintf(stderr, "DThresholder::threshImageSpecial_() does not "
		  "support multiple-channel images\n");
	  abort();
	}
	pdbl = imgSrc.dataPointer_dbl();
	pddst = imgDst.dataPointer_flt();
	len = w * h;
	for(int idx = 0; idx < len; ++idx){
	  float pxl;
	  pxl = pdbl[idx];
	  if(pxl < t1){
	    if(-1. != lowVal)
	      pddst[idx] = lowVal;
	    else
	      pddst[idx] = pxl;
	  }
	  else if(pxl > t2){
	    if(-1. != highVal)
	      pddst[idx] = highVal;
	    else
	      pddst[idx] = pxl;
	  }
	  else{
	    if(-1. != midVal)
	      pddst[idx] = midVal;
	    else
	      pddst[idx] = pxl;
	  }
	}
      }
      break;
    case DImage::DImage_cmplx:
      fprintf(stderr, "DThresholder::threshImageSpecial_() does not support "
	      "complex images\n");
      exit(1);
      break;
  } // end switch
}




/// performs Otsu's global thresholding method
/** This function just calls getOtsuThreshVal() to get the Otsu
 *  threshold value, and then calls threshImage_() with that Otsu
 *  value.  Otsu's method assumes a bimodal distribution (2 classes)
 *  and tries all threshold values, picking the one that minimizes
 *  within-class variance as the Otsu threshold.  Since this is a
 *  global threshold method, it may be best to perform background
 *  removal on the image before choosing and applying a threshold.
 */
void DThresholder::otsuThreshImage_(DImage &imgDst, const DImage &imgSrc){
  double dblOtsuVal;
  
  dblOtsuVal = DThresholder::getOtsuThreshVal(imgSrc);
  //  printf("otsuVal=%f\n", dblOtsuVal);
  DThresholder::threshImage_(imgDst, imgSrc, dblOtsuVal);
}


/// performs Doug Kennard's Connected Component global thresholding method
/** This function just calls getCCThreshVal() to get the Connected Component
 *  threshold value, and then calls threshImage_() with that thrshold
 *  value.  Since this is a global threshold
 *  method, it may be best to perform background removal on the image
 *  before choosing and applying a threshold.
 *  If tval is not NULL, then it will be filled with the threshold value
 *  used for the threshold.
 */
void DThresholder::ccThreshImage_(DImage &imgDst, const DImage &imgSrc,
				  DProgress *pProg, int numThreads,
				  int *tval){
  double tVal;

  tVal = DThresholder::getCCThreshVal(imgSrc, pProg, numThreads);
  
  if(tval != NULL)
    (*tval) = tVal;
  DThresholder::threshImage_(imgDst, imgSrc, tVal);
}


///find the otsu threshold value for an image
double DThresholder::getOtsuThreshVal(const DImage &img){
  int w, h;
  double retVal = 0.;
  w = img.width();
  h = img.height();
  
  switch(img.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint8 *pu8;
	unsigned int *pHist;

	pHist = (unsigned int *)calloc(256, sizeof(unsigned int));
	if(!pHist){
	  fprintf(stderr, "DThresholder::getOtsuThreshVal() no memory!\n");
	  exit(1);
	}
	pu8 = img.dataPointer_u8();
	for(int idx = 0, len=w*h; idx < len; ++idx){
	  ++(pHist[pu8[idx]]);
	}
	retVal = getOtsuThreshVal(pHist, 256);

	free(pHist);
// 	return retVal;
      }
      break;
    case DImage::DImage_u16:
      {
	D_uint16 *pu16;
	unsigned int *pHist;

	pHist = (unsigned int *)calloc(65536, sizeof(unsigned int));
	if(!pHist){
	  fprintf(stderr, "DThresholder::getOtsuThreshVal() no memory!\n");
	  exit(1);
	}
	pu16 = img.dataPointer_u16();
	for(int idx = 0, len=w*h; idx < len; ++idx){
	  ++(pHist[pu16[idx]]);
	}
	retVal = getOtsuThreshVal(pHist, 65536);
	free(pHist);
// 	return retVal;
      }
      break;
    case DImage::DImage_u32:
    case DImage::DImage_RGB:
    case DImage::DImage_RGB_16:
    case DImage::DImage_flt_multi:
    case DImage::DImage_dbl_multi:
    case DImage::DImage_cmplx:
      fprintf(stderr, "DThresholder::getOtsuThreshVal() currently only "
	      "supports grayscale images of type DImage_u8 or DImage_u16\n");
      abort();
      break;
  }
  return retVal;
}


///find the otsu threshold value for a histogram
/**This version is optimized better than the other version because it
   just uses means, not variances, and the means are updated iteratively
   instead of re-calculated for each possible threshold.  Pixel counts are used
   instead of probabilities since the number of pixels is constant and we are
   just maximizing the function, and only the between-class variance is
   calculated (not the total variance) since the total variance is constant
   also.  If there are not at least two levels in the histogram, then 
   the result doesn't matter, but histLen-1 is returned just in case the caller
   wants to know (histLen-1 will never be the value otherwise since thresholds
   include all values <= T).
 **/
double DThresholder::getOtsuThreshVal(unsigned int *pHist, int histLen){
  int T;
  // float V1, V2, M1, M2, Q1, Q2; // variances, means, and probs of groups 1<=t<2
  // float Vw; // within-class variance
  // float VwBest;
  int Tbest, Tmin, Tmax;

  double n1, u1, n2, u2;
  double n1_m1, n2_m1; // the n1 and n2 values for previous T
  double Vb; //between-class variance (sort of - it's scaled by a constant)
  double VbBest;

  n1= n2 = u1 = u2 = Vb = VbBest = 0.;


  Tmin = 0;
  while((pHist[Tmin] == 0) && (Tmin < (histLen-1)))
    Tmin++;

  Tmax = (histLen-1);
  while((pHist[Tmax] == 0) && (Tmax > 0))
    Tmax--;

  if(Tmin >= Tmax)//less than 2 histogram bins are used, result doesn't matter
    return histLen-1;

  //initialize n1, n2, u1, u2
  n1 = pHist[Tmin];
  u1 = pHist[Tmin]*Tmin / n1;

  for(T = Tmin+1; T <= Tmax; ++T){
    n2 += pHist[T];
    u2 += pHist[T]*T;
  }
  u2 /= n2;
  
  Vb = n1*n2*(u2-u1)*(u2-u1);
  VbBest = Vb;
  Tbest = Tmin;

  for(T = Tmin+1; T < Tmax; ++T){
    n1_m1 = n1;
    n2_m1 = n2;
    n1 += pHist[T];
    n2 -= pHist[T];
    u1 = (u1*n1_m1 + pHist[T]*T) / n1;//n1 will never be 0
    // if(n2 != 0.)//should only be 0. when T==Tmax, but we stop at Tmax-1
    u2 = (u2*n2_m1 - pHist[T]*T) / n2;
    
    Vb = n1*n2*(u1-u2)*(u1-u2);
    //if equal, choose the higher one because CC thresh uses otsu as starting
    //point and will have slightly fewer thresholds to perform CC labeling on
    if(Vb >= VbBest){
      VbBest = Vb;
      Tbest = T;
    }
  }    
  return (double)Tbest;
}



///find the otsu threshold value for a histogram
double DThresholder::getOtsuThreshVal2(unsigned int *pHist, int histLen){
  int T;
  float V1, V2, M1, M2, Q1, Q2; // variances, means, and probs of groups 1<=t<2
  float Vw; // within-class variance
  float VwBest;
  int Tbest, Tmin, Tmax;
  int i;

  Tmin = 0;
  while((pHist[Tmin] == 0) && (Tmin < (histLen-1)))
    Tmin++;

  Tmax = (histLen-1);
  while((pHist[Tmax] == 0) && (Tmax > 0))
    Tmax--;

  Tbest = Tmin;
  VwBest = 0.; /* initialize to suppress warning (this gets set when T==Tmin)*/
  //TODO: the following can probably be sped up (instead of resetting and
  // recalculating entire values each time with the inner loops, perhaps
  // we can reuse some information and just modify the values as the
  // outer loop progresses.  Although this is fast anyway for 8-bit images
  for(T = Tmin; T <= Tmax; ++T){
    Q1 = Q2 = M1 = M2 = V1 = V2 = 0.;
    for(i = 0; i <= T; i++)
      Q1 += pHist[i];
    for(i = T+1; i < histLen; i++)
      Q2 += pHist[i];
    for(i = 0; i <= T; i++)
      if(Q1 > 0.)
	M1 += i * pHist[i] / Q1;
    for(i = T+1; i < histLen; i++)
      if(Q2 > 0.)
	M2 += i * pHist[i] / Q2;

    for(i = 0; i <= T; i++)
      if(Q1 > 0.)
	V1 += (i - M1) * (i - M1) * pHist[i] / Q1;
    for(i = T+1; i < histLen; i++)
      if(Q2 > 0.)
	V2 += (i - M2) * (i - M2) * pHist[i] / Q2;
    
    Vw = Q1 * V1 + Q2 * V2;
    if((Vw <= VwBest) || (T==Tmin)){
      VwBest = Vw;
      Tbest = T;
    }
  }
  return (double)Tbest;
}


///find the otsu threshold value for a histogram - not working version.
double DThresholder::getOtsuThreshVal3(unsigned int *pHist, int histLen){
  int T;
  double V1, V2, M1, M2, Q1, Q2; // variances, means, and probs of groups 1<=t<2
  double Vw; // within-class variance
  double VwBest;
  int Tbest, Tmin, Tmax;
  int i;
  double S1,S2,SS1,SS2,QS1,QS2;

  Tmin = 0;
  while((pHist[Tmin] == 0) && (Tmin < (histLen-1)))
    Tmin++;

  Tmax = (histLen-1);
  while((pHist[Tmax] == 0) && (Tmax > 0))
    Tmax--;

  //TODO: the following can probably be sped up (instead of resetting and
  // recalculating entire values each time with the inner loops, perhaps
  // we can reuse some information and just modify the values as the
  // outer loop progresses.  Although this is fast anyway for 8-bit images

  //Q1,Q2: the sum of the values in set 1 and set 2, respectively
  //

  printf("pHist[0]=%u pHist[255]=%u\n", pHist[0], pHist[255]);

  Q1 = Q2 = M1 = M2 = V1 = V2 = 0.;
  S1=S2=SS1=SS2=QS1=QS2=0.;
  T = Tmin;
  for(i = Tmin; i <= T; i++){
    Q1 += pHist[i];
    QS1 += pHist[i]*pHist[i];
  }
  for(i = T+1; i <= Tmax; i++){
    Q2 += pHist[i];
    QS2 += pHist[i]*pHist[i];
  }

  for(i = Tmin; i <= T; i++){
    S1 += i * pHist[i];
    SS1 += (i*pHist[i])*(i*pHist[i]);
  }
  for(i = T+1; i <= Tmax; i++){
    S2 += i * pHist[i];
    SS2 += (i*pHist[i])*(i*pHist[i]);
  }


  //I think I can change this to M1 = S1/Q1
  for(i = Tmin; i <= T; i++)
    if(Q1 > 0.)
      M1 += i * pHist[i] / Q1;
  for(i = T+1; i <= Tmax; i++)
    if(Q2 > 0.)
      M2 += i * pHist[i] / Q2;

  if(Q1 > 0.){
    if(M1 != S1/Q1){
      fprintf(stderr,"getOtsuThreshVal() M1 != S1/Q1\n");
      exit(1);
    }
  }
  if(Q2 > 0.){
    if(M2 != S2/Q2){
      fprintf(stderr,"getOtsuThreshVal() M2 != S2/Q2\n");
      exit(1);
    }
  }

  //I think I can change this to V1 = (SS1/QS1) - M1*M1
  for(i = Tmin; i <= T; i++)
    if(Q1 > 0.)
      V1 += (i - M1) * (i - M1) * pHist[i] / Q1;
  for(i = T+1; i <= Tmax; i++)
    if(Q2 > 0.)
      V2 += (i - M2) * (i - M2) * pHist[i] / Q2;
  // printf("Tmin=%d Tmax=%d T=%d\n", Tmin,Tmax,T);
  // printf("Q1=%lf M1=%lf V1=%lf S1=%lf SS1=%lf QS1=%lf\n",
  // 	 Q1,M1,V1,S1,SS1,QS1);
  // printf("Q2=%lf M2=%lf V2=%lf S2=%lf SS2=%lf QS2=%lf\n",
  // 	 Q2,M2,V2,S2,SS2,QS2);

  // if(QS1 > 0.){
  //   if(V1 !=  ( (SS1/QS1) - M1 * M1)){
  //     fprintf(stderr,"getOtsuThreshVal() V1 != (SS1/QS1) - M1 * M1\n");
  //     exit(1);
  //   }
  // }
  // if(QS2 > 0.){
  //   if(V2 !=  ( (SS2/QS2) - M2 * M2)){
  //     fprintf(stderr,"getOtsuThreshVal() V2 != (SS2/QS2) - M2 * M2\n");
  //     fprintf(stderr,"(%f != %f)\n",V2,(SS2/QS2) - M2 * M2);
  //     fprintf(stderr,"SS2/Q2 - M2*M2=%lf\n",SS2/Q2 - M2*M2);
  //     exit(1);
  //   }
  // }


  Vw = Q1 * V1 + Q2 * V2;
  VwBest = Vw;
  Tbest = T;

  for(T = Tmin+1; T <= Tmax; ++T){
    //update Q1,Q2,  M1,M2,  V1,V2
    Q1 += pHist[T];
    Q2 -= pHist[T];
    
    S1 += T*pHist[T];
    S2 -= T*pHist[T];

    SS1 += (T*pHist[T])*(T*pHist[T]);
    SS2 -= (T*pHist[T])*(T*pHist[T]);

    if(Q1 > 0.)
      M1 = S1/Q1;
    if(Q2 > 0.)
      M2 = S2/Q2;

    V1 = V2 = 0.;
    for(i = 0; i <= T; i++)
      if(Q1 > 0.)
	V1 += (i - M1) * (i - M1) * pHist[i] / Q1;
    for(i = T+1; i < histLen; i++)
      if(Q2 > 0.)
	V2 += (i - M2) * (i - M2) * pHist[i] / Q2;


    // if(QS1 > 0.)
    //   V1 = SS1/QS1 - M1*M1;
    // if(QS2 > 0.)
    //   V2 = SS2/QS2 - M2*M2;

    Vw = Q1 * V1 + Q2 * V2;
    if((Vw <= VwBest) || (T==Tmin)){
      VwBest = Vw;
      Tbest = T;
    }

    {
      double tV1, tV2, tM1, tM2, tQ1, tQ2; // variances, means, and probs of groups 1<=t<2
      

      tQ1 = tQ2 = tM1 = tM2 = tV1 = tV2 = 0.;
      for(i = 0; i <= T; i++)
	tQ1 += pHist[i];
      for(i = T+1; i < histLen; i++)
	tQ2 += pHist[i];
      for(i = 0; i <= T; i++)
	if(tQ1 > 0.)
	  tM1 += i * pHist[i] / tQ1;
      for(i = T+1; i < histLen; i++)
	if(tQ2 > 0.)
	  tM2 += i * pHist[i] / tQ2;
      
      for(i = 0; i <= T; i++)
	if(tQ1 > 0.)
	  tV1 += (i - tM1) * (i - tM1) * pHist[i] / tQ1;
      for(i = T+1; i < histLen; i++)
	if(tQ2 > 0.)
	  tV2 += (i - tM2) * (i - tM2) * pHist[i] / tQ2;
      
      if((tQ1 != Q1)){
	fprintf(stderr,"getOtsuThreshVal() tQ1 != Q1\n");
	exit(1);
      }
      if((tQ2 != Q2)){
	fprintf(stderr,"getOtsuThreshVal() tQ2 != Q2\n");
	exit(1);
      }
      if((tM1 != M1)){
	fprintf(stderr,"getOtsuThreshVal() tM1 != M1\n");
	exit(1);
      }
      if((tM2 != M2)){
	fprintf(stderr,"getOtsuThreshVal() tM2 != M2\n");
	exit(1);
      }
      if((tV1 != V1)){
	fprintf(stderr,"getOtsuThreshVal() tV1 != V1\n");
	exit(1);
      }
      if((tV2 != V2)){
	fprintf(stderr,"getOtsuThreshVal() tV2 != V2\n");
	exit(1);
      }
    }

  }
  return (double)Tbest;
}




#if 0
// ///find the global threshold that results in the fewest Connected Components
// /**This method chooses the threshold value that results in the least
//  * number of connected components (CCs).  This is based on the
//  * assumption that fewer CCs will be found when words/characters are
//  * well-connected instead of disjoint (disjoint characters resulting
//  * from a threshold that is too low), and when the threshold is not so
//  * high that extra noise is introduced.  In the limit, however,
//  * everything gets connected into one big component, so heuristics are
//  * used to short-circuit turning the whole page (or a large portion of
//  * it) into one big blob.  Possible threshold values returned range
//  * from the otsu value (the algorithm's starting point) to 254.*/
// double DThresholder::getCCThreshVal(const DImage &img, DProgress *pProg,
// 				    int numThreads){
//   double tOtsu;
//   int bestT;
//   int bestNumCCs;
//   int numCCs;
//   DImage imgCCmap;
//   DImage imgThreshed;
//   D_uint32 *pMapData;
//   int otsuCCArea;
//   int *rgCCsizes;
//   int w, h; //width, height
//   bool fSingleCCExceedsOtsuArea;

//   //TODO: speed this up.  There are a few things we can do:
//   // 1-instead of using the standard DConnectedComponentLabeler::getCCimage_()
//   //   function, we can basically copy that function and add some code
//   //   that records the size of each component and as soon as one exceeds the
//   //   threshold from the otsu value, an error flag is returned to short-
//   //   circuit the processing.  Then, we wouldn't have to traverse the image
//   //   again after every CC labeling just to find out if we have gone beyond
//   //   the mark, because the flag would tell us that.
//   // 2-once #1 is in place, it would be easy to do a coarse threshold first,
//   //   then go back and fill in the gaps once we find a limit on the threshold
//   // 3-once #1 is in place, it would be easy to use multiple threads.  Just
//   //   farm the next Tvalue that hasn't been processed off to each thread that
//   //   is done with its previous Tvalue until finished.

  
//   tOtsu = getOtsuThreshVal(img);
//   threshImage_(imgThreshed, img, tOtsu);
//   DConnectedComponentLabeler::getCCimage_(imgCCmap, imgThreshed, &numCCs);

//   bestT = (int)tOtsu;
//   bestNumCCs = numCCs;

//   // record how many pixels are ink at Otsu level
//   pMapData = imgCCmap.dataPointer_u32();
//   otsuCCArea = 0;
//   w = imgCCmap.width();
//   h = imgCCmap.height();
//   for(int idx = 0, len=w*h; idx < len; ++idx){
//     otsuCCArea += (int)(pMapData[idx] != 0);
//   }

//   // iterate through all possible thresholds above the otsu threshold
//   for(int curT = 1+(int)tOtsu; curT < 255; ++curT){
//     threshImage_(imgThreshed, img, curT);
//     DConnectedComponentLabeler::getCCimage_(imgCCmap, imgThreshed, &numCCs);
//     if(numCCs < bestNumCCs){
//       // less CCs is better, but we need to make sure we haven't gone
//       // beyond the mark and started making everything into a big blob
//       // by merging all the valid CCs.  Our heuristic for checking
//       // this is to make sure no single CC has more pixels than the
//       // total amount of ink at the Otsu level.
//       rgCCsizes = (int*)calloc(numCCs+1, sizeof(int));
//       if(!rgCCsizes){
// 	fprintf(stderr, "DThresholder::getCCThreshVal() no memory!\n");
// 	exit(1);
//       }
//       w = imgCCmap.width();
//       h = imgCCmap.height();
//       pMapData = imgCCmap.dataPointer_u32();
//       for(int idx = 0, len=w*h; idx < len; ++idx){
// 	++(rgCCsizes[pMapData[idx]]);//we left slot 0 for the background count
//       }
//       fSingleCCExceedsOtsuArea = false;
//       for(int i = 1; i <= numCCs; ++i){
// 	if(rgCCsizes[i] >= otsuCCArea){
// 	  fSingleCCExceedsOtsuArea = true;
// 	  break; // break out of the inner loop: for(i=1; i <= numCCs; ...
// 	}
//       }
//       free(rgCCsizes);
//       if(fSingleCCExceedsOtsuArea)
// 	break; // break out of the outer loop: for(int curT = 1+tOtsu...
//       else{
// 	// we checked our heuristic and this thresh looks ok, so update best
// 	bestT = curT;
// 	bestNumCCs = numCCs;
//       }
//     }// end if(numCCs < bestNumCCs)
//   }// end for(int curT=1+tOtsu...
//   return (double)bestT;
// }
#endif

///performs thresholding via Milewski and Govindaraju's sine wave method
/** Most algorithm parameters are based on an estimate of pen stroke
 *  thickness.  In their research, they used 5 pixels as stroke width,
 *  which resulted in phi=5, N=5, P=3.  Technically, parameter P could
 *  fall within a range of values that is limited by the stroke
 *  thickness, but if P is unspecified (-1), then we will use the
 *  upper bound of the range.  The sineConst parameter was empirically
 *  chosen by the authors as 0.5.  It determines the frequency of the
 *  sine wave. (The equation they use is
 *  y=2*strokeThickness*sin(0.5*x).)  It makes sense that we would
 *  want to adjust that if the writing is much bigger, so I (not the
 *  paper authors) suggest maybe making it a function of the stroke
 *  thickness, (i.e. sineConst = 5./(2.*strokeThickness)).
 *
 *  This thresholding algorithm was developed for use with Carbon copy
 *  documents from medical paperwork, but it seems that noisy
 *  microfilm has many of the same characteristics, and may benefit
 *  from its use.  The paper I read to implement this is currently in
 *  press in Pattern Recognition, and pre-press manuscript is
 *  available online at ScienceDirect.com.  I have added an additional
 *  stopping criteria for implementation purposes: no sine wave will
 *  continue past three cycles, even if the 2nd condition hasn't been
 *  met yet.  Also, in this implementation I do not necessarily
 *  calculate all 8 sine wave minimums, I only trace a given wave
 *  until a mean smaller than the central mask mean is found (which
 *  guarantees that the minimum for that wave will also be smaller),
 *  and once 4 waves are guaranteed to be smaller, I don't bother
 *  computing the rest since the result will still be a foreground
 *  pixel no matter what the minimum of the rest of the waves are.
 *
 *  Citation: Robert Milewski and Venu Govindaraju, Binarization and
 *  cleanup of handwritten text from carbon copy medical form images,
 *  Pattern Recognition (2007), doi:10.1016/j.patcog.2007.08.018. */
void DThresholder::milewskiThreshImage_precollaboration(DImage &imgDst, const DImage &imgSrc,
					int strokeThickness, int P, int kappa,
					double sineConst,
					DProgress *pProg, int numThreads){
  int N; // central mask is NxN pixels where N is odd, >=3
  // outer masks are PxP pixels where P is odd, P>=3, and P<= (N+1)/2
  int Nradius;
  int Pradius;
  DKernel2D kern;
  DImage imgMeansN;
  DImage imgMeansP;
  int *rgXsSineE;
  int *rgYsSineE;
  int *rgXsSineNE;
  int *rgYsSineNE;
  double waveLen;
  int rgLen;
  double dx;
  int lastixE;
  int lastiyE;
  int lastixNE;
  int lastiyNE;
  int cycleLenE;
  int cycleLenNE;
  int rgLenE;
  int rgLenNE;
  double sinorcos45deg;
  D_uint8 *pDst8;
  D_uint8 *pMeansN;
  D_uint8 *pMeansP;
  int w, h;
  bool fContinue; // true until stopping criteria are met
  bool fCanStillVote; // this directional sine can still vote for foreground
  bool fVoteBlack; // this directional sine votes for foreground
  int meanV_plus_kappa;
  //TODO: report progress to prog
  fprintf(stderr, "milewskiThreshImage_precollaboration_() does not work.\n");
  abort();

  if(imgSrc.getImageType() != DImage::DImage_u8){
    fprintf(stderr, "DThresholder::milewskiThreshImage_() only supports 8-bit "
	    "grayscale images!\n");
    abort();
  }

  if(numThreads > 1){
    fprintf(stderr, "DThresholder::milewskiThreshImage_() doesn't support "
	    "multiple threads yet.  Using a single thread...\n");
  }
  // calculate parameters based on strokeThickness
  N = strokeThickness;
  if(0 == (N & 0x00000001)) // make sure N is an odd number
    ++N;
  if(N < 3)
    N = 3;
  if(-1 == P)
    P = (N+1)/2;
  if(0 == (P & 0x00000001)) // make sure P is an odd number
    --P;
  if(P < 3)
    P = 3;
  Nradius = N/2;
  Pradius = P/2;
  printf("N=%d P=%d Nradius=%d Pradius=%d\n", N, P, Nradius, Pradius);

  // pre-calculate all of the NxN and PxP means since they will all be used
  // more than once. (just convolve rectangular kernels)
  printf("calculating means for NxN (%dx%d) kernel...\n", N,N);
  kern.setRect(N,N);
  DConvolver::convolve_(imgMeansN, imgSrc, kern, false, false, true);
  pMeansN = imgMeansN.dataPointer_u8();
  printf("calculating means for PxP (%dx%d) kernel...\n", P,P);
  kern.setRect(P,P);
  DConvolver::convolve_(imgMeansP, imgSrc, kern, false, false, true);
  pMeansP = imgMeansP.dataPointer_u8();

  // calculate offsets from point X,Y for a sine-wave cycle (and # of offsets)
  // for each directional sine wive
  waveLen = strokeThickness * 4 * M_PI / 5.;
  dx = 1./waveLen;
  rgLen = (int)(waveLen * 3 / dx + 1);

  rgXsSineE = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgXsSineE);
  rgYsSineE = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgYsSineE);
  rgXsSineNE = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgXsSineNE);
  rgYsSineNE = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgYsSineNE);

  sinorcos45deg = sin(M_PI / 4.);
  
  lastixE = lastiyE = lastixNE = lastiyNE = -9999;
  cycleLenE = cycleLenNE = 0;
  rgLenE = rgLenNE = 0;
  for(int idx=0; idx < rgLen; ++idx){
    double xcur, ycur;
    double xcurNE, ycurNE;
    int ix, iy;
    
    xcur = dx * idx;
    ycur = 2 * strokeThickness * sin(sineConst * xcur);
    ix = (int)(xcur + 0.5);
    iy = (int)(ycur + 0.5);
    if((ix != lastixE) || (iy != lastiyE)){
      rgXsSineE[rgLenE] = ix;
      rgYsSineE[rgLenE] = iy;
      ++rgLenE;
      if(xcur < waveLen)
	cycleLenE = rgLenE;
      lastixE = ix;
      lastiyE = iy;
    }
    xcurNE = xcur * sinorcos45deg + ycur * sinorcos45deg;
    ycurNE = -1. * xcur * sinorcos45deg + ycur * sinorcos45deg;
    ix = (int)(xcurNE + 0.5);
    iy = (int)(ycurNE + 0.5);
    if((ix != lastixNE) || (iy != lastiyNE)){
      rgXsSineNE[rgLenNE] = ix;
      rgYsSineNE[rgLenNE] = iy;
      ++rgLenNE;
      if(xcur < waveLen)
	cycleLenNE = rgLenNE;
      lastixNE = ix;
      lastiyNE = iy;
    }
  }

#define MDEBUG 0

  // iterate through all pixels doing the thresholding algorithm
  w = imgSrc.width();
  h = imgSrc.height();
  imgDst.create(w, h, DImage::DImage_u8, 1, imgSrc.getAllocMethod());
  pDst8 = imgDst.dataPointer_u8();
#if MDEBUG
  imgDst = imgSrc;
  for(int y = 100; y <= 100; ++y){
    for(int x = 100; x <= 100; ++x){
	  pDst8[y*w+x] = 0xff;
#else
  for(int y = 0; y < h; ++y){
    for(int x = 0; x < w; ++x){
#endif
      int xTmp, yTmp;
      int idxTmp;
      int fgVotes; // number of votes for x,y being foreground (black)
      D_uint8 prevMwi; // previous Mwi (mean of P kernel at previous position)
      int minMwi;
      int Mwi;
      int Mv;

      fgVotes = 0;
      Mv = (int)(unsigned int)pMeansN[y*w+x];
//       if(0 == ((y+1) % 100))
// 	printf("."); fflush(stdout);
      // get vote for N direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenE)){
 	xTmp = x + rgYsSineE[idxTmp];
	yTmp = y - Nradius - Pradius - 1 - rgXsSineE[idxTmp];
// 	printf("idxTmp=%3d   xTmp=%d  yTmp=%d\n", idxTmp, xTmp, yTmp);
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) ){ // in bounds
#if MDEBUG
 	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenE) && (Mwi > prevMwi))
	    fContinue = false;
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // N direction votes for foreground


      // get vote for NE direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenNE)){
	xTmp = x + Nradius + Pradius + 1 + rgXsSineNE[idxTmp];
	yTmp = y - Nradius - Pradius - 1 + rgYsSineNE[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp < h)){//in bounds
#if MDEBUG
 	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenNE) && (Mwi > prevMwi))
	    fContinue = false;
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // NE direction votes for foreground


      // get vote for NW direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenNE)){
	xTmp = x - Nradius - Pradius - 1 - rgXsSineNE[idxTmp];
	yTmp = y - Nradius - Pradius - 1 + rgYsSineNE[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp < h)){//in bounds
#if MDEBUG
 	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenNE) && (Mwi > prevMwi))
	    fContinue = false;
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // NW direction votes for foreground



      // get vote for W direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenE)){
	xTmp = x - Nradius - Pradius - 1 - rgXsSineE[idxTmp];
	yTmp = y - rgYsSineE[idxTmp];
// 	printf("idxTmp=%3d   xTmp=%d  yTmp=%d\n", idxTmp, xTmp, yTmp);
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp<h)){ // in bounds
#if MDEBUG
//   	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenE) && (Mwi > prevMwi))
	    fContinue = false;
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // W direction votes for foreground



      // get vote for E direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenE)){
	xTmp = x + Nradius + Pradius + 1 + rgXsSineE[idxTmp];
	yTmp = y - rgYsSineE[idxTmp];
// 	printf("idxTmp=%3d   xTmp=%d  yTmp=%d\n", idxTmp, xTmp, yTmp);
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp<h)){ // in bounds
#if MDEBUG
//   	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenE) && (Mwi > prevMwi))
	    fContinue = false;
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // E direction votes for foreground


      // get vote for SW direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenNE)){
	xTmp = x - Nradius - Pradius - 1 - rgXsSineNE[idxTmp];
	yTmp = y + Nradius + Pradius + 1 - rgYsSineNE[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp < h)){//in bounds
#if MDEBUG
//  	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenNE) && (Mwi > prevMwi))
	    fContinue = false;
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // SW direction votes for foreground



      // get vote for SE direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenNE)){
	xTmp = x + Nradius + Pradius + 1 + rgXsSineNE[idxTmp];
	yTmp = y + Nradius + Pradius + 1 - rgYsSineNE[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp < h)){//in bounds
#if MDEBUG
//  	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenNE) && (Mwi > prevMwi))
	    fContinue = false;
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // SE direction votes for foreground


      // get vote for S direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenE)){
	xTmp = x + rgYsSineE[idxTmp];
	yTmp = y + Nradius + Pradius + 1 + rgXsSineE[idxTmp];
// 	printf("idxTmp=%3d   xTmp=%d  yTmp=%d\n", idxTmp, xTmp, yTmp);
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) ){ // in bounds
#if MDEBUG
//   	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenE) && (Mwi > prevMwi))
	    fContinue = false;
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // S direction votes for foreground

      if(fgVotes >= 3)
	pDst8[y*w+x] = 0x00;
      else
	pDst8[y*w+x] = 0xff;

    }// end for x
  }// end for y

  imgDst.save("/tmp/sines.pgm");

  
}


///performs Sauvola adaptive thresholding
/** The window size "window" should vary linearly from 10 to 20 pixels
    as image dpi varies from 75 to 300.  R is the dynamic range of
    standard deviation, and K gets positive values (0.5 is what
    Sauvola used). To speed up the algorithm, the paper suggests that
    a threshold for every nth pixel can be computed instead of every
    pixel and interpolation is used to pick threshold values for the
    pixels in between.  We assume n=1 (threshold is computed for every
    pixel). */
void DThresholder::sauvolaThreshImage_(DImage &imgDst, const DImage &imgSrc,
				       int window, double R,double K,
				       DProgress *pProg, int numThreads){
  float *rgTileTransDif;
  float *rgTileMean;
  int tilesX, tilesY;
  int w, h;
  int tx, ty;
  float min, max, fltScale;
  D_uint8 *pSrc;

  fprintf(stderr, "sauvolaThreshImage_() is not yet implemented.\n");
  abort();


  if(imgSrc.getImageType() != DImage::DImage_u8){
    fprintf(stderr, "DThresholder::sauvolaThreshImage_() only supports 8-bit "
	    "grayscale images!\n");
    abort();
  }

  if(numThreads > 1){
    fprintf(stderr, "DThresholder::sauvolaThreshImage_() doesn't support "
	    "multiple threads yet.  Using a single thread...\n");
  }

  w = imgSrc.width();
  h = imgSrc.height();
  pSrc = imgSrc.dataPointer_u8();

  tilesX = (w + window - 1) / window;
  tilesY = (h +  window - 1) / window;

  rgTileTransDif = (float*)calloc(tilesX * tilesY, sizeof(float));
  D_CHECKPTR(rgTileTransDif);
  rgTileMean = (float*)calloc(tilesX * tilesY, sizeof(float));
  D_CHECKPTR(rgTileMean);
  
  // calculate the mean and transient difference for each tile----------------

  // special case first row: when y=0, don't use invalid position at y-1
  rgTileMean[0] += pSrc[0]; // special case x=0
  for(int x = 1; x < w; ++x){
    tx = x / window;
    rgTileMean[tx] += pSrc[x];
    rgTileTransDif[tx] +=
      (float)(abs( ((int)(unsigned int)pSrc[x]) -
		   ((int)(unsigned int)pSrc[x-1])));
  }
  //now do the rest of the rows (y= 1...h-1)
  for(int y = 1; y < h; ++y){
    int idx;
    ty = y / window;
    //special case for x=0:
    idx = y*w;
    rgTileMean[ty*tilesX] += pSrc[idx];
    rgTileTransDif[ty*tilesX] +=
      (float)(abs( ((int)(unsigned int)pSrc[idx]) -
		   ((int)(unsigned int)pSrc[idx-w])));
    // now for the rest of the x's:
    for(int x = 1; x < w; ++x){
      idx = y*w + x;
      tx = x / window;
      rgTileMean[ty*tilesX + tx] += pSrc[idx];
      rgTileTransDif[ty*tilesX + tx] +=
	(float)(abs( 2 * ((int)(unsigned int)pSrc[idx]) -
		     ((int)(unsigned int)pSrc[idx-1] + 
		      (int)(unsigned int)pSrc[idx-w])));
    }
  }
  // now divide each tile's accumulated sum by the proper amount to get mean
  // and normalize the transDif values for number of pixels (to take care of
  // the smaller tiles on right/bottom edges)
  min = 999.;
  max = -999.;
  for(int ty=0, tidx=0; ty < tilesY; ++ty){
    for(int tx=0; tx < tilesX; ++tx, ++tidx){
      int winw, winh;
      winw = window;
      winh = window;
      if(ty == (tilesY-1)){
	if(0 != (h % window))
	  winh = h % window;
      }
      if(tx == (tilesX-1)){
	if(0 != (w % window))
	  winw = w % window;
      }
      rgTileMean[tidx] /= ((float)winw*(float)winh);
      rgTileTransDif[tidx] /= (512.*winw*winh); // this is different than
      // in the paper.  The paper says divide by (Ln)^2, which I don't
      // think makes sense (L is the number of gray-levels--I'm not
      // sure if it means the range in the image or the possible
      // range--I assume the possible range which is 256).  It seems
      // to me that the denominator should be 2*L*(n^2), not what the
      // paper says, but I could be wrong.  Since we are scaling to
      // 0..1 later, it really doesn't matter though, as all tiles in
      // the paper are normalized with the same value.  My change here
      // only affects the tiles on the right and bottom edge, which
      // the paper doesn't mention how to treat anyway.
      if(rgTileTransDif[tidx] < min)
	min = rgTileTransDif[tidx];
      if(rgTileTransDif[tidx] > max)
	max = rgTileTransDif[tidx];
    }
  }
  // scale all of the rgTileTransDif means to 0..1 based on highest and
  // lowest value.
  if( (max-min) > 0)
    fltScale = 1. / (max - min);
  else
    fltScale = 1.;
  float newMin, newMax;

  newMin = 999.;
  newMax = -999.;
  for(int ty=0, tidx=0; ty < tilesY; ++ty){
    for(int tx=0; tx < tilesX; ++tx, ++tidx){
      rgTileTransDif[tidx] = (rgTileTransDif[tidx] - min) * fltScale;
      if(rgTileTransDif[tidx] < newMin)
	newMin = rgTileTransDif[tidx];
      if(rgTileTransDif[tidx] > newMax)
	newMax = rgTileTransDif[tidx];
    }
  }


#if 1 // output a debug map of where text binarization vs non-text will be done
 {
   DImage imgTmp;
   D_uint8 *pTmp;

   imgTmp = imgSrc;
   pTmp = imgTmp.dataPointer_u8();
   for(int y = 0, idxTmp=0; y < h; ++y){
     for(int x = 0; x < w; ++x, ++idxTmp){
       int tileIdx;
       tileIdx = (y/window) * tilesX + (x/window);
       pTmp[idxTmp] = (rgTileTransDif[tileIdx] >= 0.3) ? 0x00 : 0xff;
     }
   }
   imgTmp.save("/tmp/map.pgm");
 }
#endif
  




//   fprintf(stderr, "DThresholder::sauvolaThreshImage_() not yet implemented\n");
//   exit(1);
}


///performs thresholding via Milewski and Govindaraju's sine wave method
/** Most algorithm parameters are based on an estimate of pen stroke
 *  thickness.  In their research, they used 5 pixels as stroke width,
 *  which resulted in phi=5, N=5, P=3.  Technically, parameter P could
 *  fall within a range of values that is limited by the stroke
 *  thickness, but if P is unspecified (-1), then we will use the
 *  upper bound of the range.  The sineConst parameter was empirically
 *  chosen by the authors as 0.5.  It determines the frequency of the
 *  sine wave. (The equation they use is
 *  y=2*strokeThickness*sin(0.5*x).)  It makes sense that we would
 *  want to adjust that if the writing is much bigger, so I (not the
 *  paper authors) suggest maybe making it a function of the stroke
 *  thickness, (i.e. sineConst = 5./(2.*strokeThickness)). The
 *  parameter voteThresh is how many (minimum) of the directions must
 *  have a minimum mask lighter than the central mask in order for the
 *  center pixel to be marked as foreground (black).
 *
 *  This thresholding algorithm was developed for use with Carbon copy
 *  documents from medical paperwork, but it seems that noisy
 *  microfilm has many of the same characteristics, and may benefit
 *  from its use.  The paper I read to implement this is currently in
 *  press in Pattern Recognition, and pre-press manuscript is
 *  available online at ScienceDirect.com.  I have added an additional
 *  stopping criteria for implementation purposes: no sine wave will
 *  continue past three cycles, even if the 2nd condition hasn't been
 *  met yet.  Also, in this implementation I do not necessarily
 *  calculate all 8 sine wave minimums, I only trace a given wave
 *  until a mean smaller than the central mask mean is found (which
 *  guarantees that the minimum for that wave will also be smaller),
 *  and once voteThresh waves are guaranteed to be smaller, I don't
 *  bother computing the rest since the result will still be a
 *  foreground pixel no matter what the minimum of the rest of the
 *  waves are.
 *
 *  Citation: Robert Milewski and Venu Govindaraju, Binarization and
 *  cleanup of handwritten text from carbon copy medical form images,
 *  Pattern Recognition (2007), doi:10.1016/j.patcog.2007.08.018.
 *
 *  Also, many thanks to Rob Milewski for his willingness to provide
 *  me with a sample result from his binarizer for me to compare the
 *  results of my implementation to, and who also answered several of
 *  my questions about the algorithm when I was trying to track down
 *  mistakes in my initial implementation. */
void DThresholder::milewskiThreshImage_(DImage &imgDst, const DImage &imgSrc,
					int strokeThickness, int P, int kappa,
					double sineConst, int voteThresh,
					DProgress *pProg, int numThreads){
  int N; // central mask is NxN pixels where N is odd, >=3
  // outer masks are PxP pixels where P is odd, P>=3, and P<= (N+1)/2
  int Nradius;
  int Pradius;
  DKernel2D kern;
  DImage imgMeansN;
  DImage imgMeansP;
  int *rgXsSineE;
  int *rgYsSineE;
  int *rgXsSineNE;
  int *rgYsSineNE;
  int *rgXsSineNW;
  int *rgYsSineNW;
  double waveLen;
  int rgLen;
  double dx;
  int lastixE;
  int lastiyE;
  int lastixNE;
  int lastiyNE;
  int lastixNW;
  int lastiyNW;
  int cycleLenE;
  int cycleLenNE;
  int cycleLenNW;
  int rgLenE;
  int rgLenNE;
  int rgLenNW;
  double sinorcos45deg;
  D_uint8 *pDst8;
  D_uint8 *pMeansN;
  D_uint8 *pMeansP;
  int w, h;
  bool fContinue; // true until stopping criteria are met
  bool fCanStillVote; // this directional sine can still vote for foreground
  bool fVoteBlack; // this directional sine votes for foreground
  int meanV_plus_kappa;

  fprintf(stderr, "milewskiThreshImage_() does not work yet.\n");
  abort();
  if(imgSrc.getImageType() != DImage::DImage_u8){
    fprintf(stderr, "DThresholder::milewskiThreshImage_() only supports 8-bit "
	    "grayscale images!\n");
    abort();
  }

  if(numThreads > 1){
    fprintf(stderr, "DThresholder::milewskiThreshImage_() doesn't support "
	    "multiple threads yet.  Using a single thread...\n");
  }
  // calculate parameters based on strokeThickness
  N = strokeThickness;
  if(0 == (N & 0x00000001)) // make sure N is an odd number
    ++N;
  if(N < 3)
    N = 3;
  if(-1 == P)
    P = (N+1)/2;
  if(0 == (P & 0x00000001)) // make sure P is an odd number
    --P;
  if(P < 3)
    P = 3;
  Nradius = N/2;
  Pradius = P/2;
  printf("N=%d P=%d Nradius=%d Pradius=%d\n", N, P, Nradius, Pradius);

  // pre-calculate all of the NxN and PxP means since they will all be used
  // more than once. (just convolve rectangular kernels)
  printf("calculating means for NxN (%dx%d) kernel...\n", N,N);
  kern.setRect(Nradius,Nradius);
  DConvolver::convolve_(imgMeansN, imgSrc, kern, false, false, true);
  pMeansN = imgMeansN.dataPointer_u8();
  printf("calculating means for PxP (%dx%d) kernel...\n", P,P);
  kern.setRect(Pradius,Pradius);
  DConvolver::convolve_(imgMeansP, imgSrc, kern, false, false, true);
  pMeansP = imgMeansP.dataPointer_u8();

  // calculate offsets from point X,Y for a sine-wave cycle (and # of offsets)
  // for each directional sine wive
  waveLen = 2.* M_PI / sineConst;
  dx = 1./waveLen;
  rgLen = (int)(waveLen * 3 / dx + 1);

  rgXsSineE = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgXsSineE);
  rgYsSineE = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgYsSineE);
  rgXsSineNE = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgXsSineNE);
  rgYsSineNE = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgYsSineNE);
  rgXsSineNW = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgXsSineNW);
  rgYsSineNW = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgYsSineNW);



  sinorcos45deg = sin(M_PI / 4.);
  
  lastixE = lastiyE = lastixNE = lastiyNE = lastixNW = lastiyNW = -9999;
  cycleLenE = cycleLenNE = cycleLenNW = 0;
  rgLenE = rgLenNE = rgLenNW = 0;

  D_uint8 *pSrc8;
  pSrc8 = imgSrc.dataPointer_u8();
  for(int idx=0; idx < rgLen; ++idx){
    double xcur, ycur;
    double xcurNE, ycurNE;
    double xcurNW, ycurNW;
    int ix, iy;
    
    xcur = dx * idx;
    ycur = 2 * strokeThickness * sin(sineConst * xcur);
    ix = (int)(xcur + 0.5);
    iy = (int)(ycur + 0.5);
    if((ix != lastixE) || (iy != lastiyE)){
      rgXsSineE[rgLenE] = ix;
      rgYsSineE[rgLenE] = iy;
      ++rgLenE;
      if(xcur < waveLen)
	cycleLenE = rgLenE;
      lastixE = ix;
      lastiyE = iy;
    }
    xcurNE = xcur * sinorcos45deg + ycur * sinorcos45deg;
    ycurNE = -1. * xcur * sinorcos45deg + ycur * sinorcos45deg;
    ix = (int)(xcurNE + 0.5);
    iy = (int)(ycurNE + 0.5);
    if((ix != lastixNE) || (iy != lastiyNE)){
      rgXsSineNE[rgLenNE] = ix;
      rgYsSineNE[rgLenNE] = iy;
      ++rgLenNE;
      if(xcur < waveLen)
	cycleLenNE = rgLenNE;
      lastixNE = ix;
      lastiyNE = iy;
    }
    xcurNW = xcur * (-1.*sinorcos45deg) + ycur * sinorcos45deg;
    ycurNW = -1. * xcur * sinorcos45deg + ycur * (-1.*sinorcos45deg);
    ix = (int)(xcurNW + 0.5);
    iy = (int)(ycurNW + 0.5);
    if((ix != lastixNW) || (iy != lastiyNW)){
      rgXsSineNW[rgLenNW] = ix;
      rgYsSineNW[rgLenNW] = iy;
      ++rgLenNW;
      if(xcur < waveLen)
	cycleLenNW = rgLenNW;
      lastixNW = ix;
      lastiyNW = iy;
    }
  }

#define MDEBUG 0

  // iterate through all pixels doing the thresholding algorithm
  w = imgSrc.width();
  h = imgSrc.height();
  imgDst.create(w, h, DImage::DImage_u8, 1, imgSrc.getAllocMethod());
  pDst8 = imgDst.dataPointer_u8();
#if MDEBUG
  imgDst = imgSrc;
  for(int y = 100; y <= 100; ++y){
    for(int x = 100; x <= 100; ++x){
	  pDst8[y*w+x] = 0xff;
#else
  for(int y = 0; y < h; ++y){
    for(int x = 0; x < w; ++x){
#endif
      int xTmp, yTmp;
      int idxTmp;
      int fgVotes; // number of votes for x,y being foreground (black)
      D_uint8 prevMwi; // previous Mwi (mean of P kernel at previous position)
      int minMwi;
      int Mwi;
      int Mv;
      int sum;
      int count;

      fgVotes = 0;
      Mv = (int)(unsigned int)pMeansN[y*w+x];
      //      Mv = (int)(unsigned int)pSrc8[y*w+x];
//       if(0 == ((y+1) % 100))
// 	printf("."); fflush(stdout);

      // get vote for N direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 0;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenE)){
 	xTmp = x - rgYsSineE[idxTmp];
	yTmp = y - Nradius - Pradius - 1 - rgXsSineE[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) ){ // in bounds
#if MDEBUG
 	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenE) && (Mwi > prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi > minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
// 	    if((minMwi - Mv) < kappa)
// 	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
//       if( ((255-minMwi) - (255-Mv)) >= kappa)
	++fgVotes; // N direction votes for foreground


      // get vote for NE direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 0;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenNE)){
	xTmp = x + Nradius + Pradius + 1 - rgXsSineNW[idxTmp];
	yTmp = y - Nradius - Pradius - 1 + rgYsSineNW[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp < h)){//in bounds
#if MDEBUG
 	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenNE) && (Mwi > prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi > minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
// 	    if((minMwi - Mv) < kappa)
// 	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
//       if( ((255-minMwi) - (255-Mv)) >= kappa)
	++fgVotes; // NE direction votes for foreground


      // get vote for NW direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 0;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenNE)){
	xTmp = x - Nradius - Pradius - 1 - rgXsSineNE[idxTmp];
	yTmp = y - Nradius - Pradius - 1 + rgYsSineNE[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp < h)){//in bounds
#if MDEBUG
 	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenNE) && (Mwi > prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi > minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
// 	    if((minMwi - Mv) < kappa)
// 	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
//       if( ((255-minMwi) - (255-Mv)) >= kappa)
	++fgVotes; // NW direction votes for foreground



      // get vote for W direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 0;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenE)){
	xTmp = x - Nradius - Pradius - 1 - rgXsSineE[idxTmp];
	yTmp = y + rgYsSineE[idxTmp];
// 	printf("idxTmp=%3d   xTmp=%d  yTmp=%d\n", idxTmp, xTmp, yTmp);
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp<h)){ // in bounds
#if MDEBUG
   	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenE) && (Mwi > prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi > minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
// 	    if((minMwi - Mv) < kappa)
// 	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
//       if( ((255-minMwi) - (255-Mv)) >= kappa)
	++fgVotes; // W direction votes for foreground



      // get vote for E direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 0;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenE)){
	xTmp = x + Nradius + Pradius + 1 + rgXsSineE[idxTmp];
	yTmp = y - rgYsSineE[idxTmp];
// 	printf("idxTmp=%3d   xTmp=%d  yTmp=%d\n", idxTmp, xTmp, yTmp);
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp<h)){ // in bounds
#if MDEBUG
   	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenE) && (Mwi > prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi > minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
// 	    if((minMwi - Mv) < kappa)
// 	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
//       if( ((255-minMwi) - (255-Mv)) >= kappa)
	++fgVotes; // E direction votes for foreground


      // get vote for SW direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 0;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenNE)){
	xTmp = x - Nradius - Pradius - 1 + rgXsSineNW[idxTmp];
	yTmp = y + Nradius + Pradius + 1 - rgYsSineNW[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp < h)){//in bounds
#if MDEBUG
  	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenNE) && (Mwi > prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi > minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
// 	    if((minMwi - Mv) < kappa)
// 	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
//       if( ((255-minMwi) - (255-Mv)) >= kappa)
	++fgVotes; // SW direction votes for foreground



      // get vote for SE direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 0;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenNE)){
	xTmp = x + Nradius + Pradius + 1 + rgXsSineNE[idxTmp];
	yTmp = y + Nradius + Pradius + 1 - rgYsSineNE[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp < h)){//in bounds
#if MDEBUG
  	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenNE) && (Mwi > prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi > minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
// 	    if((minMwi - Mv) < kappa)
// 	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
     if( (minMwi - Mv) >= kappa)
//       if( ((255-minMwi) - (255-Mv)) >= kappa)
	++fgVotes; // SE direction votes for foreground


      // get vote for S direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      fVoteBlack = false;
      minMwi = 0;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenE)){
	xTmp = x + rgYsSineE[idxTmp];
	yTmp = y + Nradius + Pradius + 1 + rgXsSineE[idxTmp];
// 	printf("idxTmp=%3d   xTmp=%d  yTmp=%d\n", idxTmp, xTmp, yTmp);
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) ){ // in bounds
#if MDEBUG
   	  pDst8[yTmp*w+xTmp] = 0xff;
#endif
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenE) && (Mwi > prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi > minMwi){
	    minMwi = Mwi;
#if MDEBUG
#else
// 	    if((minMwi - Mv) < kappa)
// 	      fCanStillVote = false;
#endif
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
//       if( ((255-minMwi) - (255-Mv)) >= kappa)
	++fgVotes; // S direction votes for foreground

      if(fgVotes >= voteThresh)
	pDst8[y*w+x] = 0x00;
      else
	pDst8[y*w+x] = 0xff;

    }// end for x
  }// end for y

  imgDst.save("/tmp/sines.pgm");

  free(rgXsSineE);
  free(rgYsSineE);
  free(rgXsSineNE);
  free(rgYsSineNE);
  free(rgXsSineNW);
  free(rgYsSineNW);
}




///performs thresholding via Milewski and Govindaraju's sine wave method
/** Most algorithm parameters are based on an estimate of pen stroke
 *  thickness.  In their research, they used 5 pixels as stroke width,
 *  which resulted in phi=5, N=5, P=3.  Technically, parameter P could
 *  fall within a range of values that is limited by the stroke
 *  thickness, but if P is unspecified (-1), then we will use the
 *  upper bound of the range.  The sineConst parameter was empirically
 *  chosen by the authors as 0.5.  It determines the frequency of the
 *  sine wave. (The equation they use is
 *  y=2*strokeThickness*sin(0.5*x).)  It makes sense that we would
 *  want to adjust that if the writing is much bigger, so I (not the
 *  paper authors) suggest maybe making it a function of the stroke
 *  thickness, (i.e. sineConst = 5./(2.*strokeThickness)). The
 *  parameter voteThresh is how many (minimum) of the directions must
 *  have a minimum mask lighter than the central mask in order for the
 *  center pixel to be marked as foreground (black).
 *
 *  This thresholding algorithm was developed for use with Carbon copy
 *  documents from medical paperwork, but it seems that noisy
 *  microfilm has many of the same characteristics, and may benefit
 *  from its use.  The paper I read to implement this is currently in
 *  press in Pattern Recognition, and pre-press manuscript is
 *  available online at ScienceDirect.com.  I have added an additional
 *  stopping criteria for implementation purposes: no sine wave will
 *  continue past three cycles, even if the 2nd condition hasn't been
 *  met yet.  Also, in this implementation I do not necessarily
 *  calculate all 8 sine wave minimums, I only trace a given wave
 *  until a mean smaller than the central mask mean is found (which
 *  guarantees that the minimum for that wave will also be smaller),
 *  and once voteThresh waves are guaranteed to be smaller, I don't
 *  bother computing the rest since the result will still be a
 *  foreground pixel no matter what the minimum of the rest of the
 *  waves are.
 *
 *  Citation: Robert Milewski and Venu Govindaraju, Binarization and
 *  cleanup of handwritten text from carbon copy medical form images,
 *  Pattern Recognition (2007), doi:10.1016/j.patcog.2007.08.018.
 *
 *  Also, many thanks to Rob Milewski for his willingness to provide
 *  me with a sample result from his binarizer for me to compare the
 *  results of my implementation to, and who also answered several of
 *  my questions about the algorithm when I was trying to track down
 *  mistakes in my initial implementation. */
void DThresholder::milewskiThreshImage3_(DImage &imgDst, const DImage &imgSrc,
					int strokeThickness, int P, int kappa,
					double sineConst, int voteThresh,
					DProgress *pProg, int numThreads){
  int N; // central mask is NxN pixels where N is odd, >=3
  // outer masks are PxP pixels where P is odd, P>=3, and P<= (N+1)/2
  int Nradius;
  int Pradius;
  DKernel2D kern;
  DImage imgMeansN;
  DImage imgMeansP;
  int *rgXsSineE;
  int *rgYsSineE;
  int *rgXsSineNE;
  int *rgYsSineNE;
  int *rgXsSineNW;
  int *rgYsSineNW;
  double waveLen;
  int rgLen;
  double dx;
  int lastixE;
  int lastiyE;
  int lastixNE;
  int lastiyNE;
  int lastixNW;
  int lastiyNW;
  int cycleLenE;
  int cycleLenNE;
  int cycleLenNW;
  int rgLenE;
  int rgLenNE;
  int rgLenNW;
  double sinorcos45deg;
  D_uint8 *pDst8;
  D_uint8 *pMeansN;
  D_uint8 *pMeansP;
  int w, h;
  bool fContinue; // true until stopping criteria are met
  bool fCanStillVote; // this directional sine can still vote for foreground
  bool fVoteBlack; // this directional sine votes for foreground
  int meanV_plus_kappa;

  fprintf(stderr, "milewskiThreshImage3_() does not work yet.\n");
  abort();
  if(imgSrc.getImageType() != DImage::DImage_u8){
    fprintf(stderr, "DThresholder::milewskiThreshImage_() only supports 8-bit "
	    "grayscale images!\n");
    abort();
  }

  if(numThreads > 1){
    fprintf(stderr, "DThresholder::milewskiThreshImage_() doesn't support "
	    "multiple threads yet.  Using a single thread...\n");
  }
  // calculate parameters based on strokeThickness
  N = strokeThickness;
  if(0 == (N & 0x00000001)) // make sure N is an odd number
    ++N;
  if(N < 3)
    N = 3;
  if(-1 == P)
    P = (N+1)/2;
  if(0 == (P & 0x00000001)) // make sure P is an odd number
    --P;
  if(P < 3)
    P = 3;
  Nradius = N/2;
  Pradius = P/2;
  printf("N=%d P=%d Nradius=%d Pradius=%d\n", N, P, Nradius, Pradius);

  // pre-calculate all of the NxN and PxP means since they will all be used
  // more than once. (just convolve rectangular kernels)
  printf("calculating means for NxN (%dx%d) kernel...\n", N,N);
  kern.setRect(Nradius,Nradius);
  DConvolver::convolve_(imgMeansN, imgSrc, kern, false, false, true);
  pMeansN = imgMeansN.dataPointer_u8();
  printf("calculating means for PxP (%dx%d) kernel...\n", P,P);
  kern.setRect(Pradius,Pradius);
  DConvolver::convolve_(imgMeansP, imgSrc, kern, false, false, true);
  pMeansP = imgMeansP.dataPointer_u8();

  // invert the masks so they represent "intensity" rather than actual
  // pixel value
  w = imgSrc.width();
  h = imgSrc.height();
  for(int y=0, idx=0; y < h; ++y){
    for(int x=0; x < w; ++x, ++idx){
      pMeansP[idx] = (D_uint8)(255-(unsigned int)(pMeansP[idx]));
    }
  }
  imgMeansP.save("/tmp/meansP.pgm");
  for(int y=0, idx=0; y < h; ++y){
    for(int x=0; x < w; ++x, ++idx){
      pMeansN[idx] = (D_uint8)(255-(unsigned int)(pMeansN[idx]));
    }
  }
  imgMeansN.save("/tmp/meansN.pgm");



  // calculate offsets from point X,Y for a sine-wave cycle (and # of offsets)
  // for each directional sine wive
  waveLen = 2.* M_PI / sineConst;
  dx = 1./waveLen;
  rgLen = (int)(waveLen * 3 / dx + 1);

  rgXsSineE = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgXsSineE);
  rgYsSineE = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgYsSineE);
  rgXsSineNE = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgXsSineNE);
  rgYsSineNE = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgYsSineNE);
  rgXsSineNW = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgXsSineNW);
  rgYsSineNW = (int*)malloc(sizeof(int) * rgLen);
  D_CHECKPTR(rgYsSineNW);



  sinorcos45deg = sin(M_PI / 4.);
  
  lastixE = lastiyE = lastixNE = lastiyNE = lastixNW = lastiyNW = -9999;
  cycleLenE = cycleLenNE = cycleLenNW = 0;
  rgLenE = rgLenNE = rgLenNW = 0;

  D_uint8 *pSrc8;
  pSrc8 = imgSrc.dataPointer_u8();
  for(int idx=0; idx < rgLen; ++idx){
    double xcur, ycur;
    double xcurNE, ycurNE;
    double xcurNW, ycurNW;
    int ix, iy;
    
    xcur = dx * idx;
    ycur = 2 * strokeThickness * sin(sineConst * xcur);
    ix = (int)(xcur + 0.5);
    iy = (int)(ycur + 0.5);
    if((ix != lastixE) || (iy != lastiyE)){
      rgXsSineE[rgLenE] = ix;
      rgYsSineE[rgLenE] = iy;
      ++rgLenE;
      if(xcur < waveLen)
	cycleLenE = rgLenE;
      lastixE = ix;
      lastiyE = iy;
    }
    xcurNE = xcur * sinorcos45deg + ycur * sinorcos45deg;
    ycurNE = -1. * xcur * sinorcos45deg + ycur * sinorcos45deg;
    ix = (int)(xcurNE + 0.5);
    iy = (int)(ycurNE + 0.5);
    if((ix != lastixNE) || (iy != lastiyNE)){
      rgXsSineNE[rgLenNE] = ix;
      rgYsSineNE[rgLenNE] = iy;
      ++rgLenNE;
      if(xcur < waveLen)
	cycleLenNE = rgLenNE;
      lastixNE = ix;
      lastiyNE = iy;
    }
    xcurNW = xcur * (-1.*sinorcos45deg) + ycur * sinorcos45deg;
    ycurNW = -1. * xcur * sinorcos45deg + ycur * (-1.*sinorcos45deg);
    ix = (int)(xcurNW + 0.5);
    iy = (int)(ycurNW + 0.5);
    if((ix != lastixNW) || (iy != lastiyNW)){
      rgXsSineNW[rgLenNW] = ix;
      rgYsSineNW[rgLenNW] = iy;
      ++rgLenNW;
      if(xcur < waveLen)
	cycleLenNW = rgLenNW;
      lastixNW = ix;
      lastiyNW = iy;
    }
  }


  // iterate through all pixels doing the thresholding algorithm
  w = imgSrc.width();
  h = imgSrc.height();
  imgDst.create(w, h, DImage::DImage_u8, 1, imgSrc.getAllocMethod());
  pDst8 = imgDst.dataPointer_u8();
  for(int y = 0; y < h; ++y){
    for(int x = 0; x < w; ++x){
      int xTmp, yTmp;
      int idxTmp;
      int fgVotes; // number of votes for x,y being foreground (black)
      D_uint8 prevMwi; // previous Mwi (mean of P kernel at previous position)
      int minMwi;
      int Mwi;
      int Mv;
      int sum;
      int count;

      fgVotes = 0;
      Mv = (int)(unsigned int)pMeansN[y*w+x];

      // get vote for N direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenE)){
 	xTmp = x - rgYsSineE[idxTmp];
	yTmp = y - Nradius - Pradius - 1 - rgXsSineE[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) ){ // in bounds
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenE) && (Mwi < prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // N direction votes for foreground


      // get vote for NE direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenNE)){
	xTmp = x + Nradius + Pradius + 1 - rgXsSineNW[idxTmp];
	yTmp = y - Nradius - Pradius - 1 + rgYsSineNW[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp < h)){//in bounds
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenNE) && (Mwi < prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // NE direction votes for foreground


      // get vote for NW direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenNE)){
	xTmp = x - Nradius - Pradius - 1 - rgXsSineNE[idxTmp];
	yTmp = y - Nradius - Pradius - 1 + rgYsSineNE[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp < h)){//in bounds
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenNE) && (Mwi < prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // NW direction votes for foreground



      // get vote for W direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenE)){
	xTmp = x - Nradius - Pradius - 1 - rgXsSineE[idxTmp];
	yTmp = y + rgYsSineE[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp<h)){ // in bounds
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenE) && (Mwi < prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // W direction votes for foreground



      // get vote for E direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenE)){
	xTmp = x + Nradius + Pradius + 1 + rgXsSineE[idxTmp];
	yTmp = y - rgYsSineE[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp<h)){ // in bounds
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenE) && (Mwi < prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // E direction votes for foreground


      // get vote for SW direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenNE)){
	xTmp = x - Nradius - Pradius - 1 + rgXsSineNW[idxTmp];
	yTmp = y + Nradius + Pradius + 1 - rgYsSineNW[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp < h)){//in bounds
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenNE) && (Mwi < prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // SW direction votes for foreground



      // get vote for SE direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenNE)){
	xTmp = x + Nradius + Pradius + 1 + rgXsSineNE[idxTmp];
	yTmp = y + Nradius + Pradius + 1 - rgYsSineNE[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) && (yTmp < h)){//in bounds
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenNE) && (Mwi < prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
     if( (minMwi - Mv) >= kappa)
	++fgVotes; // SE direction votes for foreground


      // get vote for S direction
      idxTmp = 0;
      fContinue = true;
      fCanStillVote = true;
      minMwi = 255;
      prevMwi = 0;
      while(fContinue && fCanStillVote && (idxTmp < rgLenE)){
	xTmp = x + rgYsSineE[idxTmp];
	yTmp = y + Nradius + Pradius + 1 + rgXsSineE[idxTmp];
	if( (xTmp >= 0) && (xTmp < w) && (yTmp >= 0) ){ // in bounds
	  Mwi = (int)(unsigned int)pMeansP[yTmp*w+xTmp];
	  if((idxTmp >= cycleLenE) && (Mwi < prevMwi)){
	    fContinue = false;
	    continue;
	  }
	  prevMwi = Mwi;
	  if(Mwi < minMwi){
	    minMwi = Mwi;
	    if((minMwi - Mv) < kappa)
	      fCanStillVote = false;
	  } // end if(Mwi < minMwi)
	} // end if in bounds
	++idxTmp;
      } // end while
      if( (minMwi - Mv) >= kappa)
	++fgVotes; // S direction votes for foreground

      if(fgVotes >= voteThresh)
	pDst8[y*w+x] = 0x00;
      else
	pDst8[y*w+x] = 0xff;

    }// end for x
  }// end for y

  imgDst.save("/tmp/sines.pgm");

  free(rgXsSineE);
  free(rgYsSineE);
  free(rgXsSineNE);
  free(rgYsSineNE);
  free(rgXsSineNW);
  free(rgYsSineNW);
}





///performs the modified Niblack portion (text part) of Sauvola's algorithm
/** The window size "window" should vary linearly from 10 to 20 pixels
    as image dpi varies from 75 to 300.  R is the dynamic range of
    standard deviation, and K gets positive values (0.5 is what
    Sauvola used). This function does NOT implement the entire Sauvola
    method.  It assumes the entire image is text and performs
    Sauvola's modified Niblack thresholding.  However, most
    comparisons I have seen in the literature ONLY implement this
    part.  I do not recall seeing any implementations of Sauvola's
    entire algorithm used for comparisons by anyone but Sauvols himself.

    This implementation uses integral images for speed, as described
    in "Efficient Implementation of Local Adaptive Thresholding
    Techniques Using Integral Images" by Shafait, Keysers, and Breuel
    that was to appear in Proc. Document Recognition and Retrieval XV,
    IST/SPIE Annual Symposium, San Jose, CA, January 2008.*/

/**TODO: double-check this to make sure I am not off by one (in the
 * meanfilter, I had to pad the left and top by radius+1 but right and
 * bottom just by radius*/
void DThresholder::sauvolaNiblackThreshImage_(DImage &imgDst,
					      const DImage &imgSrc,
					      int window, double R,double K,
					      DProgress *pProg,
					      int numThreads){
  double mean; //mean
  double stddev; // standard deviation
  double thresh; // threshold value
  double sum;
  double sumofsquares;
  int ww, ww2; // window width, windowwidthSquared
  int w, h;
  DImage imgPadded;
  DImage imgDstPadded;
  D_uint8 *pu8Dst;
  D_uint8 *pu8Src;
  D_uint64 *pIntegralImgOfValues;
  D_uint64 *pIntegralImgOfSquaredValues;
  int radius1, radius2; // in case window is not odd, we use two radii here
  int progMax; // for reporting progress
  int progCur = 0; // for reporting progress

  if(imgSrc.getImageType() != DImage::DImage_u8){
    fprintf(stderr, "DThresholder::sauvolaNiblackThreshImage_() "
	    "only supports 8-bit grayscale images!\n");
    abort();
  }

  if(numThreads > 1){
    fprintf(stderr, "DThresholder::sauvolaNiblackThreshImage_() doesn't "
	    "support multiple threads yet.  Using a single thread...\n");
  }


  // first, pad the source image (mirror the borders) and a temporary dst image
  radius1 = (window-1)/2;
  radius2 = window/2;
  imgSrc.padEdges_(imgPadded, radius1+1, radius2, radius1+1, radius2,
		   DImage::DImagePadReplicate);

  w = imgPadded.width();
  h = imgPadded.height();
  imgDstPadded.create(w, h, DImage::DImage_u8, 1);
  imgDstPadded.fill(128.);

  pu8Src = imgPadded.dataPointer_u8();
  pu8Dst = imgDstPadded.dataPointer_u8();

  progMax = h+1;


  // now perform the modified Niblack filtering that Sauvola used for text
  ww = radius1 + radius2 + 1;
#ifdef DEBUG
  if(ww != window){
    fprintf(stderr, "sauvolaNiblackThreshImage_ radius calculation is off\n");
    exit(1);
  }
#endif
  ww2 = ww * ww;
  memset(pu8Dst, 0, w*h*sizeof(unsigned char));//clear the destination image

  // create the two integral images (one for values and one for squared vals)
  pIntegralImgOfValues = (D_uint64*)malloc(sizeof(D_uint64) * w * h);
  D_CHECKPTR(pIntegralImgOfValues);
  pIntegralImgOfSquaredValues = (D_uint64*)malloc(sizeof(D_uint64) * w * h);
  D_CHECKPTR(pIntegralImgOfSquaredValues);
  pIntegralImgOfValues[0] = (D_uint64)(pu8Src[0]);
  pIntegralImgOfSquaredValues[0] =
    ((D_uint64)(pu8Src[0])) * ((D_uint64)(pu8Src[0]));
  for(int x = 1; x < w; ++x){ // first row
    pIntegralImgOfValues[x] = pIntegralImgOfValues[x-1]+(D_uint64)(pu8Src[x]);
    pIntegralImgOfSquaredValues[x] = pIntegralImgOfSquaredValues[x-1] +
      ((D_uint64)(pu8Src[x])) * ((D_uint64)(pu8Src[x]));
  }
  for(int y = 1, idx=w; y < h; ++y){
    pIntegralImgOfValues[idx] = pIntegralImgOfValues[idx-w] +
      (D_uint64)(pu8Src[idx]);
    pIntegralImgOfSquaredValues[idx] =pIntegralImgOfSquaredValues[idx-w] +
      ((D_uint64)(pu8Src[idx])) * ((D_uint64)(pu8Src[idx]));
    ++idx;
    for(int x = 1; x < w; ++x, ++idx){
      pIntegralImgOfValues[idx] = pIntegralImgOfValues[idx-1]+
	pIntegralImgOfValues[idx-w] + (D_uint64)(pu8Src[idx]) -
	pIntegralImgOfValues[idx-w-1];

      pIntegralImgOfSquaredValues[idx] = pIntegralImgOfSquaredValues[idx-1] +
	pIntegralImgOfSquaredValues[idx-w] +
	((D_uint64)(pu8Src[idx])) * ((D_uint64)(pu8Src[idx])) -
	pIntegralImgOfSquaredValues[idx-w-1];
    }
  }

  // now do the thresholding using the integral images
  for(int y = radius1+1; y < (h-radius2); ++y){
    int idxA, idxB, idxC, idxD;

    // update progress report and check if user cancelled the operation
    if(NULL != pProg){
      if(0 == (progCur & 0x000000ff)){
	if(0 != pProg->reportStatus(progCur, 0, progMax)){
	  // the operation has been cancelled
	  pProg->reportStatus(-1, 0, progMax); // report cancel acknowledged
	  return;
	}
      }
      ++progCur;//current progress
    }

    idxB = (y-radius1-1) * w;
    idxC = idxB + ww ;
    idxD = (y+radius2) * w;
    idxA = idxD + ww ;
    
    sum = 0.;
    sumofsquares = 0.;
    for(int x = radius1+1; x < (w-radius2); ++x){
      // calculate result for the new pixel using updates sum and sumofsquares
      mean =
	(pIntegralImgOfValues[idxA]+ pIntegralImgOfValues[idxB]-
	 pIntegralImgOfValues[idxC]- pIntegralImgOfValues[idxD]) / (double)ww2;
      stddev = (sqrt(  ((pIntegralImgOfSquaredValues[idxA] +
			 pIntegralImgOfSquaredValues[idxB] -
			 pIntegralImgOfSquaredValues[idxC] -
			 pIntegralImgOfSquaredValues[idxD])/((double)ww2))
		       - (mean*mean) ));


      thresh = mean * (1.+ K * (stddev / R - 1));
      //      thresh = mean + -0.2*stddev;
      
      // now set the pixel to black if less than thresh, white otherwise
      pu8Dst[y*w+x] = (pu8Src[y*w+x] < thresh) ? 0x00 : 0xff;
      ++idxA;
      ++idxB;
      ++idxC;
      ++idxD;
    }// end for x
  }//end for y
  free(pIntegralImgOfValues);
  free(pIntegralImgOfSquaredValues);

  // now copy the result from the padded temporary dst image into imgDst
  imgDstPadded.copy_(imgDst, radius1+1, radius1+1,
		     imgSrc.width(), imgSrc.height(),
		     imgSrc.getAllocMethod());
  if(NULL != pProg){ // report progress (complete)
    pProg->reportStatus(progMax, 0, progMax);
  }
}


/// performs Niblack locally adaptive thresholding
/** The window size is (2*radius+1).  K is the constant.

    This implementation uses integral images for speed, as described
    in "Efficient Implementation of Local Adaptive Thresholding
    Techniques Using Integral Images" by Shafait, Keysers, and Breuel
    to appear in Document Recognition and Retrieval (DRR) 2008.*/
void DThresholder::niblackThreshImage_(DImage &imgDst, const DImage &imgSrc,
				       int radius, double K,
				       DProgress *pProg, int numThreads){
  double mean; //mean
  double stddev; // standard deviation
  double thresh; // threshold value
  double sum;
  double sumofsquares;
  int ww, ww2; // window width, windowwidthSquared
  int w, h;
  DImage imgPadded;
  DImage imgDstPadded;
  D_uint8 *pu8Dst;
  D_uint8 *pu8Src;
  D_uint64 *pIntegralImgOfValues;
  D_uint64 *pIntegralImgOfSquaredValues;
  int progMax; // for reporting progress
  int progCur = 0; // for reporting progress

  if(imgSrc.getImageType() != DImage::DImage_u8){
    fprintf(stderr, "DThresholder::niblackThreshImage_() "
	    "only supports 8-bit grayscale images!\n");
    abort();
  }

  if(numThreads > 1){
    fprintf(stderr, "DThresholder::niblackThreshImage_() doesn't "
	    "support multiple threads yet.  Using a single thread...\n");
  }


  // first, pad the source image (mirror the borders) and a temporary dst image
  imgSrc.padEdges_(imgPadded, radius+1, radius, radius+1, radius,
		   DImage::DImagePadReplicate);

  w = imgPadded.width();
  h = imgPadded.height();
  imgDstPadded.create(w, h, DImage::DImage_u8, 1);
  imgDstPadded.fill(128.);

  pu8Src = imgPadded.dataPointer_u8();
  pu8Dst = imgDstPadded.dataPointer_u8();

  progMax = h+1;


  // now perform the modified Niblack filtering that Sauvola used for text
  ww = radius + radius + 1;
  ww2 = ww * ww;
  memset(pu8Dst, 0, w*h*sizeof(unsigned char));//clear the destination image

  // create the two integral images (one for values and one for squared vals)
  pIntegralImgOfValues = (D_uint64*)malloc(sizeof(D_uint64) * w * h);
  D_CHECKPTR(pIntegralImgOfValues);
  pIntegralImgOfSquaredValues = (D_uint64*)malloc(sizeof(D_uint64) * w * h);
  D_CHECKPTR(pIntegralImgOfSquaredValues);
  pIntegralImgOfValues[0] = (D_uint64)(pu8Src[0]);
  pIntegralImgOfSquaredValues[0] =
    ((D_uint64)(pu8Src[0])) * ((D_uint64)(pu8Src[0]));
  for(int x = 1; x < w; ++x){ // first row
    pIntegralImgOfValues[x] = pIntegralImgOfValues[x-1]+(D_uint64)(pu8Src[x]);
    pIntegralImgOfSquaredValues[x] = pIntegralImgOfSquaredValues[x-1] +
      ((D_uint64)(pu8Src[x])) * ((D_uint64)(pu8Src[x]));
  }
  for(int y = 1, idx=w; y < h; ++y){
    pIntegralImgOfValues[idx] = pIntegralImgOfValues[idx-w] +
      (D_uint64)(pu8Src[idx]);
    pIntegralImgOfSquaredValues[idx] =pIntegralImgOfSquaredValues[idx-w] +
      ((D_uint64)(pu8Src[idx])) * ((D_uint64)(pu8Src[idx]));
    ++idx;
    for(int x = 1; x < w; ++x, ++idx){
      pIntegralImgOfValues[idx] = pIntegralImgOfValues[idx-1]+
	pIntegralImgOfValues[idx-w] + (D_uint64)(pu8Src[idx]) -
	pIntegralImgOfValues[idx-w-1];

      pIntegralImgOfSquaredValues[idx] = pIntegralImgOfSquaredValues[idx-1] +
	pIntegralImgOfSquaredValues[idx-w] +
	((D_uint64)(pu8Src[idx])) * ((D_uint64)(pu8Src[idx])) -
	pIntegralImgOfSquaredValues[idx-w-1];
    }
  }

  // now do the thresholding using the integral images
  for(int y = radius+1; y < (h-radius); ++y){
    int idxA, idxB, idxC, idxD;

    // update progress report and check if user cancelled the operation
    if(NULL != pProg){
      if(0 == (progCur & 0x000000ff)){
	if(0 != pProg->reportStatus(progCur, 0, progMax)){
	  // the operation has been cancelled
	  pProg->reportStatus(-1, 0, progMax); // report cancel acknowledged
	  return;
	}
      }
      ++progCur;//current progress
    }

    idxB = (y-radius-1) * w;
    idxC = idxB + ww ;
    idxD = (y+radius) * w;
    idxA = idxD + ww ;
    
    sum = 0.;
    sumofsquares = 0.;
    for(int x = radius+1; x < (w-radius); ++x){
      // calculate result for the new pixel using updates sum and sumofsquares
      mean =
	(pIntegralImgOfValues[idxA]+ pIntegralImgOfValues[idxB]-
	 pIntegralImgOfValues[idxC]- pIntegralImgOfValues[idxD]) / (double)ww2;
      stddev = (sqrt(  ((pIntegralImgOfSquaredValues[idxA] +
			 pIntegralImgOfSquaredValues[idxB] -
			 pIntegralImgOfSquaredValues[idxC] -
			 pIntegralImgOfSquaredValues[idxD])/((double)ww2))
		       - (mean*mean) ));


      thresh = mean + K*stddev;
      
      // now set the pixel to black if less than thresh, white otherwise
      pu8Dst[y*w+x] = (pu8Src[y*w+x] < thresh) ? 0x00 : 0xff;
      ++idxA;
      ++idxB;
      ++idxC;
      ++idxD;
    }// end for x
  }//end for y
  free(pIntegralImgOfValues);
  free(pIntegralImgOfSquaredValues);

  // now copy the result from the padded temporary dst image into imgDst
  imgDstPadded.copy_(imgDst, radius+1, radius+1,
		     imgSrc.width(), imgSrc.height(),
		     imgSrc.getAllocMethod());
  if(NULL != pProg){ // report progress (complete)
    pProg->reportStatus(progMax, 0, progMax);
  }
}
 

 
 









// typedef struct{
//   DImage *pImgSrc;//pointer to the image (shared by all threads)
//   int numThreads;//how many threads are processing
//   int threadNum;//which thread number this is (0..numThreads-1)
//   int *rgNumCCs;//all threads point to same array
//   int tOtsu;
//   int otsuCCArea;
//   int otsuNumCCs;
// } CCTHRESHVAL_THREAD_PARMS;

void* DThresholder::getCCThreshVal_thread_func(void *params){
  CCTHRESHVAL_THREAD_PARMS *pparms;
  int numCCs;
  int bestT;
  int bestNumCCs;
  DImage imgCCmap;
  DImage imgThreshed;
  D_uint32 *pMapData;
  int otsuCCArea;
  int *rgCCsizes;
  int w, h; //width, height
  bool fSingleCCExceedsOtsuArea = false;

  pparms = (CCTHRESHVAL_THREAD_PARMS*)params;
  otsuCCArea = pparms->otsuCCArea;
  bestT = pparms->tOtsu;
  bestNumCCs = pparms->otsuNumCCs;

  // printf("Thread%d: otsuCCArea=%d tOtsu=%d numCCsOtsu=%d numThreads=%d\n",
  // 	 pparms->threadNum, otsuCCArea, bestT, bestNumCCs, pparms->numThreads);

  w = pparms->pImgSrc->width();
  h = pparms->pImgSrc->height();

  // iterate through all possible thresholds above the otsu threshold
  for(int curT = 1+(pparms->tOtsu)+(pparms->threadNum); curT < 255;
      curT += (pparms->numThreads)){
    // printf("curT=%d",curT);
    threshImage_(imgThreshed, *(pparms->pImgSrc), curT);
    DConnectedComponentLabeler::getCCimage_(imgCCmap, imgThreshed, &numCCs);
    // printf(" numCCs=%d\n", numCCs);
    if(numCCs < bestNumCCs){
      // less CCs is better, but we need to make sure we haven't gone
      // beyond the mark and started making everything into a big blob
      // by merging all the valid CCs.  Our heuristic for checking
      // this is to make sure no single CC has more pixels than the
      // total amount of ink at the Otsu level.
      rgCCsizes = (int*)calloc(numCCs+1, sizeof(int));
      if(!rgCCsizes){
	fprintf(stderr, "DThresholder::getCCThreshVal() no memory!\n");
	exit(1);
      }
      w = imgCCmap.width();
      h = imgCCmap.height();
      pMapData = imgCCmap.dataPointer_u32();
      for(int idx = 0, len=w*h; idx < len; ++idx){
	++(rgCCsizes[pMapData[idx]]);//we left slot 0 for the background count
      }
      fSingleCCExceedsOtsuArea = false;
      for(int i = 1; i <= numCCs; ++i){
	if(rgCCsizes[i] >= otsuCCArea){
	  // printf("  thread%d found CC#%d with area=%d\n",pparms->threadNum,
	  // 	 i,rgCCsizes[i]);
	  fSingleCCExceedsOtsuArea = true;
	  break; // break out of the inner loop: for(i=1; i <= numCCs; ...
	}
      }
      free(rgCCsizes);
      if(fSingleCCExceedsOtsuArea){
	// printf("  thread%d setting numCCs for curT=%d to -1\n",
	//        pparms->threadNum,curT);
	pparms->rgNumCCs[curT] = -1;
	break; // break out of the outer loop: for(int curT = 1+tOtsu...
      }
      else{
	pparms->rgNumCCs[curT] = numCCs;
	// we checked our heuristic and this thresh looks ok, so update best
	bestT = curT;
	bestNumCCs = numCCs;
      }
    }// end if(numCCs < bestNumCCs)
    else{
      // pparms->rgNumCCs[curT] = 1998877;
      pparms->rgNumCCs[curT] = numCCs;
    }
  }// end for(int curT=1+tOtsu...
  return NULL;
}

///find the global threshold that results in the fewest Connected Components
/**This method chooses the threshold value that results in the least
 * number of connected components (CCs).  This is based on the
 * assumption that fewer CCs will be found when words/characters are
 * well-connected instead of disjoint (disjoint characters resulting
 * from a threshold that is too low), and when the threshold is not so
 * high that extra noise is introduced.  In the limit, however,
 * everything gets connected into one big component, so heuristics are
 * used to short-circuit turning the whole page (or a large portion of
 * it) into one big blob.  Possible threshold values returned range
 * from the otsu value (the algorithm's starting point) to 254.*/
double DThresholder::getCCThreshVal(const DImage &img, DProgress *pProg,
				    int numThreads){
  double tOtsu;
  int bestT;
  int bestNumCCs;
  int numCCs;
  DImage imgCCmap;
  DImage imgThreshed;
  D_uint32 *pMapData;
  int otsuCCArea;
  int *rgCCsizes;
  int w, h; //width, height
  int *rgNumCCs;
  CCTHRESHVAL_THREAD_PARMS *rgThreadParms;
#ifndef D_NOTHREADS
  pthread_t rgThreadID[MAX_CCTHRESHVAL_THREADS];
#else
  numThreads = 1;
#endif

  rgNumCCs = new int[256];
  D_CHECKPTR(rgNumCCs);

  // printf("getCCThreshVal2 numThreads=%d\n", numThreads);

  rgThreadParms =
    (CCTHRESHVAL_THREAD_PARMS*)malloc(sizeof(CCTHRESHVAL_THREAD_PARMS) *
				      numThreads);
  D_CHECKPTR(rgThreadParms);

  memset(rgNumCCs, 99, sizeof(int)*256);
  tOtsu = getOtsuThreshVal(img);
  threshImage_(imgThreshed, img, tOtsu);
  DConnectedComponentLabeler::getCCimage_(imgCCmap, imgThreshed, &numCCs);

  bestT = (int)tOtsu;
  bestNumCCs = numCCs;

  // record how many pixels are ink at Otsu level
  pMapData = imgCCmap.dataPointer_u32();
  otsuCCArea = 0;
  w = imgCCmap.width();
  h = imgCCmap.height();
  for(int idx = 0, len=w*h; idx < len; ++idx){
    otsuCCArea += (int)(pMapData[idx] != 0);
  }

  // for(int tnum = 0; tnum < numThreads; ++tnum){
  for(int tnum = numThreads-1; tnum >= 0; --tnum){//so other threads launch
    //before thread 0 calls the function so they don't have to wait for it
    rgThreadParms[tnum].pImgSrc = &img;
    rgThreadParms[tnum].numThreads = numThreads;
    rgThreadParms[tnum].threadNum = tnum;
    rgThreadParms[tnum].rgNumCCs = rgNumCCs;
    rgThreadParms[tnum].tOtsu = (int)tOtsu;
    rgThreadParms[tnum].otsuCCArea = otsuCCArea;
    rgThreadParms[tnum].otsuNumCCs = numCCs;
#ifdef D_NOTHREADS
    getCCThreshVal_thread_func(rgThreadParms);
#else
    if(0 == tnum){//don't spawn thread zero. Use the current thread
      getCCThreshVal_thread_func(rgThreadParms);
    }
    else{//spawn all other threads besides thread zero
      if(0 != pthread_create(&rgThreadID[tnum], NULL,
			     DThresholder::getCCThreshVal_thread_func,
			     &rgThreadParms[tnum])){
	fprintf(stderr,"DThresholder::getCCThreshVal() failed to spawn thread "
		"#%d. Exiting.\n",tnum);
	exit(1);
      }
    }
#endif

  }
#ifndef D_NOTHREADS
  // wait for all threads to finish
  for(int tnum = 1; tnum < numThreads; ++tnum){
    if(pthread_join(rgThreadID[tnum],NULL)){
      fprintf(stderr, "DThresholder::getCCThreshVal() failed to join "
	      "thread %d. Exiting.\n", tnum);
      exit(1);
    }
  }
#endif

  // iterate through all possible thresholds above the otsu threshold
  // and see which is best (until -1 is hit because that means a CC is too big)
  for(int curT = 1+(int)tOtsu; curT < 255; ++curT){
    if(-1 == rgNumCCs[curT])
      break;
    if(rgNumCCs[curT] < bestNumCCs){
      bestNumCCs = rgNumCCs[curT];
      bestT = curT;
    }
  }
  // printf("otsu=%d:%d\n",(int)tOtsu,rgNumCCs[(int)tOtsu]);
  // for(int i=(int)tOtsu; i < 256; ++i){
  //   printf("%d:%d   ",i,rgNumCCs[i]);
  // }
  // printf("\n");
  // printf("bestT=%d\n",bestT);
  delete [] rgNumCCs;
  return (double)bestT;
}
