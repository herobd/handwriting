#include "dhough.h"
#include "dimage.h"
#include "dprogress.h"
#include "dinstancecounter.h"
#include "dedgedetector.h"
#include "dmath.h"

#include "dconvolver.h"
#include "dkernel2d.h"

///Put the line Hough transform accumulator for imgSrc into imgAcum
/**angleResDegs is the resolution (in degrees) of the Theta angles in
   the accumulator (smaller value means finer resolution).  pixelRes
   is the resolution of the R values in the accumulator.  
   
   If not NULL, imgEdgeDirs should have the gradient direction for
   each pixel (in radians), stored as doubles.  angleWindowDegs (if
   not 360.) is the window of angles (relative to the edge direction
   at any given pixel) for which accumulation is performed.  The
   window is actually twice that angle since the angle is for either
   side of the edge.  fWeightByGradient is currently not
   used. Accumulation only happens for pixels where the edge strength
   (gradient magnitude) is not less than ignoreGradientBelowThresh.
 */
void DHough::houghEdgeImage_(DImage &imgAccum, const DImage &imgEdge,
			     double angleResDegs,
			     double pixelRes,
			     D_uint8 ignoreGradientBelowThresh,
			     bool fWeightByGradient,
			     DImage *imgEdgeDirs,
			     double angleWindowDegs){
  int accumW_r;
  int accumH_theta;
  int numSinesCosines; // how many angles will be used to vote at each position
  double *rgCosines;// cosines used to vote at each position
  double *rgSines;
  int w, h;
  double angleWindowRads; // angleWindowDegs converted to degrees
  double dAngRads; // delta for incrementing angle in loop (radians)
  D_uint8 *pSrc;
  double *pdblEdgeDirs;
  D_uint32 *pAccum;

  DImage imgTmp = imgEdge;
  imgTmp.save("/tmp/edges1.pgm");

  if((angleWindowDegs < 360.) && (imgEdgeDirs==NULL)){
    fprintf(stderr, "DHough::houghImage_() if angleWindowDegs < 360, then "
	    "imgEdgeDirs must be provided (and already calculated)\n");
    abort();
  }
  if(imgEdge.getImageType() != DImage::DImage_u8){
    fprintf(stderr, "DHough::houghImage_() currently only supports 8-bit "
	    "grayscale edge images\n");
    abort();
  }
  if(angleResDegs <= 0.){
    fprintf(stderr, "DHough::houghImage_() angleResDegs <= 0.\n");
    abort();
  }
  if(pixelRes <= 0.){
    fprintf(stderr, "DHough::houghImage_() pixelRes <= 0.\n");
    abort();
  }

  w = imgEdge.width();
  h = imgEdge.height();

  // create the accumulator
  accumW_r = (int)(0.9999 + (sqrt(w*w+h*h)) / pixelRes);
  accumH_theta = (int)(0.9999 + 180. / angleResDegs);
  imgAccum.create(accumW_r, accumH_theta, DImage::DImage_u32, 1);
  imgAccum.clear();
  pAccum = imgAccum.dataPointer_u32();
  dAngRads = DMath::degreesToRadians(angleResDegs);

//   printf("w=%d h=%d  accumW_r=%d accumH_theta=%d\n",w,h,accumW_r,accumH_theta);
  
  // run through the image recording votes in the accumulator
  if(angleWindowDegs < 360.){ // need to limit the angles
    int numAngles;
    int numAngles_2;
    int thetaIdxCur;
    pSrc = imgEdge.dataPointer_u8();
    pdblEdgeDirs = imgEdgeDirs->dataPointer_dbl();

    numSinesCosines = (int)(0.9 + 1 * M_PI / dAngRads);
//     printf("numSinesCosines=%d\n", numSinesCosines);
    rgCosines = (double*)malloc(numSinesCosines * sizeof(double));
    D_CHECKPTR(rgCosines);
    rgSines = (double*)malloc(numSinesCosines * sizeof(double));
    D_CHECKPTR(rgSines);
    for(int i = 0; i < numSinesCosines; ++i){
      rgSines[i] = sin(i * dAngRads);
      rgCosines[i] = cos(i * dAngRads);
    }
    numAngles = (int)(0.5+ angleWindowDegs / angleResDegs) * 2 + 1;
    numAngles_2 = numAngles / 2;
//     printf("numAngles=%d numAngles_2=%d\n", numAngles, numAngles_2);
    for(int y = 0, idx = 0; y < h; ++y){
      for(int x = 0; x < w; ++x, ++idx){
	if(pSrc[idx] < ignoreGradientBelowThresh)
	  continue;
	thetaIdxCur = (int)(0.5 + pdblEdgeDirs[idx] / dAngRads) - numAngles_2;
	if(thetaIdxCur < 0)
	  thetaIdxCur += numSinesCosines;
	else if(thetaIdxCur >= numSinesCosines)
	  thetaIdxCur -= numSinesCosines;
	for(int angNum = 0; angNum <= numAngles; ++angNum){
	  double rCur;
	  double thetaCur;
	  int rIdx, thetaIdx;

	  thetaCur = thetaIdxCur * dAngRads;
	  rCur = y * rgCosines[thetaIdxCur] + x * rgSines[thetaIdxCur];
// 	  if(rCur < 0)
// 	    rCur += accumW_r;
	  thetaIdx = (int)(thetaIdxCur);
	  rIdx = (int)(rCur / pixelRes);
// 	  if((x==201) && (y==150)){
// 	    printf("rCur=%.2f thetaCur=%.2f(%.2fdeg) rIdx=%d thetaIdxCur=%d "
// 		   "thetaIdx=%d edgeDir=%f(%.2fdeg)\n", rCur, thetaCur,
// 		   DMath::radiansToDegrees(thetaCur), rIdx, thetaIdxCur,
// 		   thetaIdx, pdblEdgeDirs[idx], 
// 		   DMath::radiansToDegrees(pdblEdgeDirs[idx]));
// 	  }
	  if( (rIdx >= 0) && (rIdx < accumW_r) &&
	      (thetaIdx >= 0) && (thetaIdx < accumH_theta)){
	    pAccum[thetaIdx*accumW_r + rIdx] += 1;
	  }
// 	  else{
//   	    printf("x=%d y=%d can't accumulate! rIdx=%d thetaIdx=%d rCur=%.01f"
// 		   " thetaCur=%.01f(%.01fdeg) angNum=%d "
// 		   "edgeDir=%f(%.01fdeg)\n",
// 		   x, y, rIdx, thetaIdx, rCur, thetaCur,
// 		   DMath::radiansToDegrees(thetaCur), angNum,
// 		   DMath::radiansToDegrees(pdblEdgeDirs[idx]));
// 	  }
	  
	  ++thetaIdxCur;
	  if(thetaIdxCur >= numSinesCosines)
	    thetaIdxCur -= numSinesCosines;
	}
      }
    }
  }
  else{
#if 0
    pSrc = imgEdge.dataPointer_u8();
    for(int y = 0, idx = 0; y < h; ++y){
      for(int x = 0; x < w; ++x, ++idx){
	if(pSrc[idx] < ignoreGradientBelowThresh)
	  continue;
	for(double thetaCur = 0.; thetaCur < M_PI; thetaCur += dAngRads){
	  double rCur;
	  rCur = x * cos(thetaCur) + y * sin(thetaCur);
	  if( (rCur >= 0) && (rCur < accumW_r) &&
	      (thetaCur >= 0) && (thetaCur < accumH_theta))
	    pAccum[((int)thetaCur)*accumW_r + (int)rCur] +=
	      (fWeightByGradient ? (pSrc[idx]) : 1);
	}
      }
    }
#else //--------------------------------------Try this:
    pSrc = imgEdge.dataPointer_u8();
    numSinesCosines = M_PI / dAngRads;
//     printf("numSinesCosines=%d\n", numSinesCosines);
    rgCosines = (double*)malloc(numSinesCosines * sizeof(double));
    D_CHECKPTR(rgCosines);
    rgSines = (double*)malloc(numSinesCosines * sizeof(double));
    D_CHECKPTR(rgSines);
    for(int i = 0; i < numSinesCosines; ++i){
      rgSines[i] = sin(i * dAngRads);
      rgCosines[i] = cos(i * dAngRads);
    }

    for(int y = 0, idx = 0; y < h; ++y){
      for(int x = 0; x < w; ++x, ++idx){
	if(pSrc[idx] < ignoreGradientBelowThresh)
	  continue;
	for(int thetaIdxCur = 0; thetaIdxCur < numSinesCosines; ++thetaIdxCur){
	  double rCur;
	  double thetaCur;
	  int rIdx, thetaIdx;
	  thetaCur = thetaIdxCur * dAngRads;
	  rCur = y * rgCosines[thetaIdxCur] + x * rgSines[thetaIdxCur];
// 	  if(rCur < 0)
// 	    rCur += accumW_r;
	  thetaIdx = (int)(thetaIdxCur /** angleResDegs*/);
	  rIdx = (int)(rCur / pixelRes);
// 	  if((x==201) && (y==150)){
// 	    printf("rCur=%.2f thetaCur=%.2f rIdx=%d thetaIdxCur=%d thetaIdx=%d\n", rCur, thetaCur, rIdx, thetaIdxCur, thetaIdx);
// 	  }
	  if( (rIdx >= 0) && (rIdx < accumW_r) &&
	      (thetaIdx >= 0) && (thetaIdx < accumH_theta))
	    pAccum[thetaIdx*accumW_r + rIdx] += 1;
	}
      }
    }

#endif
  }
#if 0
 {
   DImage imgA2;
   DKernel2D kern;
   DImage imgTmp3;

   imgA2 = imgAccum;
   kern.setGauss(1,1);
   DConvolver::convolve_(imgTmp3, imgA2, kern, false, false, true);
   imgA2 = imgTmp3;
   imgA2.setDataRange(0, 255);
   imgA2 = imgA2.convertedImgType(DImage::DImage_u8);
   imgA2.save("/tmp/accum.pgm");
 }
#endif
}

///Calculates edge image and calls houghEdgeImage_()
/**This function is provided for convenience.  It calculates the edge
   image and the edge directions, then calls
   houghEdgeImage_(). Currently, the Sobel filter is used for the edge
   detector.
*/
void DHough::houghImageCalcEdges_(DImage &imgAccum, const DImage &imgSrc,
				  double angleResDegs,
				  double pixelRes,
				  D_uint8 ignoreGradientBelowThresh,
				  bool fWeightByGradient,
				  double angleWindowDegs){
  DImage imgEdges;
  DImage imgEdgeDirs;

  if(angleWindowDegs < 360.){
    DEdgeDetector::sobel_(imgEdges, imgSrc, &imgEdgeDirs);
    DEdgeDetector::convertGradDirsToEdgeDirs(imgEdgeDirs, true);

    imgEdges.save("/tmp/edges1.pgm");
    houghEdgeImage_(imgAccum, imgEdges, angleResDegs, pixelRes,
		    ignoreGradientBelowThresh, fWeightByGradient,
		    &imgEdgeDirs, angleWindowDegs);
  }
  else{
    DEdgeDetector::sobel_(imgEdges, imgSrc, NULL);
    houghEdgeImage_(imgAccum, imgEdges, angleResDegs, pixelRes,
		    ignoreGradientBelowThresh, fWeightByGradient,
		    NULL, angleWindowDegs);
  }
}
