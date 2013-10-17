#include "dglobalskew.h"
#include "dprofile.h"
#include "dedgedetector.h"
#include "dthresholder.h"
#include "dmath.h"

double DGlobalSkew::getSkewAng_postl(const DImage &img,
				     double angFirstDeg,
				     double angLastDeg,
				     double angDeltaDeg,
				     double angDeltaDegLowRes,
				     bool fUseEdges){
  fprintf(stderr, "Not Yet Implemented!  (see www.leptonica.com/papers/skew-measurement.pdf for implementation details)\n");
  exit(1);
}


///Calculate the skew angle as the angle for which profile variance is maximum
double DGlobalSkew::getSkewAng_var(const DImage &img, double angFirstDeg,
				   double angLastDeg, double angDeltaDeg){
  DProfile prof;
  double bestAng, bestVar, var;
  double ang;
  double *rgProf;
  int profLen;
  bestAng = angFirstDeg;
  switch(img.getImageType()){
    case DImage::DImage_u8:
      bestVar = 0.; /* will get fixed first time in loop as long as var>0. */
      bestAng = angFirstDeg;
      for(ang = angFirstDeg; ang <= angLastDeg; ang += angDeltaDeg){
	prof.getAngledVertProfile(img, ang, true);
	rgProf = prof.dataPointer();
	profLen = prof.dataLen();
	if(profLen > 3){// ignore end conditions on profile
	  var = DMath::variance(&(rgProf[1]), profLen-2);
	}
	else
	  var = DMath::variance(rgProf, profLen);
	if(var > bestVar){
	  bestVar = var;
	  bestAng = ang;
	}
      }// end for(ang...
      break;
    default:
      fprintf(stderr,
	      "DGlobalSkew::getSkewAng_var() unsupported image type\n");
      break;
  }
  return bestAng;
}




///Calculate the skew angle as the angle for which profile variance is maximum
/**This version of the function shrinks the image down to a reasonably
 * small size (if it is large to begin with), and finds the angle for
 * which variance of the angled vertical profile is maximal with the
 * small version of the image by calling getSkewAng_var(), then
 * refines that angle on the full-resolution image.  The user can
 * specify the interval of angles checked for both the regular and low
 * resolution images by setting the parameters angDeltaDeg and
 * angDeltaDegLowRes, respectively. The starting and ending angles to
 * check at the low resolution are specified by angFirstDeg and
 * angLastDeg.
 */
double DGlobalSkew::getSkewAng_fast(const DImage &img, double angFirstDeg,
				    double angLastDeg, double angDeltaDeg,
				    double angDeltaDegLowRes, bool fUseEdges){
  double bestAng;
  double bestAng2;
  int wLowRes, hLowRes;
  DImage imgLowRes;
  DImage imgEdges;
  int numHalves;
  const int bgMinW = 400;//was 100x100
  const int bgMinH = 400;
  const int bgMaxW = 1000;//was 800x800
  const int bgMaxH = 1000;
  double otsuVal;

  numHalves = 0;
  wLowRes = img.width();
  hLowRes = img.height();
  while(((wLowRes > bgMinW*2) && (hLowRes > bgMinH*2)) &&
	((wLowRes > bgMaxW) || (hLowRes > bgMaxH))){
    wLowRes /= 2;
    hLowRes /= 2;
    ++numHalves;
  }
  if(0 == numHalves){
    if(fUseEdges){
      DEdgeDetector::prewittHorizEdges_(imgEdges, img);
      otsuVal = DThresholder::getOtsuThreshVal(imgEdges);
      DThresholder::threshImageSpecial_(imgEdges, imgEdges, otsuVal, otsuVal,
					0., -1., -1.);
      return getSkewAng_var(imgEdges, angFirstDeg, angLastDeg, angDeltaDeg);
    }
    return getSkewAng_var(img, angFirstDeg, angLastDeg, angDeltaDeg);
  }
  // scale down
  img.scaledDownPow2_(imgLowRes, numHalves);
  if(fUseEdges){
    DEdgeDetector::prewittHorizEdges_(imgEdges, imgLowRes);
    otsuVal = DThresholder::getOtsuThreshVal(imgEdges);
    DThresholder::threshImageSpecial_(imgEdges, imgEdges, otsuVal, otsuVal,
				      0., -1., -1.);
    bestAng =
      getSkewAng_var(imgEdges, angFirstDeg, angLastDeg, angDeltaDegLowRes);
  }
  else{
    bestAng =
      getSkewAng_var(imgLowRes, angFirstDeg, angLastDeg, angDeltaDegLowRes);
  }
//   printf("on lowRes image, bestAng=%.2f\n", bestAng);



  if(fUseEdges){
    DEdgeDetector::prewittHorizEdges_(imgEdges, img);
    otsuVal = DThresholder::getOtsuThreshVal(imgEdges);
    DThresholder::threshImageSpecial_(imgEdges, imgEdges, otsuVal, otsuVal,
				      0., -1., -1.);
    bestAng2 =getSkewAng_var(imgEdges, bestAng-angDeltaDegLowRes+angDeltaDeg,
			      bestAng+angDeltaDegLowRes-angDeltaDeg,
			      angDeltaDeg);
  }
  else{
    bestAng2 =
      getSkewAng_var(img, bestAng-angDeltaDegLowRes+angDeltaDeg,
		     bestAng+angDeltaDegLowRes-angDeltaDeg,
		     angDeltaDeg);
  }
//   printf("on regular image refinement, bestAng=%.2f\n", bestAng2);

  return bestAng2;
}
