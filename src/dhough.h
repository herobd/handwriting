#ifndef DHOUGH_H
#define DHOUGH_H

#include "dimage.h"
#include "dmath.h"
///This class provides functionality for finding lines with hough transform


class DHough{
public:
  static void houghEdgeImage_(DImage &imgAccum, const DImage &imgEdge,
			      double angleResDegs = 1.0,
			      double pixelRes = 1.0,
			      D_uint8 ignoreGradientBelowThresh = 10,
			      bool fWeightByGradient = true,
			      DImage *imgEdgeDirs = NULL,
			      double angleWindowDegs = 360.);
  static void houghImageCalcEdges_(DImage &imgAccum, const DImage &imgSrc,
				   double angleResDegs = 1.0,
				   double pixelRes = 1.0,
				   D_uint8 ignoreGradientBelowThresh = 10,
				   bool fWeightByGradient = true,
				   double angleWindowDegs = 360.);
  static void getRTheta(int w, int h, double angleResDegs,
			double pixelRes, int x_R, int y_Theta,
			double *R, double *Theta);
  static void getRTheta(DImage &imgAccum, double angleResDegs,
			double pixelRes, int x_R, int y_Theta,
			double *R, double *Theta);
  static void getXY(int w, int h, double angleResDegs,
		    double pixelRes, double R, double Theta,
		    int *x_R, int *y_Theta);
  static void getXY(DImage &imgAccum, double angleResDegs,
		    double pixelRes, double R, double Theta,
		    int *x_R, int *y_Theta);
  
//   static void getLinePoints(int w, int h, double R, double Theta,
// 			    double *x0, double *y0, double *x1, double *y1);
  static void getLinePointsDouble(int w, int h, double R, double Theta,
				  bool fThetaInDegrees,
				  double *x0, double *y0,
				  double *x1, double *y1);
  static void getLinePointsDouble(DImage &imgAccum, double R, double Theta,
				  bool fThetaInDegrees,
				  double *x0, double *y0,
				  double *x1, double *y1);
};

///Calculate the R and Theta value for an x,y position in an accumulator
/** w, h are the width and height of the accumulator image,
    angleResDegs is the angular resolution in degrees (i.e, if the
    accumulator handles angles of 0.1 degrees, angleResDegs==0.1),
    pixelRes is the pixel resolution of the accumulator R values, x_R
    and y_Theta are the x,y (R_index, Theta_index) within the
    accumulator for which we want to know the corresponding values of
    R and Theta.*/
inline void DHough::getRTheta(int w, int h, double angleResDegs,
			      double pixelRes, int x_R, int y_Theta,
			      double *R, double *Theta){
  abort();
}

///Calculate the R and Theta value for an x,y position in an accumulator
/** This function is provided for convenience.  It just calls the
    other version of this function using the width and height of
    imgAccum. */
inline void DHough::getRTheta(DImage &imgAccum, double angleResDegs,
			      double pixelRes, int x_R, int y_Theta,
			      double *R, double *Theta){
  getRTheta(imgAccum.width(), imgAccum.height(), angleResDegs, pixelRes,
	    x_R, y_Theta, R, Theta);
}

///Calculate the x,y position corresponding to R,Theta in an accumulator
/** w, h are the width and height of the accumulator image,
    angleResDegs is the angular resolution in degrees (i.e, if the
    accumulator handles angles of 0.1 degrees, angleResDegs==0.1),
    pixelRes is the pixel resolution of the accumulator R values, x_R
    and y_Theta are the x,y (R_index, Theta_index) coordinate pair
    within the accumulator that corresponds to R,Theta.  Please note,
    the x,y returned are the location within the accumulator array
    (since the resolution may not be 1), NOT an x,y position in the
    source image that the Hough transform is performed on.  (In the
    source image, an R,Theta pair represent a line, not a point).*/
inline void DHough::getXY(int w, int h, double angleResDegs,
			  double pixelRes, double R, double Theta,
			  int *x_R, int *y_Theta){
  abort();
}

///Calculate the x,y position corresponding to R,Theta in an accumulator
/** This function is provided for convenience.  It just calls the
    other version of this function using the width and height of
    imgAccum. */
inline void DHough::getXY(DImage &imgAccum, double angleResDegs,
			  double pixelRes, double R, double Theta,
			  int *x_R, int *y_Theta){
  getXY(imgAccum.width(), imgAccum.height(), angleResDegs, pixelRes,
	R, Theta, x_R, y_Theta);
}


///Find the endpoints of the line defined by R and Theta within the image
/**If the line defined by R,Theta passes through this image, the
   endpoints of where the line crosses the borders of the image are
   returned in the locations of x0,y0, x1,y1.  If the line does not
   pass through the image, then all values will be -1. */
inline void DHough::getLinePointsDouble(int w, int h, double R, double Theta,
					bool fThetaInDegrees,
					double *x0, double *y0,
					double *x1, double *y1){
  double thetaRad;
  double m;
  double b; // y-position where line crosses x=0
  double c; // y-pos where line crosses x=w-1
  double d; // x-pos where line crosses y=0
  double e; // x-pos where line crosses y=h-1

  (*x0) = -1;
  (*y0) = -1;
  (*x1) = -1;
  (*y1) = -1;

  if(fThetaInDegrees){ // degrees
    thetaRad = DMath::degreesToRadians(Theta);
  }
  else{ // radians
    thetaRad = Theta;
  }
  if((Theta < (-50 * M_PI)) || (Theta >= (50 * M_PI))){
    fprintf(stderr, "DHough::getLinePointsDouble() Theta=%f should be in "
	    "range [0..2PI)\n", Theta);
    abort();
  }
  while(Theta < 0.){
    Theta += 2*M_PI;
  }
  while(Theta > 2*M_PI){
    Theta -= 2*M_PI;
  }
  if( (fabs(thetaRad) < 0.00001) ||
      (fabs(M_PI - thetaRad) < 0.00001) ||
      (fabs(2*M_PI - thetaRad) < 0.00001)){
    (*x0) = 0;
    (*y0) = R;
    (*x1) = w-1;
    (*y1) = (*y0);
    return;
  }
  if( (fabs(M_PI_2 - thetaRad) < 0.00001) ||
      (fabs(M_PI+M_PI_2 - thetaRad) < 0.00001) ){
    (*x0) = R;
    (*y0) = 0;
    (*x1) = (*x0);
    (*y1) = h-1;
    return;
  }

  //  m = -1. / tan(thetaRad);
//   b = R / sin(thetaRad);
  m = -1. * tan(thetaRad);
  b = R / sin(M_PI_2 -thetaRad);
  c = m * (w-1) + b;
  d = (0-b) / m;
  e = ( (h-1) - b) / m;

  printf("m=%.2f b=%.1f c=%.1f d=%.1f e=%.1f thetaRad=%.1f(%.1fdeg) "
	 "w=%d h=%d R=%.1f\n\n",m,b,c,d,e,thetaRad,
	 DMath::radiansToDegrees(thetaRad), w, h, R);
  if( (b>= 0.) && (b <= (h-1))){
    (*x0) = 0;
    (*y0) = b;
    if( (c>= 0.) && (c <= (h-1))){
      (*x1) = w-1;
      (*y1) = c;
      return;
    }
    if( (d >= 0.) && (d <= (w-1))){
      (*x1) = d;
      (*y1) = 0;
      return;
    }
    (*x1) = e;
    (*y1) = h-1;
    return;
  }
  if( b < 0.){
    if( (d >=0) && (d <= (w-1))){
      (*x0) = d;
      (*y0) = 0;
      if( (e >=0) && (e <= (w-1))){
	(*x1) = e;
	(*y1) = h-1;
	return;
      }
      (*x1) = w-1;
      (*y1) = c;
      return;
    }
    return; // all -1
  }
  if( b > (h-1) ){
    if( (e >=0) && (e <= (w-1))){
      (*x0) = e;
      (*y0) = h-1;
      if( (d >=0) && (d <= (w-1))){
	(*x1) = d;
	(*y1) = 0;
	return;
      }
      (*x1) = w-1;
      (*y1) = c;
      return;
    }
    return; // all -1
  }
  printf("huh?\n");
  abort();
  return;
}

///Provided for convenience.  Calls getLinePointsDouble with imAccum's w and h
inline void DHough::getLinePointsDouble(DImage &imgAccum,
					double R, double Theta,
					bool fThetaInDegrees,
					double *x0, double *y0,
					double *x1, double *y1){
  getLinePointsDouble(imgAccum.width(), imgAccum.height(),
		      R, Theta, fThetaInDegrees, x0, y0, x1, y1);
}

#endif


