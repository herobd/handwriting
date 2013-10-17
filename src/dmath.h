#ifndef DMATH_H
#define DMATH_H

#include <math.h>

#if 1
#include <stdio.h>
#include <stdlib.h>
#endif

/** This class has some miscellaneous math functions for use throughout the
    project

This code includes a function pnpoly() that tests whether a point is in
a polygon.  It is available at 
http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html

We include it using the following License:

Copyright (c) 1970-2003, Wm. Randolph Franklin

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
   2. Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution.
   3. The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


class DMath{
public:

  static double degreesToRadians(const double d);
  static double radiansToDegrees(const double r);
  static double variance(double *rgVals, int len);
  static double variance(int *rgVals, int len);
  static double distPointLine(double x, double y, double x0, double y0,
			      double x1, double y1);
  static double distPointLineSegment3D(double x, double y, double z,
				       double x0, double y0, double z0,
				       double x1, double y1, double z1);
  static double whichSideOfLine(double x, double y,
				double x0, double y0, double x1, double y1);
  static bool pnpoly(int nvert, double *vertx, double *verty,
		     double x, double y);
  static bool getBarycentricCoords(double *L1, double *L2, double *L3,
				   double x, double y,
				   double Ax, double Ay, double Bx, double By,
				   double Cx, double Cy);
  static bool getCoordsFromBarycentric(double *x, double *y,
				       double L1, double L2, double L3,
				       double Ax, double Ay,
				       double Bx, double By,
				       double Cx, double Cy);
};

///returns the number of radians in d degrees
inline double DMath::degreesToRadians(const double d){
  return d * (M_PI / 180.);
}

///returns the number of degrees in r radians
inline double DMath::radiansToDegrees(const double r){
  return r * (180./M_PI);
}

//TODO: speed this up!
// inline double DMath::variance(double *rgVals, int len){
//   double sum = 0.;
//   double avg = 0.;
//   int i;

//   for(i = 0; i < len; ++i)
//     sum += rg[i];
//   avg = sum / len;
//   sum = 0.;
//   for(i = 0; i < len; ++i)
//     sum += ((rg[i] - avg) * (rg[i] - avg));
//   return (sum / len);
// }

///Returns the variance of the values in an array of doubles that is len long
inline double DMath::variance(double *rgVals, int len){
  double sumX = 0.;
  double sumX2 = 0.;

  if(0 == len)
    return 0.;
  for(int i = 0; i < len; ++i){
    register double val;
    val = rgVals[i];
    sumX += val;
    sumX2 += (val * val);
  }
  return (sumX2 - sumX*sumX/len)/len;
}

///Returns the variance of the values in an array of ints that is len long
inline double DMath::variance(int *rgVals, int len){
  double sumX = 0.;
  double sumX2 = 0.;

  if(0 == len)
    return 0.;
  for(int i = 0; i < len; ++i){
    register double val;
    val = (double)(rgVals[i]);
    sumX += val;
    sumX2 += (val * val);
  }
  return (sumX2 - sumX*sumX/len)/len;
}

///Return the distance between a point and a line defined by 2 pts on the line
/**Returns the Euclidean distance between point (x,y) and the infinite line defined by two points, (x0,y0) and (x1,y1), that are on the line.
 */
inline double DMath::distPointLine(double x, double y, double x0, double y0,
				   double x1, double y1){
  register double x1mx0, y1my0;
  //The method used to calculate this is described at the following URLs:
  //http://www.allegro.cc/forums/thread/589720
  //http://www.softsurfer.com/Archive/algorithm_0102/algorithm_0102.htm
  //(the two equations are equivalent, as seen through algebraic manipulation)
  x1mx0 = x1-x0;
  y1my0 = y1-y0;
  return fabs((x-x0)*y1my0-x1mx0*(y-y0)) / sqrt(x1mx0*x1mx0+y1my0*y1my0);
  //The directed distance could be calculated by using the signed value in the
  // numberator instead of taking the absolute value (the sign of the result
  // would indicate which side of the line the point is on)

}

///Return the distance between a point and a segment defined by 2 pts in 3D
/**Returns the Euclidean distance between point (x,y,z) and the line
   segment defined by the two endpoints, (x0,y0,z0) and (x1,y1,z1).
   (This may be different than the distance between the point and the
   infinite line, as only the segment between the endpoints is used).
 */
inline double DMath::distPointLineSegment3D(double x, double y, double z,
					    double x0, double y0, double z0,
					    double x1, double y1, double z1){
  double dblU;

  
  if( (x0 == x1) && (y0 == y1) && (z0 == z1)){
    // the segment is really a point.  just return dist to point
    return sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0));
  }

  dblU = ( (x-x0)*(x1-x0) + (y-y0)*(y1-y0) + (z-z0)*(z1-z0) ) /
    ( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0) );

  if(dblU < 0.){
#if 1
    if( sqrt( (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1)) <
	sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0)) ){
      fprintf(stderr, "DMath::distPointLineSegment3D() logic error (1)!\n");
      abort();
    }
#endif
    return sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0));
  }
  else if(dblU > 1.){
#if 1
    if( sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0)) <
	sqrt( (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1)) ){
      fprintf(stderr, "DMath::distPointLineSegment3D() logic error! (2)\n");
      abort();
    }
#endif
    return sqrt( (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1));
  }
  else{
    double xp, yp, zp;
    xp = x0 + dblU * (x1-x0);
    yp = y0 + dblU * (y1-y0);
    zp = z0 + dblU * (z1-z0);

    return sqrt( (x-xp)*(x-xp) + (y-yp)*(y-yp) + (z-zp)*(z-zp) );
  }


  return 0.;
}

///determine which side of a (directed) line point x,y is on
/**assumes directional line from x0,y0 to x1,y1.  If the point x,y is
   left of the line, then the return value is negative, if on the
   line, then zero, if right of the line then positive.**/
inline double DMath::whichSideOfLine(double x, double y,
				     double x0, double y0, double x1,double y1){
  return (x1-x0)*(y-y0)-(y1-y0)*(x-x0);
}


///Return true if point x, y is within a polygon having nvert vertices
//I changed the parameter types from float to double and the return type
//(as well as type of variable "c") from int to bool, and the name of the
//point variables from testx, testy to x, y.  -djk
inline bool DMath::pnpoly(int nvert, double *vertx, double *verty,
			  double x, double y){
  int i, j;
  bool c = false;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>y) != (verty[j]>y)) &&
	 (x < (vertx[j]-vertx[i]) * (y-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}


///given 2-D triangle vertices A,B,C calculate bary coords (L1,L2,L3) of pt x,y
inline bool DMath::getBarycentricCoords(double *L1, double *L2, double *L3,
					double x, double y,
					double Ax, double Ay,
					double Bx, double By,
					double Cx, double Cy){
  double denom;
  denom = (By-Cy)*(Ax-Cx)+(Cx-Bx)*(Ay-Cy);
  if(denom==0.)
    denom = 0.00000001;
  (*L1) = (By-Cy)*(x-Cx)+(Cx-Bx)*(y-Cy)/denom;
  denom = (Cy-Ay)*(x-Cx)+(Ax-Cx)*(y-Cy);
  if(denom==0.)
    denom = 0.00000001;
  (*L2) = (Cy-Ay)*(Bx-Cx)+(Ax-Cx)*(By-Cy)/denom;
  (*L3) = 1. - (*L1) - (*L2);
  return true;
}

inline bool DMath::getCoordsFromBarycentric(double *x, double *y,
					    double L1, double L2, double L3,
					    double Ax, double Ay,
					    double Bx, double By,
					    double Cx, double Cy){
  (*x) = L1*Ax+L2*Bx+L3*Cx;
  (*y) = L1*Ay+L2*By+L3*Cy;
  return true;
}


#endif
