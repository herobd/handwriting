#ifndef DCOLORSPACE_H
#define DCOLORSPACE_H

#include "dimage.h"
#include <math.h>

class DColorSpace{
public:
  // convert image from RGB
  static void convertRGBImageToCIELab(DImage &imgDst, const DImage &imgSrc);
  static void convertRGBImageToHSV(DImage &imgDst, const DImage &imgSrc);

  // convert image back to RGB
  static void convertCIELabImageToRGB(DImage &imgDst, const DImage &imgSrc);
  static void convertHSVImageToRGB(DImage &imgDst, const DImage &imgSrc);


  // convert a single pixel value to/from RGB
  static void getCIELabFromRGB(D_uint8 R, D_uint8 G, D_uint8 B,
			       float *L, float *a, float *b);
  static void getHSVFromRGB(D_uint8 R, D_uint8 G, D_uint8 B,
			    float *H, float *S, float *V);
  static void getRGBFromHSV(float H, float S, float V,
			    D_uint8 *R, D_uint8 *G, D_uint8 *B);
  static void getRGBFromCIELab(float L, float a, float b,
			       D_uint8 *R, D_uint8 *G, D_uint8 *B);
  
};

///Returns the L,a,B values given R,G,B (where R,G,B are [0..255])
inline void DColorSpace::getCIELabFromRGB(D_uint8 R, D_uint8 G, D_uint8 B,
					  float *L, float *a, float *b){
  float fltR, fltG, fltB;
  float X, Y, Z;
  const float XYZn = 1440.954; // since Xn==Yn==Zn, I just use a single const

  fltR = (float)R;
  fltG = (float)G;
  fltB = (float)B;
  X = 2.769*fltR + 1.7518*fltG + 1.13*fltB;
  Y = fltR + 4.5907*fltG + 0.0601*fltB;
  Z = 0.0565*fltG + 5.5943*fltB;
  
  (*L) = (float)(116. * cbrt(Y / XYZn) - 16.);
  (*a) = (float)(500. * (cbrt(X / XYZn) - cbrt(Y / XYZn)));
  (*b) = (float)(200. * (cbrt(Y / XYZn) - cbrt(Z / XYZn)));
}

///Returns the R,G,B values [0..255] corresponding to the L,a,b values
inline void DColorSpace::getRGBFromCIELab(float L, float a, float b,
					  D_uint8 *R, D_uint8 *G, D_uint8 *B){
  int iR, iG, iB;
  float Xreverse, Yreverse, Zreverse;
  float fltTmpX;
  float fltTmpY;
  float fltTmpZ;
  const float XYZn = 1440.954f; // since Xn==Yn==Zn, I just use a single const

  fltTmpY = (L+16.f)/116.f;
  Yreverse = XYZn * fltTmpY * fltTmpY * fltTmpY;
  fltTmpX = a / 500.f + cbrt(Yreverse / XYZn);
  Xreverse = XYZn * fltTmpX * fltTmpX * fltTmpX;
  fltTmpZ = cbrt(Yreverse/XYZn) - b / 200.f;
  Zreverse = XYZn * fltTmpZ * fltTmpZ * fltTmpZ;

  iR = (int)nearbyintf(0.4184383528f * Xreverse + -0.158655799f * Yreverse +
		       -0.0828164605f * Zreverse);
  iG = (int)nearbyintf(-0.0911611925f * Xreverse + 0.2524253419f * Yreverse +
		       0.0157019438f * Zreverse);
  iB = (int)nearbyintf(0.0009206884461f * Xreverse + -0.0025493863f * Yreverse
		       + 0.1785947912f * Zreverse);
  // cap the range between 0 and 255 (inclusive)
  if(iR < 0)
    iR = 0;
  if(iR > 255)
    iR = 255;

  if(iG < 0)
    iG = 0;
  if(iG > 255)
    iG = 255;

  if(iB < 0)
    iB = 0;
  if(iB > 255)
    iB = 255;

  (*R) = (D_uint8)iR;
  (*G) = (D_uint8)iG;
  (*B) = (D_uint8)iB;
}

///Computes the H,S,V values corresponding to the R,G,B values
/** h in [0,360), s and v in [0,1]. <b>Note: if s==0 then h is undefined (but this function sets it to 0 for convenience).</b>
 **/
inline void DColorSpace::getHSVFromRGB(D_uint8 R, D_uint8 G, D_uint8 B,
				       float *H, float *S, float *V){
  float max, min, delta;
  max = min = (float)R;
  if(G > max)
    max = G;
  else if(G < min)
    min = G;
  if(B > max)
    max = B;
  else if(B < min)
    min = B;
  
  (*V) = max / 255.;
  (*S) = (0.f != max) ? ((max-min) / max) : 0.f;
  if(0.f == (*S)){
    (*H) = 0.f; // hue is undefined when S==0, but set to 0 for convenience
  }
  else{
    delta = max - min;
    if(max == (float)R){
      (*H) = (G - B) / delta;
    }
    else if(max == (float)G){
      (*H) = 2.f + (B-R) / delta;
    }
    else{ // max is B
      (*H) = 4.f + (R-G) / delta;
    }
    (*H) *= 60.f;
    if((*H) < 0.f)
      (*H) += 360.f;
  }
}
///Computes the R,G,B values for H,S,V
/**0<=H<360, S and V should be 0 to 1 (inclusive).  <b>Note: If S is zero, then Hue is undefined (doesn't matter) and this is a grayscale value equal to V<b>.
 */
inline void DColorSpace::getRGBFromHSV(float H, float S, float V,
				       D_uint8 *R, D_uint8 *G, D_uint8 *B){
  if(0.f == S){ // gray value
    (*R) = (*G) = (*B) = (D_uint8)(255.f*V);
#ifdef DEBUG
//     if(H != 361.361f){
//       fprintf(stderr, "DColorSpace::getRGBFromHSV() bad hue value for gray\n");
//     }
#endif
  }
  else{
    float f,p,q,t;
    int i;
    if(360.f == H)
      H = 0.f;
    H /= 60.f;
    i = (int)H;
    f = H-i;
    V *= 255.f;
    p = V * (1.f-S);
    q = V * (1.f-(S*f));
    t = V * (1.f-(S*(1.f-f)));
    switch(i){
      case 0: (*R) = (D_uint8)V; (*G) = (D_uint8)t; (*B) = (D_uint8)p; break;
      case 1: (*R) = (D_uint8)q; (*G) = (D_uint8)V; (*B) = (D_uint8)p; break;
      case 2: (*R) = (D_uint8)p; (*G) = (D_uint8)V; (*B) = (D_uint8)t; break;
      case 3: (*R) = (D_uint8)p; (*G) = (D_uint8)q; (*B) = (D_uint8)V; break;
      case 4: (*R) = (D_uint8)t; (*G) = (D_uint8)p; (*B) = (D_uint8)V; break;
      case 5: (*R) = (D_uint8)V; (*G) = (D_uint8)p; (*B) = (D_uint8)q; break;
      default:
	fprintf(stderr,
		"DColorSpace::getRGBFromHSV() bad logic in switch %d\n",i);
	break;
    }
  }
}


#endif
