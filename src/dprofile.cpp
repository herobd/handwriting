#include "dprofile.h"
#include "dmath.h"
#include "dinstancecounter.h"
#include <string.h>

///default constructor
/**The constructor initializes the length to 0 and data pointer to NULL*/
DProfile::DProfile(){
  DInstanceCounter::addInstance("DProfile");
  rgProf = NULL;
  len = 0;
}


DProfile::DProfile(const DProfile &src){
  DInstanceCounter::addInstance("DProfile");
  rgProf = NULL;
  len = 0;
  if(NULL != src.rgProf){
    len = src.len;
    rgProf = (double*)malloc(sizeof(double)*len);
    D_CHECKPTR(rgProf);
    memcpy(rgProf, src.rgProf, len * sizeof(double));
  }
}

DProfile::~DProfile(){
  DInstanceCounter::removeInstance("DProfile");
  if(rgProf)
    free(rgProf);
  rgProf = NULL;
  len = 0;
}



///computes the profile of an image, projected onto an axis of arbitray angle
/**Computes the projection profile of img onto an axis with angle
axisAngDeg.  Horizontal and vertical profiles (when axisAngDeg is 0 or
90 degress, respectively) are special-cased for speed by calling
getImageVerticalProfile() or getImageHorizontalProfile().  I have seen
inconsistent usage of the terms "vertical profile" and "horizontal
profile" since some people describe the direction of projection
instead of the direction of the axis onto which the image is
projected.  It seems more common, however, to use the direction of the
axis.  Therefore, when I say "vertical profile," I mean that the
profile length is equal to the height of the image, and the Rth value
in the profile array is the sum of the values in the Rth row of the
image.  Likewise, a horizontal profile has length equal to the image
width, and each value is the projection of the corresponding image
column onto the horizontal axis.  For a horizontal profile, set
axisAngDeg to 0.  For a vertical profile, set axisAngDeg to 90.  Note
that angles increase clockwise since the y-coordinate of images
increases from top to bottom.  For RGB images and multi-channel float
or double images, the sum of all channels is used (a single profile is
calculated).  If this is not what is desired, you could create
separate images for each channel and take the profiles
seperately. Complex images are not supported directly, so they must be
converted into a different type before a profile can be taken.  If
fInterp is false (default), then the nearest image pixel value is used
when the position the profile passes through is between pixels.  If
fInterp is true, bilinear interpolation is used to estimate the value
that should go in the profile.*/
void DProfile::getImageProfile(const DImage &img, double axisAngDeg, 
			       bool fNormalize, bool fInterp){
  int w, h;
  double wm1, hm1;
  int numPixels;
  double dblSum;
  double theta;
  double dblTmp;
  int alen;//length of anchor segment.  The anchor segment is the
	   //segment that passes through the center of the image and
	   //is oriented perpendicular to axisAngDeg.  The length of
	   //the profile will be equal to alen.
  int olen;//length of offset segments.  The offset segments are the 
           //integration paths parallel to axisAngDeg (one segment per
           //profile value.
  double asx, asy;//start (x,y) of the anchor segment
  double adx, ady;//deltaX and deltaY of the anchor segment
  double ax, ay;  // current anchor point (x,y along the anchor segment)
  double osx, osy;//starty x,y offset from ax,ay for an offset segment
  double odx, ody;//deltaX and deltaY of offset segments
  double ox, oy;  // current offset point (x,y along the offset segment)

  if(DImage::DImage_cmplx == img.getImageType()){
    fprintf(stderr,
	    "DProfile::getImageProfile() does not support complex images\n");
    abort();
  }
//   if(((axisAngDeg > -0.00001) && (axisAngDeg < 0.00001)) ||
//      ((axisAngDeg > 179.00001) && (axisAngDeg < 180.00001)) ||
//      ((axisAngDeg > 359.00001) && (axisAngDeg < 360.00001))){
//     getImageHorizontalProfile(img, fNormalize);
//     return;
//   }
//   if(((axisAngDeg > 89.00001) && (axisAngDeg < 90.00001)) ||
//      ((axisAngDeg > 269.00001) && (axisAngDeg < 270.00001)) ||
//      ((axisAngDeg < -269.00001) && (axisAngDeg > -270.00001)) ||
//      ((axisAngDeg < -89.00001) && (axisAngDeg > -90.00001))){
//     getImageVerticalProfile(img, fNormalize);
//     return;
//   }
  w = img.width();
  h = img.height();
  wm1 = w-1;
  hm1 = h-1;
  theta = DMath::degreesToRadians(axisAngDeg);
  //figure out how long rgProf should be and how wide the integration should be
  odx = cos(theta);
  ody = sin(theta);
  adx = -1. * ody;
  ady = odx;
  olen = (int)(2*DMath::distPointLine(0, 0,  w/2., h/2.,  w/2.+adx, h/2.+ady));
  dblTmp = 2*DMath::distPointLine(w, 0,  w/2., h/2.,  w/2.+adx, h/2.+ady);
  if(dblTmp > olen)
    olen = (int)dblTmp;

  alen = (int)(2*DMath::distPointLine(0, 0,  w/2., h/2.,  w/2.+odx, h/2.+ody));
  dblTmp = 2*DMath::distPointLine(w, 0,  w/2., h/2.,  w/2.+odx, h/2.+ody);
  if(dblTmp > alen)
    alen = (int)dblTmp;

//   printf("w=%d h=%d\n", w, h);
//   printf("adx=%.2f ady=%.2f  alen=%d\n", adx, ady, alen);
//   printf("odx=%.2f ody=%.2f  olen=%d\n", odx, ody, olen);

  // allocate the rgProf array
  if(NULL == rgProf){
    rgProf = (double*)malloc(alen * sizeof(double));
    D_CHECKPTR(rgProf);
    len = alen;
  }
  else{
    if(len != alen){
      rgProf = (double*)realloc(rgProf,alen*sizeof(double));
      D_CHECKPTR(rgProf);
      len = alen;
    }
  }

  ax = asx = w/2.-(alen/2*adx);
  ay = asy = h/2.-(alen/2*ady);
  osx = -(olen/2*odx);
  osy = -(olen/2*ody);

  switch(img.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint8 *pu8;
	pu8=img.dataPointer_u8();
	if(fInterp){
	  fprintf(stderr, "WARNING! DProfile::getImageProfile() may be buggy "
		  "when fInterp is true\n");
	  for(int aa=0; aa < alen; ++aa){
	    int ix, iy;
	    double w1, w2; //weight of left vs right pixels
	    double  w3, w4;//weight of top vs bottom pixels
	    double left, right;
	    int oidx;
	    ox = ax+osx;
	    oy = ay+osy;
	    numPixels = 0;
	    dblSum = 0.;
	    for(int oo=0; oo < olen; ++oo){
	      if( (ox>=0) && (ox < wm1) && (oy>=0) && (oy<hm1)){
		ix = (int)ox;
		iy = (int)oy;
		oidx = iy*w + ix;
		w2 = ox-ix;
		w1 = 1.-w2;
		w4 = oy-iy;
		w3 = 1.-w4;
		left = (pu8[oidx] * w3 + pu8[oidx+w] * w4);
		right = (pu8[oidx+1] * w3 + pu8[oidx+w+1] * w4);
		dblSum += left*w1 + right*w2;
		++numPixels;
	      }
	      else if((ox==wm1) || (oy == hm1)){
		// don't interpolate on the edges, just approximate
		ix = (int)ox;
		iy = (int)oy;
		oidx = iy*w + ix;
		dblSum += pu8[oidx];
		++numPixels;
	      }
	      ox += odx;
	      oy += ody;
	    }
	    ax += adx;
	    ay += ady;
	    if(fNormalize && (numPixels>0))
	      rgProf[aa] = dblSum / numPixels;
	    else
	      rgProf[aa] = dblSum;
	  }
	}
	else{
	  for(int aa=0; aa < alen; ++aa){
	    int ix, iy;
	    int oidx;
	    ox = ax+osx;
	    oy = ay+osy;
	    numPixels = 0;
	    dblSum = 0.;
	    for(int oo=0; oo < olen; ++oo){
	      if( (ox>=0) && (ox <w) && (oy>=0) && (oy<h)){
		ix = (int)ox;
		iy = (int)oy;
		oidx = iy*w + ix;
		dblSum += pu8[oidx];
		++numPixels;
	      }
	      ox += odx;
	      oy += ody;
	    }
	    ax += adx;
	    ay += ady;
	    if(fNormalize && (numPixels>0))
	      rgProf[aa] = dblSum / numPixels;
	    else
	      rgProf[aa] = dblSum;
	  }
	}
      }
      break;
    default:
      fprintf(stderr, "Not yet implemented!\n");
      abort();
  }//end switch(img.getImageType())

}

///Set this to the vertical profile of img (this is a grayscale profile)
/**Each row of img is projected onto the vertical axis.  Resulting
   data length will be equal to the height of img.  The profile is a summation
   of the grayscale values in each row.  If fNormalize is true, then each value
   is divided by img.width() so it is the average grayscale value for the row
   instead of the sum.

   If fNormalize is true, the resulting profile values are divided by the image
   width.
 */
void DProfile::getImageVerticalProfile(const DImage &img, bool fNormalize){
  int w, h;

  w = img.width();
  h = img.height();
  // allocate the rgProf array
  if(NULL == rgProf){
    rgProf = (double*)malloc(h * sizeof(double));
    D_CHECKPTR(rgProf);
    len = h;
  }
  else{
    if(len != h){
      rgProf = (double*)realloc(rgProf,h*sizeof(double));
      D_CHECKPTR(rgProf);
      len = h;
    }
  }
  switch(img.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint8 *pu8;
	pu8=img.dataPointer_u8();
	for(int y = 0, idx=0; y < h; ++y){
	  rgProf[y] = 0.;
	  for(int x = 0; x < w; ++x, ++idx){
	    rgProf[y] += pu8[idx];
	  }
	  if(fNormalize)
	    rgProf[y] /= w;
	}
      }
      break;
    case DImage::DImage_flt_multi:
      {
	float *pflt;
	if(img.numChannels() > 1){
	  fprintf(stderr,"DProfile::getImageVerticalProfile() floats only "
		  "supported with a single channel\n");
	  abort();
	}
	pflt=img.dataPointer_flt(0);
	for(int y = 0, idx=0; y < h; ++y){
	  rgProf[y] = 0.;
	  for(int x = 0; x < w; ++x, ++idx){
	    rgProf[y] += pflt[idx];
	  }
	  if(fNormalize)
	    rgProf[y] /= w;
	}
      }
      break;
    default:
      fprintf(stderr, "Not yet implemented!\n");
      abort();
  }//end switch(img.getImageType())
}

///Set this to the horizontal profile of img
/**Each column of img is projected onto the horizontal axis.
   Resulting data length will be equal to the width of img.

   If fNormalize is true, the resulting profile values are divided by the image
   height.
 */
void DProfile::getImageHorizontalProfile(const DImage &img, bool fNormalize){
    int w, h;

  w = img.width();
  h = img.height();
  // allocate the rgProf array
  if(NULL == rgProf){
    rgProf = (double*)malloc(w * sizeof(double));
    D_CHECKPTR(rgProf);
    len = w;
  }
  else{
    if(len != w){
      rgProf = (double*)realloc(rgProf,w*sizeof(double));
      D_CHECKPTR(rgProf);
      len = w;
    }
  }
  memset(rgProf, 0, sizeof(double) * w);
  switch(img.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint8 *pu8;
	pu8=img.dataPointer_u8();
	for(int y = 0, idx=0; y < h; ++y){
	  for(int x = 0; x < w; ++x, ++idx){
	    rgProf[x] += pu8[idx];
	  }
	}
	if(fNormalize){
	  for(int x = 0; x < w; ++x)
	    rgProf[x] /= h;
	}
      }
      break;
    default:
      fprintf(stderr, "Not yet implemented!\n");
      abort();
  }//end switch(img.getImageType())
	
	
}

///Take an angled vertical profile down the center of the image
/**This function calculates the vertical profile (the length of
 * profile will be the same as the height of the image).  However,
 * instead of projecting pixels straight across to the vertical axis,
 * the projection is taken using lines angled through the middle
 * (x=width/2) of the image.  Linear interpolation is used for
 * sampling the image pixel values since the y-position on the
 * projection lines will normally be "between" two pixels for any
 * given x-position.  The function is intended only for angles between
 * -45 and 45 degrees, since otherwise the slope will be too steep for
 * the assumption I am making (I am steping through the x-values and
 * calculating the y-values.  If the slope is steeper, I should do the
 * opposite.  As a "to-do," maybe we should just check the angle and
 * handle the two cases individually). */
void DProfile::getAngledVertProfile(const DImage &img, double ang,
				    int fNormalize){
  double m;
  int x;
  double y, val;
  int hm1;
  double xc;
  int numPixels;
  D_uint8 *pimg;
  int w, h;
  int initialOffset=0;
  int *rgYoffsets;
  double *rgBotWeights;
  int yTop;
  int yTopPrev=0;
  int imglen;

  if(DImage::DImage_u8 != img.getImageType()){
    fprintf(stderr,
	    "DProfile::getAngledVertProfile() currently only supports 8-bit "
	    "grayscale images\n");
    abort();
  }
  w=img.width();
  h=img.height();
  pimg = img.dataPointer_u8();

  if(NULL == rgProf){
    rgProf = (double*)malloc(h * sizeof(double));
    D_CHECKPTR(rgProf);
    len = h;
  }
  else{
    if(len != h){
      rgProf = (double*)realloc(rgProf,h*sizeof(double));
      D_CHECKPTR(rgProf);
      len = h;
    }
  }
//   memset(rgProf, 0, sizeof(double) * h);

  rgYoffsets = (int*)malloc(sizeof(int)*w);
  D_CHECKPTR(rgYoffsets);
  rgBotWeights = (double*)malloc(sizeof(double)*w);
  D_CHECKPTR(rgBotWeights);

  m = 1 * tan(DMath::degreesToRadians(ang)); /* dy per dx=1 */
  xc = w / 2.;

  // initialOffset is the y-offset of first pixel, rgYoffsets[i] for i=1...w-1
  // are the relative offsets from rgYoffsets[i-1] of the top pixel
  // rgBotWeights[i] is the interpolation weight of the bottom pixel, while
  // 1.-rgBotWeights[i] is the interpolation weight of the top pixel.
  for(x = 0; x < w; ++x){
    y = m * ((double)x - xc);
    yTop = (int)(floor(y));
    rgBotWeights[x] = y - (double)(yTop);
    if(0 == x){
      rgYoffsets[0] = 0;
      initialOffset = yTop;
    }
    else
      rgYoffsets[x] = (yTop - yTopPrev);
    yTopPrev = yTop;
  }

  hm1 = h - 1;
  imglen = w * hm1;

  for(int i = 0; i < h; ++i){ /* profile index */
    int idx;
    rgProf[i] = 0.;
    numPixels = 0;
    idx = w * (i + initialOffset);
    for(x = 0; x < w; ++x, ++idx){
      idx += (w*rgYoffsets[x]);// go to next row if needed
      if((idx <0) || (idx >= imglen)) // out of bounds
	continue;
      ++numPixels;
      val = rgBotWeights[x] * pimg[idx] + (1.-rgBotWeights[x]) * pimg[idx+w];
      rgProf[i] += val;
    }
    if((numPixels > 0) && fNormalize)
      rgProf[i] /= numPixels;
  }

  free(rgYoffsets);
  free(rgBotWeights);
}

///save the profile data in a file named stPath in gnuplot format
void DProfile::saveGnuplot(const char *stPath){
  FILE *fout;

  fout = fopen(stPath, "wb");
  if(!fout){
    fprintf(stderr, "DProfile::saveGnuplot() can't write file '%s'\n", stPath);
    abort();
  }
  for(int i = 0; i < len; ++i){
    fprintf(fout,"%d %.24f\n", i, rgProf[i]);
  }
  fclose(fout);
}

///smooth the profile using a rectangular convolution kernel
/**The profile is treated as if it were padded on each end with
   zeros. */
void DProfile::smoothAvg(int windowRadius){
  double sum;
  double *rgTmp;
  if((NULL == rgProf) || (0 == len)){
    fprintf(stderr, "DProfile::smoothAvg() called on uninitialized profile\n");
    abort();
  }

  rgTmp = (double*)malloc(sizeof(double) * len);
  D_CHECKPTR(rgTmp);
  // I copy the source and copy result back (instead of just creating
  // a new result and swapping the pointer) so I don't invalidate the
  // pointer, just in case a user of this function has grabbed the
  // pointer using dataPointer()
  memcpy(rgTmp, rgProf, sizeof(double)*len);

  for(int i=0; i < len; ++i){
    sum = 0.;
    for(int j=-windowRadius; j <=windowRadius; ++j){
      int ii = i+j;
      if((ii >= 0) && (ii < len))
	sum += rgTmp[ii];
    }
    rgProf[i] = sum / (windowRadius * 2 + 1);
  }
}


///convert the profile into an image
/**The resulting image will be 8-bit grayscale and if fVertical is
 * true, then the image will be numPixels wide and dataLen() high (if
 * fVertical is false, then the image will be dataLen() wide and
 * numPixels high).  The grayscale values for the profile foreground
 * and background can be specified with fg and bg, respectively. (fg
 * specifies the pixels that are part of the profile, bg specifies the
 * pixels that are NOT part of the profile)
 */
DImage DProfile::toDImage(int numPixels, bool fVertical,
			  D_uint8 fg, D_uint8 bg){
  DImage img;
  double max;
  D_uint8 *pdata;

  img.create(numPixels, len, DImage::DImage_u8);
  pdata = img.dataPointer_u8();
  max = 0.;
  for(int i = 0; i < len; ++i){
    if(rgProf[i] > max)
      max = rgProf[i];
  }
  // printf("DProfile::toDImage() max=%f\n",max);
  for(int y = 0, idx = 0; y < len; ++y){
    double fillTo;
    fillTo = rgProf[y] * numPixels / max;
    for(int x = 0; x < numPixels; ++x, ++idx){
      if(x <= fillTo)
	pdata[idx] = fg;
      else
	pdata[idx] = bg;
    }
  }
  if(!fVertical){
    DImage imgRot;
    img.rotate90_(imgRot, -90);
    return imgRot;
  }
  return img;
}

///return the maximum data value in the profile
double DProfile::max(){
  double max = 0.;
  for(int i = 0; i < len; ++i){
    if(rgProf[i] > max)
      max = rgProf[i];
  }
  return max;
}
///return the minimum data value in the profile
double DProfile::min(){
  double min = 0.;
  if(len>0)
    min = rgProf[0];
  for(int i = 0; i < len; ++i){
    if(rgProf[i] < min)
      min = rgProf[i];
  }
  return min;
}


///Profile of max runlength of some pixel value (8-bit RGB or grayscale)
/**Max runlength in each column of img is projected onto the horizontal axis.
   Resulting data length will be equal to the width of img.
   If fNormalize is true, each profile value will be divided by image height,
   so the value is a fraction of the image height instead of a number of pixels.
 */
void DProfile::getHorizMaxRunlengthProfile(const DImage &img, D_uint32 rgbVal,
					  bool fNormalize){
  int w, h;
  unsigned int *rgRunlengths;

  w = img.width();
  h = img.height();
  // allocate the rgProf array
  if(NULL == rgProf){
    rgProf = (double*)malloc(w * sizeof(double));
    D_CHECKPTR(rgProf);
    len = w;
  }
  else{
    if(len != w){
      rgProf = (double*)realloc(rgProf,w*sizeof(double));
      D_CHECKPTR(rgProf);
      len = w;
    }
  }
  rgRunlengths = (unsigned int*)malloc(sizeof(unsigned int)*w);
  D_CHECKPTR(rgRunlengths);
  memset(rgRunlengths, 0, sizeof(unsigned int)*w);
  memset(rgProf, 0, sizeof(double) * w);
  switch(img.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint8 *pu8;
	pu8=img.dataPointer_u8();
	for(int y = 0, idx=0; y < h; ++y){
	  for(int x = 0; x < w; ++x, ++idx){
	    if((D_uint8)rgbVal == pu8[idx]){//increment run length for this col
	      ++(rgRunlengths[x]);
	      if(rgRunlengths[x] > rgProf[x])
		rgProf[x] = (double)rgRunlengths[x];
	    }
	    else{
	      rgRunlengths[x] = 0;
	    }
	  }
	}
	if(fNormalize){
	  for(int x = 0; x < w; ++x)
	    rgProf[x] /= h;
	}
      }
      break;
    default:
      fprintf(stderr, "Not yet implemented!\n");
      abort();
  }//end switch(img.getImageType())
	
  free(rgRunlengths);
}


///Profile of avg runlength of some pixel value (8-bit RGB or grayscale)
/**Avg runlength in each column of img is projected onto the horizontal axis.
   Resulting data length will be equal to the width of img.
   If fNormalize is true, each profile value will be divided by image height,
   so the value is a fraction of the image height instead of a number of pixels.
 */
void DProfile::getHorizAvgRunlengthProfile(const DImage &img, D_uint32 rgbVal,
					   bool fNormalize){
  int w, h;
  unsigned int *rgRunlengths;
  unsigned int *rgNumRuns;

  w = img.width();
  h = img.height();
  // allocate the rgProf array
  if(NULL == rgProf){
    rgProf = (double*)malloc(w * sizeof(double));
    D_CHECKPTR(rgProf);
    len = w;
  }
  else{
    if(len != w){
      rgProf = (double*)realloc(rgProf,w*sizeof(double));
      D_CHECKPTR(rgProf);
      len = w;
    }
  }
  rgRunlengths = (unsigned int*)malloc(sizeof(unsigned int)*w);
  D_CHECKPTR(rgRunlengths);
  rgNumRuns = (unsigned int*)malloc(sizeof(unsigned int)*w);
  D_CHECKPTR(rgNumRuns);
  memset(rgRunlengths, 0, sizeof(unsigned int)*w);
  memset(rgNumRuns, 0, sizeof(unsigned int)*w);
  memset(rgProf, 0, sizeof(double) * w);
  switch(img.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint8 *pu8;
	pu8=img.dataPointer_u8();
	for(int y = 0, idx=0; y < h; ++y){
	  for(int x = 0; x < w; ++x, ++idx){
	    if((D_uint8)rgbVal == pu8[idx]){//increment run length for this col
	      if(0 == rgRunlengths[x])
		++(rgNumRuns[x]);
	      ++(rgRunlengths[x]);
	      ++(rgProf[x]);
	    }
	    else{
	      rgRunlengths[x] = 0;
	    }
	  }
	}
	for(int x = 0; x < w; ++x){
	  //divide sum by number of runs to get avg
	  if(rgNumRuns[x] > 0)
	    rgProf[x] /= rgNumRuns[x];
	}
	if(fNormalize){
	  for(int x = 0; x < w; ++x){
	    rgProf[x] /= h;
	  }
	}
      }
      break;
    default:
      fprintf(stderr, "Not yet implemented!\n");
      abort();
  }//end switch(img.getImageType())
	
  free(rgRunlengths);
  free(rgNumRuns);
}




///Profile of max runlength of some pixel value (8-bit RGB or grayscale)
/**Max runlength in each column of img is projected onto the vertical axis.
   Resulting data length will be equal to the height of img.
   If fNormalize is true, each profile value will be divided by image width,
   so the value is a fraction of the image width instead of a number of pixels.
 */
void DProfile::getVertMaxRunlengthProfile(const DImage &img, D_uint32 rgbVal,
					  bool fNormalize){
  int w, h;
  unsigned int runLength;

  w = img.width();
  h = img.height();
  // allocate the rgProf array
  if(NULL == rgProf){
    rgProf = (double*)malloc(h * sizeof(double));
    D_CHECKPTR(rgProf);
    len = h;
  }
  else{
    if(len != h){
      rgProf = (double*)realloc(rgProf,h*sizeof(double));
      D_CHECKPTR(rgProf);
      len = h;
    }
  }
  switch(img.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint8 *pu8;
	pu8=img.dataPointer_u8();
	for(int y = 0, idx=0; y < h; ++y){
	  rgProf[y] = 0.;
	  runLength = 0;
	  for(int x = 0; x < w; ++x, ++idx){
	    if((D_uint8)rgbVal == pu8[idx]){//increment run length for this row
	      ++runLength;
	      if(runLength > rgProf[y])
		rgProf[y] = runLength;
	    }
	    else{
	      runLength = 0;
	    }
	  }
	  if(fNormalize)
	    rgProf[y] /= w;
	}
      }
      break;
    default:
      fprintf(stderr, "Not yet implemented!\n");
      abort();
  }//end switch(img.getImageType())
}

///Profile of avg runlength of some pixel value (8-bit RGB or grayscale)
/**Avg runlength in each column of img is projected onto the vertical axis.
   Resulting data length will be equal to the height of img.
   If fNormalize is true, each profile value will be divided by image width,
   so the value is a fraction of the image width instead of a number of pixels.
 */
void DProfile::getVertAvgRunlengthProfile(const DImage &img, D_uint32 rgbVal,
					  bool fNormalize){
  int w, h;
  unsigned int runLength, numRuns;

  w = img.width();
  h = img.height();
  // allocate the rgProf array
  if(NULL == rgProf){
    rgProf = (double*)malloc(h * sizeof(double));
    D_CHECKPTR(rgProf);
    len = h;
  }
  else{
    if(len != h){
      rgProf = (double*)realloc(rgProf,h*sizeof(double));
      D_CHECKPTR(rgProf);
      len = h;
    }
  }
  switch(img.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint8 *pu8;
	pu8=img.dataPointer_u8();
	for(int y = 0, idx=0; y < h; ++y){
	  rgProf[y] = 0.;
	  runLength = 0;
	  numRuns = 0;
	  for(int x = 0; x < w; ++x, ++idx){
	    if((D_uint8)rgbVal == pu8[idx]){//increment run length for this row
	      if(0==runLength)
		++numRuns;
	      ++runLength;
	      ++(rgProf[y]);
	    }
	    else{
	      runLength = 0;
	    }
	  }
	  if(numRuns > 0)
	    rgProf[y] /= numRuns; //(we have sum and need to divide for avg)
	  if(fNormalize)
	    rgProf[y] /= w;
	}
      }
      break;
    default:
      fprintf(stderr, "Not yet implemented!\n");
      abort();
  }//end switch(img.getImageType())
}

void DProfile::copyFromData(double *rgData, int dataLen){
  if(rgProf)
    free(rgProf);
  rgProf = (double*)malloc(sizeof(double)*dataLen);
  D_CHECKPTR(rgData);
  len = dataLen;
  for(int i=0; i < dataLen; ++i)
    rgProf[i] = rgData[i];
}

void DProfile::copyFromData(float *rgData, int dataLen){
  if(rgProf)
    free(rgProf);
  rgProf = (double*)malloc(sizeof(double)*dataLen);
  D_CHECKPTR(rgData);
  len = dataLen;
  for(int i=0; i < dataLen; ++i)
    rgProf[i] = (double)rgData[i];
}

void DProfile::copyFromData(int *rgData, int dataLen){
  if(rgProf)
    free(rgProf);
  rgProf = (double*)malloc(sizeof(double)*dataLen);
  D_CHECKPTR(rgData);
  len = dataLen;
  for(int i=0; i < dataLen; ++i)
    rgProf[i] = (double)rgData[i];
}
