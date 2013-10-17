#include <stdio.h>
#include "dedgedetector.h"
#include "dimage.h"
#include "dkernel2d.h"
#include "dconvolver.h"


/**TODO: I think there may be a bug in the gradient direction code.
   (sobel_) and may also need to change other code that depends on it.
   When testing the directions for nonMaximalSuppression, I found that
   the directions aren't quite what I thought they were.  I am just
   going to change the nonMaximalSuppression code for now to do what I
   want it to until I have time to verify the directions.

*/


///Similar to sobel_() but for RGB images (result images have 3 channels)
/** The gradient magnitude and direction of each channel (R,G,B) is
    calculated separately and stored in the corresponding channel of the 
    result image.  This function just splits the R,G,B image into 3 channels,
    calls sobel_() on each channel, and combines the results back into an
    R,G,B image for the magnitude and 3-channel DImage_dbl_multi for imgDir.
**/
void DEdgeDetector::sobel_RGB_(DImage &imgGradMag, const DImage &imgSrc,
			       DImage *imgDir){
  DImage imgR, imgG, imgB;
  DImage magR, magG, magB;
  DImage dirR, dirG, dirB;
  DImage *pdirR, *pdirG, *pdirB;

  if(DImage::DImage_RGB != imgSrc.getImageType()){
    fprintf(stderr, "DEdgeDetector::sobel_RGB_() only handles 8-bit RGB\n");
    abort();
  }

  imgSrc.splitRGB(imgR, imgG, imgB);
  if(NULL == imgDir){
    pdirR = pdirG = pdirB = NULL;
  }
  else{
    pdirR = &dirR;
    pdirG = &dirG;
    pdirB = &dirB;
  }
  DEdgeDetector::sobel_(magR, imgR, pdirR);
  DEdgeDetector::sobel_(magG, imgG, pdirG);
  DEdgeDetector::sobel_(magB, imgB, pdirB);

  imgGradMag.combineRGB(magR, magG, magB);
  if(NULL != imgDir){
    imgDir->create(imgSrc.width(), imgSrc.height(), DImage::DImage_dbl_multi,3,
		  imgSrc.getAllocMethod());
    int w = imgSrc.width();
    int h = imgSrc.height();
    double *pdblDir, *pdbl;
    pdblDir = imgDir->dataPointer_dbl(0);
    pdbl = pdirR->dataPointer_dbl(0);
    for(int y = 0; y < h; ++y){
      for(int x = 0; x < w; ++x){
	(*pdblDir) = (*pdbl);
	++pdblDir;
	++pdbl;
      }
    }
    pdblDir = imgDir->dataPointer_dbl(1);
    pdbl = pdirG->dataPointer_dbl(0);
    for(int y = 0; y < h; ++y){
      for(int x = 0; x < w; ++x){
	(*pdblDir) = (*pdbl);
	++pdblDir;
	++pdbl;
      }
    }
    pdblDir = imgDir->dataPointer_dbl(2);
    pdbl = pdirB->dataPointer_dbl(0);
    for(int y = 0; y < h; ++y){
      for(int x = 0; x < w; ++x){
	(*pdblDir) = (*pdbl);
	++pdblDir;
	++pdbl;
      }
    }
  }
}


///Similar to sobel_() but for RGB images (result image has SINGLE channel)
/** The gradient magnitude and direction of each channel (R,G,B) is
    calculated separately and then COMBINED into a single-channel
    8-bit result.  This function splits the R,G,B image into 3
    channels, calls sobel_() on each channel, and then combines the
    results by using the magnitude of vector of the result values
    (i.e., if rr,rg,rb are the result magnitudes for the r,g,b
    channels at some pixel, then imgGradMag will have the value
    255./441.67295593 * SQRT(rr^2 + rg^2 + rb^2) at that pixel.
    (the max value if rr==rg==rb==255 is 441.67295593)
    The gradient direction
    will be based on grayscale values (however, it is not yet
    implemented, so if you pass in anything other than NULL for imgDir,
    abort() will be called.

**/
void DEdgeDetector::sobel_RGB_combined_(DImage &imgGradMag,
					const DImage &imgSrc, DImage *imgDir){
  DImage imgR, imgG, imgB;
  DImage magR, magG, magB;
  DImage dirR, dirG, dirB;
  D_uint8* p8, *pR, *pG, *pB;
  const double dblSCALE = 255. / sqrt(3*255.*255.);

  if(DImage::DImage_RGB != imgSrc.getImageType()){
    fprintf(stderr, "DEdgeDetector::sobel_RGB_combined_() only handles "
	    "8-bit RGB images\n");
    abort();
  }
  if(NULL != imgDir){
    fprintf(stderr, "DEdgeDetector::sobel_RGB_combined_() does not yet support"
	    "gradient/edge directions. NULL must be used until implemented.\n");
    abort();
  }

  imgSrc.splitRGB(imgR, imgG, imgB);
  DEdgeDetector::sobel_(magR, imgR, NULL);
  DEdgeDetector::sobel_(magG, imgG, NULL);
  DEdgeDetector::sobel_(magB, imgB, NULL);

  int w,h;
  w = imgSrc.width();
  h = imgSrc.height();
  imgGradMag.create(w, h, DImage::DImage_u8, 3, imgSrc.getAllocMethod());
  p8 = imgGradMag.dataPointer_u8();
  pR = magR.dataPointer_u8();
  pG = magG.dataPointer_u8();
  pB = magB.dataPointer_u8();
  for(int y=0; y < h; ++y){
    for(int x = 0; x < w; ++x, ++p8, ++pR, ++pG, ++pB){
      double dbl;
      dbl = dblSCALE * sqrt(((*pR)*(*pR)) + ((*pG)*(*pG)) + ((*pB)*(*pB)));
      if(dbl > 255.)
	dbl = 255.;
      (*p8) = (D_uint8)dbl;
    }
  }

}



///Compute gradient magnitude (and direction) using Sobel kernels
/** The gradient magnitude is defined as sqrt(dx^2 + dy^2), where dx
    and dy are the gradient in the x and y directions using the sobel
    kernels, and the gradient direction is arctan(dy/dx), which is
    perpendicular to the edge.  The gradient direction will only be
    calculated if imgDir is not NULL.  The resulting imgDir will be of
    type DImage_dbl_multi with a single channel, and the direction
    angles are in radians (not degrees).  For 8-bit grayscale images,
    the resulting imgGradMag will also be 8-bit grayscale, clipped at
    255.  This function currently only handles 8-bit grayscale images
    It should be faster for 8-bit grayscale than using the DConvolver
    class with the appropriate kernels and using a second pass to
    figure out gradient direction.  The Sobel kernels are:
\code
            -1  0  1                     -1 -2 -1
        dx= -2  0  2       and      dy=   0  0  0
            -1  0  1                      1  2  1
\endcode
*/
void DEdgeDetector::sobel_(DImage &imgGradMag, const DImage &imgSrc,
			   DImage *imgDir){
  int w, h, wm1, hm1, wp1;
  int idx;
  int dx, dy;
  D_uint8 *pData;
  D_uint8 *pDataMag;
  double *pDataDir = NULL;
  double mag;
  if((&imgGradMag) == (&imgSrc)){
    fprintf(stderr, "DEdgeDetector::sobel_() imgGradMag==imgSrc!\n");
    abort();
  }
  if(DImage::DImage_u8 != imgSrc.getImageType()){
    fprintf(stderr, "DEdgeDetector::sobel_() only handles 8-bit grayscale\n");
    abort();
  }
  w = imgSrc.width();
  h = imgSrc.height();
  wm1 = w-1;
  wp1 = w+1;
  hm1 = h-1;

  imgGradMag.create(w, h, DImage::DImage_u8);
  pData = imgSrc.dataPointer_u8();
  pDataMag = imgGradMag.dataPointer_u8();
  if(imgDir != NULL){
    imgDir->create(w, h, DImage::DImage_dbl_multi, 1);
    pDataDir = imgDir->dataPointer_dbl();
  }
  // handle first row:
  //   first column of first row
  dx = ((int)((3*((unsigned int)pData[1])) +
	      ((unsigned int)pData[wp1]))) -
    ((int)((3*((unsigned int)pData[0])) +
	   ((unsigned int)pData[w])));
  dy = ((int)((3*((unsigned int)pData[w])) +
	      ((unsigned int)pData[wp1]))) -
    ((int)((3*((unsigned int)pData[0])) +
	   ((unsigned int)pData[1])));
  mag = sqrt(dx*dx+dy*dy);
  pDataMag[0] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  if(imgDir != NULL){
    pDataDir[0] = atan2(dy,dx);
  }
    
  //   main body of first row
  for(int x=1; x < wm1; ++x){
    dx = ((int)((3*((unsigned int)pData[x+1])) +
		((unsigned int)pData[x+wp1]))) -
      ((int)((3*((unsigned int)pData[x-1])) +
	     ((unsigned int)pData[x+wm1])));
    dy = ((int)(((unsigned int)pData[x+wm1]) +
		(2*((unsigned int)pData[x+w])) +
		((unsigned int)pData[x+wp1]))) -
      ((int)(((unsigned int)pData[x-1]) +
	     (2*((unsigned int)pData[x])) +
	     ((unsigned int)pData[x+1])));
    mag = sqrt(dx*dx+dy*dy);
    pDataMag[x] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    if(imgDir != NULL){
      pDataDir[x] = atan2(dy,dx);
    }
  }
  //   last column of first row
  dx = ((int)((3*((unsigned int)pData[wm1])) +
	      ((unsigned int)pData[wm1+w]))) -
    ((int)((3*((unsigned int)pData[wm1-1])) +
	   ((unsigned int)pData[wm1+wm1])));
  dy = ((int)(((unsigned int)pData[wm1+wm1]) +
	      (3*((unsigned int)pData[wm1+w])))) -
    ((int)(((unsigned int)pData[wm1-1]) +
	   (3*((unsigned int)pData[wm1]))));
  mag = sqrt(dx*dx+dy*dy);
  pDataMag[wm1] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  if(imgDir != NULL){
    pDataDir[wm1] = atan2(dy,dx);
  }

  idx = w;
  // handle main body
  for(int y=1; y < hm1; ++y){
    int idxNW, idxSW;
    //   first col of current row
    dx = ((int)(((unsigned int)pData[idx-wm1]) +
		(2*((unsigned int)pData[idx+1])) +
		((unsigned int)pData[idx+wp1]))) -
      ((int)(((unsigned int)pData[idx-w]) +
	     (2*((unsigned int)pData[idx])) +
	     ((unsigned int)pData[idx+w])));
    dy = ((int)((3*((unsigned int)pData[idx+w])) +
		((unsigned int)pData[idx+wp1]))) -
      ((int)((3*((unsigned int)pData[idx-w])) +
	     ((unsigned int)pData[idx-wm1])));
    mag = sqrt(dx*dx+dy*dy);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    if(imgDir != NULL){
      pDataDir[idx] = atan2(dy,dx);
    }
    ++idx;
    //   main body of current row
    for(int x=1; x < wm1; ++x){
      idxNW = idx-wp1;
      idxSW = idx+wm1;
      dx = ((int)(((unsigned int)pData[idxNW+2]) +
		  (2*((unsigned int)pData[idx+1])) +
		  ((unsigned int)pData[idxSW+2]))) -
	((int)(((unsigned int)pData[idxNW]) +
	       (2*((unsigned int)pData[idx-1])) +
	       ((unsigned int)pData[idxSW])));
      dy = ((int)(((unsigned int)pData[idxSW]) +
		  (2*((unsigned int)pData[idxSW+1])) +
		  ((unsigned int)pData[idxSW+2]))) -
	((int)(((unsigned int)pData[idxNW]) +
	       (2*((unsigned int)pData[idxNW+1])) +
	       ((unsigned int)pData[idxNW+2])));
      mag = sqrt(dx*dx+dy*dy);
      pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
      if(imgDir != NULL){
	pDataDir[idx] = atan2(dy,dx);
      }
      ++idx;
    }
    //   last col of current row
    idxNW = idx-wp1;
    idxSW = idx+wm1;
    dx = ((int)(((unsigned int)pData[idxNW+1]) +
		(2*((unsigned int)pData[idx])) +
		((unsigned int)pData[idxSW+1]))) -
      ((int)(((unsigned int)pData[idxNW]) +
	     (2*((unsigned int)pData[idx-1])) +
	     ((unsigned int)pData[idxSW])));
    dy = ((int)(((unsigned int)pData[idxSW]) +
		(3*((unsigned int)pData[idxSW+1])))) -
      ((int)(((unsigned int)pData[idxNW]) +
	     (3*((unsigned int)pData[idxNW+1]))));
    mag = sqrt(dx*dx+dy*dy);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    if(imgDir != NULL){
      pDataDir[idx] = atan2(dy,dx);
    }
    ++idx;
  }

  // handle last row
  //   first col last row
  dx = ((int)(((unsigned int)pData[idx-wm1]) +
	      (3*((unsigned int)pData[idx+1])))) -
    ((int)(((unsigned int)pData[idx-w]) +
	   (3*((unsigned int)pData[idx]))));
  dy = ((int)((3*((unsigned int)pData[idx])) +
	      ((unsigned int)pData[idx+1]))) -
    ((int)((3*((unsigned int)pData[idx-w])) +
	   ((unsigned int)pData[idx-wm1])));
  mag = sqrt(dx*dx+dy*dy);
  pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  if(imgDir != NULL){
    pDataDir[idx] = atan2(dy,dx);
  }
  ++idx;
  
  //   main body of last row
  for(int x=1; x < wm1; ++x){
    dx = ((int)((3*((unsigned int)pData[idx+1])) +
		((unsigned int)pData[idx-wm1]))) -
      ((int)((3*((unsigned int)pData[idx-1])) +
	     ((unsigned int)pData[idx-wp1])));
    dy = ((int)(((unsigned int)pData[idx-1]) +
		(2*((unsigned int)pData[idx])) +
		((unsigned int)pData[idx+1]))) -
      ((int)(((unsigned int)pData[idx-wp1]) +
	     (2*((unsigned int)pData[idx-w])) +
	     ((unsigned int)pData[idx-wm1])));
    mag = sqrt(dx*dx+dy*dy);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    if(imgDir != NULL){
      pDataDir[idx] = atan2(dy,dx);
    }
    ++idx;
  }
  //   last col of last row
  dx = ((int)(((unsigned int)pData[idx-w]) +
	      (3*((unsigned int)pData[idx])))) -
    ((int)(((unsigned int)pData[idx-wp1]) +
	   (3*((unsigned int)pData[idx-1]))));
  dy = ((int)((3*((unsigned int)pData[idx])) +
	      ((unsigned int)pData[idx-1]))) -
    ((int)((3*((unsigned int)pData[idx-w])) +
	   ((unsigned int)pData[idx-wp1])));
  mag = sqrt(dx*dx+dy*dy);
  pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  if(imgDir != NULL){
    pDataDir[idx] = atan2(dy,dx);
  }
  ++idx;
  
#if 0
  // this is to make sure it worked as expected
 {
   DImage imgDx;
   DImage imgDy;
   DImage imgMag2;
   DKernel2D kernDx;
   DKernel2D kernDy;
   int count = 0;
   float *pdataDx;
   float *pdataDy;
   D_uint8 *pdata2;
   float rgfltKernDx[] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
   float rgfltKernDy[] = {-1, -2, -1,   0, 0,0,    1, 2, 1};
   
   kernDx.setData_flt(rgfltKernDx, 3, 3, false);
   kernDy.setData_flt(rgfltKernDy, 3, 3, false);

   printf("kernDx:\n");
   kernDx.print();
   printf("\nkernDy:\n");
   kernDy.print();
   imgMag2.create(w,h,DImage::DImage_u8);
   pdata2 = imgMag2.dataPointer_u8();
   DConvolver::convolve_(imgDx, imgSrc, kernDx, false, false, false);
   DConvolver::convolve_(imgDy, imgSrc, kernDy, false, false, false);

   pdataDx = imgDx.dataPointer_flt();
   pdataDy = imgDy.dataPointer_flt();
   
   for(int y = 0, idx3 = 0; y < h; ++y){
     for(int x = 0; x < w; ++x, ++idx3){
       mag = 0.5+
	 sqrt( pdataDx[idx3]*pdataDx[idx3] + pdataDy[idx3]*pdataDy[idx3]);
       if(mag > 255.)
	 mag = 255.;
       pdata2[idx3] = (D_uint8)mag;
       if(count <= 100){
	 if( fabs(((double)(unsigned int)pDataMag[idx3]) - mag) > 1.0){
	   printf("x=%d y=%d mag=%d  mag2=%.2f\n", x,y,
		  (unsigned int)pDataMag[idx3],mag);
	   ++count;
	 }
       }
     }
   }
   imgMag2.save("/tmp/sobel2.pgm");
 }

#endif
}




///Compute gradient magnitude of horizontal edges using the dy Sobel kernel
/** The gradient magnitude is defined as abs(dy), where dy is the
    gradient in the y direction using the sobel kernels.  For 8-bit
    grayscale images, the resulting imgGradMag will also be 8-bit
    grayscale, clipped at 255.  This function should be slightly
    faster for 8-bit grayscale, than using the DConvolver class with
    the same kernel.  The dy Sobel kernel is:
\code
                             -1 -2 -1
                        dy=   0  0  0
                              1  2  1
\endcode
*/
void DEdgeDetector::sobelHorizEdges_(DImage &imgGradMag, const DImage &imgSrc){
  int w, h, wm1, hm1, wp1;
  int idx;
  int dy;
  D_uint8 *pData;
  D_uint8 *pDataMag;
  double mag;
  if((&imgGradMag) == (&imgSrc)){
    fprintf(stderr, "DEdgeDetector::sobelHorizEdges_() imgGradMag==imgSrc!\n");
    abort();
  }
  if(DImage::DImage_u8 != imgSrc.getImageType()){
    fprintf(stderr, "DEdgeDetector::sobelHorizEdges_() only handles "
	    "8-bit grayscale\n");
    abort();
  }
  w = imgSrc.width();
  h = imgSrc.height();
  wm1 = w-1;
  wp1 = w+1;
  hm1 = h-1;

  imgGradMag.create(w, h, DImage::DImage_u8);
  pData = imgSrc.dataPointer_u8();
  pDataMag = imgGradMag.dataPointer_u8();
  // handle first row:
  //   first column of first row
  dy = ((int)((3*((unsigned int)pData[w])) +
	      ((unsigned int)pData[wp1]))) -
    ((int)((3*((unsigned int)pData[0])) +
	   ((unsigned int)pData[1])));
  mag = fabs(dy);
  pDataMag[0] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    
  //   main body of first row
  for(int x=1; x < wm1; ++x){
    dy = ((int)(((unsigned int)pData[x+wm1]) +
		(2*((unsigned int)pData[x+w])) +
		((unsigned int)pData[x+wp1]))) -
      ((int)(((unsigned int)pData[x-1]) +
	     (2*((unsigned int)pData[x])) +
	     ((unsigned int)pData[x+1])));
    mag = fabs(dy);
    pDataMag[x] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  }
  //   last column of first row
  dy = ((int)(((unsigned int)pData[wm1+wm1]) +
	      (3*((unsigned int)pData[wm1+w])))) -
    ((int)(((unsigned int)pData[wm1-1]) +
	   (3*((unsigned int)pData[wm1]))));
  mag = fabs(dy);
  pDataMag[wm1] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;

  idx = w;
  // handle main body
  for(int y=1; y < hm1; ++y){
    int idxNW, idxSW;
    //   first col of current row
    dy = ((int)((3*((unsigned int)pData[idx+w])) +
		((unsigned int)pData[idx+wp1]))) -
      ((int)((3*((unsigned int)pData[idx-w])) +
	     ((unsigned int)pData[idx-wm1])));
    mag = fabs(dy);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    ++idx;
    //   main body of current row
    for(int x=1; x < wm1; ++x){
      idxNW = idx-wp1;
      idxSW = idx+wm1;
      dy = ((int)(((unsigned int)pData[idxSW]) +
		  (2*((unsigned int)pData[idxSW+1])) +
		  ((unsigned int)pData[idxSW+2]))) -
	((int)(((unsigned int)pData[idxNW]) +
	       (2*((unsigned int)pData[idxNW+1])) +
	       ((unsigned int)pData[idxNW+2])));
      mag = fabs(dy);
      pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
      ++idx;
    }
    //   last col of current row
    idxNW = idx-wp1;
    idxSW = idx+wm1;
    dy = ((int)(((unsigned int)pData[idxSW]) +
		(3*((unsigned int)pData[idxSW+1])))) -
      ((int)(((unsigned int)pData[idxNW]) +
	     (3*((unsigned int)pData[idxNW+1]))));
    mag = fabs(dy);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    ++idx;
  }

  // handle last row
  //   first col last row
  dy = ((int)((3*((unsigned int)pData[idx])) +
	      ((unsigned int)pData[idx+1]))) -
    ((int)((3*((unsigned int)pData[idx-w])) +
	   ((unsigned int)pData[idx-wm1])));
  mag = fabs(dy);
  pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  ++idx;
  
  //   main body of last row
  for(int x=1; x < wm1; ++x){
    dy = ((int)(((unsigned int)pData[idx-1]) +
		(2*((unsigned int)pData[idx])) +
		((unsigned int)pData[idx+1]))) -
      ((int)(((unsigned int)pData[idx-wp1]) +
	     (2*((unsigned int)pData[idx-w])) +
	     ((unsigned int)pData[idx-wm1])));
    mag = fabs(dy);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    ++idx;
  }
  //   last col of last row
  dy = ((int)((3*((unsigned int)pData[idx])) +
	      ((unsigned int)pData[idx-1]))) -
    ((int)((3*((unsigned int)pData[idx-w])) +
	   ((unsigned int)pData[idx-wp1])));
  mag = fabs(dy);
  pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  ++idx;
  
#if 0
  // this is to make sure it worked as expected
 {
   DImage imgDy;
   DImage imgMag2;
   DKernel2D kernDy;
   int count = 0;
   float *pdataDy;
   D_uint8 *pdata2;
   float rgfltKernDy[] = {-1, -2, -1,   0, 0,0,    1, 2, 1};
   
   kernDy.setData_flt(rgfltKernDy, 3, 3, false);

   printf("\nkernDy:\n");
   kernDy.print();
   imgMag2.create(w,h,DImage::DImage_u8);
   pdata2 = imgMag2.dataPointer_u8();
   DConvolver::convolve_(imgDy, imgSrc, kernDy, false, false, false);

   pdataDy = imgDy.dataPointer_flt();
   
   for(int y = 0, idx3 = 0; y < h; ++y){
     for(int x = 0; x < w; ++x, ++idx3){
       mag = 0.5+ fabs(pdataDy[idx3]);
       if(mag > 255.)
	 mag = 255.;
       pdata2[idx3] = (D_uint8)mag;
       if(count <= 100){
	 if( fabs(((double)(unsigned int)pDataMag[idx3]) - mag) > 1.0){
	   printf("x=%d y=%d mag=%d  mag2=%.2f\n", x,y,
		  (unsigned int)pDataMag[idx3],mag);
	   ++count;
	 }
       }
     }
   }
   imgMag2.save("/tmp/sobel2h.pgm");
 }

#endif
}



///Compute gradient magnitude of vertical edges using the dx Sobel kernel
/** The gradient magnitude is defined as abs(dx), where dx is the
    gradient in the x direction using the sobel kernels.  For 8-bit
    grayscale images, the resulting imgGradMag will also be 8-bit
    grayscale, clipped at 255.  This function should be slightly
    faster for 8-bit grayscale, than using the DConvolver class with
    the same kernel.  The dx Sobel kernel is:
\code
                      -1  0  1 
                  dx= -2  0  2 
                      -1  0  1 
\endcode
*/
void DEdgeDetector::sobelVertEdges_(DImage &imgGradMag, const DImage &imgSrc){
  int w, h, wm1, hm1, wp1;
  int idx;
  int dx;
  D_uint8 *pData;
  D_uint8 *pDataMag;
  double mag;
  if((&imgGradMag) == (&imgSrc)){
    fprintf(stderr, "DEdgeDetector::sobelVertEdges_() imgGradMag==imgSrc!\n");
    abort();
  }
  if(DImage::DImage_u8 != imgSrc.getImageType()){
    fprintf(stderr, "DEdgeDetector::sobelVertEdges_() only handles "
	    "8-bit grayscale\n");
    abort();
  }
  w = imgSrc.width();
  h = imgSrc.height();
  wm1 = w-1;
  wp1 = w+1;
  hm1 = h-1;

  imgGradMag.create(w, h, DImage::DImage_u8);
  pData = imgSrc.dataPointer_u8();
  pDataMag = imgGradMag.dataPointer_u8();
  // handle first row:
  //   first column of first row
  dx = ((int)((3*((unsigned int)pData[1])) +
	      ((unsigned int)pData[wp1]))) -
    ((int)((3*((unsigned int)pData[0])) +
	   ((unsigned int)pData[w])));
  mag = fabs(dx);
  pDataMag[0] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    
  //   main body of first row
  for(int x=1; x < wm1; ++x){
    dx = ((int)((3*((unsigned int)pData[x+1])) +
		((unsigned int)pData[x+wp1]))) -
      ((int)((3*((unsigned int)pData[x-1])) +
	     ((unsigned int)pData[x+wm1])));
    mag = fabs(dx);
    pDataMag[x] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  }
  //   last column of first row
  dx = ((int)((3*((unsigned int)pData[wm1])) +
	      ((unsigned int)pData[wm1+w]))) -
    ((int)((3*((unsigned int)pData[wm1-1])) +
	   ((unsigned int)pData[wm1+wm1])));
  mag = fabs(dx);
  pDataMag[wm1] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;

  idx = w;
  // handle main body
  for(int y=1; y < hm1; ++y){
    int idxNW, idxSW;
    //   first col of current row
    dx = ((int)(((unsigned int)pData[idx-wm1]) +
		(2*((unsigned int)pData[idx+1])) +
		((unsigned int)pData[idx+wp1]))) -
      ((int)(((unsigned int)pData[idx-w]) +
	     (2*((unsigned int)pData[idx])) +
	     ((unsigned int)pData[idx+w])));
    mag = fabs(dx);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    ++idx;
    //   main body of current row
    for(int x=1; x < wm1; ++x){
      idxNW = idx-wp1;
      idxSW = idx+wm1;
      dx = ((int)(((unsigned int)pData[idxNW+2]) +
		  (2*((unsigned int)pData[idx+1])) +
		  ((unsigned int)pData[idxSW+2]))) -
	((int)(((unsigned int)pData[idxNW]) +
	       (2*((unsigned int)pData[idx-1])) +
	       ((unsigned int)pData[idxSW])));
      mag = fabs(dx);
      pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
      ++idx;
    }
    //   last col of current row
    idxNW = idx-wp1;
    idxSW = idx+wm1;
    dx = ((int)(((unsigned int)pData[idxNW+1]) +
		(2*((unsigned int)pData[idx])) +
		((unsigned int)pData[idxSW+1]))) -
      ((int)(((unsigned int)pData[idxNW]) +
	     (2*((unsigned int)pData[idx-1])) +
	     ((unsigned int)pData[idxSW])));
    mag = fabs(dx);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    ++idx;
  }

  // handle last row
  //   first col last row
  dx = ((int)(((unsigned int)pData[idx-wm1]) +
	      (3*((unsigned int)pData[idx+1])))) -
    ((int)(((unsigned int)pData[idx-w]) +
	   (3*((unsigned int)pData[idx]))));
  mag = fabs(dx);
  pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  ++idx;
  
  //   main body of last row
  for(int x=1; x < wm1; ++x){
    dx = ((int)((3*((unsigned int)pData[idx+1])) +
		((unsigned int)pData[idx-wm1]))) -
      ((int)((3*((unsigned int)pData[idx-1])) +
	     ((unsigned int)pData[idx-wp1])));
    mag = fabs(dx);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    ++idx;
  }
  //   last col of last row
  dx = ((int)(((unsigned int)pData[idx-w]) +
	      (3*((unsigned int)pData[idx])))) -
    ((int)(((unsigned int)pData[idx-wp1]) +
	   (3*((unsigned int)pData[idx-1]))));
  mag = fabs(dx);
  pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  ++idx;
  
#if 0
  // this is to make sure it worked as expected
 {
   DImage imgDx;
   DImage imgMag2;
   DKernel2D kernDx;
   int count = 0;
   float *pdataDx;
   D_uint8 *pdata2;
   float rgfltKernDx[] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
   
   kernDx.setData_flt(rgfltKernDx, 3, 3, false);

   printf("kernDx:\n");
   kernDx.print();
   imgMag2.create(w,h,DImage::DImage_u8);
   pdata2 = imgMag2.dataPointer_u8();
   DConvolver::convolve_(imgDx, imgSrc, kernDx, false, false, false);

   pdataDx = imgDx.dataPointer_flt();
   
   for(int y = 0, idx3 = 0; y < h; ++y){
     for(int x = 0; x < w; ++x, ++idx3){
       mag = 0.5+ fabs(pdataDx[idx3]);
       if(mag > 255.)
	 mag = 255.;
       pdata2[idx3] = (D_uint8)mag;
       if(count <= 100){
	 if( fabs(((double)(unsigned int)pDataMag[idx3]) - mag) > 1.0){
	   printf("x=%d y=%d mag=%d  mag2=%.2f\n", x,y,
		  (unsigned int)pDataMag[idx3],mag);
	   ++count;
	 }
       }
     }
   }
   imgMag2.save("/tmp/sobel2v.pgm");
 }

#endif
}




void DEdgeDetector::sobelHorizEdgesKeepSign_(DImage &imgGradMag,
					     const DImage &imgSrc){
  int w, h, wm1, hm1, wp1;
  int idx;
  int dy;
  D_uint8 *pData;
  D_uint8 *pDataMag;
  if((&imgGradMag) == (&imgSrc)){
    fprintf(stderr, "DEdgeDetector::sobelHorizEdgesKeepSign_() "
	    "imgGradMag==imgSrc!\n");
    abort();
  }
  if(DImage::DImage_u8 != imgSrc.getImageType()){
    fprintf(stderr, "DEdgeDetector::sobelHorizEdgesKeepSign_() only handles "
	    "8-bit grayscale\n");
    abort();
  }
  w = imgSrc.width();
  h = imgSrc.height();
  wm1 = w-1;
  wp1 = w+1;
  hm1 = h-1;

  imgGradMag.create(w, h, DImage::DImage_u8);
  pData = imgSrc.dataPointer_u8();
  pDataMag = imgGradMag.dataPointer_u8();
  // handle first row:
  //   first column of first row
  dy = 128 + ((int)((3*((unsigned int)pData[w])) +
		    ((unsigned int)pData[wp1]))) -
    ((int)((3*((unsigned int)pData[0])) +
	   ((unsigned int)pData[1])));
  if(dy < 0)
    dy = 0;
  if(dy > 255)
    dy = 255;
  pDataMag[0] = (D_uint8)dy;
    
  //   main body of first row
  for(int x=1; x < wm1; ++x){
    dy = 128 + ((int)(((unsigned int)pData[x+wm1]) +
		      (2*((unsigned int)pData[x+w])) +
		      ((unsigned int)pData[x+wp1]))) -
      ((int)(((unsigned int)pData[x-1]) +
	     (2*((unsigned int)pData[x])) +
	     ((unsigned int)pData[x+1])));
    if(dy < 0)
      dy = 0;
    if(dy > 255)
      dy = 255;
    pDataMag[x] = (D_uint8)dy;
  }
  //   last column of first row
  dy = 128 + ((int)(((unsigned int)pData[wm1+wm1]) +
		    (3*((unsigned int)pData[wm1+w])))) -
    ((int)(((unsigned int)pData[wm1-1]) +
	   (3*((unsigned int)pData[wm1]))));
  if(dy < 0)
    dy = 0;
  if(dy > 255)
    dy = 255;
  pDataMag[wm1] = (D_uint8)dy;

  idx = w;
  // handle main body
  for(int y=1; y < hm1; ++y){
    int idxNW, idxSW;
    //   first col of current row
    dy = 128 + ((int)((3*((unsigned int)pData[idx+w])) +
		      ((unsigned int)pData[idx+wp1]))) -
      ((int)((3*((unsigned int)pData[idx-w])) +
	     ((unsigned int)pData[idx-wm1])));
    if(dy < 0)
      dy = 0;
    if(dy > 255)
      dy = 255;
    pDataMag[idx] = (D_uint8)dy;
    ++idx;
    //   main body of current row
    for(int x=1; x < wm1; ++x){
      idxNW = idx-wp1;
      idxSW = idx+wm1;
      dy = 128 + ((int)(((unsigned int)pData[idxSW]) +
			(2*((unsigned int)pData[idxSW+1])) +
		  ((unsigned int)pData[idxSW+2]))) -
	((int)(((unsigned int)pData[idxNW]) +
	       (2*((unsigned int)pData[idxNW+1])) +
	       ((unsigned int)pData[idxNW+2])));
      if(dy < 0)
	dy = 0;
      if(dy > 255)
	dy = 255;
      pDataMag[idx] = (D_uint8)dy;
      ++idx;
    }
    //   last col of current row
    idxNW = idx-wp1;
    idxSW = idx+wm1;
    dy = 128 + ((int)(((unsigned int)pData[idxSW]) +
		      (3*((unsigned int)pData[idxSW+1])))) -
      ((int)(((unsigned int)pData[idxNW]) +
	     (3*((unsigned int)pData[idxNW+1]))));
    if(dy < 0)
      dy = 0;
    if(dy > 255)
      dy = 255;
    pDataMag[idx] = (D_uint8)dy;
    ++idx;
  }

  // handle last row
  //   first col last row
  dy = 128 + ((int)((3*((unsigned int)pData[idx])) +
		    ((unsigned int)pData[idx+1]))) -
    ((int)((3*((unsigned int)pData[idx-w])) +
	   ((unsigned int)pData[idx-wm1])));
  if(dy < 0)
    dy = 0;
  if(dy > 255)
    dy = 255;
  pDataMag[idx] = (D_uint8)dy;
  ++idx;
  
  //   main body of last row
  for(int x=1; x < wm1; ++x){
    dy = 128 + ((int)(((unsigned int)pData[idx-1]) +
		      (2*((unsigned int)pData[idx])) +
		      ((unsigned int)pData[idx+1]))) -
      ((int)(((unsigned int)pData[idx-wp1]) +
	     (2*((unsigned int)pData[idx-w])) +
	     ((unsigned int)pData[idx-wm1])));
    if(dy < 0)
      dy = 0;
    if(dy > 255)
      dy = 255;
    pDataMag[idx] = (D_uint8)dy;
    ++idx;
  }
  //   last col of last row
  dy = 128 + ((int)((3*((unsigned int)pData[idx])) +
		    ((unsigned int)pData[idx-1]))) -
    ((int)((3*((unsigned int)pData[idx-w])) +
	   ((unsigned int)pData[idx-wp1])));
  if(dy < 0)
    dy = 0;
  if(dy > 255)
    dy = 255;
  pDataMag[idx] = (D_uint8)dy;
  ++idx;
  
#if 0
  // this is to make sure it worked as expected
 {
   DImage imgDy;
   DImage imgMag2;
   DKernel2D kernDy;
   int count = 0;
   float *pdataDy;
   D_uint8 *pdata2;
   float rgfltKernDy[] = {-1, -2, -1,   0, 0,0,    1, 2, 1};
   double mag;
   
   kernDy.setData_flt(rgfltKernDy, 3, 3, false);

   printf("\nkernDy:\n");
   kernDy.print();
   imgMag2.create(w,h,DImage::DImage_u8);
   pdata2 = imgMag2.dataPointer_u8();
   DConvolver::convolve_(imgDy, imgSrc, kernDy, false, false, false);

   pdataDy = imgDy.dataPointer_flt();
   
   for(int y = 0, idx3 = 0; y < h; ++y){
     for(int x = 0; x < w; ++x, ++idx3){
       mag = 128. + pdataDy[idx3];
       if(mag > 255.)
	 mag = 255.;
       if(mag < 0.)
	 mag = 0.;
       pdata2[idx3] = (D_uint8)mag;
       if(count <= 100){
	 if( fabs(((double)(unsigned int)pDataMag[idx3]) - mag) > 1.0){
	   printf("x=%d y=%d mag=%d  mag2=%.2f\n", x,y,
		  (unsigned int)pDataMag[idx3],mag);
	   ++count;
	 }
       }
     }
   }
   imgMag2.save("/tmp/sobel2hs.pgm");
 }

#endif
}



void DEdgeDetector::sobelVertEdgesKeepSign_(DImage &imgGradMag,
					    const DImage &imgSrc){
  int w, h, wm1, hm1, wp1;
  int idx;
  int dx;
  D_uint8 *pData;
  D_uint8 *pDataMag;
  double mag;
  if((&imgGradMag) == (&imgSrc)){
    fprintf(stderr, "DEdgeDetector::sobelVertEdgesKeepSign_() "
	    "imgGradMag==imgSrc!\n");
    abort();
  }
  if(DImage::DImage_u8 != imgSrc.getImageType()){
    fprintf(stderr, "DEdgeDetector::sobelVertEdgesKeepSign_() only handles "
	    "8-bit grayscale\n");
    abort();
  }
  w = imgSrc.width();
  h = imgSrc.height();
  wm1 = w-1;
  wp1 = w+1;
  hm1 = h-1;

  imgGradMag.create(w, h, DImage::DImage_u8);
  pData = imgSrc.dataPointer_u8();
  pDataMag = imgGradMag.dataPointer_u8();
  // handle first row:
  //   first column of first row
  dx = 128 + ((int)((3*((unsigned int)pData[1])) +
		    ((unsigned int)pData[wp1]))) -
    ((int)((3*((unsigned int)pData[0])) +
	   ((unsigned int)pData[w])));
  if(dx < 0)
    dx = 0;
  if(dx > 255)
    dx = 255;
  pDataMag[0] = (D_uint8)dx;
    
  //   main body of first row
  for(int x=1; x < wm1; ++x){
    dx = 128 + ((int)((3*((unsigned int)pData[x+1])) +
		      ((unsigned int)pData[x+wp1]))) -
      ((int)((3*((unsigned int)pData[x-1])) +
	     ((unsigned int)pData[x+wm1])));
    if(dx < 0)
      dx = 0;
    if(dx > 255)
      dx = 255;
    pDataMag[x] = (D_uint8)dx;
  }
  //   last column of first row
  dx = 128 + ((int)((3*((unsigned int)pData[wm1])) +
		    ((unsigned int)pData[wm1+w]))) -
    ((int)((3*((unsigned int)pData[wm1-1])) +
	   ((unsigned int)pData[wm1+wm1])));
  if(dx < 0)
    dx = 0;
  if(dx > 255)
    dx = 255;
  pDataMag[wm1] = (D_uint8)dx;

  idx = w;
  // handle main body
  for(int y=1; y < hm1; ++y){
    int idxNW, idxSW;
    //   first col of current row
    dx = 128 + ((int)(((unsigned int)pData[idx-wm1]) +
		      (2*((unsigned int)pData[idx+1])) +
		      ((unsigned int)pData[idx+wp1]))) -
      ((int)(((unsigned int)pData[idx-w]) +
	     (2*((unsigned int)pData[idx])) +
	     ((unsigned int)pData[idx+w])));
    if(dx < 0)
      dx = 0;
    if(dx > 255)
      dx = 255;
    pDataMag[idx] = (D_uint8)dx;
    ++idx;
    //   main body of current row
    for(int x=1; x < wm1; ++x){
      idxNW = idx-wp1;
      idxSW = idx+wm1;
      dx = 128 + ((int)(((unsigned int)pData[idxNW+2]) +
			(2*((unsigned int)pData[idx+1])) +
			((unsigned int)pData[idxSW+2]))) -
	((int)(((unsigned int)pData[idxNW]) +
	       (2*((unsigned int)pData[idx-1])) +
	       ((unsigned int)pData[idxSW])));
      if(dx < 0)
	dx = 0;
      if(dx > 255)
	dx = 255;
      pDataMag[idx] = (D_uint8)dx;
      ++idx;
    }
    //   last col of current row
    idxNW = idx-wp1;
    idxSW = idx+wm1;
    dx = 128 + ((int)(((unsigned int)pData[idxNW+1]) +
		      (2*((unsigned int)pData[idx])) +
		      ((unsigned int)pData[idxSW+1]))) -
      ((int)(((unsigned int)pData[idxNW]) +
	     (2*((unsigned int)pData[idx-1])) +
	     ((unsigned int)pData[idxSW])));
    if(dx < 0)
      dx = 0;
    if(dx > 255)
      dx = 255;
    pDataMag[idx] = (D_uint8)dx;
    ++idx;
  }

  // handle last row
  //   first col last row
  dx = 128 + ((int)(((unsigned int)pData[idx-wm1]) +
		    (3*((unsigned int)pData[idx+1])))) -
    ((int)(((unsigned int)pData[idx-w]) +
	   (3*((unsigned int)pData[idx]))));
  if(dx < 0)
    dx = 0;
  if(dx > 255)
    dx = 255;
  pDataMag[idx] = (D_uint8)dx;
  ++idx;
  
  //   main body of last row
  for(int x=1; x < wm1; ++x){
    dx = 128 + ((int)((3*((unsigned int)pData[idx+1])) +
		      ((unsigned int)pData[idx-wm1]))) -
      ((int)((3*((unsigned int)pData[idx-1])) +
	     ((unsigned int)pData[idx-wp1])));
    if(dx < 0)
      dx = 0;
    if(dx > 255)
      dx = 255;
    pDataMag[idx] = (D_uint8)dx;
    ++idx;
  }
  //   last col of last row
  dx = 128 + ((int)(((unsigned int)pData[idx-w]) +
		    (3*((unsigned int)pData[idx])))) -
    ((int)(((unsigned int)pData[idx-wp1]) +
	   (3*((unsigned int)pData[idx-1]))));
  if(dx < 0)
    dx = 0;
  if(dx > 255)
    dx = 255;
  pDataMag[idx] = (D_uint8)dx;
  ++idx;
  
#if 0
  // this is to make sure it worked as expected
 {
   DImage imgDx;
   DImage imgMag2;
   DKernel2D kernDx;
   int count = 0;
   float *pdataDx;
   D_uint8 *pdata2;
   float rgfltKernDx[] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
   
   kernDx.setData_flt(rgfltKernDx, 3, 3, false);

   printf("kernDx:\n");
   kernDx.print();
   imgMag2.create(w,h,DImage::DImage_u8);
   pdata2 = imgMag2.dataPointer_u8();
   DConvolver::convolve_(imgDx, imgSrc, kernDx, false, false, false);

   pdataDx = imgDx.dataPointer_flt();
   
   for(int y = 0, idx3 = 0; y < h; ++y){
     for(int x = 0; x < w; ++x, ++idx3){
       mag = 128.+ (pdataDx[idx3]);
       if(mag > 255.)
	 mag = 255.;
       if(mag < 0.)
	 mag = 0.;
       pdata2[idx3] = (D_uint8)mag;
       if(count <= 100){
	 if( fabs(((double)(unsigned int)pDataMag[idx3]) - mag) > 1.0){
	   printf("x=%d y=%d mag=%d  mag2=%.2f\n", x,y,
		  (unsigned int)pDataMag[idx3],mag);
	   ++count;
	 }
       }
     }
   }
   imgMag2.save("/tmp/sobel2vs.pgm");
 }

#endif
}

///Convert the gradient directions (in radians) to edge directions (in radians)
/**imgGradMag is assumed to be a single-channel double DImage filled
   with gradient directions.  This function is provided as a
   convenience since edge directions are perpendicular to the gradient
   direction.  This function returns the edge directions in radians
   with angles from [0..2*PI) if fQuads1and2 is false or angles from
   [0..PI) if fQuads1and2 is true. (To limit angles to quadrants 1 and
   2, set fQuads1and2 to true).
**/
void DEdgeDetector::convertGradDirsToEdgeDirs(DImage &imgGradMag, bool
					      fQuads1and2){
  double *pdbl; int w, h;

  if((DImage::DImage_dbl_multi != imgGradMag.getImageType()) ||
     (1 != imgGradMag.numChannels())){
    fprintf(stderr, "DEdgeDetector::convertGradDirsToEdgeDirs() only works "
	    "for single-channel images of type DImage_dbl_multi\n");
    abort();
  }
  w = imgGradMag.width();
  h = imgGradMag.height();
  if((w < 1) || (h < 1)){
    fprintf(stderr, "DEdgeDetector::convertGradDirsToEdgeDirs() image size "
	    "less than 1\n");
    abort();
  }
  pdbl = imgGradMag.dataPointer_dbl();
  if(fQuads1and2){
    for(int y = 0, idx = 0; y < h; ++y){
      for(int x = 0; x < w; ++x, ++idx){
	pdbl[idx] -= M_PI_2;
	if(pdbl[idx] >= M_PI)
	  pdbl[idx] -= M_PI;
	else if(pdbl[idx] < 0.){
	  pdbl[idx] += M_PI;
	  if(pdbl[idx] < 0.)// need this second time due to the -90 degs
	    pdbl[idx] += M_PI;
	}
      }
    }
  }
  else{ // allow all 4 quadrants, but make all angles positive
    for(int y = 0, idx = 0; y < h; ++y){
      for(int x = 0; x < w; ++x, ++idx){
	pdbl[idx] -= M_PI_2;
	if(pdbl[idx] >= (2 * M_PI)){
	  pdbl[idx] -= (2 * M_PI);
	}
	else if(pdbl[idx] < 0.){
	  pdbl[idx] += (2 * M_PI);
	  if(pdbl[idx] < 0.) // need this second time due to the -90 degs
	    pdbl[idx] += (2 * M_PI);
	}
      }
    }
  }
}


























///Compute gradient magnitude of horizontal edges using the dy Prewitt kernel
/** The gradient magnitude is defined as abs(dy), where dy is the
    gradient in the y direction using the prewitt kernels.  For 8-bit
    grayscale images, the resulting imgGradMag will also be 8-bit
    grayscale, clipped at 255.  This function should be slightly
    faster for 8-bit grayscale, than using the DConvolver class with
    the same kernel.  The dy Prewitt kernel is:
\code
                             -1 -1 -1
                        dy=   0  0  0
                              1  1  1
\endcode
*/
void DEdgeDetector::prewittHorizEdges_(DImage &imgGradMag, const DImage &imgSrc){
  int w, h, wm1, hm1, wp1;
  int idx;
  int dy;
  D_uint8 *pData;
  D_uint8 *pDataMag;
  double mag;
  if((&imgGradMag) == (&imgSrc)){
    fprintf(stderr,"DEdgeDetector::prewittHorizEdges_() "
	    "imgGradMag==imgSrc!\n");
    abort();
  }
  if(DImage::DImage_u8 != imgSrc.getImageType()){
    fprintf(stderr, "DEdgeDetector::prewittHorizEdges_() only handles "
	    "8-bit grayscale\n");
    abort();
  }
  w = imgSrc.width();
  h = imgSrc.height();
  wm1 = w-1;
  wp1 = w+1;
  hm1 = h-1;

  imgGradMag.create(w, h, DImage::DImage_u8);
  pData = imgSrc.dataPointer_u8();
  pDataMag = imgGradMag.dataPointer_u8();
  // handle first row:
  //   first column of first row
  dy = ((int)((2*((unsigned int)pData[w])) +
	      ((unsigned int)pData[wp1]))) -
    ((int)((2*((unsigned int)pData[0])) +
	   ((unsigned int)pData[1])));
  mag = fabs(dy);
  pDataMag[0] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    
  //   main body of first row
  for(int x=1; x < wm1; ++x){
    dy = ((int)(((unsigned int)pData[x+wm1]) +
		(((unsigned int)pData[x+w])) +
		((unsigned int)pData[x+wp1]))) -
      ((int)(((unsigned int)pData[x-1]) +
	     (((unsigned int)pData[x])) +
	     ((unsigned int)pData[x+1])));
    mag = fabs(dy);
    pDataMag[x] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  }
  //   last column of first row
  dy = ((int)(((unsigned int)pData[wm1+wm1]) +
	      (2*((unsigned int)pData[wm1+w])))) -
    ((int)(((unsigned int)pData[wm1-1]) +
	   (2*((unsigned int)pData[wm1]))));
  mag = fabs(dy);
  pDataMag[wm1] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;

  idx = w;
  // handle main body
  for(int y=1; y < hm1; ++y){
    int idxNW, idxSW;
    //   first col of current row
    dy = ((int)((2*((unsigned int)pData[idx+w])) +
		((unsigned int)pData[idx+wp1]))) -
      ((int)((2*((unsigned int)pData[idx-w])) +
	     ((unsigned int)pData[idx-wm1])));
    mag = fabs(dy);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    ++idx;
    //   main body of current row
    for(int x=1; x < wm1; ++x){
      idxNW = idx-wp1;
      idxSW = idx+wm1;
      dy = ((int)(((unsigned int)pData[idxSW]) +
		  (((unsigned int)pData[idxSW+1])) +
		  ((unsigned int)pData[idxSW+2]))) -
	((int)(((unsigned int)pData[idxNW]) +
	       (((unsigned int)pData[idxNW+1])) +
	       ((unsigned int)pData[idxNW+2])));
      mag = fabs(dy);
      pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
      ++idx;
    }
    //   last col of current row
    idxNW = idx-wp1;
    idxSW = idx+wm1;
    dy = ((int)(((unsigned int)pData[idxSW]) +
		(2*((unsigned int)pData[idxSW+1])))) -
      ((int)(((unsigned int)pData[idxNW]) +
	     (2*((unsigned int)pData[idxNW+1]))));
    mag = fabs(dy);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    ++idx;
  }

  // handle last row
  //   first col last row
  dy = ((int)((2*((unsigned int)pData[idx])) +
	      ((unsigned int)pData[idx+1]))) -
    ((int)((2*((unsigned int)pData[idx-w])) +
	   ((unsigned int)pData[idx-wm1])));
  mag = fabs(dy);
  pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  ++idx;
  
  //   main body of last row
  for(int x=1; x < wm1; ++x){
    dy = ((int)(((unsigned int)pData[idx-1]) +
		(((unsigned int)pData[idx])) +
		((unsigned int)pData[idx+1]))) -
      ((int)(((unsigned int)pData[idx-wp1]) +
	     (((unsigned int)pData[idx-w])) +
	     ((unsigned int)pData[idx-wm1])));
    mag = fabs(dy);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    ++idx;
  }
  //   last col of last row
  dy = ((int)((2*((unsigned int)pData[idx])) +
	      ((unsigned int)pData[idx-1]))) -
    ((int)((2*((unsigned int)pData[idx-w])) +
	   ((unsigned int)pData[idx-wp1])));
  mag = fabs(dy);
  pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  ++idx;
  
#if 0
  // this is to make sure it worked as expected
 {
   DImage imgDy;
   DImage imgMag2;
   DKernel2D kernDy;
   int count = 0;
   float *pdataDy;
   D_uint8 *pdata2;
   float rgfltKernDy[] = {-1, -1, -1,   0, 0,0,    1, 1, 1};
   
   kernDy.setData_flt(rgfltKernDy, 3, 3, false);

   printf("\nkernDy:\n");
   kernDy.print();
   imgMag2.create(w,h,DImage::DImage_u8);
   pdata2 = imgMag2.dataPointer_u8();
   DConvolver::convolve_(imgDy, imgSrc, kernDy, false, false, false);

   pdataDy = imgDy.dataPointer_flt();
   
   for(int y = 0, idx3 = 0; y < h; ++y){
     for(int x = 0; x < w; ++x, ++idx3){
       mag = 0.5+ fabs(pdataDy[idx3]);
       if(mag > 255.)
	 mag = 255.;
       pdata2[idx3] = (D_uint8)mag;
       if(count <= 100){
	 if( fabs(((double)(unsigned int)pDataMag[idx3]) - mag) > 1.0){
	   printf("x=%d y=%d mag=%d  mag2=%.2f\n", x,y,
		  (unsigned int)pDataMag[idx3],mag);
	   ++count;
	 }
       }
     }
   }
   imgMag2.save("/tmp/prewitt2h.pgm");
 }

#endif
}



///Compute gradient magnitude of vertical edges using the dx Prewitt kernel
/** The gradient magnitude is defined as abs(dx), where dx is the
    gradient in the x direction using the prewitt kernels.  For 8-bit
    grayscale images, the resulting imgGradMag will also be 8-bit
    grayscale, clipped at 255.  This function should be slightly
    faster for 8-bit grayscale, than using the DConvolver class with
    the same kernel.  The dx Prewitt kernel is:
\code
                      -1  0  1 
                  dx= -1  0  1 
                      -1  0  1 
\endcode
*/
void DEdgeDetector::prewittVertEdges_(DImage &imgGradMag, const DImage &imgSrc){
  int w, h, wm1, hm1, wp1;
  int idx;
  int dx;
  D_uint8 *pData;
  D_uint8 *pDataMag;
  double mag;
  if((&imgGradMag) == (&imgSrc)){
    fprintf(stderr,"DEdgeDetector::prewittVertEdges_() imgGradMag==imgSrc!\n");
    abort();
  }
  if(DImage::DImage_u8 != imgSrc.getImageType()){
    fprintf(stderr, "DEdgeDetector::prewittVertEdges_() only handles "
	    "8-bit grayscale\n");
    abort();
  }
  w = imgSrc.width();
  h = imgSrc.height();
  wm1 = w-1;
  wp1 = w+1;
  hm1 = h-1;

  imgGradMag.create(w, h, DImage::DImage_u8);
  pData = imgSrc.dataPointer_u8();
  pDataMag = imgGradMag.dataPointer_u8();
  // handle first row:
  //   first column of first row
  dx = ((int)((2*((unsigned int)pData[1])) +
	      ((unsigned int)pData[wp1]))) -
    ((int)((2*((unsigned int)pData[0])) +
	   ((unsigned int)pData[w])));
  mag = fabs(dx);
  pDataMag[0] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    
  //   main body of first row
  for(int x=1; x < wm1; ++x){
    dx = ((int)((2*((unsigned int)pData[x+1])) +
		((unsigned int)pData[x+wp1]))) -
      ((int)((2*((unsigned int)pData[x-1])) +
	     ((unsigned int)pData[x+wm1])));
    mag = fabs(dx);
    pDataMag[x] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  }
  //   last column of first row
  dx = ((int)((2*((unsigned int)pData[wm1])) +
	      ((unsigned int)pData[wm1+w]))) -
    ((int)((2*((unsigned int)pData[wm1-1])) +
	   ((unsigned int)pData[wm1+wm1])));
  mag = fabs(dx);
  pDataMag[wm1] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;

  idx = w;
  // handle main body
  for(int y=1; y < hm1; ++y){
    int idxNW, idxSW;
    //   first col of current row
    dx = ((int)(((unsigned int)pData[idx-wm1]) +
		(((unsigned int)pData[idx+1])) +
		((unsigned int)pData[idx+wp1]))) -
      ((int)(((unsigned int)pData[idx-w]) +
	     (((unsigned int)pData[idx])) +
	     ((unsigned int)pData[idx+w])));
    mag = fabs(dx);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    ++idx;
    //   main body of current row
    for(int x=1; x < wm1; ++x){
      idxNW = idx-wp1;
      idxSW = idx+wm1;
      dx = ((int)(((unsigned int)pData[idxNW+2]) +
		  (((unsigned int)pData[idx+1])) +
		  ((unsigned int)pData[idxSW+2]))) -
	((int)(((unsigned int)pData[idxNW]) +
	       (((unsigned int)pData[idx-1])) +
	       ((unsigned int)pData[idxSW])));
      mag = fabs(dx);
      pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
      ++idx;
    }
    //   last col of current row
    idxNW = idx-wp1;
    idxSW = idx+wm1;
    dx = ((int)(((unsigned int)pData[idxNW+1]) +
		(((unsigned int)pData[idx])) +
		((unsigned int)pData[idxSW+1]))) -
      ((int)(((unsigned int)pData[idxNW]) +
	     (((unsigned int)pData[idx-1])) +
	     ((unsigned int)pData[idxSW])));
    mag = fabs(dx);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    ++idx;
  }

  // handle last row
  //   first col last row
  dx = ((int)(((unsigned int)pData[idx-wm1]) +
	      (2*((unsigned int)pData[idx+1])))) -
    ((int)(((unsigned int)pData[idx-w]) +
	   (2*((unsigned int)pData[idx]))));
  mag = fabs(dx);
  pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  ++idx;
  
  //   main body of last row
  for(int x=1; x < wm1; ++x){
    dx = ((int)((2*((unsigned int)pData[idx+1])) +
		((unsigned int)pData[idx-wm1]))) -
      ((int)((2*((unsigned int)pData[idx-1])) +
	     ((unsigned int)pData[idx-wp1])));
    mag = fabs(dx);
    pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
    ++idx;
  }
  //   last col of last row
  dx = ((int)(((unsigned int)pData[idx-w]) +
	      (2*((unsigned int)pData[idx])))) -
    ((int)(((unsigned int)pData[idx-wp1]) +
	   (2*((unsigned int)pData[idx-1]))));
  mag = fabs(dx);
  pDataMag[idx] = (mag <= 254.5) ? ((D_uint8)(mag+0.5)) : 0xff;
  ++idx;
  
#if 1
  // this is to make sure it worked as expected
 {
   DImage imgDx;
   DImage imgMag2;
   DKernel2D kernDx;
   int count = 0;
   float *pdataDx;
   D_uint8 *pdata2;
   float rgfltKernDx[] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
   
   kernDx.setData_flt(rgfltKernDx, 3, 3, false);

   printf("kernDx:\n");
   kernDx.print();
   imgMag2.create(w,h,DImage::DImage_u8);
   pdata2 = imgMag2.dataPointer_u8();
   DConvolver::convolve_(imgDx, imgSrc, kernDx, false, false, false);

   pdataDx = imgDx.dataPointer_flt();
   
   for(int y = 0, idx3 = 0; y < h; ++y){
     for(int x = 0; x < w; ++x, ++idx3){
       mag = 0.5+ fabs(pdataDx[idx3]);
       if(mag > 255.)
	 mag = 255.;
       pdata2[idx3] = (D_uint8)mag;
       if(count <= 100){
	 if( fabs(((double)(unsigned int)pDataMag[idx3]) - mag) > 1.0){
	   printf("x=%d y=%d mag=%d  mag2=%.2f\n", x,y,
		  (unsigned int)pDataMag[idx3],mag);
	   ++count;
	 }
       }
     }
   }
   imgMag2.save("/tmp/prewitt2v.pgm");
 }

#endif
}




void DEdgeDetector::prewittHorizEdgesKeepSign_(DImage &imgGradMag,
					     const DImage &imgSrc){
  int w, h, wm1, hm1, wp1;
  int idx;
  int dy;
  D_uint8 *pData;
  D_uint8 *pDataMag;
  if((&imgGradMag) == (&imgSrc)){
    fprintf(stderr, "DEdgeDetector::prewittHorizEdgesKeepSign_() "
	    "imgGradMag==imgSrc!\n");
    abort();
  }
  if(DImage::DImage_u8 != imgSrc.getImageType()){
    fprintf(stderr, "DEdgeDetector::prewittHorizEdgesKeepSign_() only handles "
	    "8-bit grayscale\n");
    abort();
  }
  w = imgSrc.width();
  h = imgSrc.height();
  wm1 = w-1;
  wp1 = w+1;
  hm1 = h-1;

  imgGradMag.create(w, h, DImage::DImage_u8);
  pData = imgSrc.dataPointer_u8();
  pDataMag = imgGradMag.dataPointer_u8();
  // handle first row:
  //   first column of first row
  dy = 128 + ((int)((2*((unsigned int)pData[w])) +
		    ((unsigned int)pData[wp1]))) -
    ((int)((2*((unsigned int)pData[0])) +
	   ((unsigned int)pData[1])));
  if(dy < 0)
    dy = 0;
  if(dy > 255)
    dy = 255;
  pDataMag[0] = (D_uint8)dy;
    
  //   main body of first row
  for(int x=1; x < wm1; ++x){
    dy = 128 + ((int)(((unsigned int)pData[x+wm1]) +
		      (((unsigned int)pData[x+w])) +
		      ((unsigned int)pData[x+wp1]))) -
      ((int)(((unsigned int)pData[x-1]) +
	     (((unsigned int)pData[x])) +
	     ((unsigned int)pData[x+1])));
    if(dy < 0)
      dy = 0;
    if(dy > 255)
      dy = 255;
    pDataMag[x] = (D_uint8)dy;
  }
  //   last column of first row
  dy = 128 + ((int)(((unsigned int)pData[wm1+wm1]) +
		    (2*((unsigned int)pData[wm1+w])))) -
    ((int)(((unsigned int)pData[wm1-1]) +
	   (2*((unsigned int)pData[wm1]))));
  if(dy < 0)
    dy = 0;
  if(dy > 255)
    dy = 255;
  pDataMag[wm1] = (D_uint8)dy;

  idx = w;
  // handle main body
  for(int y=1; y < hm1; ++y){
    int idxNW, idxSW;
    //   first col of current row
    dy = 128 + ((int)((2*((unsigned int)pData[idx+w])) +
		      ((unsigned int)pData[idx+wp1]))) -
      ((int)((2*((unsigned int)pData[idx-w])) +
	     ((unsigned int)pData[idx-wm1])));
    if(dy < 0)
      dy = 0;
    if(dy > 255)
      dy = 255;
    pDataMag[idx] = (D_uint8)dy;
    ++idx;
    //   main body of current row
    for(int x=1; x < wm1; ++x){
      idxNW = idx-wp1;
      idxSW = idx+wm1;
      dy = 128 + ((int)(((unsigned int)pData[idxSW]) +
			(((unsigned int)pData[idxSW+1])) +
		  ((unsigned int)pData[idxSW+2]))) -
	((int)(((unsigned int)pData[idxNW]) +
	       (((unsigned int)pData[idxNW+1])) +
	       ((unsigned int)pData[idxNW+2])));
      if(dy < 0)
	dy = 0;
      if(dy > 255)
	dy = 255;
      pDataMag[idx] = (D_uint8)dy;
      ++idx;
    }
    //   last col of current row
    idxNW = idx-wp1;
    idxSW = idx+wm1;
    dy = 128 + ((int)(((unsigned int)pData[idxSW]) +
		      (2*((unsigned int)pData[idxSW+1])))) -
      ((int)(((unsigned int)pData[idxNW]) +
	     (2*((unsigned int)pData[idxNW+1]))));
    if(dy < 0)
      dy = 0;
    if(dy > 255)
      dy = 255;
    pDataMag[idx] = (D_uint8)dy;
    ++idx;
  }

  // handle last row
  //   first col last row
  dy = 128 + ((int)((2*((unsigned int)pData[idx])) +
		    ((unsigned int)pData[idx+1]))) -
    ((int)((2*((unsigned int)pData[idx-w])) +
	   ((unsigned int)pData[idx-wm1])));
  if(dy < 0)
    dy = 0;
  if(dy > 255)
    dy = 255;
  pDataMag[idx] = (D_uint8)dy;
  ++idx;
  
  //   main body of last row
  for(int x=1; x < wm1; ++x){
    dy = 128 + ((int)(((unsigned int)pData[idx-1]) +
		      (((unsigned int)pData[idx])) +
		      ((unsigned int)pData[idx+1]))) -
      ((int)(((unsigned int)pData[idx-wp1]) +
	     (((unsigned int)pData[idx-w])) +
	     ((unsigned int)pData[idx-wm1])));
    if(dy < 0)
      dy = 0;
    if(dy > 255)
      dy = 255;
    pDataMag[idx] = (D_uint8)dy;
    ++idx;
  }
  //   last col of last row
  dy = 128 + ((int)((2*((unsigned int)pData[idx])) +
		    ((unsigned int)pData[idx-1]))) -
    ((int)((2*((unsigned int)pData[idx-w])) +
	   ((unsigned int)pData[idx-wp1])));
  if(dy < 0)
    dy = 0;
  if(dy > 255)
    dy = 255;
  pDataMag[idx] = (D_uint8)dy;
  ++idx;
  
#if 1
  // this is to make sure it worked as expected
 {
   DImage imgDy;
   DImage imgMag2;
   DKernel2D kernDy;
   int count = 0;
   float *pdataDy;
   D_uint8 *pdata2;
   float rgfltKernDy[] = {-1, -1, -1,   0, 0,0,    1, 1, 1};
   double mag;
   
   kernDy.setData_flt(rgfltKernDy, 3, 3, false);

   printf("\nkernDy:\n");
   kernDy.print();
   imgMag2.create(w,h,DImage::DImage_u8);
   pdata2 = imgMag2.dataPointer_u8();
   DConvolver::convolve_(imgDy, imgSrc, kernDy, false, false, false);

   pdataDy = imgDy.dataPointer_flt();
   
   for(int y = 0, idx3 = 0; y < h; ++y){
     for(int x = 0; x < w; ++x, ++idx3){
       mag = 128. + pdataDy[idx3];
       if(mag > 255.)
	 mag = 255.;
       if(mag < 0.)
	 mag = 0.;
       pdata2[idx3] = (D_uint8)mag;
       if(count <= 100){
	 if( fabs(((double)(unsigned int)pDataMag[idx3]) - mag) > 1.0){
	   printf("x=%d y=%d mag=%d  mag2=%.2f\n", x,y,
		  (unsigned int)pDataMag[idx3],mag);
	   ++count;
	 }
       }
     }
   }
   imgMag2.save("/tmp/prewitt2hs.pgm");
 }

#endif
}



void DEdgeDetector::prewittVertEdgesKeepSign_(DImage &imgGradMag,
					    const DImage &imgSrc){
  int w, h, wm1, hm1, wp1;
  int idx;
  int dx;
  D_uint8 *pData;
  D_uint8 *pDataMag;
  double mag;
  if((&imgGradMag) == (&imgSrc)){
    fprintf(stderr, "DEdgeDetector::prewittVertEdgesKeepSign_() "
	    "imgGradMag==imgSrc!\n");
    abort();
  }
  if(DImage::DImage_u8 != imgSrc.getImageType()){
    fprintf(stderr, "DEdgeDetector::prewittVertEdgesKeepSign_() only handles "
	    "8-bit grayscale\n");
    abort();
  }
  w = imgSrc.width();
  h = imgSrc.height();
  wm1 = w-1;
  wp1 = w+1;
  hm1 = h-1;

  imgGradMag.create(w, h, DImage::DImage_u8);
  pData = imgSrc.dataPointer_u8();
  pDataMag = imgGradMag.dataPointer_u8();
  // handle first row:
  //   first column of first row
  dx = 128 + ((int)((2*((unsigned int)pData[1])) +
		    ((unsigned int)pData[wp1]))) -
    ((int)((2*((unsigned int)pData[0])) +
	   ((unsigned int)pData[w])));
  if(dx < 0)
    dx = 0;
  if(dx > 255)
    dx = 255;
  pDataMag[0] = (D_uint8)dx;
    
  //   main body of first row
  for(int x=1; x < wm1; ++x){
    dx = 128 + ((int)((2*((unsigned int)pData[x+1])) +
		      ((unsigned int)pData[x+wp1]))) -
      ((int)((2*((unsigned int)pData[x-1])) +
	     ((unsigned int)pData[x+wm1])));
    if(dx < 0)
      dx = 0;
    if(dx > 255)
      dx = 255;
    pDataMag[x] = (D_uint8)dx;
  }
  //   last column of first row
  dx = 128 + ((int)((2*((unsigned int)pData[wm1])) +
		    ((unsigned int)pData[wm1+w]))) -
    ((int)((2*((unsigned int)pData[wm1-1])) +
	   ((unsigned int)pData[wm1+wm1])));
  if(dx < 0)
    dx = 0;
  if(dx > 255)
    dx = 255;
  pDataMag[wm1] = (D_uint8)dx;

  idx = w;
  // handle main body
  for(int y=1; y < hm1; ++y){
    int idxNW, idxSW;
    //   first col of current row
    dx = 128 + ((int)(((unsigned int)pData[idx-wm1]) +
		      (((unsigned int)pData[idx+1])) +
		      ((unsigned int)pData[idx+wp1]))) -
      ((int)(((unsigned int)pData[idx-w]) +
	     (((unsigned int)pData[idx])) +
	     ((unsigned int)pData[idx+w])));
    if(dx < 0)
      dx = 0;
    if(dx > 255)
      dx = 255;
    pDataMag[idx] = (D_uint8)dx;
    ++idx;
    //   main body of current row
    for(int x=1; x < wm1; ++x){
      idxNW = idx-wp1;
      idxSW = idx+wm1;
      dx = 128 + ((int)(((unsigned int)pData[idxNW+2]) +
			(((unsigned int)pData[idx+1])) +
			((unsigned int)pData[idxSW+2]))) -
	((int)(((unsigned int)pData[idxNW]) +
	       (((unsigned int)pData[idx-1])) +
	       ((unsigned int)pData[idxSW])));
      if(dx < 0)
	dx = 0;
      if(dx > 255)
	dx = 255;
      pDataMag[idx] = (D_uint8)dx;
      ++idx;
    }
    //   last col of current row
    idxNW = idx-wp1;
    idxSW = idx+wm1;
    dx = 128 + ((int)(((unsigned int)pData[idxNW+1]) +
		      (((unsigned int)pData[idx])) +
		      ((unsigned int)pData[idxSW+1]))) -
      ((int)(((unsigned int)pData[idxNW]) +
	     (((unsigned int)pData[idx-1])) +
	     ((unsigned int)pData[idxSW])));
    if(dx < 0)
      dx = 0;
    if(dx > 255)
      dx = 255;
    pDataMag[idx] = (D_uint8)dx;
    ++idx;
  }

  // handle last row
  //   first col last row
  dx = 128 + ((int)(((unsigned int)pData[idx-wm1]) +
		    (2*((unsigned int)pData[idx+1])))) -
    ((int)(((unsigned int)pData[idx-w]) +
	   (2*((unsigned int)pData[idx]))));
  if(dx < 0)
    dx = 0;
  if(dx > 255)
    dx = 255;
  pDataMag[idx] = (D_uint8)dx;
  ++idx;
  
  //   main body of last row
  for(int x=1; x < wm1; ++x){
    dx = 128 + ((int)((2*((unsigned int)pData[idx+1])) +
		      ((unsigned int)pData[idx-wm1]))) -
      ((int)((2*((unsigned int)pData[idx-1])) +
	     ((unsigned int)pData[idx-wp1])));
    if(dx < 0)
      dx = 0;
    if(dx > 255)
      dx = 255;
    pDataMag[idx] = (D_uint8)dx;
    ++idx;
  }
  //   last col of last row
  dx = 128 + ((int)(((unsigned int)pData[idx-w]) +
		    (2*((unsigned int)pData[idx])))) -
    ((int)(((unsigned int)pData[idx-wp1]) +
	   (2*((unsigned int)pData[idx-1]))));
  if(dx < 0)
    dx = 0;
  if(dx > 255)
    dx = 255;
  pDataMag[idx] = (D_uint8)dx;
  ++idx;
  
#if 1
  // this is to make sure it worked as expected
 {
   DImage imgDx;
   DImage imgMag2;
   DKernel2D kernDx;
   int count = 0;
   float *pdataDx;
   D_uint8 *pdata2;
   float rgfltKernDx[] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
   
   kernDx.setData_flt(rgfltKernDx, 3, 3, false);

   printf("kernDx:\n");
   kernDx.print();
   imgMag2.create(w,h,DImage::DImage_u8);
   pdata2 = imgMag2.dataPointer_u8();
   DConvolver::convolve_(imgDx, imgSrc, kernDx, false, false, false);

   pdataDx = imgDx.dataPointer_flt();
   
   for(int y = 0, idx3 = 0; y < h; ++y){
     for(int x = 0; x < w; ++x, ++idx3){
       mag = 128.+ (pdataDx[idx3]);
       if(mag > 255.)
	 mag = 255.;
       if(mag < 0.)
	 mag = 0.;
       pdata2[idx3] = (D_uint8)mag;
       if(count <= 100){
	 if( fabs(((double)(unsigned int)pDataMag[idx3]) - mag) > 1.0){
	   printf("x=%d y=%d mag=%d  mag2=%.2f\n", x,y,
		  (unsigned int)pDataMag[idx3],mag);
	   ++count;
	 }
       }
     }
   }
   imgMag2.save("/tmp/prewitt2vs.pgm");
 }

#endif
}


///return an 8-bit grayscale image with Laplacian zeros set to 0xff
/**Zero crossings are determined by convolving a 3x3 Laplace kernel
   (-8 in the center, 1 everywhere else) and then setting the pixels
   nearest to zero crossings to 0xff (everything else gets set to
   0x00).  For each pixel (except row and column 0), the west, north,
   and northwest neighbors' values are compared against the current
   pixel's value.  For each of the three comparisons, if the signs are
   opposite, then whichever of the two pixels being compared is closer
   to zero is set to 0xff.
*/
void DEdgeDetector::laplace3x3Zeros_(DImage &imgZeros, DImage &imgSrc,
				     bool fNearestOnly){
  DImage imgLaplace;
  DKernel2D kern;
  int w, h;
  float *pLaplace;
  D_uint8 *pZeros;

  if(DImage::DImage_u8 != imgSrc.getImageType()){
    fprintf(stderr, "DEdgeDetector::laplace3x3Zeros_() currently only "
	    "supports 8-bit grayscale images\n");
    abort();
  }

  w = imgSrc.width();
  h = imgSrc.height();
  kern.setLaplace();
  DConvolver::convolve_(imgLaplace, imgSrc, kern);
  pLaplace = imgLaplace.dataPointer_flt();
  imgZeros.create(w, h, DImage::DImage_u8);
  pZeros = imgZeros.dataPointer_u8();
  // set first row to 0
  for(int x = 0; x < w; ++x)
    pZeros[x] = 0x00;
  for(int y = 1, idx = w; y < h; ++y){
    pZeros[idx] = 0x00;
    ++idx; // since we start the inner loop at x=1 instead of x=0
    for(int x = 1; x < w; ++x, ++idx){
      if((( (pLaplace[idx] > 0.000001) && (pLaplace[idx-1] < -0.000001)) ||
	  ( (pLaplace[idx] < -0.000001) && (pLaplace[idx-1] > 0.000001)))){
	if(fNearestOnly){
	  // set the pixel that is closest to zero to 0xff
	  if( fabs(pLaplace[idx]) <= fabs(pLaplace[idx-1]))
	    pZeros[idx] = 0xff;
	  else
	    pZeros[idx-1] = 0xff;
	}
	else{
	  // set both to 0xff
	  pZeros[idx] = 0xff;
	  pZeros[idx-1] = 0xff;
	}
      }	 
      if((( (pLaplace[idx] > 0.000001) && (pLaplace[idx-w] < -0.000001)) ||
	  ( (pLaplace[idx] < -0.000001) && (pLaplace[idx-w] > 0.000001)))){
	if(fNearestOnly){
	  // set the pixel that is closest to zero to 0xff
	  if( fabs(pLaplace[idx]) <= fabs(pLaplace[idx-w]))
	    pZeros[idx] = 0xff;
	  else
	    pZeros[idx-w] = 0xff;
	}
	else{
	  // set both to 0xff
	  pZeros[idx] = 0xff;
	  pZeros[idx-w] = 0xff;
	}
      }
      if((( (pLaplace[idx] > 0.000001) && (pLaplace[idx-w-1] < -0.000001)) ||
	  ( (pLaplace[idx] < -0.000001) && (pLaplace[idx-w-1] > 0.000001)))){
	if(fNearestOnly){
	  // set the pixel that is closest to zero to 0xff
	  if( fabs(pLaplace[idx]) <= fabs(pLaplace[idx-w-1]))
	    pZeros[idx] = 0xff;
	  else
	    pZeros[idx-w-1] = 0xff;
	}
	else{
	  // set both to 0xff
	  pZeros[idx] = 0xff;
	  pZeros[idx-w-1] = 0xff;
	}
      }
    }//end for(x...
  }//end for(y...
}


// /** imgSrc is the source gradient magnitude image. imgDirs is the
//     gradient (not edge) direction image (in radians)
//  */
// void DEdgeDetector::nonMaximalSuppression_(DImage &imgDst, DImage &imgSrc,
// 					   DImage &imgDirs){
//   int w, h;
//   double *pDirs;
//   D_uint8 *pDst;
//   int dirIdx;
//   int rgDx[17]={1,1,1,0,0,-1,-1,-1,//x-offset for directions 0..16 (clockwise)
// 		-1,-1,-1,0,0,1,1,1,1};
//   int rgDy[17]={0,1,1,1,1,1,1,0,//y-offset for directions 0..16 (clockwise)
// 		0,-1,-1,-1,-1,-1,-1,0,0};
//   w = imgSrc.width();
//   h = imgSrc.height();

//   if( (w!=imgDirs.width()) || (h != imgDirs.height())){
//     fprintf(stderr, "DEdgeDetector::nonMaximalSuppression_() direction image "
// 	    "size (%dx%d) != source image size(%dx%d)\n",
// 	    imgDirs.width(), imgDirs.height(), w, h);
//     abort();
//   }
//   imgDst = imgSrc;
//   pDst = imgDst.dataPointer_u8();
//   pDirs = imgDirs.dataPointer_dbl();
//   for(int y = 0, idx = 0; y < h; ++y){
//     for(int x = 0; x < w; ++x, ++idx){
//       int testX, testY;
//       int testIdx;
//       dirIdx = 8 + (int)(pDirs[idx] / (M_PI/8.));
//       if((dirIdx < 0) || (dirIdx > 16)){
// 	fprintf(stderr,
// 		"DEdgeDetector::nonMaximalSuppression_() bad dirIdx=%d at "
// 		"x=%d,y=%d (dir=%.3f)\n",
// 		dirIdx, x, y, pDirs[idx]);
// 	abort();
//       }
//       testX = x + rgDx[dirIdx];
//       testY = y + rgDy[dirIdx];
//       testIdx = testY*w+testX;
//       if( (testX >= 0) && (testX < w) && (testY >= 0) && (testY < h) ){
// 	if(pDst[testIdx] >= pDst[idx])
// 	  pDst[idx] = 0x00;
//       }
//     }
//   }
// }


/** imgSrc is the source gradient magnitude image. imgDirs is the
    gradient (not edge) direction image (in radians)
 */
// void DEdgeDetector::nonMaximalSuppression_(DImage &imgDst, DImage &imgSrc,
// 					   DImage &imgDirs){
//   int w, h;
//   double *pDirs;
//   D_uint8 *pDst;
//     DImage imgDirs2;
//   int dirIdx;
//   int rgDx[9]={-1,0,0,1,//x-offset for directions 0..16 (clockwise)
// 		1,0,0,-1,-1};
//   int rgDy[9]={0,1,1,0,//y-offset for directions 0..16 (clockwise)
// 		0,-1,-1,0,0};
//   w = imgSrc.width();
//   h = imgSrc.height();

//   if( (w!=imgDirs.width()) || (h != imgDirs.height())){
//     fprintf(stderr, "DEdgeDetector::nonMaximalSuppression_() direction image "
// 	    "size (%dx%d) != source image size(%dx%d)\n",
// 	    imgDirs.width(), imgDirs.height(), w, h);
//     abort();
//   }
//   imgDst = imgSrc;
//   pDst = imgDst.dataPointer_u8();
//   pDirs = imgDirs.dataPointer_dbl();


//   /// for debug, fill the directions image with color-coded directions
//   {
//     D_uint8 *pDirs2;
//     imgDirs2.create(w,h,DImage::DImage_RGB);
//     pDirs2 = imgDirs2.dataPointer_u8();
//     for(int y = 0, idx = 0; y < h; ++y){
//       for(int x = 0; x < w; ++x, ++idx){
// 	int testX, testY;
// 	dirIdx = 4 + (int)(pDirs[idx] / (M_PI/4.));
// 	testX =  rgDx[dirIdx];
// 	testY =  rgDy[dirIdx];
// 	if(pDst[idx] < 5){
// 	  pDirs2[idx*3] = 0x00;
// 	  pDirs2[idx*3+1] = 0x00;
// 	  pDirs2[idx*3+2] = 0x00;
// 	}
// 	else if((testX == 1) && (testY == 0)){
// 	  pDirs2[idx*3] = 0xff;
// 	  pDirs2[idx*3+1] = 0x00;
// 	  pDirs2[idx*3+2] = 0x00;
// 	}
// 	else if((testX == 0) && (testY == -1)){
// 	  pDirs2[idx*3] = 0x00;
// 	  pDirs2[idx*3+1] = 0xff;
// 	  pDirs2[idx*3+2] = 0x00;
// 	}
// 	else if((testX == -1) && (testY == 0)){
// 	  pDirs2[idx*3] = 0x00;
// 	  pDirs2[idx*3+1] = 0x00;
// 	  pDirs2[idx*3+2] = 0xff;
// 	}
// 	else if((testX == 0) && (testY == 1)){
// 	  pDirs2[idx*3] = 0xff;
// 	  pDirs2[idx*3+1] = 0xff;
// 	  pDirs2[idx*3+2] = 0x00;
// 	}
// 	else{
// 	  pDirs2[idx*3] = 0xff;
// 	  pDirs2[idx*3+1] = 0xff;
// 	  pDirs2[idx*3+2] = 0xff;
// 	}
//       }
//     }
//   }	  


//   for(int y = 0, idx = 0; y < h; ++y){
//     for(int x = 0; x < w; ++x, ++idx){
//       int testX, testY;
//       int testIdx;
//       dirIdx = 4 + (int)(pDirs[idx] / (M_PI/4.));
//       if((dirIdx < 0) || (dirIdx > 9)){
// 	fprintf(stderr,
// 		"DEdgeDetector::nonMaximalSuppression_() bad dirIdx=%d at "
// 		"x=%d,y=%d (dir=%.3f)\n",
// 		dirIdx, x, y, pDirs[idx]);
// 	abort();
//       }
//       testX = x + rgDx[dirIdx];
//       testY = y + rgDy[dirIdx];
//       testIdx = testY*w+testX;
//       if( (testX >= 0) && (testX < w) && (testY >= 0) && (testY < h) ){
// 	if(pDst[testIdx] >= pDst[idx])
// 	  pDst[idx] = 0x00;
//       }
//     }
//   }
	

//     imgDirs = imgDirs2;
// }

void DEdgeDetector::showDirsInRGBImage_(DImage &imgDst, DImage &imgGradMag,
					DImage &imgDirs, int thresh){
  int w, h;
  double *pDirs;
  D_uint8 *pDst;
  D_uint8 *pGradMag;
  w = imgGradMag.width();
  h = imgGradMag.height();

  if( (w!=imgDirs.width()) || (h != imgDirs.height())){
    fprintf(stderr, "DEdgeDetector::nonMaximalSuppression_() direction image "
	    "size (%dx%d) != gradient magnitude image size(%dx%d)\n",
	    imgDirs.width(), imgDirs.height(), w, h);
    abort();
  }
  imgDst.create(w,h,DImage::DImage_RGB);
  pDst = imgDst.dataPointer_u8();
  pDirs = imgDirs.dataPointer_dbl();
  pGradMag = imgGradMag.dataPointer_u8();
  for(int y = 0, idx = 0; y < h; ++y){
    for(int x = 0; x < w; ++x, ++idx){
      if(pGradMag[idx] > thresh){
	double dir;
	double w0, w1;
	int R0, G0, B0;
	int R1, G1, B1;
	int dstR, dstG, dstB;
	dir = pDirs[idx];

	//NW quadrant 
	if( (dir >= -M_PI) && (dir <= -M_PI_2)){ // 
	  w1 = fabs(-M_PI - dir) / (M_PI_2);
	  R0 = 0; G0 = 255; B0 = 0;
	  R1 = 255; G1 = 255; B1 = 0;
	}
	//NE
	else if( (dir >= -M_PI_2) && (dir <= 0.)){
	  w1 = fabs(-M_PI_2 - dir) / (M_PI_2);
	  R0 = 255; G0 = 255; B0 = 0;
	  R1 = 255; G1 = 0; B1 = 0;
	}
	//SE
	else if( (dir >= 0.) && (dir <= M_PI_2)){
	  w1 = fabs(0. - dir) / (M_PI_2);
	  R0 = 255; G0 = 0; B0 = 0;
	  R1 = 0; G1 = 0; B1 = 255;
	}
	//SW
	else if( (dir >= M_PI_2) && (dir <= M_PI)){
	  w1 = fabs(M_PI_2 - dir) / (M_PI_2);
	  R0 = 0; G0 = 0; B0 = 255;
	  R1 = 0; G1 = 255; B1 = 0;
	}
	else{
	  w1 = 1.;
	  R0 = 255; G0 = 255; B0 = 255;
	  R1 = 255; G1 = 255; B1 = 255;
	}

	w0 = 1. - w1;
	dstR = (int)(w0 * R0 + w1 * R1);
	if(dstR > 255)
	  dstR = 255;
	dstG = (int)(w0 * G0 + w1 * G1);
	if(dstG > 255)
	  dstG = 255;
	dstB = (int)(w0 * B0 + w1 * B1);
	if(dstB > 255)
	  dstB = 255;
	pDst[idx*3] = (D_uint8)dstR;
	pDst[idx*3+1] = (D_uint8)dstG;
	pDst[idx*3+2] = (D_uint8)dstB;
      }
      else{
	pDst[idx*3] = 0x00;
	pDst[idx*3+1] = 0x00;
	pDst[idx*3+2] = 0x00;
      }
    }
  }

}



/** imgSrc is the source gradient magnitude image. imgDirs is the
    gradient (not edge) direction image (in radians)
 */
void DEdgeDetector::nonMaximalSuppression_(DImage &imgDst, DImage &imgSrc,
					   DImage &imgDirs){
  int w, h;
  int hm1, wm1; //h and w minus 1
  double *pDirs;
  D_uint8 *pDst;
  D_uint8 *pSrc;
  int dirIdx;
  w = imgSrc.width();
  h = imgSrc.height();
  wm1 = w-1;
  hm1 = h-1;
  
  if( (w!=imgDirs.width()) || (h != imgDirs.height())){
    fprintf(stderr, "DEdgeDetector::nonMaximalSuppression_() direction image "
	    "size (%dx%d) != source image size(%dx%d)\n",
	    imgDirs.width(), imgDirs.height(), w, h);
    abort();
  }
  imgDst = imgSrc;
  pSrc = imgSrc.dataPointer_u8();
  pDst = imgDst.dataPointer_u8();
  pDirs = imgDirs.dataPointer_dbl();



  for(int y = 0, idx = 0; y < h; ++y){
    for(int x = 0; x < w; ++x, ++idx){
      double dblDir;
      dblDir = pDirs[idx];
      if( ((dblDir >= (M_PI/4.)) && (dblDir <= (M_PI * 3. / 4.))) ||
	  ((dblDir <= (-M_PI/4.)) && (dblDir >= (-M_PI * 3. / 4.)))){
	// check above and below the pixel
	if(y > 0){
	  if( (pSrc[idx-w] > pSrc[idx]) )
	    pDst[idx] = 0x00;
	}
	if(y < hm1){
	  if( (pSrc[idx+w] > pSrc[idx]) )
	    pDst[idx] = 0x00;
	}
      }
      else{
	// check left and right of the pixel
	if(x > 0){
	  if( (pSrc[idx-1] > pSrc[idx]) )
	    pDst[idx] = 0x00;
	}
	if(x < wm1){
	  if( (pSrc[idx+1] > pSrc[idx]) )
	    pDst[idx] = 0x00;
	}
      }
    }
  }
}





