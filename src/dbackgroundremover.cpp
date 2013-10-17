#include "dbackgroundremover.h"
#include "dmedianfilter.h"
#include "dmaxfilter.h"
#include "dkernel2d.h"
#include "dconvolver.h"
#include "dedgedetector.h"
#include <string.h>

#include "dtimer.h"
#include "dvariancefilter.h"

/*comparison function for qsort()ing unsigned chars in nondecreasing order */
int bgr_compareChars(const void *p1, const void *p2){
  if ( (*((unsigned char*)(p1))) < (*((unsigned char*)(p2))) )
    return -1;
  else if ( (*((unsigned char*)(p1))) > (*((unsigned char*)(p2))) )
    return 1;
  return 0;
}



///Remove document background using Luke Hutchison's method
/**If pImgBG is non-NULL, the DImage object pointed to will be filled
 * with the background image.  If radiusX is -1, then the default
 * window size will be used depending on the size of imgSrc.  If
 * radiusY is -1, it will be set equal to radiusX.
 */
void DBackgroundRemover::removeBackground(DImage &imgDst, const DImage &imgSrc,
					  DImage *pImgBG,
					  int radiusX, int radiusY,
					  DProgress *pProg, int numThreads){
  DImage imgBGtmp;
  DImage *pBG;
  int len; // width * height * (times number of channels)

  if(NULL == pImgBG)
    pBG = &imgBGtmp;
  else
    pBG = pImgBG;

  if(-1 == radiusX){
    radiusX =
      (imgSrc.height() > imgSrc.width()) ?
      (imgSrc.height() * 3 / 200) : (imgSrc.width() * 3 / 200);
  }
  if(-1 == radiusY){
    radiusY = radiusX;
  }

  imgDst.create(imgSrc.width(), imgSrc.height(), imgSrc.getImageType(),
		imgSrc.numChannels(), imgSrc.getAllocMethod());
  len = imgSrc.width() * imgSrc.height() * imgSrc.numChannels();

  DMedianFilter::medianFilterImage(*pBG, imgSrc, false, radiusX,radiusY,
				   DMedianFilter::DMedFilt_default,pProg,
				   numThreads);
  switch(imgSrc.getImageType()){
    case DImage::DImage_u8:
    case DImage::DImage_RGB:
      {
	D_uint8 *pDataSrc;
	D_uint8 *pDataBG;
	D_uint8 *pDataDst;
	int rgHist[256];
	int pxlVal;
	unsigned char rgLUT[256];
	int pxlCount;
	int top, bot;
	
	memset(rgHist, 0, sizeof(rgHist));
	
	pDataSrc = imgSrc.dataPointer_u8();
	pDataBG = pBG->dataPointer_u8();
	pDataDst = imgDst.dataPointer_u8();
	// do the background subtraction to leave a foreground image
	for(int idx = 0; idx < len; ++idx){
// 	  pxlVal = 255 - ((int)(unsigned int)pDataBG[idx] -
// 			  (int)(unsigned int)pDataSrc[idx]);
// 	  if(pxlVal < 0)
// 	    pxlVal = 0;
// 	  if(pxlVal > 255)
// 	    pxlVal = 255;
// 	  pDataDst[idx] = (unsigned char)pxlVal;
// 	  ++(rgHist[pxlVal]);



	  float f;

#if 1
	  pxlVal =  ((int)(unsigned int)pDataBG[idx] -
		     (int)(unsigned int)pDataSrc[idx]);
	  if(pxlVal < 0)
	    pxlVal = 0;
	  f = pxlVal / 255.;
	  pxlVal = (int)((f * (int)(unsigned int)pDataSrc[idx]) +
			 (1.-f)*255.);
#else // try division instead of subtraction
	  pxlVal = (0.0001+(float)(unsigned int)(pDataSrc[idx])) /
	    (0.0001+(float)(unsigned int)(pDataBG[idx])) * 255.;

#endif
	  if(pxlVal < 0)
	    pxlVal = 0;
	  if(pxlVal > 255)
	    pxlVal = 255;
	  pDataDst[idx] = (unsigned char)pxlVal;
	  ++(rgHist[pxlVal]);
	}
	// TODO: when RGB, may want to convert to HSV, histogram
	// stretch the V channel, then convert back to RGB using the
	// stretched V instead of just stretching based on max/min of
	// any channel


	// find middle 96% of histogram to stretch over complete 0..255 range
	pxlCount = 0;
	bot = 0;
	while(pxlCount <= 0.02*len){
	  pxlCount += rgHist[bot];
	  ++bot;
	}
	--bot;
	pxlCount = 0;
	top = 255;
	while(pxlCount <= 0.02*len){
	  pxlCount += rgHist[top];
	  --top;
	}
	++top;
	if(top > bot){
	  // create the LUT to map the source value in histogram to dest value
	  for(int i = 0; i <= bot; ++i)
	    rgLUT[i] = 0;
	  for(int i = bot; i <= top; ++i)
	    rgLUT[i] = (i-bot)*255/(top-bot);
	  for(int i = top+1; i < 256; ++i)
	    rgLUT[i] = 255;
	  // map pixels to the stretched-histogram values
	  for(int i = 0; i < len; ++i){
	    pDataDst[i] = rgLUT[pDataDst[i]];
	  }
	}
      }
      break;
//     case DImage::DImage_u16:
//     case DImage::DImage_RGB_16:
//       {
	
//       }
//       break;
    default:
      fprintf(stderr, "DBackgroundRemover::removeBackground() unsupported "
	      "image type\n");
      exit(1);
      return;
  }
}






///Remove document background using Luke Hutchison's method at low resolution
/**If pImgBG is non-NULL, the DImage object pointed to will be filled
 * with the background image.  If radiusX is -1, then the default
 * window size will be used depending on the size of imgSrc.  If
 * radiusY is -1, it will be set equal to radiusX.  Only the
 * background is computed at low-res (shrink the image so that width
 * and height are between 100 and 800 before doing median filter, then
 * scale up the low-res background before difference).  If radiusX and
 * radiusY are specified, they are the actual radii used, even though
 * the image will be scaled before processing.  If no scaling is
 * required (i.e. the image is already small enough), then
 * removeBackground() will be called.
 */
void DBackgroundRemover::removeLowResBackground(DImage &imgDst,
						const DImage &imgSrc,
						DImage *pImgBG,
						int radiusX, int radiusY,
						DProgress *pProg,
						int numThreads,
						int bgMinW, int bgMinH,
						int bgMaxW, int bgMaxH){
  DImage imgBGtmp;
  DImage imgBGjaggy;
  DImage *pBG;
  DImage imgLowRes;
  DImage imgLowResBG;
  int len; // width * height * (times number of channels)
  int numHalves;
  int wLowRes, hLowRes;
  DKernel2D kern;


  numHalves = 0;
  wLowRes = imgSrc.width();
  hLowRes = imgSrc.height();
  while(((wLowRes > bgMinW) && (hLowRes > bgMinH)) &&
	((wLowRes > bgMaxW) || (hLowRes > bgMaxH))){
    wLowRes /= 2;
    hLowRes /= 2;
    ++numHalves;
  }
  if(0 == numHalves){
    removeBackground(imgDst, imgSrc, pImgBG, radiusX, radiusY, pProg,
		     numThreads);
    return;
  }

  printf("scaling down... numHalves=%d\n", numHalves);
  imgSrc.scaledDownPow2_(imgLowRes, numHalves);
  //  imgLowRes.save("/tmp/_imgLowRes.ppm");

  if(-1 == radiusX){
    radiusX =
      (imgLowRes.height() > imgLowRes.width()) ?
      (imgLowRes.height() * 3 / 200) : (imgLowRes.width() * 3 / 200);
  }
  if(-1 == radiusY){
    radiusY = radiusX;
  }

  DMedianFilter::medianFilterImage(imgLowResBG, imgLowRes,
				   false, radiusX,radiusY,
				   DMedianFilter::DMedFilt_default,pProg,
				   numThreads);
  if(NULL == pImgBG)
    pBG = &imgBGtmp;
  else
    pBG = pImgBG;
  //  imgLowResBG.save("/tmp/_imgLowResBG.pnm");

  
  imgLowResBG.scaledUpPow2_(imgBGjaggy, numHalves,
			    imgSrc.width(), imgSrc.height(),
			    DImage::DImageTransSample);
  //  imgBGjaggy.save("/tmp/_imgBGjaggy.pnm");
  
  // smooth the scaled bg before subtraction to get pid of pixelation
//   kern.setGauss(2,2);
  printf("smoothing background...\n");
  kern.setRect(numHalves, numHalves);
  DConvolver::convolve_(*pBG, imgBGjaggy, kern,
			false /*fAlreadyPadded*/,
			false /*fResize*/,
			true /*fConvertBack*/,
			false /*fUseDoublePrec*/,
			NULL);
  pBG->save("/tmp/_imgBGsmooth.pnm");
  printf("done smoothing background...\n");
  

  imgDst.create(imgSrc.width(), imgSrc.height(), imgSrc.getImageType(),
		imgSrc.numChannels(), imgSrc.getAllocMethod());
  len = imgSrc.width() * imgSrc.height() * imgSrc.numChannels();

  switch(imgSrc.getImageType()){
    case DImage::DImage_u8:
    case DImage::DImage_RGB:
      {
	D_uint8 *pDataSrc;
	D_uint8 *pDataBG;
	D_uint8 *pDataDst;
	int rgHist[256];
	int pxlVal;
	unsigned char rgLUT[256];
	int pxlCount;
	int top, bot;
	
	memset(rgHist, 0, sizeof(rgHist));
	
	pDataSrc = imgSrc.dataPointer_u8();
	pDataBG = pBG->dataPointer_u8();
	pDataDst = imgDst.dataPointer_u8();
	// do the background subtraction to leave a foreground image
	for(int idx = 0; idx < len; ++idx){
// 	  pxlVal = 255 - ((int)(unsigned int)pDataBG[idx] -
// 			  (int)(unsigned int)pDataSrc[idx]);
// 	  if(pxlVal < 0)
// 	    pxlVal = 0;
// 	  if(pxlVal > 255)
// 	    pxlVal = 255;
// 	  pDataDst[idx] = (unsigned char)pxlVal;
// 	  ++(rgHist[pxlVal]);



	  float f;

#if 1
	  pxlVal =  ((int)(unsigned int)pDataBG[idx] -
		     (int)(unsigned int)pDataSrc[idx]);
	  if(pxlVal < 0)
	    pxlVal = 0;
	  f = pxlVal / 255.;
	  pxlVal = (int)((f * (int)(unsigned int)pDataSrc[idx]) +
			 (1.-f)*255.);
#else // try division instead of subtraction
	  pxlVal = (0.0001+(float)(unsigned int)(pDataSrc[idx])) /
	    (0.0001+(float)(unsigned int)(pDataBG[idx])) * 255.;

#endif
	  if(pxlVal < 0)
	    pxlVal = 0;
	  if(pxlVal > 255)
	    pxlVal = 255;
	  pDataDst[idx] = (unsigned char)pxlVal;
	  ++(rgHist[pxlVal]);
	}
	// TODO: when RGB, may want to convert to HSV, histogram
	// stretch the V channel, then convert back to RGB using the
	// stretched V instead of just stretching based on max/min of
	// any channel


	// find middle 96% of histogram to stretch over complete 0..255 range
	pxlCount = 0;
	bot = 0;
	while(pxlCount <= 0.02*len){
	  pxlCount += rgHist[bot];
	  ++bot;
	}
	--bot;
	pxlCount = 0;
	top = 255;
	while(pxlCount <= 0.02*len){
	  pxlCount += rgHist[top];
	  --top;
	}
	++top;
	if(top > bot){
	  // create the LUT to map the source value in histogram to dest value
	  for(int i = 0; i <= bot; ++i)
	    rgLUT[i] = 0;
	  for(int i = bot; i <= top; ++i)
	    rgLUT[i] = (i-bot)*255/(top-bot);
	  for(int i = top+1; i < 256; ++i)
	    rgLUT[i] = 255;
	  // map pixels to the stretched-histogram values
	  for(int i = 0; i < len; ++i){
	    pDataDst[i] = rgLUT[pDataDst[i]];
	  }
	}
      }
      break;
//     case DImage::DImage_u16:
//     case DImage::DImage_RGB_16:
//       {
	
//       }
//       break;
    default:
      fprintf(stderr, "DBackgroundRemover::removeBackground() unsupported "
	      "image type\n");
      exit(1);
      return;
  }
}

