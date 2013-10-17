#include "dvariancefilter.h"

void DVarianceFilter::varianceImage_(DImage &imgDst, const DImage &imgSrc,
				     bool fAlreadyPadded,
				     int radiusX, int radiusY){
  DImage *pimgSrcPad;
  int ww, wh, wa; // window width, window height, window area
  int w, h;
  double *pDst;
  
  w = imgSrc.width();
  h = imgSrc.height();

  if(radiusX == -1){
    radiusX = (imgSrc.height() > imgSrc.width()) ?
      (imgSrc.height() * 3 / 200) : (imgSrc.width() * 3 / 200);
  }
  if(radiusY == -1){
    radiusY = (imgSrc.height() > imgSrc.width()) ?
      (imgSrc.height() * 3 / 200) : (imgSrc.width() * 3 / 200);
  }
  
  if(fAlreadyPadded){
    pimgSrcPad = (DImage*)&imgSrc;
  }
  else{
    pimgSrcPad = new DImage;
    D_CHECKPTR(pimgSrcPad);
    imgSrc.padEdges_(*pimgSrcPad, radiusX+1, radiusX, radiusY+1, radiusY,
		     DImage::DImagePadReplicate);
  }

  w = pimgSrcPad->width();
  h = pimgSrcPad->height();
  imgDst.create(imgSrc.width(), imgSrc.height(),
		DImage::DImage_dbl_multi, imgSrc.numChannels());
  
  ww = 2*radiusX+1;
  wh = 2*radiusY+1;
  wa = ww * wh;

  switch(imgSrc.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint64 *pIntegralImg64;//integral image of values
	D_uint64 *pIntegral2Img64;//integral image of squared values
	D_uint8 *pu8SrcPad;
	pIntegralImg64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h);
	D_CHECKPTR(pIntegralImg64);
	pIntegral2Img64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h);
	D_CHECKPTR(pIntegral2Img64);
	pDst = imgDst.dataPointer_dbl();
	pu8SrcPad = pimgSrcPad->dataPointer_u8();
	
	// fill in the integral image values
	pIntegralImg64[0] = (D_uint64)(pu8SrcPad[0]);
	pIntegral2Img64[0] = (D_uint64)(pu8SrcPad[0])*(D_uint64)(pu8SrcPad[0]);
	for(int x = 1; x < w; ++x){ // first row
	  pIntegralImg64[x] = pIntegralImg64[x-1]+(D_uint64)(pu8SrcPad[x]);
	  pIntegral2Img64[x] = pIntegral2Img64[x-1]+
	    (D_uint64)(pu8SrcPad[x])*(D_uint64)(pu8SrcPad[x]);
	}
	for(int y = 1, idx=w; y < h; ++y){
	  pIntegralImg64[idx] = pIntegralImg64[idx-w] +
	    (D_uint64)(pu8SrcPad[idx]);
	  pIntegral2Img64[idx] = pIntegral2Img64[idx-w] +
	    (D_uint64)(pu8SrcPad[idx])*(D_uint64)(pu8SrcPad[idx]);
	  ++idx;
	  for(int x = 1; x < w; ++x, ++idx){
	    pIntegralImg64[idx] = pIntegralImg64[idx-1]+
	      pIntegralImg64[idx-w] + (D_uint64)(pu8SrcPad[idx]) -
	      pIntegralImg64[idx-w-1];
	    pIntegral2Img64[idx] = pIntegral2Img64[idx-1]+
	      pIntegral2Img64[idx-w] +
	      (D_uint64)(pu8SrcPad[idx])*(D_uint64)(pu8SrcPad[idx]) -
	      pIntegral2Img64[idx-w-1];
	  }
	}
	// now calculate the variance at each point using the integral images
	for(int y = radiusY+1, idxDst = 0; y < (h-radiusY); ++y){
	  int idxA, idxB, idxC, idxD;
	  idxB = (y-radiusY-1) * w;
	  idxC = idxB + ww ;
	  idxD = (y+radiusY) * w;
	  idxA = idxD + ww ;
	  
	  for(int x = radiusX+1; x < (w-radiusX); ++x, ++idxDst){
	    double mean;
	    mean = (double)
	      ((pIntegralImg64[idxA] + pIntegralImg64[idxB]-
		pIntegralImg64[idxC] - pIntegralImg64[idxD]) / (double)wa);
	    pDst[idxDst] = (double)
	      (((pIntegral2Img64[idxA] + pIntegral2Img64[idxB]-
		 pIntegral2Img64[idxC] - pIntegral2Img64[idxD]) / (double)wa) -
	       mean*mean);
	    ++idxA;
	    ++idxB;
	    ++idxC;
	    ++idxD;
	  }//end for(x...
	}//end for(y...
	free(pIntegralImg64);
	free(pIntegral2Img64);
      }
      break;
    case DImage::DImage_u16:
      {
	D_uint64 *pIntegralImg64;
	D_uint64 *pIntegral2Img64;
	D_uint16 *pu16SrcPad;
	pIntegralImg64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h);
	D_CHECKPTR(pIntegralImg64);
	pIntegral2Img64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h);
	D_CHECKPTR(pIntegral2Img64);
	pDst = imgDst.dataPointer_dbl();
	pu16SrcPad = pimgSrcPad->dataPointer_u16();
	
	// fill in the integral image values
	pIntegralImg64[0] = (D_uint64)(pu16SrcPad[0]);
	pIntegral2Img64[0] =
	  (D_uint64)(pu16SrcPad[0])*(D_uint64)(pu16SrcPad[0]);
	for(int x = 1; x < w; ++x){ // first row
	  pIntegralImg64[x] = pIntegralImg64[x-1]+(D_uint64)(pu16SrcPad[x]);
	  pIntegral2Img64[x] = pIntegral2Img64[x-1]+
	    (D_uint64)(pu16SrcPad[x])*(D_uint64)(pu16SrcPad[x]);
	}
	for(int y = 1, idx=w; y < h; ++y){
	  pIntegralImg64[idx] = pIntegralImg64[idx-w] +
	    (D_uint64)(pu16SrcPad[idx]);
	  pIntegral2Img64[idx] = pIntegral2Img64[idx-w] +
	    (D_uint64)(pu16SrcPad[idx])*(D_uint64)(pu16SrcPad[idx]);
	  ++idx;
	  for(int x = 1; x < w; ++x, ++idx){
	    pIntegralImg64[idx] = pIntegralImg64[idx-1]+
	      pIntegralImg64[idx-w] + (D_uint64)(pu16SrcPad[idx]) -
	      pIntegralImg64[idx-w-1];
	    pIntegral2Img64[idx] = pIntegral2Img64[idx-1]+
	      pIntegral2Img64[idx-w] +
	      (D_uint64)(pu16SrcPad[idx])*(D_uint64)(pu16SrcPad[idx]) -
	      pIntegral2Img64[idx-w-1];
	  }
	}
	// now calculate the mean at each point using the integral image
	for(int y = radiusY+1, idxDst = 0; y < (h-radiusY); ++y){
	  int idxA, idxB, idxC, idxD;
	  idxB = (y-radiusY-1) * w;
	  idxC = idxB + ww ;
	  idxD = (y+radiusY) * w;
	  idxA = idxD + ww ;
	  
	  for(int x = radiusX+1; x < (w-radiusX); ++x, ++idxDst){
	    double mean;
	    mean = (double)
	      ((pIntegralImg64[idxA] + pIntegralImg64[idxB]-
		pIntegralImg64[idxC] - pIntegralImg64[idxD]) / (double)wa);
	    pDst[idxDst] = (double)
	      (((pIntegral2Img64[idxA] + pIntegral2Img64[idxB]-
		 pIntegral2Img64[idxC] - pIntegral2Img64[idxD]) / (double)wa) -
	       mean*mean);
	    ++idxA;
	    ++idxB;
	    ++idxC;
	    ++idxD;
	  }//end for(x...
	}//end for(y...
	free(pIntegralImg64);
	free(pIntegral2Img64);
      }
      break;
    case DImage::DImage_u32:
      {
	D_uint64 *pIntegralImg64;
	D_uint64 *pIntegral2Img64;
	D_uint32 *pu32SrcPad;
	pIntegralImg64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h);
	D_CHECKPTR(pIntegralImg64);
	pIntegral2Img64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h);
	D_CHECKPTR(pIntegral2Img64);
	pDst = imgDst.dataPointer_dbl();
	pu32SrcPad = pimgSrcPad->dataPointer_u32();
	
	// fill in the integral image values
	pIntegralImg64[0] = (D_uint64)(pu32SrcPad[0]);
	pIntegral2Img64[0] =
	  (D_uint64)(pu32SrcPad[0])*(D_uint64)(pu32SrcPad[0]);
	for(int x = 1; x < w; ++x){ // first row
	  pIntegralImg64[x] = pIntegralImg64[x-1]+(D_uint64)(pu32SrcPad[x]);
	  pIntegral2Img64[x] = pIntegral2Img64[x-1]+
	    (D_uint64)(pu32SrcPad[x])*(D_uint64)(pu32SrcPad[x]);
	}
	for(int y = 1, idx=w; y < h; ++y){
	  pIntegralImg64[idx] = pIntegralImg64[idx-w] +
	    (D_uint64)(pu32SrcPad[idx]);
	  pIntegral2Img64[idx] = pIntegral2Img64[idx-w] +
	    (D_uint64)(pu32SrcPad[idx])*(D_uint64)(pu32SrcPad[idx]);
	  ++idx;
	  for(int x = 1; x < w; ++x, ++idx){
	    pIntegralImg64[idx] = pIntegralImg64[idx-1]+
	      pIntegralImg64[idx-w] + (D_uint64)(pu32SrcPad[idx]) -
	      pIntegralImg64[idx-w-1];
	    pIntegral2Img64[idx] = pIntegral2Img64[idx-1]+
	      pIntegral2Img64[idx-w] +
	      (D_uint64)(pu32SrcPad[idx])*(D_uint64)(pu32SrcPad[idx]) -
	      pIntegral2Img64[idx-w-1];
	  }
	}
	// now calculate the mean at each point using the integral image
	for(int y = radiusY+1, idxDst = 0; y < (h-radiusY); ++y){
	  int idxA, idxB, idxC, idxD;
	  idxB = (y-radiusY-1) * w;
	  idxC = idxB + ww ;
	  idxD = (y+radiusY) * w;
	  idxA = idxD + ww ;
	  
	  for(int x = radiusX+1; x < (w-radiusX); ++x, ++idxDst){
	    double mean;
	    mean = (double)
	      ((pIntegralImg64[idxA] + pIntegralImg64[idxB]-
		pIntegralImg64[idxC] - pIntegralImg64[idxD]) / (double)wa);
	    pDst[idxDst] = (double)
	      (((pIntegral2Img64[idxA] + pIntegral2Img64[idxB]-
		 pIntegral2Img64[idxC] - pIntegral2Img64[idxD]) / (double)wa) -
	       mean*mean);
	    ++idxA;
	    ++idxB;
	    ++idxC;
	    ++idxD;
	  }//end for(x...
	}//end for(y...
	free(pIntegralImg64);
	free(pIntegral2Img64);
      }
      break;
    case DImage::DImage_RGB:
      {
	D_uint64 *pIntegralImg64;
	D_uint64 *pIntegral2Img64;
	D_uint8 *pu8SrcPad;
	DImage *rgImgSrcChannels;

	rgImgSrcChannels = new DImage[imgSrc.numChannels()];
	D_CHECKPTR(rgImgSrcChannels);
	pIntegralImg64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h) ;
	D_CHECKPTR(pIntegralImg64);
	pIntegral2Img64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h) ;
	D_CHECKPTR(pIntegral2Img64);

	pimgSrcPad->splitRGB(rgImgSrcChannels[0], rgImgSrcChannels[1],
			     rgImgSrcChannels[2]);

	for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
	  pu8SrcPad = rgImgSrcChannels[chan].dataPointer_u8();
	  pDst = imgDst.dataPointer_dbl(chan);
	
	  // fill in the integral image values
	  pIntegralImg64[0] = (D_uint64)(pu8SrcPad[0]);
	  pIntegral2Img64[0] =
	    (D_uint64)(pu8SrcPad[0])*(D_uint64)(pu8SrcPad[0]);
	  for(int x = 1; x < w; ++x){ // first row
	    pIntegralImg64[x] = pIntegralImg64[x-1]+
	      (D_uint64)(pu8SrcPad[x]);
	    pIntegral2Img64[x] = pIntegral2Img64[x-1]+
	      (D_uint64)(pu8SrcPad[x])*(D_uint64)(pu8SrcPad[x]);
	  }
	  for(int y = 1, idx=w; y < h; ++y){
	    pIntegralImg64[idx] = pIntegralImg64[idx-w] +
	      (D_uint64)(pu8SrcPad[idx]);
	    pIntegral2Img64[idx] = pIntegral2Img64[idx-w] +
	      (D_uint64)(pu8SrcPad[idx])*(D_uint64)(pu8SrcPad[idx]);
	    ++idx;
	    for(int x = 1; x < w; ++x, ++idx){
	      pIntegralImg64[idx] = pIntegralImg64[idx-1]+
		pIntegralImg64[idx-w] + (D_uint64)(pu8SrcPad[idx]) -
		pIntegralImg64[idx-w-1];
	      pIntegral2Img64[idx] = pIntegral2Img64[idx-1]+
		pIntegral2Img64[idx-w] +
		(D_uint64)(pu8SrcPad[idx])*(D_uint64)(pu8SrcPad[idx]) -
		pIntegral2Img64[idx-w-1];
	    }
	  }
	  // now calculate the mean at each point using the integral image
	  for(int y = radiusY+1, idxDst = 0; y < (h-radiusY); ++y){
	    int idxA, idxB, idxC, idxD;
	    idxB = (y-radiusY-1) * w;
	    idxC = idxB + ww ;
	    idxD = (y+radiusY) * w;
	    idxA = idxD + ww ;
	    
	    for(int x = radiusX+1; x < (w-radiusX); ++x, ++idxDst){
	      double mean;
	      mean = (double)
		((pIntegralImg64[idxA] + pIntegralImg64[idxB]-
		  pIntegralImg64[idxC] - pIntegralImg64[idxD]) / (double)wa);
	      pDst[idxDst] = (double)
		(((pIntegral2Img64[idxA]+pIntegral2Img64[idxB]-
		   pIntegral2Img64[idxC]-pIntegral2Img64[idxD]) / (double)wa) -
		 mean*mean);
	      ++idxA;
	      ++idxB;
	      ++idxC;
	      ++idxD;
	    }//end for(x...
	  }//end for(y...
	}//end for(chan...
	free(pIntegralImg64);
	free(pIntegral2Img64);
	delete [] rgImgSrcChannels;
      }
      break;
    case DImage::DImage_RGB_16:
      {
	D_uint64 *pIntegralImg64;
	D_uint64 *pIntegral2Img64;
	D_uint16 *pu16SrcPad;
	DImage *rgImgSrcChannels;

	rgImgSrcChannels = new DImage[imgSrc.numChannels()];
	D_CHECKPTR(rgImgSrcChannels);
	pIntegralImg64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h) ;
	D_CHECKPTR(pIntegralImg64);
	pIntegral2Img64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h) ;
	D_CHECKPTR(pIntegral2Img64);

	pimgSrcPad->splitRGB(rgImgSrcChannels[0], rgImgSrcChannels[1],
			     rgImgSrcChannels[2]);

	for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
	  pu16SrcPad = rgImgSrcChannels[chan].dataPointer_u16();
	  pDst = imgDst.dataPointer_dbl(chan);
	
	  // fill in the integral image values
	  pIntegralImg64[0] = (D_uint64)(pu16SrcPad[0]);
	  pIntegral2Img64[0] =
	    (D_uint64)(pu16SrcPad[0])*(D_uint64)(pu16SrcPad[0]);
	  for(int x = 1; x < w; ++x){ // first row
	    pIntegralImg64[x] = pIntegralImg64[x-1]+
	      (D_uint64)(pu16SrcPad[x]);
	    pIntegral2Img64[x] = pIntegral2Img64[x-1]+
	      (D_uint64)(pu16SrcPad[x])*(D_uint64)(pu16SrcPad[x]);
	  }
	  for(int y = 1, idx=w; y < h; ++y){
	    pIntegralImg64[idx] = pIntegralImg64[idx-w] +
	      (D_uint64)(pu16SrcPad[idx]);
	    pIntegral2Img64[idx] = pIntegral2Img64[idx-w] +
	      (D_uint64)(pu16SrcPad[idx])*(D_uint64)(pu16SrcPad[idx]);
	    ++idx;
	    for(int x = 1; x < w; ++x, ++idx){
	      pIntegralImg64[idx] = pIntegralImg64[idx-1]+
		pIntegralImg64[idx-w] + (D_uint64)(pu16SrcPad[idx]) -
		pIntegralImg64[idx-w-1];
	      pIntegral2Img64[idx] = pIntegral2Img64[idx-1]+
		pIntegral2Img64[idx-w] +
		(D_uint64)(pu16SrcPad[idx])*(D_uint64)(pu16SrcPad[idx]) -
		pIntegral2Img64[idx-w-1];
	    }
	  }
	  // now calculate the mean at each point using the integral image
	  for(int y = radiusY+1, idxDst = 0; y < (h-radiusY); ++y){
	    int idxA, idxB, idxC, idxD;
	    idxB = (y-radiusY-1) * w;
	    idxC = idxB + ww ;
	    idxD = (y+radiusY) * w;
	    idxA = idxD + ww ;
	    
	    for(int x = radiusX+1; x < (w-radiusX); ++x, ++idxDst){
	      double mean;
	      mean = (double)
		((pIntegralImg64[idxA] + pIntegralImg64[idxB]-
		  pIntegralImg64[idxC] - pIntegralImg64[idxD]) / (double)wa);
	      pDst[idxDst] = (double)
		(((pIntegral2Img64[idxA]+pIntegral2Img64[idxB]-
		   pIntegral2Img64[idxC]-pIntegral2Img64[idxD]) / (double)wa) -
		 mean*mean);
	      ++idxA;
	      ++idxB;
	      ++idxC;
	      ++idxD;
	    }//end for(x...
	  }//end for(y...
	}//end for(chan...
	free(pIntegralImg64);
	free(pIntegral2Img64);
	delete [] rgImgSrcChannels;
      }
      break;
    case DImage::DImage_flt_multi:
      {
	long double *pIntegralImgDbl;
	long double *pIntegral2ImgDbl;
	float *pFltDst;
	float *pFltSrcPad;

	pIntegralImgDbl = (long double*)malloc(sizeof(long double) * w * h) ;
	D_CHECKPTR(pIntegralImgDbl);
	pIntegral2ImgDbl = (long double*)malloc(sizeof(long double) * w * h) ;
	D_CHECKPTR(pIntegral2ImgDbl);

	for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
	  pFltDst = imgDst.dataPointer_flt(chan);
	  pFltSrcPad = pimgSrcPad->dataPointer_flt(chan);
	  
	  // fill in the integral image values
	  pIntegralImgDbl[0] = (long double)(pFltSrcPad[0]);
	  pIntegral2ImgDbl[0] =
	    (long double)(pFltSrcPad[0])*(long double)(pFltSrcPad[0]);
	  for(int x = 1; x < w; ++x){ // first row
	    pIntegralImgDbl[x] = pIntegralImgDbl[x-1]+
	      (long double)(pFltSrcPad[x]);
	    pIntegral2ImgDbl[x] = pIntegral2ImgDbl[x-1]+
	      (long double)(pFltSrcPad[x])*(long double)(pFltSrcPad[x]);
	  }
	  for(int y = 1, idx=w; y < h; ++y){
	    pIntegralImgDbl[idx] = pIntegralImgDbl[idx-w] +
	      (long double)(pFltSrcPad[idx]);
	    pIntegral2ImgDbl[idx] = pIntegral2ImgDbl[idx-w] +
	      (long double)(pFltSrcPad[idx])*(long double)(pFltSrcPad[idx]);
	    ++idx;
	    for(int x = 1; x < w; ++x, ++idx){
	      pIntegralImgDbl[idx] = pIntegralImgDbl[idx-1]+
		pIntegralImgDbl[idx-w] + (long double)(pFltSrcPad[idx]) -
		pIntegralImgDbl[idx-w-1];
	      pIntegral2ImgDbl[idx] = pIntegral2ImgDbl[idx-1]+
		pIntegral2ImgDbl[idx-w] +
		(long double)(pFltSrcPad[idx])*(long double)(pFltSrcPad[idx]) -
		pIntegral2ImgDbl[idx-w-1];
	    }
	  }
	  // now calculate the mean at each point using the integral image
	  for(int y = radiusY+1, idxDst = 0; y < (h-radiusY); ++y){
	    int idxA, idxB, idxC, idxD;
	    idxB = (y-radiusY-1) * w;
	    idxC = idxB + ww ;
	    idxD = (y+radiusY) * w;
	    idxA = idxD + ww ;
	    
	    for(int x = radiusX+1; x < (w-radiusX); ++x, ++idxDst){
	      long double mean;
	      mean = (long double)
		((pIntegralImgDbl[idxA] + pIntegralImgDbl[idxB]-
		  pIntegralImgDbl[idxC] - pIntegralImgDbl[idxD]) /
		 (long double)wa);
	      pFltDst[idxDst] = (double)
		(((pIntegral2ImgDbl[idxA]+pIntegral2ImgDbl[idxB]-
		   pIntegral2ImgDbl[idxC]-pIntegral2ImgDbl[idxD]) /
		  (long double)wa) -
		 mean*mean);
	      ++idxA;
	      ++idxB;
	      ++idxC;
	      ++idxD;
	    }//end for(x...
	  }//end for(y...
	}//end for(chan...
	free(pIntegralImgDbl);
	free(pIntegral2ImgDbl);
      }
      break;
    case DImage::DImage_dbl_multi:
      {
	long double *pIntegralImgDbl;
	long double *pIntegral2ImgDbl;
	double *pDblDst;
	double *pDblSrcPad;

	pIntegralImgDbl = (long double*)malloc(sizeof(long double) * w * h) ;
	D_CHECKPTR(pIntegralImgDbl);
	pIntegral2ImgDbl = (long double*)malloc(sizeof(long double) * w * h) ;
	D_CHECKPTR(pIntegral2ImgDbl);

	for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
	  pDblDst = imgDst.dataPointer_dbl(chan);
	  pDblSrcPad = pimgSrcPad->dataPointer_dbl(chan);
	  
	  // fill in the integral image values
	  pIntegralImgDbl[0] = (long double)(pDblSrcPad[0]);
	  pIntegral2ImgDbl[0] =
	    (long double)(pDblSrcPad[0])*(long double)(pDblSrcPad[0]);
	  for(int x = 1; x < w; ++x){ // first row
	    pIntegralImgDbl[x] = pIntegralImgDbl[x-1]+
	      (long double)(pDblSrcPad[x]);
	    pIntegral2ImgDbl[x] = pIntegral2ImgDbl[x-1]+
	      (long double)(pDblSrcPad[x])*(long double)(pDblSrcPad[x]);
	  }
	  for(int y = 1, idx=w; y < h; ++y){
	    pIntegralImgDbl[idx] = pIntegralImgDbl[idx-w] +
	      (long double)(pDblSrcPad[idx]);
	    pIntegral2ImgDbl[idx] = pIntegral2ImgDbl[idx-w] +
	      (long double)(pDblSrcPad[idx])*(long double)(pDblSrcPad[idx]);
	    ++idx;
	    for(int x = 1; x < w; ++x, ++idx){
	      pIntegralImgDbl[idx] = pIntegralImgDbl[idx-1]+
		pIntegralImgDbl[idx-w] + (long double)(pDblSrcPad[idx]) -
		pIntegralImgDbl[idx-w-1];
	      pIntegral2ImgDbl[idx] = pIntegral2ImgDbl[idx-1]+
		pIntegral2ImgDbl[idx-w] +
		(long double)(pDblSrcPad[idx])*(long double)(pDblSrcPad[idx]) -
		pIntegral2ImgDbl[idx-w-1];
	    }
	  }
	  // now calculate the mean at each point using the integral image
	  for(int y = radiusY+1, idxDst = 0; y < (h-radiusY); ++y){
	    int idxA, idxB, idxC, idxD;
	    idxB = (y-radiusY-1) * w;
	    idxC = idxB + ww ;
	    idxD = (y+radiusY) * w;
	    idxA = idxD + ww ;
	    
	    for(int x = radiusX+1; x < (w-radiusX); ++x, ++idxDst){
	      long double mean;
	      mean = (long double)
		((pIntegralImgDbl[idxA] + pIntegralImgDbl[idxB]-
		  pIntegralImgDbl[idxC] - pIntegralImgDbl[idxD]) /
		 (long double)wa);
	      pDblDst[idxDst] = (double)
		(((pIntegral2ImgDbl[idxA]+pIntegral2ImgDbl[idxB]-
		   pIntegral2ImgDbl[idxC]-pIntegral2ImgDbl[idxD]) /
		  (long double)wa) -
		 mean*mean);
	      ++idxA;
	      ++idxB;
	      ++idxC;
	      ++idxD;
	    }//end for(x...
	  }//end for(y...
	}//end for(chan...
	free(pIntegralImgDbl);
	free(pIntegral2ImgDbl);
      }
      break;
    default:
      fprintf(stderr, "DVarianceFilter::varianceImage_() not "
	      "implemented for complex images\n");
      abort();
  }

  if(!fAlreadyPadded){
    delete pimgSrcPad;
  }
}


///Convenience function that takes square roots after calling varianceImage_()
/**The standard deviation is just the square rott of the
 * variance. This function calls varianceImage_() and then takes the
 * square root at each pixel.  The result is stored in imgDst.*/
void DVarianceFilter::standardDeviationImage_(DImage &imgDst,
					      const DImage &imgSrc,
					      bool fAlreadyPadded,
					      int radiusX, int radiusY){
  int w, h;
  double *pDst;
  varianceImage_(imgDst, imgSrc, fAlreadyPadded, radiusX, radiusY);
  w = imgDst.width();
  h = imgDst.height();

  for(int chan = 0; chan < imgDst.numChannels(); ++chan){
    pDst = imgDst.dataPointer_dbl(chan);
    for(int y = 0, idx = 0; y < h; ++y){
      for(int x = 0; x < w; ++x, ++idx){
	pDst[idx] = sqrt(pDst[idx]);
      }
    }
  }
}












