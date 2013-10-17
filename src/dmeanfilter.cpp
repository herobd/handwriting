#include "dmeanfilter.h"
#include "dimage.h"
#include "dprogress.h"
#include "dinstancecounter.h"
#include <string.h>


///Default constructor
DMeanFilter::DMeanFilter(){
  DInstanceCounter::addInstance("DMeanFilter");
  _meanFiltType = DMeanFilt_circle;
  rgKern = NULL;
  rgRightEdge = NULL;
  _radiusX = 1;
  _radiusY = 1;
}


///Destructor
DMeanFilter::~DMeanFilter(){
  DInstanceCounter::removeInstance("DMeanFilter");
  if(rgKern){
    free(rgKern);
    rgKern = NULL;
  }
  if(rgRightEdge){
    free(rgRightEdge);
    rgRightEdge = NULL;
  }
}

///Set the kernel (window) radii for the mean filter
/**The actual width of the window will be 2*radiusX+1*/
void DMeanFilter::setRadii(int radiusX, int radiusY){
  if((radiusX < 0) || (radiusY < 0)){
    fprintf(stderr, "DMeanFilter::setRadii() negative value not allowed\n");
    return;
  }
  _radiusX = radiusX;
  _radiusY = radiusY;
  if(rgKern){
    free(rgKern);
    rgKern = NULL;
    free(rgRightEdge);
    rgRightEdge = NULL;
  }
}

///set the type of mean filter algorithm
/**By default, RGB and _u8 images will use Huang_circle and all others
 * (including RGB_16 and _u16) will use the slow method */
void DMeanFilter::setType(DMeanFiltType filtType){
  if(rgKern){
    free(rgKern);
    rgKern = NULL;
    free(rgRightEdge);
    rgRightEdge = NULL;
  }
  _meanFiltType = filtType;
}

///return the type of the mean filter
DMeanFilter::DMeanFiltType DMeanFilter::getType(){
  return _meanFiltType;
}

///private function that actually performs mean filter on 8-bit images
/** imgDst will be 2*radiusX pixels less wide and 2*radiusY pixels
 * less high than imgSrc because of the padding that is added before
 * calling this function.  This function requires that imgDst.create()
 * has already been called with the proper w,h,imgType,etc.
 */
void DMeanFilter::meanFilt_u8(DImage &imgDst, const DImage &imgSrc,
			      int radiusX, int radiusY,
			      D_uint8 *rgKern,
			      int numKernPxls, int *rgRightEdge,
			      DProgress *pProg,
			      int progStart, int progMax){
  int mean;
  int idxDst;
  int idx3;
  D_uint8 *pTmp; // pointer to padded image data
  int wTmp, hTmp; // width, height of imgSrc
  int w, h; // width, height of imgDst
  D_uint8 *pDst;
  int wKern, hKern;
  unsigned int sum;

  wKern = radiusX * 2 + 1;
  hKern = radiusY * 2 + 1;

  wTmp = imgSrc.width();
  hTmp = imgSrc.height();
  w = wTmp - radiusX*2;
  h = hTmp - radiusY*2;
  pDst = imgDst.dataPointer_u8();
  pTmp = imgSrc.dataPointer_u8();

  for(int y = 0; y < h; ++y){
    // update progress report and check if user cancelled the operation
    if((NULL != pProg) && (0 == (y & 0x0000003f))){
      if(0 != pProg->reportStatus(progStart + y, 0, progMax)){
	// the operation has been cancelled
	pProg->reportStatus(-1, 0, progMax); // report cancel acknowledged
	return;
      }
    }

    sum = 0;

    // position window at the beginning of a new row and fill the kernel
    for(int kr = 0, kidx =0; kr < hKern; ++kr){
      for(int kc = 0; kc < wKern; ++kc, ++kidx){
	if(rgKern[kidx]){ // pixel is part of the kernel mask
	  sum += pTmp[(y+kr)*wTmp+kc];//add pixel val to sum
	}
      }
    }
    // calculate mean for first spot
    mean = sum / numKernPxls;

    // put the mean in the spot we're at
    idxDst = y*w;
    pDst[idxDst] = (unsigned char)mean;
    
    // remove pixels from leftmost column
    idx3 = y*wTmp+radiusX;
    for(int ky = 0; ky < hKern; ++ky){
      sum -= pTmp[idx3 - rgRightEdge[ky]];
      idx3 += wTmp;
    }

    for(int x=1;  x < w;  ++x){
      ++idxDst;
      // add pixels from the right-hand side of kernel (after moving over one)
      idx3 = y*wTmp+x+radiusX;
      for(int ky = 0; ky < hKern; ++ky){
	sum += pTmp[idx3 + rgRightEdge[ky]];
	idx3 += wTmp;
      }

      // find mean
      mean = sum / numKernPxls;
      
      // put the mean value in the destination pixel
      pDst[idxDst] = (unsigned char)mean;

      // remove pixels from leftmost column for next time through loop
      if(x < (w-1)){//don't need to remove left edge if going to a new row
	idx3 = y*wTmp+x+radiusX;
	for(int ky = 0; ky < hKern; ++ky){
	  sum -= pTmp[idx3 - rgRightEdge[ky]];
	  idx3 += wTmp;
	} // end for(ky...
      } // end if

    } // end for (x=1; ...

  }// end for(y=0...
  // report progress
  if(NULL != pProg){
    pProg->reportStatus(progStart + h, 0, progMax);
  }
}


///Mean filter imgSrc with current settings and return the result
DImage DMeanFilter::filterImage(const DImage &imgSrc,
				  bool fAlreadyPadded, DProgress *pProg){
  DImage dst;
  filterImage_(dst, imgSrc, fAlreadyPadded, pProg);
  return dst;
}

///Mean filter imgSrc with current settings and store result in imgDst
void DMeanFilter::filterImage_(DImage &imgDst, const DImage &imgSrc,
				 bool fAlreadyPadded, DProgress *pProg){
  DMeanFiltType filtType;
  DImage *pImgPad;
  int wKern, hKern;
  int numKernPxls;
  int wUnpad, hUnpad;


  filtType = _meanFiltType;

  pImgPad = (DImage*)&imgSrc;
  if(!fAlreadyPadded){
    pImgPad = new DImage();
    imgSrc.padEdges_(*pImgPad, _radiusX, _radiusX, _radiusY, _radiusY,
		     DImage::DImagePadReplicate);
  }

  wUnpad = pImgPad->width()-(_radiusX*2);
  hUnpad = pImgPad->height()-(_radiusY*2);

  wKern = _radiusX * 2 + 1;
  hKern = _radiusY * 2 + 1;
  if(NULL == rgKern){
    rgKern = (unsigned char*)malloc(sizeof(unsigned char) * wKern * hKern);
    if(!rgKern){
      fprintf(stderr, "DMeanFilter::filterImage_() out of memory\n");
      exit(1);
    }
    rgRightEdge = (int*)malloc(sizeof(int)*hKern);
    if(!rgRightEdge){
      fprintf(stderr, "DMeanFilter::filterImage_() out of memory\n");
      exit(1);
    }
    if(DMeanFilt_circle == filtType){
      fill_circle_kern_offsets(_radiusX, _radiusY, rgKern,
			       rgRightEdge, &numKernPxls);
    }
    else{
      fill_square_kern_offsets(_radiusX, _radiusY, rgKern,
			       rgRightEdge, &numKernPxls);
    }
  }
  
  switch(imgSrc.getImageType()){
    case DImage::DImage_u8:
      {
	imgDst.create(wUnpad, hUnpad, DImage::DImage_u8, 1,
		      imgSrc.getAllocMethod());
	meanFilt_u8(imgDst, *pImgPad, _radiusX, _radiusY, rgKern,
		    numKernPxls, rgRightEdge, pProg, 0, hUnpad+1);
	if(NULL != pProg){
	  pProg->reportStatus(hUnpad+1, 0, hUnpad+1);//report progress complete
	}	  
      }
      break;
    case DImage::DImage_RGB:
      {
	DImage imgR, imgG, imgB;
	DImage imgRDst, imgGDst, imgBDst;

	imgRDst.create(wUnpad, hUnpad, DImage::DImage_u8, 1,
		       imgSrc.getAllocMethod());
	imgGDst.create(wUnpad, hUnpad, DImage::DImage_u8, 1,
		       imgSrc.getAllocMethod());
	imgBDst.create(wUnpad, hUnpad, DImage::DImage_u8, 1,
		       imgSrc.getAllocMethod());

	pImgPad->splitRGB(imgR, imgG, imgB);

	meanFilt_u8(imgRDst, imgR, _radiusX,_radiusY, rgKern, numKernPxls,
		    rgRightEdge, pProg, 0, 3 * hUnpad);
	meanFilt_u8(imgGDst, imgG, _radiusX,_radiusY, rgKern, numKernPxls,
		    rgRightEdge, pProg, hUnpad, 3 * hUnpad);
	meanFilt_u8(imgBDst, imgB, _radiusX,_radiusY, rgKern, numKernPxls,
		    rgRightEdge, pProg, 2 * hUnpad, 3 * hUnpad+1);
	if(NULL != pProg){
	  pProg->reportStatus(3*hUnpad+1,0,3*hUnpad+1);//report complete
	}	  
	imgDst.combineRGB(imgRDst, imgGDst, imgBDst);
      }
      break;
    default:
      //TODO: finish implementing this
      fprintf(stderr, "DMeanFilter::filterImage_() not yet implemented for "
	      "some image types\n");
      exit(1);
      break;
  }
  if(!fAlreadyPadded){
    delete pImgPad;
    pImgPad = NULL;
  }

}


///Fills rgKernFull with a circle and sets rgRightEdge and numKernPxls
/**Pixels that are part of the circle will be 0xff and the rest will
 * be 0x00.  rgRightEdge is an array that says how many pixels to the
 * right the leading edge of the kernel is (the right-most pixel), and
 * it is assumed that the trailing edge is the same number of pixels
 * to the left. numKernPxls is a count of how many pixels in the
 * kernel are set to 0xff.
 */
void DMeanFilter::fill_circle_kern_offsets(int radiusX, int radiusY,
					   unsigned char *rgKernFull,
					   int *rgRightEdge,
					   int *numKernPxls){
  int kernWidth;
  int radiusSquared;
  int kySquared;

  (*numKernPxls)= 0;
  kernWidth = radiusX*2+1;
  radiusSquared = (int)((radiusX+0.5)*(radiusX+0.5));

  for(int ky = 0-radiusY, idx = 0; ky <= radiusY; ++ky){
    kySquared = ky*ky;
    for(int kx = 0-radiusX; kx <= radiusX; ++kx, ++idx){
      if( (kx*kx+kySquared) <= radiusSquared){
	rgKernFull[idx] = 0xff;
	rgRightEdge[ky+radiusX] = kx;
	++(*numKernPxls);
      }
      else
	rgKernFull[idx] = 0;
    }
  }
}

///Fills rgKernFull with a square and sets rgRightEdge and numKernPxls
/**Similar to fill_circle_kern_offsets() but all of the pixels are set
 * to 0xff, in effect, forming a square kernel that is (radius*2+1)
 * wide and high.
 */
void DMeanFilter::fill_square_kern_offsets(int radiusX, int radiusY,
					     unsigned char *rgKernFull,
					     int *rgRightEdge,
					     int *numKernPxls){
  int kernWidth;
  int kernHeight;
  kernWidth = radiusX*2+1;
  kernHeight = radiusY*2+1;
  
  memset(rgKernFull, 1, kernWidth*kernHeight);
  for(int ky = 0; ky < kernHeight; ++ky){
    rgRightEdge[ky] = radiusX;
  }
  (*numKernPxls) = kernWidth*kernHeight;
}

///Static function provided for convenience
void DMeanFilter::meanFilterImage(DImage &imgDst, const DImage &imgSrc,
				  bool fAlreadyPadded, int radiusX,int radiusY,
				  DMeanFiltType filtType, DProgress *pProg){
  DMeanFilter mfilt;
  mfilt.setRadii(radiusX, radiusY);
  mfilt.setType(filtType);
  mfilt.filterImage_(imgDst, imgSrc, fAlreadyPadded, pProg);
}

///use fast integral image method for performing mean filter
/**Compute the mean filtered image using a square kernel of size
   (2*radiusX+1) by (2*radiusY+1) pixels.  For RGB images (or
   multi-channel float or double images), each channel is filtered
   independently.
   
   This implementation uses integral images for speed, as described in
   "Efficient Implementation of Local Adaptive Thresholding Techniques
   Using Integral Images" by Shafait, Keysers, and Breuel that was to
   appear in Proc. Document Recognition and Retrieval XV, IST/SPIE
   Annual Symposium, San Jose, CA, January 2008. It uses more memory
   than a convolution approach since the integral image must be
   created and uses long longs or long doubles to prevent overflow.*/
void DMeanFilter::meanFilterUseIntegralImage_(DImage &imgDst,
					      const DImage &imgSrc,
					      int radiusX, int radiusY,
					      bool fAlreadyPadded){
  // pad the image if necessary
  DImage *pimgSrcPad;
  int ww, wh, wa; // window width, window height, window area
  int w, h;
  
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
		imgSrc.getImageType(), imgSrc.numChannels());
  
  ww = 2*radiusX+1;
  wh = 2*radiusY+1;
  wa = ww * wh;
  
  switch(imgSrc.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint64 *pIntegralImg64;
	D_uint8 *pu8Dst;
	D_uint8 *pu8SrcPad;
	pIntegralImg64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h);
	D_CHECKPTR(pIntegralImg64);
	pu8Dst = imgDst.dataPointer_u8();
	pu8SrcPad = pimgSrcPad->dataPointer_u8();
	
	// fill in the integral image values
	pIntegralImg64[0] = (D_uint64)(pu8SrcPad[0]);
	for(int x = 1; x < w; ++x){ // first row
	  pIntegralImg64[x] = pIntegralImg64[x-1]+(D_uint64)(pu8SrcPad[x]);
	}
	for(int y = 1, idx=w; y < h; ++y){
	  pIntegralImg64[idx] = pIntegralImg64[idx-w] +
	    (D_uint64)(pu8SrcPad[idx]);
	  ++idx;
	  for(int x = 1; x < w; ++x, ++idx){
	    pIntegralImg64[idx] = pIntegralImg64[idx-1]+
	      pIntegralImg64[idx-w] + (D_uint64)(pu8SrcPad[idx]) -
	      pIntegralImg64[idx-w-1];
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
	    pu8Dst[idxDst] = (D_uint8)
	      ((pIntegralImg64[idxA] + pIntegralImg64[idxB]-
		pIntegralImg64[idxC] - pIntegralImg64[idxD]) / (double)wa);
	    ++idxA;
	    ++idxB;
	    ++idxC;
	    ++idxD;
	  }//end for(x...
	}//end for(y...
	free(pIntegralImg64);
      }
      break;
    case DImage::DImage_u16:
      {
	D_uint64 *pIntegralImg64;
	D_uint16 *pu16Dst;
	D_uint16 *pu16SrcPad;
	pIntegralImg64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h);
	D_CHECKPTR(pIntegralImg64);
	pu16Dst = imgDst.dataPointer_u16();
	pu16SrcPad = pimgSrcPad->dataPointer_u16();
	
	// fill in the integral image values
	pIntegralImg64[0] = (D_uint64)(pu16SrcPad[0]);
	for(int x = 1; x < w; ++x){ // first row
	  pIntegralImg64[x] = pIntegralImg64[x-1]+(D_uint64)(pu16SrcPad[x]);
	}
	for(int y = 1, idx=w; y < h; ++y){
	  pIntegralImg64[idx] = pIntegralImg64[idx-w] +
	    (D_uint64)(pu16SrcPad[idx]);
	  ++idx;
	  for(int x = 1; x < w; ++x, ++idx){
	    pIntegralImg64[idx] = pIntegralImg64[idx-1]+
	      pIntegralImg64[idx-w] + (D_uint64)(pu16SrcPad[idx]) -
	      pIntegralImg64[idx-w-1];
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
	    pu16Dst[idxDst] = (D_uint16)
	      ((pIntegralImg64[idxA] + pIntegralImg64[idxB]-
		pIntegralImg64[idxC] - pIntegralImg64[idxD]) / (double)wa);
	    ++idxA;
	    ++idxB;
	    ++idxC;
	    ++idxD;
	  }//end for(x...
	}//end for(y...
	free(pIntegralImg64);
      }
      break;
    case DImage::DImage_u32:
      {
	D_uint64 *pIntegralImg64;
	D_uint32 *pu32Dst;
	D_uint32 *pu32SrcPad;
	pIntegralImg64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h);
	D_CHECKPTR(pIntegralImg64);
	pu32Dst = imgDst.dataPointer_u32();
	pu32SrcPad = pimgSrcPad->dataPointer_u32();
	
	// fill in the integral image values
	pIntegralImg64[0] = (D_uint64)(pu32SrcPad[0]);
	for(int x = 1; x < w; ++x){ // first row
	  pIntegralImg64[x] = pIntegralImg64[x-1]+(D_uint64)(pu32SrcPad[x]);
	}
	for(int y = 1, idx=w; y < h; ++y){
	  pIntegralImg64[idx] = pIntegralImg64[idx-w] +
	    (D_uint64)(pu32SrcPad[idx]);
	  ++idx;
	  for(int x = 1; x < w; ++x, ++idx){
	    pIntegralImg64[idx] = pIntegralImg64[idx-1]+
	      pIntegralImg64[idx-w] + (D_uint64)(pu32SrcPad[idx]) -
	      pIntegralImg64[idx-w-1];
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
	    pu32Dst[idxDst] = (D_uint32)
	      ((pIntegralImg64[idxA] + pIntegralImg64[idxB]-
		pIntegralImg64[idxC] - pIntegralImg64[idxD]) / (double)wa);
	    ++idxA;
	    ++idxB;
	    ++idxC;
	    ++idxD;
	  }//end for(x...
	}//end for(y...
	free(pIntegralImg64);
      }
      break;
    case DImage::DImage_RGB:
      {
	D_uint64 *pIntegralImg64;
	D_uint8 *pu8Dst;
	D_uint8 *pu8SrcPad;
	DImage *rgImgSrcChannels;

	rgImgSrcChannels = new DImage[imgSrc.numChannels()];
	D_CHECKPTR(rgImgSrcChannels);
	pIntegralImg64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h) ;
	D_CHECKPTR(pIntegralImg64);
	pu8Dst = imgDst.dataPointer_u8();

	pimgSrcPad->splitRGB(rgImgSrcChannels[0], rgImgSrcChannels[1],
			     rgImgSrcChannels[2]);

	for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
	  pu8SrcPad = rgImgSrcChannels[chan].dataPointer_u8();
	
	  // fill in the integral image values
	  pIntegralImg64[0] = (D_uint64)(pu8SrcPad[0]);
	  for(int x = 1; x < w; ++x){ // first row
	    pIntegralImg64[x] = pIntegralImg64[x-1]+
	      (D_uint64)(pu8SrcPad[x]);
	  }
	  for(int y = 1, idx=w; y < h; ++y){
	    pIntegralImg64[idx] = pIntegralImg64[idx-w] +
	      (D_uint64)(pu8SrcPad[idx]);
	    ++idx;
	    for(int x = 1; x < w; ++x, ++idx){
	      pIntegralImg64[idx] = pIntegralImg64[idx-1]+
		pIntegralImg64[idx-w] + (D_uint64)(pu8SrcPad[idx]) -
		pIntegralImg64[idx-w-1];
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
	      pu8Dst[idxDst*3+chan] = (D_uint8)
		((pIntegralImg64[idxA] + pIntegralImg64[idxB]-
		  pIntegralImg64[idxC] - pIntegralImg64[idxD]) / (double)wa);
	      ++idxA;
	      ++idxB;
	      ++idxC;
	      ++idxD;
	    }//end for(x...
	  }//end for(y...
	}//end for(chan...
	free(pIntegralImg64);
	delete [] rgImgSrcChannels;
      }
      
      break;
    case DImage::DImage_RGB_16:
      {
	D_uint64 *pIntegralImg64;
	D_uint16 *pu16Dst;
	D_uint16 *pu16SrcPad;
	DImage *rgImgSrcChannels;

	rgImgSrcChannels = new DImage[imgSrc.numChannels()];
	D_CHECKPTR(rgImgSrcChannels);
	pIntegralImg64 = (D_uint64*)malloc(sizeof(D_uint64) * w * h) ;
	D_CHECKPTR(pIntegralImg64);
	pu16Dst = imgDst.dataPointer_u16();

	pimgSrcPad->splitRGB(rgImgSrcChannels[0], rgImgSrcChannels[1],
			     rgImgSrcChannels[2]);

	for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
	  pu16SrcPad = rgImgSrcChannels[chan].dataPointer_u16();
	
	  // fill in the integral image values
	  pIntegralImg64[0] = (D_uint64)(pu16SrcPad[0]);
	  for(int x = 1; x < w; ++x){ // first row
	    pIntegralImg64[x] = pIntegralImg64[x-1]+
	      (D_uint64)(pu16SrcPad[x]);
	  }
	  for(int y = 1, idx=w; y < h; ++y){
	    pIntegralImg64[idx] = pIntegralImg64[idx-w] +
	      (D_uint64)(pu16SrcPad[idx]);
	    ++idx;
	    for(int x = 1; x < w; ++x, ++idx){
	      pIntegralImg64[idx] = pIntegralImg64[idx-1]+
		pIntegralImg64[idx-w] + (D_uint64)(pu16SrcPad[idx]) -
		pIntegralImg64[idx-w-1];
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
	      pu16Dst[idxDst*3+chan] = (D_uint16)
		((pIntegralImg64[idxA] + pIntegralImg64[idxB]-
		  pIntegralImg64[idxC] - pIntegralImg64[idxD]) / (double)wa);
	      ++idxA;
	      ++idxB;
	      ++idxC;
	      ++idxD;
	    }//end for(x...
	  }//end for(y...
	}//end for(chan...
	free(pIntegralImg64);
	delete [] rgImgSrcChannels;
      }
      
      break;
    case DImage::DImage_flt_multi:
      {
	long double *pIntegralImgDbl;
	float *pFltDst;
	float *pFltSrcPad;

	pIntegralImgDbl = (long double*)malloc(sizeof(long double) * w * h) ;
	D_CHECKPTR(pIntegralImgDbl);

	for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
	  pFltDst = imgDst.dataPointer_flt(chan);
	  pFltSrcPad = pimgSrcPad->dataPointer_flt(chan);
	  
	  // fill in the integral image values
	  pIntegralImgDbl[0] = (long double)(pFltSrcPad[0]);
	  for(int x = 1; x < w; ++x){ // first row
	    pIntegralImgDbl[x] = pIntegralImgDbl[x-1]+
	      (long double)(pFltSrcPad[x]);
	  }
	  for(int y = 1, idx=w; y < h; ++y){
	    pIntegralImgDbl[idx] = pIntegralImgDbl[idx-w] +
	      (long double)(pFltSrcPad[idx]);
	    ++idx;
	    for(int x = 1; x < w; ++x, ++idx){
	      pIntegralImgDbl[idx] = pIntegralImgDbl[idx-1]+
		pIntegralImgDbl[idx-w] + (long double)(pFltSrcPad[idx]) -
		pIntegralImgDbl[idx-w-1];
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
	      pFltDst[idxDst] = (float)
		((pIntegralImgDbl[idxA] + pIntegralImgDbl[idxB]-
		  pIntegralImgDbl[idxC] - pIntegralImgDbl[idxD]) / (float)wa);
	      ++idxA;
	      ++idxB;
	      ++idxC;
	      ++idxD;
	    }//end for(x...
	  }//end for(y...
	}//end for(chan...
	free(pIntegralImgDbl);
      }
      break;
    case DImage::DImage_dbl_multi:
      {
	long double *pIntegralImgDbl;
	double *pDblDst;
	double *pDblSrcPad;

	pIntegralImgDbl = (long double*)malloc(sizeof(long double) * w * h) ;
	D_CHECKPTR(pIntegralImgDbl);

	for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
	  pDblDst = imgDst.dataPointer_dbl(chan);
	  pDblSrcPad = pimgSrcPad->dataPointer_dbl(chan);
	  
	  // fill in the integral image values
	  pIntegralImgDbl[0] = (long double)(pDblSrcPad[0]);
	  for(int x = 1; x < w; ++x){ // first row
	    pIntegralImgDbl[x] = pIntegralImgDbl[x-1]+
	      (long double)(pDblSrcPad[x]);
	  }
	  for(int y = 1, idx=w; y < h; ++y){
	    pIntegralImgDbl[idx] = pIntegralImgDbl[idx-w] +
	      (long double)(pDblSrcPad[idx]);
	    ++idx;
	    for(int x = 1; x < w; ++x, ++idx){
	      pIntegralImgDbl[idx] = pIntegralImgDbl[idx-1]+
		pIntegralImgDbl[idx-w] + (long double)(pDblSrcPad[idx]) -
		pIntegralImgDbl[idx-w-1];
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
	      pDblDst[idxDst] = (double)
		((pIntegralImgDbl[idxA] + pIntegralImgDbl[idxB]-
		  pIntegralImgDbl[idxC] - pIntegralImgDbl[idxD]) / (double)wa);
	      ++idxA;
	      ++idxB;
	      ++idxC;
	      ++idxD;
	    }//end for(x...
	  }//end for(y...
	}//end for(chan...
	free(pIntegralImgDbl);
      }
      break;
    default:
      fprintf(stderr, "DMeanFilter::meanFilterUseIntegralImage_() not "
	      "implemented for complex images\n");
      abort();
  }

  if(!fAlreadyPadded){
    delete pimgSrcPad;
  }
}
