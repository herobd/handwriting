#include "dconvolver.h"
#include "ddefs.h"
#include "dimage.h"
#include "dkernel2d.h"
#include "dprogress.h"


///Convolve imgSrc with kern and put the result into imgDst
/**If fSrcAlreadyPadded is false, a padded version of the imgSrc will
 * be created (by replicating the edge pixels) so that edge conditions
 * can be ignored during convolution, and the pad will be removed in
 * the final results (imgDst) if fResize is false or left in imgDst if
 * fResize is true.  If fSrcAlreadyPadded is true, then it is assumed
 * that imgSrc has already been padded with (kern.getRadiusX()*2)
 * pixels on each side and (kern.getRadiusY()*2) pixels on top and
 * bottom. If fResize is true, then imgDst will be the size of the
 * original image plus 2*radiusX wide and 2*radiusY high (half the
 * kernel radius on each side).  If fConvertBack is true, imgDst will
 * be converted back to the same type as imgSrc (an 8-bit grayscale or
 * RGB image, for instance) instead of being returned as a
 * DImage_flt_multi or DImage_dbl_multi image.  If fUseDoublePrec is
 * true, double-precision floating-point numbers (double) will be used
 * instead of single-precision floating-point numbers (float) during
 * convolution.  If an image has multiple channels (an RGB image or
 * flt_multi, for example), then each channel will be convolved
 * separately.  If pProg is NULL, it is ignored, otherwise the
 * progress is updated periodically.
 */
void DConvolver::convolve_(DImage &imgDst, const DImage &imgSrc,
			   DKernel2D &kern, bool fSrcAlreadyPadded,
			   bool fResize, bool fConvertBack,
			   bool fUseDoublePrec, DProgress *pProg){
  DImage imgTmpPad;
  DImage *pPad;
  DImage imgTmpDst;
  DImage *pDst;
  DImage imgTmpSrc;
  DImage *pSrc;

  if(!kern.isInitialized()){
    fprintf(stderr,
	    "DConvolver::convolve_() called with uninitialized kernel\n");
    exit(1);
  }


  pPad = (DImage*)&imgSrc;
  if(!fSrcAlreadyPadded){
    // pad the image edges first (by replicating the edges) to account for
    // convolution filter size
    imgSrc.padEdges_(imgTmpPad, 2*kern.getRadiusX(),2*kern.getRadiusX(),
		     2*kern.getRadiusY(),2*kern.getRadiusY(),
		     DImage::DImagePadReplicate);
    pPad = &imgTmpPad;
  }
  // convert to either float or double image (if it isn't already one)
  pSrc = pPad;
  if((fUseDoublePrec) && (DImage::DImage_dbl_multi != pSrc->getImageType())){
    pPad->convertedImgType_(imgTmpSrc, DImage::DImage_dbl_multi,
			    -1, 0xffffffff, AllocationMethod_malloc);
    pSrc = &imgTmpSrc;
  }
  else if((!fUseDoublePrec) &&
	  (DImage::DImage_flt_multi != pSrc->getImageType())){
    pPad->convertedImgType_(imgTmpSrc, DImage::DImage_flt_multi,
			    -1, 0xffffffff, AllocationMethod_malloc);
    pSrc = &imgTmpSrc;
  }
  // point pDst to either imgTmpDst or imgDst, depending on fConvertBack
  pDst = &imgDst;
  if(fConvertBack){
    pDst = &imgTmpDst;
  }

  //call correct private method for double/float separable/non-separable kernel
  if(fUseDoublePrec){
    if(kern.isSeparable())
      convSepDbl(*pDst, *pSrc, kern, fResize, pProg);
    else
      convRegDbl(*pDst, *pSrc, kern, fResize, pProg);
  }
  else{
    if(kern.isSeparable()){
      convSepFlt(*pDst, *pSrc, kern, fResize, pProg);
    }
    else
      convRegFlt(*pDst, *pSrc, kern, fResize, pProg);
  }
  
  // now convert back to the appropriate type (if requested)
  if(fConvertBack){
    pDst->convertedImgType_(imgDst, imgSrc.getImageType(),
			    -1, 0xffffffff, imgSrc.getAllocMethod());
  }
}

void DConvolver::convSepFlt(DImage &imgDst, const DImage &imgSrc,
			    DKernel2D &kern,
			    bool fResize, DProgress *pProg){
  DImage imgInter; // intermediate image to store vertical results
  int sw,sh; // width, height of imgSrc (source image)
  int dw,dh; // width, height of imgDst (destination image)
  int kw,kh; // width, height of kern
  int krx,kry; // radiusX and radiusY of kern
  int soffsX,soffsY; // x and y offset for source image pixels
  float *pDataSrc;
  float *pDataInter;
  float *pDataDst;
  float *pDst;
  float *pKern;
  float *pSrc;
  float *pInter;
  float fltSum;
  int progMax;
  int progCur = 0;
  float *colBuf;
  float *pColBuf;
  
  sw = imgSrc.width();
  sh = imgSrc.height();
  krx = kern.getRadiusX();
  kry = kern.getRadiusY();
  kw = kern.getWidth();
  kh = kern.getHeight();
  // sw and sh are already paded to 2*radius on each side

  if(fResize){
    dw = sw - 2*krx;
    dh = sh - 2*kry;
    soffsX = 0;
    soffsY = 0;
  }
  else{
    dw = sw - 4*krx;
    dh = sh - 4*kry;
    soffsX = krx;
    soffsY = kry;
  }

  imgDst.create(dw, dh, DImage::DImage_flt_multi, imgSrc.numChannels(),
		imgSrc.getAllocMethod());
  imgInter.create(sw, sh, DImage::DImage_flt_multi, imgSrc.numChannels(),
		  AllocationMethod_malloc);

  progMax = imgSrc.numChannels() * (dh+sh-2*kry);
  colBuf = (float*)malloc(sizeof(float)*sh);
  if(!colBuf){
    fprintf(stderr, "DConvolver::convSepFlt() out of memory!\n");
    exit(1);
  }
  for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
    pDataSrc = imgSrc.dataPointer_flt(chan);
    pDataDst = imgDst.dataPointer_flt(chan);
    pDataInter = imgInter.dataPointer_flt(chan);

    // do vertical direction first
    pKern = &((kern.getData_flt())[kern.getWidth()]);//point to vertical kern
    pInter = &pDataInter[kry*sw];
    pSrc = pDataSrc;
    for(int x = 0; x < sw; ++x){ // each intermediate column
      // update progress report and check if user cancelled the operation
      if(NULL != pProg){
	if(0 == (progCur & 0x000000ff)){
	  if(0 != pProg->reportStatus(progCur, 0, progMax)){
	    // the operation has been cancelled
	    pProg->reportStatus(-1, 0, progMax); // report cancel acknowledged
	    return;
	  }
	}
	++progCur;//current progress
      }
      pColBuf = colBuf;
      //To reduce cache misses, we copy whichever column we are
      //dealing with into a contiguous buffer and treat it the same as
      //we would a row.  This cut the total convolution time down to
      //less than half of what it was on a 3K x 4K RGB image with
      //35x35 kernel.
      //
      //TODO: instead of copying column by column into the
      //contiguous buffer, it might be faster to flip the whole image
      //sideways all at once, do the vertical conv as if it were
      //horizontal, flip the intermediate back to the original
      //orientation and then do the horizontal.  It may save a little
      //time, but we are already almost as fast as the GIMP filter, so
      //I don't think it will save more than 3 or 4 seconds (at most).
      //Whatever we change here should be done in both the float and double
      //version of this method.

      // copy column into contiguous buffer
      for(int y = kry, limy = sh-kry; y < limy; ++y){ // each intermediate row
	colBuf[y] = pSrc[y*sw];
      }
      // perform the convolution using the contiguous buffer
      for(int y = kry, limy = sh-kry; y < limy; ++y){ // each intermediate row
	fltSum = 0.;
	for(int ky = 0; ky < kh; ++ky){ // each kernel position
	  fltSum += pColBuf[ky] * pKern[ky];
	}
	(*pColBuf) = fltSum;
	++pColBuf;
      }

      // copy the "vertical" results from contiguous buffer to intermediate img
      for(int y = kry, limy = sh-kry; y < limy; ++y){ // each intermediate row
	pInter[y*sw] = colBuf[y];
      }
      ++pInter;
      ++pSrc;
    }
      
    // do horizontal direction now ------------------------------------------
    pKern = kern.getData_flt();//point to horizontal kern
    pDst = pDataDst;
    for(int y = 0; y < dh; ++y){ // each dst row
      // update progress report and check if user cancelled the operation
      if(NULL != pProg){
	if(0 == (progCur & 0x000000ff)){
	  if(0 != pProg->reportStatus(progCur, 0, progMax)){
	    // the operation has been cancelled
	    pProg->reportStatus(-1, 0, progMax); // report cancel acknowledged
	    return;
	  }
	}
	++progCur;//current progress
      }

      pInter = &pDataInter[(kry+soffsY+y)*sw+ soffsX];//go to start of next row
      for(int x = 0; x < dw; ++x){ // each dst column
	fltSum = 0;
	for(int kx = 0; kx < kw; ++kx){ // each kernel position
	  fltSum += pInter[kx] * pKern[kx];
	}
	pDst[x] = fltSum;
	++pInter;// move pSrc to left end of kernel for next col
      }
      pDst += dw;
    }
  }//end for chan

  if(NULL != pProg){ // report progress (complete)
    pProg->reportStatus(progMax, 0, progMax);
  }
}

void DConvolver::convSepDbl(DImage &imgDst, const DImage &imgSrc,
			    DKernel2D &kern,
			    bool fResize, DProgress *pProg){
  DImage imgInter; // intermediate image to store vertical results
  int sw,sh; // width, height of imgSrc (source image)
  int dw,dh; // width, height of imgDst (destination image)
  int kw,kh; // width, height of kern
  int krx,kry; // radiusX and radiusY of kern
  int soffsX,soffsY; // x and y offset for source image pixels
  double *pDataSrc;
  double *pDataInter;
  double *pDataDst;
  double *pDst;
  double *pKern;
  double *pSrc;
  double *pInter;
  double dblSum;
  int progMax;
  int progCur = 0;
  double *colBuf;
  double *pColBuf;
  
  sw = imgSrc.width();
  sh = imgSrc.height();
  krx = kern.getRadiusX();
  kry = kern.getRadiusY();
  kw = kern.getWidth();
  kh = kern.getHeight();
  // sw and sh are already paded to 2*radius on each side

  if(fResize){
    dw = sw - 2*krx;
    dh = sh - 2*kry;
    soffsX = 0;
    soffsY = 0;
  }
  else{
    dw = sw - 4*krx;
    dh = sh - 4*kry;
    soffsX = krx;
    soffsY = kry;
  }

  imgDst.create(dw, dh, DImage::DImage_dbl_multi, imgSrc.numChannels(),
		imgSrc.getAllocMethod());
  imgInter.create(sw, sh, DImage::DImage_dbl_multi, imgSrc.numChannels(),
		  AllocationMethod_malloc);

  progMax = imgSrc.numChannels() * (dh+sh-2*kry);
  colBuf = (double*)malloc(sizeof(double)*sh);
  if(!colBuf){
    fprintf(stderr, "DConvolver::convSepDbl() out of memory!\n");
    exit(1);
  }
  for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
    pDataSrc = imgSrc.dataPointer_dbl(chan);
    pDataDst = imgDst.dataPointer_dbl(chan);
    pDataInter = imgInter.dataPointer_dbl(chan);

    // do vertical direction first
    pKern = &((kern.getData_dbl())[kern.getWidth()]);//point to vertical kern
    pInter = &pDataInter[kry*sw];
    pSrc = pDataSrc;
    for(int x = 0; x < sw; ++x){ // each intermediate column
      // update progress report and check if user cancelled the operation
      if(NULL != pProg){
	if(0 == (progCur & 0x000000ff)){
	  if(0 != pProg->reportStatus(progCur, 0, progMax)){
	    // the operation has been cancelled
	    pProg->reportStatus(-1, 0, progMax); // report cancel acknowledged
	    return;
	  }
	}
	++progCur;//current progress
      }
      pColBuf = colBuf;
      //To reduce cache misses, we copy whichever column we are
      //dealing with into a contiguous buffer and treat it the same as
      //we would a row.  This cut the total convolution time down to
      //less than half of what it was on a 3K x 4K RGB image with
      //35x35 kernel.
      //
      //TODO: instead of copying column by column into the
      //contiguous buffer, it might be faster to flip the whole image
      //sideways all at once, do the vertical conv as if it were
      //horizontal, flip the intermediate back to the original
      //orientation and then do the horizontal.  It may save a little
      //time, but we are already almost as fast as the GIMP filter, so
      //I don't think it will save more than 3 or 4 seconds (at most).
      //Whatever we change here should be done in both the float and double
      //version of this method.

      // copy column into contiguous buffer
      for(int y = kry, limy = sh-kry; y < limy; ++y){ // each intermediate row
	colBuf[y] = pSrc[y*sw];
      }
      // perform the convolution using the contiguous buffer
      for(int y = kry, limy = sh-kry; y < limy; ++y){ // each intermediate row
	dblSum = 0.;
	for(int ky = 0; ky < kh; ++ky){ // each kernel position
	  dblSum += pColBuf[ky] * pKern[ky];
	}
	(*pColBuf) = dblSum;
	++pColBuf;
      }

      // copy the "vertical" results from contiguous buffer to intermediate img
      for(int y = kry, limy = sh-kry; y < limy; ++y){ // each intermediate row
	pInter[y*sw] = colBuf[y];
      }
      ++pInter;
      ++pSrc;
    }
      
    // do horizontal direction now ------------------------------------------
    pKern = kern.getData_dbl();//point to horizontal kern
    pDst = pDataDst;
    for(int y = 0; y < dh; ++y){ // each dst row
      // update progress report and check if user cancelled the operation
      if(NULL != pProg){
	if(0 == (progCur & 0x000000ff)){
	  if(0 != pProg->reportStatus(progCur, 0, progMax)){
	    // the operation has been cancelled
	    pProg->reportStatus(-1, 0, progMax); // report cancel acknowledged
	    return;
	  }
	}
	++progCur;//current progress
      }

      pInter = &pDataInter[(kry+soffsY+y)*sw+ soffsX];//go to start of next row
      for(int x = 0; x < dw; ++x){ // each dst column
	dblSum = 0;
	for(int kx = 0; kx < kw; ++kx){ // each kernel position
	  dblSum += pInter[kx] * pKern[kx];
	}
	pDst[x] = dblSum;
	++pInter;// move pSrc to left end of kernel for next col
      }
      pDst += dw;
    }
  }//end for chan

  if(NULL != pProg){ // report progress (complete)
    pProg->reportStatus(progMax, 0, progMax);
  }
}
void DConvolver::convRegFlt(DImage &imgDst, const DImage &imgSrc,
			    DKernel2D &kern,
			    bool fResize, DProgress *pProg){
  int sw,sh; // width, height of imgSrc (source image)
  int dw,dh; // width, height of imgDst (destination image)
  int kw,kh; // width, height of kern
  int krx,kry; // radiusX and radiusY of kern
  int soffsX,soffsY; // x and y offset for source image pixels
  float *pDataSrc;
  float *pDataDst;
  float *pDst;
  float *pKern;
  float *pSrc;
  float fltSum;
  int progMax;
  int progCur = 0;
  
  sw = imgSrc.width();
  sh = imgSrc.height();
  krx = kern.getRadiusX();
  kry = kern.getRadiusY();
  kw = kern.getWidth();
  kh = kern.getHeight();
  // sw and sh are already paded to 2*radius on each side
  if(fResize){
    dw = sw - 2*krx;
    dh = sh - 2*kry;
    soffsX = 0;
    soffsY = 0;
  }
  else{
    dw = sw - 4*krx;
    dh = sh - 4*kry;
    soffsX = krx;
    soffsY = kry;
  }

  imgDst.create(dw, dh, DImage::DImage_flt_multi, imgSrc.numChannels(),
		imgSrc.getAllocMethod());
  progMax = imgSrc.numChannels() * dh;
  for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
    pDataSrc = imgSrc.dataPointer_flt(chan);
    pDataDst = imgDst.dataPointer_flt(chan);


    pKern = kern.getData_flt();//point to kernel floating-point data
    
    // pSrc will start where top-left of kernel is, then proceed right/down
    // as kernel multiplied, then will be reset to top-left for next
    // dst position
    pSrc = &pDataSrc[soffsY*sw + soffsX];
    pDst = pDataDst;
    for(int y = 0; y < dh; ++y){ // each dst row
      // update progress report and check if user cancelled the operation
      if(NULL != pProg){
	if(0 == (progCur & 0x0000003f)){
	  if(0 != pProg->reportStatus(progCur, 0, progMax)){
	    // the operation has been cancelled
	    pProg->reportStatus(-1, 0, progMax); // report cancel acknowledged
	    return;
	  }
	}
	++progCur;//current progress
      }
      for(int x = 0; x < dw; ++x){ // each dst column
	fltSum = 0;
	for(int ky = 0, kidx=0; ky < kh; ++ky){ // each kernel row
	  for(int kx = 0; kx < kw; ++kx, ++kidx){ // each kernel column
	    fltSum += pSrc[kx] * pKern[kidx];
	  }
	  pSrc += sw; // move pSrc down to next row of kernel
	}
	pDst[x] = fltSum;
	pSrc -= (kh*sw-1);// move pSrc to top of kernel for next column
      }
      pSrc += (sw - dw); // position to start of next row
      pDst += dw;
    }
    
  }      
  if(NULL != pProg){ // report progress (complete)
    pProg->reportStatus(progMax, 0, progMax);
  }
}
void DConvolver::convRegDbl(DImage &imgDst, const DImage &imgSrc,
			    DKernel2D &kern,
			    bool fResize, DProgress *pProg){
  int sw,sh; // width, height of imgSrc (source image)
  int dw,dh; // width, height of imgDst (destination image)
  int kw,kh; // width, height of kern
  int krx,kry; // radiusX and radiusY of kern
  int soffsX,soffsY; // x and y offset for source image pixels
  double *pDataSrc;
  double *pDataDst;
  double *pDst;
  double *pKern;
  double *pSrc;
  double dblSum;
  int progMax;
  int progCur = 0;
  
  sw = imgSrc.width();
  sh = imgSrc.height();
  krx = kern.getRadiusX();
  kry = kern.getRadiusY();
  kw = kern.getWidth();
  kh = kern.getHeight();
  // sw and sh are already paded to 2*radius on each side
  if(fResize){
    dw = sw - 2*krx;
    dh = sh - 2*kry;
    soffsX = 0;
    soffsY = 0;
  }
  else{
    dw = sw - 4*krx;
    dh = sh - 4*kry;
    soffsX = krx;
    soffsY = kry;
  }

  imgDst.create(dw, dh, DImage::DImage_dbl_multi, imgSrc.numChannels(),
		imgSrc.getAllocMethod());
  progMax = imgSrc.numChannels() * dh;
  for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
    pDataSrc = imgSrc.dataPointer_dbl(chan);
    pDataDst = imgDst.dataPointer_dbl(chan);


    pKern = kern.getData_dbl();//point to kernel double data
    
    // pSrc will start where top-left of kernel is, then proceed right/down
    // as kernel multiplied, then will be reset to top-left for next
    // dst position
    pSrc = &pDataSrc[soffsY*sw + soffsX];
    pDst = pDataDst;
    for(int y = 0; y < dh; ++y){ // each dst row
      // update progress report and check if user cancelled the operation
      if(NULL != pProg){
	if(0 == (progCur & 0x0000003f)){
	  if(0 != pProg->reportStatus(progCur, 0, progMax)){
	    // the operation has been cancelled
	    pProg->reportStatus(-1, 0, progMax); // report cancel acknowledged
	    return;
	  }
	}
	++progCur;//current progress
      }
      for(int x = 0; x < dw; ++x){ // each dst column
	dblSum = 0;
	for(int ky = 0, kidx=0; ky < kh; ++ky){ // each kernel row
	  for(int kx = 0; kx < kw; ++kx, ++kidx){ // each kernel column
	    dblSum += pSrc[kx] * pKern[kidx];
	  }
	  pSrc += sw; // move pSrc down to next row of kernel
	}
	pDst[x] = dblSum;
	pSrc -= (kh*sw-1);// move pSrc to top of kernel for next column
      }
      pSrc += (sw - dw); // position to start of next row
      pDst += dw;
    }
  }      
  if(NULL != pProg){ // report progress (complete)
    pProg->reportStatus(progMax, 0, progMax);
  }
}
