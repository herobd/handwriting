#include "dcolorspace.h"

///Creates a 3-channel float CIELab image with channels: 0=L, 1=a, 2=b
void DColorSpace::convertRGBImageToCIELab(DImage &imgDst,
					  const DImage &imgSrc){
  int w, h;
  D_uint8 *pSrc;
  float *pL;
  float *pA;
  float *pB;
  if(DImage::DImage_RGB != imgSrc.getImageType()){
    fprintf(stderr, "convertRGBImageToCIELab() called with non-RGB image\n");
    exit(1);
  }
  w = imgSrc.width();
  h = imgSrc.height();
  imgDst.create(w,h,DImage::DImage_flt_multi,3,imgSrc.getAllocMethod());
  pSrc = imgSrc.dataPointer_u8();
  pL = imgDst.dataPointer_flt(0);
  pA = imgDst.dataPointer_flt(1);
  pB = imgDst.dataPointer_flt(2);
  for(int y = 0, idx = 0, idxD = 0; y < h; ++y){
    for(int x = 0; x < w; ++x, idx+=3, ++idxD){
      getCIELabFromRGB(pSrc[idx], pSrc[idx+1], pSrc[idx+2],
		       &pL[idxD], &pA[idxD], &pB[idxD]);
    }
  }
}

void DColorSpace::convertRGBImageToHSV(DImage &imgDst,
				       const DImage &imgSrc){
  int w, h;
  D_uint8 *pSrc;
  float *pH;
  float *pS;
  float *pV;
  if(DImage::DImage_RGB != imgSrc.getImageType()){
    fprintf(stderr, "convertRGBImageToHSV() called with non-RGB image\n");
    exit(1);
  }
  w = imgSrc.width();
  h = imgSrc.height();
  imgDst.create(w,h,DImage::DImage_flt_multi,3,imgSrc.getAllocMethod());
  pSrc = imgSrc.dataPointer_u8();
  pH = imgDst.dataPointer_flt(0);
  pS = imgDst.dataPointer_flt(1);
  pV = imgDst.dataPointer_flt(2);
  for(int y = 0, idx = 0, idxD = 0; y < h; ++y){
    for(int x = 0; x < w; ++x, idx+=3, ++idxD){
      getHSVFromRGB(pSrc[idx], pSrc[idx+1], pSrc[idx+2],
		       &pH[idxD], &pS[idxD], &pV[idxD]);
    }
  }
}


void DColorSpace::convertCIELabImageToRGB(DImage &imgDst,
					  const DImage &imgSrc){
  int w, h;
  D_uint8 *pDst;
  float *pL;
  float *pA;
  float *pB;
  if(DImage::DImage_flt_multi != imgSrc.getImageType()){
    fprintf(stderr, "convertCIELabImageToRGB() called with non-float image\n");
    exit(1);
  }
  if(3 != imgSrc.numChannels()){
    fprintf(stderr, "convertCIELabImageToRGB() needs 3-channel float image\n");
    exit(1);
  }
  w = imgSrc.width();
  h = imgSrc.height();
  imgDst.create(w,h,DImage::DImage_RGB,3,imgSrc.getAllocMethod());
  pL = imgSrc.dataPointer_flt(0);
  pA = imgSrc.dataPointer_flt(1);
  pB = imgSrc.dataPointer_flt(2);
  pDst = imgDst.dataPointer_u8();
  for(int y = 0, idx = 0, idxD = 0; y < h; ++y){
    for(int x = 0; x < w; ++x, ++idx, idxD+=3){
      getRGBFromCIELab(pL[idx], pA[idx], pB[idx],
		       &pDst[idxD], &pDst[idxD+1], &pDst[idxD+2]);
    }
  }
}


void DColorSpace::convertHSVImageToRGB(DImage &imgDst,
					  const DImage &imgSrc){
  int w, h;
  D_uint8 *pDst;
  float *pH;
  float *pS;
  float *pV;
  if(DImage::DImage_flt_multi != imgSrc.getImageType()){
    fprintf(stderr, "convertHSVImageToRGB() called with non-float image\n");
    exit(1);
  }
  if(3 != imgSrc.numChannels()){
    fprintf(stderr, "convertHSVImageToRGB() needs 3-channel float image\n");
    exit(1);
  }
  w = imgSrc.width();
  h = imgSrc.height();
  imgDst.create(w,h,DImage::DImage_RGB,3,imgSrc.getAllocMethod());
  pH = imgSrc.dataPointer_flt(0);
  pS = imgSrc.dataPointer_flt(1);
  pV = imgSrc.dataPointer_flt(2);
  pDst = imgDst.dataPointer_u8();
  for(int y = 0, idx = 0, idxD = 0; y < h; ++y){
    for(int x = 0; x < w; ++x, ++idx, idxD+=3){
      getRGBFromHSV(pH[idx], pS[idx], pV[idx],
		       &pDst[idxD], &pDst[idxD+1], &pDst[idxD+2]);
    }
  }
}
