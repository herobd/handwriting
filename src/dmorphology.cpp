#include "dmorphology.h"
#include "dinstancecounter.h"




void DMorphology::dilate3x3_(DImage &imgDst, DImage &imgSrc, bool fBlackIsFG){
  int w, h;
  D_uint8 *p8src, *p8dst;
  if(DImage::DImage_u8 != imgSrc.getImageType()){
    fprintf(stderr,"DMorphology::dilate3x3_() currently only supports 8-bit GS images\n");
    exit(1);
  }
  if(!fBlackIsFG){
    fprintf(stderr,"DMorphology::dilate3x3_() currently only supports fBlackIsFG=TRUE\n");
    exit(1);
  }
  if((&imgDst)==(&imgSrc)){
    fprintf(stderr,"DMorphology::dilate3x3_() currently currently does not allow imgDst to be the same as imgSrc\n");
    exit(1);
  }
    
  w = imgSrc.width();
  h = imgSrc.height();
  p8src = imgSrc.dataPointer_u8();

  imgDst = imgSrc;
  p8dst = imgDst.dataPointer_u8();
  
  for(int y=0; y < h; ++y){
    for(int x = 0; x < w; ++x){
      bool fFoundFG;
      fFoundFG = false;
      if((p8src[y*w+x] != 0) && (p8src[y*w+x] != 0xff)){
	fprintf(stderr,"DMorphology::dilate3x3_() found a pixel that is not 0 or 255!\n");
	exit(1);
      }
      for(int selY = -1; selY <= 1; ++selY){
	for(int selX = -1; selX <= 1; ++selX){
	  if(((x+selX)>=0) && ((x+selX)<w) &&
	     ((y+selY)>=0) && ((y+selY)<h) &&
	     (0x00 == p8src[(y+selY)*w+(x+selX)])){
	    fFoundFG = true;
	    break;
	  }
	}
      }
      if(fFoundFG){
	p8dst[y*w+x] = 0x00;
      }
      else{
	p8dst[y*w+x] = 0xff;
      }
    }
  }
}


void DMorphology::erode3x3_(DImage &imgDst, DImage &imgSrc, bool fBlackIsFG){
  int w, h;
  D_uint8 *p8src, *p8dst;
  if(DImage::DImage_u8 != imgSrc.getImageType()){
    fprintf(stderr,"DMorphology::erode3x3_() currently only supports 8-bit GS images\n");
    exit(1);
  }
  if(!fBlackIsFG){
    fprintf(stderr,"DMorphology::erode3x3_() currently only supports fBlackIsFG=TRUE\n");
    exit(1);
  }
  if((&imgDst)==(&imgSrc)){
    fprintf(stderr,"DMorphology::erode3x3_() currently currently does not allow imgDst to be the same as imgSrc\n");
    exit(1);
  }
    
  w = imgSrc.width();
  h = imgSrc.height();
  p8src = imgSrc.dataPointer_u8();

  imgDst = imgSrc;
  p8dst = imgDst.dataPointer_u8();
  
  for(int y=0; y < h; ++y){
    for(int x = 0; x < w; ++x){
      bool fFoundBG;
      fFoundBG = false;
      if((p8src[y*w+x] != 0) && (p8src[y*w+x] != 0xff)){
	fprintf(stderr,"DMorphology::erode3x3_() found a pixel that is not 0 or 255!\n");
	exit(1);
      }
      for(int selY = -1; selY <= 1; ++selY){
	for(int selX = -1; selX <= 1; ++selX){
	  if(((x+selX)>=0) && ((x+selX)<w) &&
	     ((y+selY)>=0) && ((y+selY)<h) &&
	     (0xff == p8src[(y+selY)*w+(x+selX)])){
	    fFoundBG = true;
	    break;
	  }
	}
      }
      if(fFoundBG){
	p8dst[y*w+x] = 0xff;
      }
      else{
	p8dst[y*w+x] = 0x00;
      }
    }
  }
}
