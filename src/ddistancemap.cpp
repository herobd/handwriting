#include "ddistancemap.h"


/**Assumes that ink is black (0x00) and any non-zero pixel is background.  Positive distances are distance from ink.  Negative distances are within the ink (becoming more negative the deeper into the ink they get).  Distance is Manhattan distance, so it is always an integer.  Values are clamped to maxDist and MaxNegDist.  The image type is unsigned, so when using them, the data pointer should be cast to signed 32-bit integer. **/
void DDistanceMap::getDistFromInkBitonal_(DImage &dst, DImage &src,
					  int maxDist, int maxNegDist){
  //  DImage imgResult;
  int w, h;
  D_uint8 *pu8;
  D_sint32 *ps32;
  const int MAXDIST = 9999;
  if(src.getImageType() != DImage::DImage_u8){
    fprintf(stderr, "DDistanceMap::getDistFromInkBitonal() only works on "
	    "grayscale images\n");
    return;
  }
  w = src.width();
  h = src.height();
  pu8 = src.dataPointer_u8();
  dst.create(w,h,DImage::DImage_u32,1);
  //dst.fill(0.); //don't fill - we are going to overwrite it all anyway
  ps32 = (D_sint32*)dst.dataPointer_u32();

  //first pass - forward -----------------------------------
  //first pixel (x=0, y=0)
  if(0x00 == pu8[0]){
    //&&((0x00 != pu8[1])||(0x00!=pu8[w]))){
    ps32[0] = -1*MAXDIST;
  }
  else{
    ps32[0] = MAXDIST;
  }
  //first row (y=0)
  for(int x=1; x < w; ++x){
    if(0x00 == pu8[x]){//black (ink)
      if(0x00 == pu8[x-1])
	ps32[x] = ps32[x-1]-1;
      else
	ps32[x] = 0;
    }
    else{
      ps32[x] = ps32[x-1]+1;
      if(ps32[x] < 1)
	ps32[x] = 1;
      // if(ps32[x] > MAXDIST)
      // 	ps32[x] = MAXDIST;
    }
  }
  for(int y=1,idx=w; y < h; ++y){
    //first col of row y
    if(0x00 == pu8[idx]){//black (ink)
      if(0x00 == pu8[idx-w])
	ps32[idx] = ps32[idx-w]-1;
      else
	ps32[idx] = 0;
    }
    else{
      ps32[idx] = ps32[idx-w]+1;
      if(ps32[idx] < 1)
	ps32[idx] = 1;
      // if(ps32[idx] > MAXDIST)
      // 	ps32[idx] = MAXDIST;
    }
    ++idx;
    //rest of cols in row y
    for(int x=1; x < w; ++x,++idx){
      if(0x00 == pu8[idx]){//black (ink)
	D_sint32 min;
	if((ps32[idx-1] > 0) || (ps32[idx-w] > 0))
	  ps32[idx] = 0;
	else{//min is actually min abs value here for internal distance
	  min = ((ps32[idx-1]>ps32[idx-w]) ? (ps32[idx-1]) : (ps32[idx-w]));
	  ps32[idx] = min - 1;
	}
      }
      else{//not ink
	D_sint32 min;
	min = ((ps32[idx-1]<ps32[idx-w]) ? (ps32[idx-1]) : (ps32[idx-w]));
	if(min < 0)
	  min = 0;
	ps32[idx] = min + 1;
      }

    }
  }


  //second pass - backward -----------------------------------
  //adjust dist for bottom-right pixel (x=w-1, y=h-1) if ink
  // if(0x00 == pu8[w*h-1]){
  //   ps32[w*h-1] = 0; // only change if setting to zero
  // }
  //bottom row: (y=h-1) so idx=y*w-2 to start at next to last pixel
  for(int x=w-2, idx=w*h-2; x >= 0; --x, --idx){
    if(ps32[idx] < 0){ // black (ink) -- see if EAST is closer than current
      if(pu8[idx+1] != 0x00)
	ps32[idx] = 0;
      else if((ps32[idx+1]-1) > (ps32[idx]))
	ps32[idx] = (ps32[idx+1]-1);
    }
    else if(ps32[idx] > 0){
      if(pu8[idx+1] == 0x00)
	ps32[idx] = 1;
      else{
	if((ps32[idx+1]+1) < (ps32[idx])){
	  ps32[idx] = (ps32[idx+1]+1);
	  if(ps32[idx] < 1)
	    ps32[idx] = 1;
	}
	//else don't make it bigger than already is
      }
    }
    //else if already zero don't change it
  }
  //rest of rows
  for(int y=h-2,idx=w*(h-1)-1; y >= 0; --y){
    //right col of row y
    if(ps32[idx] < 0){ // black (ink) - see if SOUTH is closer than current
      if(pu8[idx+w] != 0x00)
	ps32[idx] = 0;
      else{
	if((ps32[idx+w]-1) > ps32[idx])
	  ps32[idx] = (ps32[idx+w]-1);
      }
    }
    else{
      if((ps32[idx+w]+1) < ps32[idx]){
	ps32[idx] = (ps32[idx+w]+1);
	if(ps32[idx] < 1)
	  ps32[idx] = 1;
      }
      // if(ps32[idx] > MAXDIST)
      // 	ps32[idx] = MAXDIST;
    }
    --idx;
    //rest of cols in row y
    for(int x=w-2; x >= 0; --x,--idx){
      if(ps32[idx]<0){ // black (ink) - see if EAST/SOUTH is closer than current
	if((pu8[idx+1] != 0x00) || (pu8[idx+w] != 0x00))
	  ps32[idx] = 0;
	else{
	  D_sint32 min;//min is actually min abs value here for internal dist
	  min = ((ps32[idx+1]>ps32[idx+w]) ? (ps32[idx+1]) : (ps32[idx+w]));
	  if((min-1) > ps32[idx])
	    ps32[idx] = min-1;
	}
      }
      else{
	D_sint32 min;
	min = ((ps32[idx+1]<ps32[idx+w]) ? (ps32[idx+1]) : (ps32[idx+w]));
	if(min < 0)
	  min = 0;
	if((min+1) < ps32[idx])
	  ps32[idx] = min + 1;
	// if(ps32[idx] > MAXDIST)
	//   ps32[idx] = MAXDIST;
      }
    }
  }

  //fix distance range from maxNegDist to maxDist
  for(int idx = 0, len = w*h; idx < len; ++idx){
    if(ps32[idx] > maxDist)
      ps32[idx] = maxDist;
    if(ps32[idx] < maxNegDist)
      ps32[idx] = maxNegDist;
  }


#if 0
  //debug image
  DImage imgDeb;
  imgDeb.create(w,h,DImage::DImage_RGB);
  D_uint8 *pdeb;
  pdeb=imgDeb.dataPointer_u8();
  for(int y=0, idx=0; y < h; ++y){
    for(int x=0; x < w; ++x, ++idx){
      if((int)(ps32[idx]) < 0){
	// pdeb[idx*3]=(unsigned char)(ps32[idx]);
	pdeb[idx*3]=(unsigned char)(-16*ps32[idx]);
	pdeb[idx*3+1]=0;
      }
      else if((int)(ps32[idx]) > 0){
	pdeb[idx*3]=0;
	pdeb[idx*3+1]=(unsigned char)(ps32[idx]*8);
      }
      else{
	pdeb[idx*3]=255;
	pdeb[idx*3+1]=255;
      }
      pdeb[idx*3+2]=0;
    }
  }
  imgDeb.save("/tmp/dist.ppm");
#endif // debug image
}



