#include "dmedialaxis.h"
#include "ddistancemap.h"
#include "dthresholder.h"


/**
	As an alternative to Doug's method using the medail axis, I am 
	implementing additional thinning algorithms to see if they produce
	better skeletons.
	
	This first is a modification of the Zhang-Suen method proposed by Chen
	and Hsu in 1988, which supposedly fixes some of the originals flaws.
	(http://ac.els-cdn.com/0167865588901249/1-s2.0-0167865588901249-main.pdf?_tid=994be7c0-393f-11e4-9f70-00000aacb360&acdnat=1410390784_863de2b56be91cad8f478ed2dd1b7fae)
	
	I will try to make this as similar to Doug's as possible, although
	fThin is currently ignored.
	Assumes that ink is black (0x00) and any non-zero pixel is background 
	for source image. Returns an 8-bit grayscale and has 0x00 everywhere 
	except for the points that have medial axis (on), which will be 0xff.
**/

//Can be parallelized (with GPU I think)
DImage DMedialAxis::getMZhangSkeletonFromBinaryImage(DImage &src,
						  bool fThin, int *numPoints){
	int w, h;
	//D_uint8 *p8_src;
	//D_uint8 *p8_imgs[2];
	//DImage img1;
	//DImage img2;
	DImage imgs[2];
	if(DImage::DImage_u8 != src.getImageType()){
		fprintf(stderr, "DMedialAxis::getMZhangSkeletonFromImage() expects an "
		    "image of type DImage_u8, which is filled with binary data, 0x00 "
		    "being foreground, all else being background.\n");
		exit(1);
	}
	w = src.width();
	h = src.height();
	imgs[0].create(w,h,DImage::DImage_u8);
	imgs[1].create(w,h,DImage::DImage_u8);
	//p8_imgs[0] = img1.dataPointer_u8();
	//p8_imgs[1] = img1.dataPointer_u8();
	//p8_src = src.dataPointer_u8();
	
	//set up first image
	for (int x=0; x<w; x++){
		for (int y=0; y<h; y++){
			setBPixel(imgs[0],x,y,!isOn(src,x,y));//invert as input data is different
		}
	}
	
	int k = 0;
	int s = 1;
	bool flag=true;
	int count;
	while (flag){
		count=0;
		flag=false;
		for (int x=0; x<w; x++){
			for (int y=0; y<h; y++){
				if (isOn(imgs[k],x,y)){
					int B_P1 = neighborOnCount(imgs[k],x,y);
					int A_P1 = neighborOffToOnTransitionCount(imgs[k],x,y);
					if (k==0){
						
						bool d=(2 <= B_P1 && B_P1 <=7 &&			//a
							  A_P1 == 1 &&					//b
							  !(isOn(imgs[k],x,y-1) && isOn(imgs[k],x+1,y) && 
							  isOn(imgs[k],x,y+1)) && 			//c
							  !(isOn(imgs[k],x-1,y) && isOn(imgs[k],x+1,y) && 
							  isOn(imgs[k],x,y+1)))				//d
							 ||
							 (A_P1 == 2 &&					//e
							  isOn(imgs[k],x,y-1) && isOn(imgs[k],x+1,y) && 
							  !isOn(imgs[k],x,y+1) && !isOn(imgs[k],x-1,y+1) && 
							  !isOn(imgs[k],x-1,y) &&			//f
							  isOn(imgs[k],x+1,y) && isOn(imgs[k],x,y+1) && 
							  !isOn(imgs[k],x,y-1) && !isOn(imgs[k],x-1,y) && 
							  !isOn(imgs[k],x-1,y-1));			//g
						setBPixel(imgs[s],x,y,d);
						if (d)
							count++;
						else
							flag=true;
							
					}
					else{
						
						bool d=(2 <= B_P1 && B_P1 <=7 &&			//a
							  A_P1 == 1 &&					//b
							  !(isOn(imgs[k],x,y-1) && isOn(imgs[k],x+1,y) && 
							  isOn(imgs[k],x-1,y)) && 			//c'
							  !(isOn(imgs[k],x,y-1) && isOn(imgs[k],x-1,y) && 
							  isOn(imgs[k],x,y+1)))				//d'
							 ||
							 (A_P1 == 2 &&					//e
							  isOn(imgs[k],x,y-1) && isOn(imgs[k],x-1,y) && 
							  !isOn(imgs[k],x+1,y) && !isOn(imgs[k],x+1,y+1) && 
							  !isOn(imgs[k],x,y+1) &&			//f'
							  isOn(imgs[k],x-1,y) && isOn(imgs[k],x,y+1) && 
							  !isOn(imgs[k],x,y-1) && !isOn(imgs[k],x+1,y-1) && 
							  !isOn(imgs[k],x+1,y));			//g'
						setBPixel(imgs[s],x,y,d);
						if (d)
							count++;
						else
							flag=true;
							
					}
				}
				else{
					setBPixel(imgs[s],x,y,false);
				}
			}
		}
		
		k = (k+1)%2;
		s = (s+1)%2;
	}
	*numPoints=count;
	
	imgs[k].save("./tmp/debug_images/zhang.ppm");
	
	return imgs[k];
}

int DMedialAxis::neighborOffToOnTransitionCount(const DImage &img, int x, int y){
	int count=0;
	
	if (!isOn(img,x,y-1) && isOn(img,x+1,y-1)) count++;
	else if (!isOn(img,x+1,y-1) && isOn(img,x+1,y)) count++;
	if (!isOn(img,x+1,y) && isOn(img,x+1,y+1)) count++;
	else if (!isOn(img,x+1,y+1) && isOn(img,x,y+1)) count++;
	if (!isOn(img,x,y+1) && isOn(img,x-1,y+1)) count++;
	else if (!isOn(img,x-1,y+1) && isOn(img,x-1,y)) count++;
	if (!isOn(img,x-1,y) && isOn(img,x-1,y-1)) count++;
	else if (!isOn(img,x-1,y-1) && isOn(img,x,y-1)) count++;
	
	return count;
}

int DMedialAxis::neighborOnCount(const DImage &img, int x, int y){
	int count=0;
	
	if (isOn(img,x,y-1)) count++;
	if (isOn(img,x+1,y-1)) count++;
	if (isOn(img,x+1,y)) count++;
	if (isOn(img,x+1,y+1)) count++;
	if (isOn(img,x,y+1)) count++;
	if (isOn(img,x-1,y+1)) count++;
	if (isOn(img,x-1,y)) count++;
	if (isOn(img,x-1,y-1)) count++;
	
	return count;
}

bool DMedialAxis::isOn(const DImage &img, int x, int y){
	const D_uint8 *p8_img = img.dataPointer_u8();
	if (x>=0 && y>=0 && x<img.width() && y<img.height()){
		int idx=x+y*img.width();
		return 0x00 != p8_img[idx];
	}
	else{
		return false;
	}
}

void DMedialAxis::setBPixel(DImage &img, int x, int y, bool on){
	D_uint8 *p8_img = img.dataPointer_u8();
	if (x>=0 && y>=0 && x<img.width() && y<img.height()){
		int idx=x+y*img.width();
		p8_img[idx] = on ? 0xFF : 0x00;
	}
}


/**Assumes the input is a distance map found with DDistanceMap.
   Returns an 8-bit grayscale and has 0x00 everywhere except for the
   points that have medial axis, which will be 0xff.  The medial axis
   is found by finding local minima in the distance map (negative
   points that have no neighbors that are "more negative" than the
   points themselves).  If fThin is true, another pass of the image will be
   made and wherever superpixel composed of 4 pixels is set, one will be
   eliminated.

   If I need a reference for this, the URL 
   http://portal.acm.org/citation.cfm?id=321357 is a good one.  It is
   Azriel Rosenfeld and John L. Pfaltz, "Sequential Operations in Digital
   Picture Processing" Journal of the ACM (JACM) Vol 13 Issue 4, Oct 1996,
   pp.471-494.  Cited elsewhere as "Comm. ACM" instead of JACM.  Perhaps the
   name changed?  Anyway, my algorithm isn't exactly the same, but pretty
   close.
*/
DImage DMedialAxis::getMedialAxisImageFromDistMap(DImage &imgDistMap,
						  bool fThin, int *numPoints){
  int w, h;
  signed int *ps32;
  D_uint8 *p8;
  DImage imgMA;
  int numPxls = 0;
  if(DImage::DImage_u32 != imgDistMap.getImageType()){
    fprintf(stderr, "DMedialAxis::getMedialAxisImageFromDistMap() expects an "
	    "image of type DImage_u32, which is actually filled with signed "
	    "32-bit integers representing a distance map. (use DDistanceMap to "
	    "compute it from a bitonal image).\n");
    exit(1);
  }
  w = imgDistMap.width();
  h = imgDistMap.height();
  imgMA.create(w,h,DImage::DImage_u8);
  p8 = imgMA.dataPointer_u8();
  ps32 = (signed int*)imgDistMap.dataPointer_u32();

  for(int y=0, idx=0; y < h; ++y){
    int ym1, yp1;
    ym1 = (y>0) ? (0-w) : 0;
    yp1 = (y<(h-1)) ? w: 0;
    for(int x=0; x < w; ++x, ++idx){
      //The following code was only checking 4-neighbors (leaving incorrect
      // points as median pixels) and was also not allowing image edges as MA
      // if( (ps32[idx]<=0) && // only non-positive (ink) can be medial axis
      // 	  ((x>0)&&(ps32[idx] <= ps32[idx-1])) &&
      // 	  ((x<(w-1))&&(ps32[idx] <= ps32[idx+1])) &&
      // 	  ((y>0)&&(ps32[idx] <= ps32[idx-w])) &&
      // 	  ((y<(h-1))&&(ps32[idx] <= ps32[idx+w])) ){
      // 	p8[idx] = 0xff;
      // 	++numPxls;
      // 	continue;
      // }
      if(ps32[idx]<=0){ // only non-positive (ink) can be medial axis
      	int xm1, xp1;
      	xm1 = (x>0) ? -1 : 0;
      	xp1 = (x<(w-1)) ? 1 : 0;
      	if(
      	    (ps32[idx] <= ps32[idx+0+ym1]) && //N
      	    (ps32[idx] <= ps32[idx+xm1+0]) && //W
      	    (ps32[idx] <= ps32[idx+xp1+0]) && //E
      	    (ps32[idx] <= ps32[idx+0+yp1])){ //S
      	  p8[idx] = 0xff;
      	  ++numPxls;
      	  continue;
      	}
      }
      // changed to the following:
      // if(ps32[idx]<=0){ // only non-positive (ink) can be medial axis
      // 	int xm1, xp1;
      // 	xm1 = (x>0) ? -1 : 0;
      // 	xp1 = (x<(w-1)) ? 1 : 0;
      // 	if( (ps32[idx] <= ps32[idx+xm1+ym1]) && //NW
      // 	    (ps32[idx] <= ps32[idx+0+ym1]) && //N
      // 	    (ps32[idx] <= ps32[idx+xp1+ym1]) && //NE
      // 	    (ps32[idx] <= ps32[idx+xm1+0]) && //W
      // 	    (ps32[idx] <= ps32[idx+xp1+0]) && //E
      // 	    (ps32[idx] <= ps32[idx+xm1+yp1]) && //SW
      // 	    (ps32[idx] <= ps32[idx+0+yp1]) && //S
      // 	    (ps32[idx] <= ps32[idx+xp1+yp1])){ //SE
      // 	  p8[idx] = 0xff;
      // 	  ++numPxls;
      // 	  continue;
      // 	}
      // }
      p8[idx] = 0;
    }
  }

  //do morphological close
  if(0)
  {
    DImage imgClose;
    D_uint8 *p8c;
    imgClose.create(w,h,DImage::DImage_u8);
    p8c = imgClose.dataPointer_u8();
    numPxls = 0;
    for(int y=0,idx=0; y<h; ++y){
      int ym1, yp1;
      ym1 = (y>0) ? (0-w) : 0;
      yp1 = (y<(h-1)) ? w: 0;
      for(int x=0;x<w;++x,++idx){
      	int xm1, xp1;
      	xm1 = (x>0) ? -1 : 0;
      	xp1 = (x<(w-1)) ? 1 : 0;
	if((p8[idx+ym1+xm1]>0)||(p8[idx+ym1]>0)||(p8[idx+ym1+xp1]>0)||
	   (p8[idx+xm1]>0)||(p8[idx]>0)||(p8[idx+xp1]>0)||
	   (p8[idx+yp1+xm1]>0)||(p8[idx+yp1]>0)||(p8[idx+yp1+xp1]>0)){
	  p8c[idx] = 0xff;
	  ++numPxls;
	}
	else
	  p8c[idx] = 0x00;
      }
    }
   #if 1
    numPxls = 0;
    for(int y=0,idx=0; y<h; ++y){
      int ym1, yp1;
      ym1 = (y>0) ? (0-w) : 0;
      yp1 = (y<(h-1)) ? w: 0;
      for(int x=0;x<w;++x,++idx){
      	int xm1, xp1;
      	xm1 = (x>0) ? -1 : 0;
      	xp1 = (x<(w-1)) ? 1 : (w-1);
	if((p8c[idx+ym1+xm1]==0)||(p8c[idx+ym1]==0)||(p8c[idx+ym1+xp1]==0)||
	   (p8c[idx+xm1]==0)||(p8c[idx]==0)||(p8c[idx+xp1]==0)||
	   (p8c[idx+yp1+xm1]==0)||(p8c[idx+yp1]==0)||(p8c[idx+yp1+xp1]==0))
	  p8[idx] = 0x00;
	else{
	  p8[idx] = 0xff;
	  ++numPxls;
	}
      }
    }
   #else
    imgMA = imgClose;
    p8 = imgMA.dataPointer_u8();
   #endif
  }


  if(fThin){
    // if(0){
    DImage imgTmp;
    D_uint8 *p8t;
    imgTmp = imgMA;
    p8t = imgTmp.dataPointer_u8();
    for(int y=1, idx=w; y < h; ++y){
      ++idx;
      for(int x=1; x < w; ++x, ++idx){
	if((0==p8t[idx])||(0==p8t[idx-1])||(0==p8t[idx-w])||(0==p8t[idx-w-1]))
	  continue;
	p8[idx] = 0;
	--numPxls;
      }
    }
  }
  if(NULL != numPoints)
    (*numPoints) = numPxls;
    
  imgMA.save("./tmp/debug_images/MA.ppm");
    
  return imgMA;
}

/**Given an 8-bit grayscale image, finds the medial axis.  Result
   image is 8-bit grayscale and has 0x00 everywhere except for the
   points that have medial axis, which will be 0xff.  The medial axis
   is found by computing a 4-connected (Manhattan Distance) distance
   map of img (thresholded at thresholdVal) (see DDistanceMap class)
   and then finding local minima within the ink areas.*/
DImage DMedialAxis::getMedialAxisImage(DImage &img, bool fZeroIsInk,
				       int thresholdVal){
  DImage imgDst;
  DImage imgDistanceMap;
  int lenMA = 0;
  imgDst = img;
  if(!fZeroIsInk)
    imgDst.invertGrayscale();
  DThresholder::threshImage_(imgDst,imgDst,thresholdVal);
  //  printf("w=%d h=%d\n", img.width(), img.height());
  // imgDistanceMap = DDistanceMap::getDistFromInkBitonal(imgDst,10000,-10000);
  DDistanceMap::getDistFromInkBitonal_(imgDistanceMap,imgDst,10000,-10000);
  imgDst = DMedialAxis::getMedialAxisImageFromDistMap(imgDistanceMap,
						      true,&lenMA);

  // fprintf(stderr, "DMedialAxis::getMedialAxisImageFromBitonal() NYI!\n");
  // exit(1);
  return imgDst;
}

DImage DMedialAxis::colorizeMedialAxisImage(DImage &imgMA, int R, int G, int B,
					    int bgR, int bgG, int bgB){
  int w, h;
  D_uint8 *p8;
  D_uint8 *pdst8;
  DImage imgResult;

  w = imgMA.width();
  h = imgMA.height();
  imgResult.create(w,h,DImage::DImage_RGB);

  p8 = imgMA.dataPointer_u8();
  pdst8 = imgResult.dataPointer_u8();
  for(int y=0, idx=0; y < h; ++y){
    for(int x=0; x < w; ++x, ++idx){
      if(0xff == p8[idx]){
	pdst8[idx*3]=(D_uint8)R;
	pdst8[idx*3+1]=(D_uint8)G;
	pdst8[idx*3+2]=(D_uint8)B;
      }
      else{
	pdst8[idx*3]=(D_uint8)bgR;
	pdst8[idx*3+1]=(D_uint8)bgG;
	pdst8[idx*3+2]=(D_uint8)bgB;
      }
    }
  }
  return imgResult;
}
