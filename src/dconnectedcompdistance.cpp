#include "dconnectedcompdistance.h"

///for now, this just calls getCCDistAndNearestCC_() and discards imgNearestCC
void DConnectedComponentDistance::getCCDist_(DImage &imgCCDist,
					     const DImage &imgCCs){
  DImage imgTmp;
  getCCDistAndNearestCC_(imgCCDist, imgTmp, imgCCs);
}

///for now, this just calls getCCDistAndNearestCC_() and discards imgCCDist
void DConnectedComponentDistance::getNearestCC_(DImage &imgNearestCC,
						const DImage &imgCCs){
  DImage imgTmp;
  getCCDistAndNearestCC_(imgTmp, imgNearestCC, imgCCs);
}


///determine distance to (and ID of) nearest CC to each pixel
/**This function uses a two-pass algorithm to determine the manhattan
   distance (diagonals are ignored) from any given pixel to the
   connected component "nearest" to it, as well as the CC label of the
   connected component that it is "nearest" to.  This may not actually
   be the nearest component due to the fact that diagonal directions
   are ignored, the pixels are always checked in the same order of
   directions, and distance/CClabel are only updated if the distance
   from the direction being checked is LESS THAN the distance reorded
   so far.  So this is a pretty good approximation, but not
   necessarily completely accurate.  The input image (imgCCs) should
   be a connected component map image (of type DImage_u32) with CC
   label 0 being the background pixels that are ignored by this
   algorithm.
 */
void DConnectedComponentDistance::getCCDistAndNearestCC_(DImage &imgCCDist,
							 DImage &imgNearestCC,
							 const DImage &imgCCs){
  D_uint32 *pCCs; // source connected component map image
  D_uint32 *pCCDist;
  D_uint32 *pNearestCC;
  int w, h; // image width and height
  int wm1, hm1; // width minus 1 and height minus 1
  
  if(imgCCs.getImageType() != DImage::DImage_u32){
    fprintf(stderr, "DConnectedComponentDistance::getCCDistAndNearestCC_() "
	    "expects imgCCs to be a CC map (type DImage_u32 with "
	    "0=background\n");
    abort();
  }
  w = imgCCs.width();
  h = imgCCs.height();
  wm1 = w - 1;
  hm1 = h - 1;
  imgCCDist.create(w, h, DImage::DImage_u32);
  imgNearestCC.create(w, h, DImage::DImage_u32);
  pCCs = imgCCs.dataPointer_u32();
  pCCDist = imgCCDist.dataPointer_u32();
  pNearestCC = imgNearestCC.dataPointer_u32();

  // STEP 1:
  // scan from top-left and record the shorter distance N,W along with the
  // ID of the component from that direction
  for(int y = 0, idx=0; y < h; ++y){
    for(int x = 0; x < w; ++x, ++idx){
      if(0 != pCCs[idx]){
	pCCDist[idx] = 0;
	pNearestCC[idx] = pCCs[idx];
      }
      else{
	pCCDist[idx] = 0xfffffffe;//(we need to add 1 without wrapping to 0)
	pNearestCC[idx] = 0;
      }
      // West
      if((x>0) && ((pCCDist[idx-1]+1) < pCCDist[idx])){
	pCCDist[idx] = pCCDist[idx-1]+1;
	pNearestCC[idx] = pNearestCC[idx-1];
      }
      // North
      if((y>0) && ((pCCDist[idx-w]+1) < pCCDist[idx])){
	pCCDist[idx] = pCCDist[idx-w]+1;
	pNearestCC[idx] = pNearestCC[idx-w];
      }
    }
  }

  // STEP 2:
  // scan from bottom right and update for S,E along with the
  // ID of the component from that shortest direction if shorter than current
  for(int y = hm1, idx=w*h-1; y >= 0; --y){
    for(int x = wm1; x >= 0; --x, --idx){
      // East
      if((x<wm1) && ((pCCDist[idx+1]+1) < pCCDist[idx])){
	pCCDist[idx] = pCCDist[idx+1]+1;
	pNearestCC[idx] = pNearestCC[idx+1];
      }
      // South
      if((y<hm1) && ((pCCDist[idx+w]+1) < pCCDist[idx])){
	pCCDist[idx] = pCCDist[idx+w]+1;
	pNearestCC[idx] = pNearestCC[idx+w];
      }
    }
  }
}
