#include "dconnectedcomplabeler.h"
#include "dconnectedcompinfo.h"
#include <stdio.h>
#include <stdlib.h>

// increment size of array for the equivalence table
// (this many slots allocated at a time, as needed)
#define CC_SIZE_INCR 10000

///returns a DImage that is a connected component map (type=DImage_u32)
/**pixels with value BGval will be labeled with component ID 0
 * (background), whether they are connected or not.  For 8-bit images,
 * the value of BGval is cast to an 8-bit value.  Also, any adjacent
 * non-background pixels are considered to be connected.  This may not
 * be what you want.  Only 8-connected neighborhoods are currently
 * implemented, and only for images of type DImage_u8. The number of
 * non-background components is stored in *numCCs. (background pixels
 * are labeled as 0, and component 0 does not count as a component in
 * the number of components)
 */
DImage DConnectedComponentLabeler::getCCimage(const DImage &imgSrc,
					      int *numCCs,
					      D_uint32 BGval,
					      bool f8connected,
					      bool fSequentialIDs){
  DImage dst;
  getCCimage_(dst, imgSrc, numCCs, BGval, f8connected, fSequentialIDs);
  return dst;
}




///sets imgDst to a DImage that is a connected component map (type=DImage_u32)
/**pixels with value BGval will be labeled with component ID 0
 * (background), whether they are connected or not.  For 8-bit images,
 * the value of BGval is cast to an 8-bit value.  Also, any adjacent
 * non-background pixels are considered to be connected.  This may not
 * be what you want.  Only 8-connected neighborhoods are currently
 * implemented, and only for images of type DImage_u8. The number of
 * non-background components is stored in *numCCs. (background pixels
 * are labeled as 0, and component 0 does not count as a component in
 * the number of components)
 */
void DConnectedComponentLabeler::getCCimage_(DImage &imgDst,
					     const DImage &imgSrc,
					     int *numCCs,
					     D_uint32 BGval,
					     bool f8connected,
					     bool fSequentialIDs){
  int w, h, wm1;
  D_uint32 *pDst;
  D_uint16 *pSrc16;
  D_uint32 *pSrc32;
  int numBlobIds = 0; // blob IDs used so far (before merging equivalent IDs)
  D_uint32 lblW, lblNW, lblN, lblNE;
  D_uint32 rgCurBlobs[4]; //distinct blob IDs for current top/left neighborhood
  int numCurBlobs;// num of distinct blob IDs in current top/left neighborhood
  D_uint32 minLabel;
  D_uint32 newLbl;
  D_uint32 oldLbl;
  D_uint32 *rgEquivTable;
  D_uint32 *rgBlobIDs;
  int equivTableSize = 0;
  int cc_size_incr = CC_SIZE_INCR;

  if(!f8connected){
    fprintf(stderr, "DConnectedComponentLabeler::getCCimage_() "
	    "only 8-connected is implemented\n");
    exit(1);
  }

  w = imgSrc.width();
  wm1 = w-1;
  h = imgSrc.height();

  if(cc_size_incr < w){
    cc_size_incr = w;
  }

  if((w != imgDst.width()) || (h != imgDst.height()) ||
     (imgDst.getImageType() != DImage::DImage_u32))
    imgDst.create(w, h, DImage::DImage_u32);
  pDst = imgDst.dataPointer_u32();
  
  rgEquivTable = (D_uint32*)malloc(sizeof(D_uint32)*cc_size_incr);
  D_CHECKPTR(rgEquivTable);
  equivTableSize = cc_size_incr;
  for(int i = 0; i < cc_size_incr; ++i)
    rgEquivTable[i] = i;

  switch(imgSrc.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint8 bgVal8;
	D_uint8 *pSrc8;
	D_uint32 *pPrevDst;
	pSrc8 = imgSrc.dataPointer_u8();
	bgVal8 = (D_uint8)BGval;
	// do first row (any time the current pixel isn't equal to the previous
	// and not background, then new blobID
	if(pSrc8[0] == bgVal8){ // first pixel (on first row)
	  pDst[0] = 0;
	}
	else{
	  ++numBlobIds;
	  pDst[0] = numBlobIds;
	}
	// the rest of the pixels on first row
	for(int x = 1; x < w; ++x){
	  if(pSrc8[x] == bgVal8){
	    pDst[x] = 0;
	    continue;
	  }
	  if(pSrc8[x-1] == bgVal8){
	    ++numBlobIds;
	  }
	  pDst[x] = numBlobIds;
	}

	// now do the rest of the rows
	pPrevDst = pDst;
	pSrc8 += w;
	pDst += w;
	for(int y = 1; y < h; ++y){
	  if((w+numBlobIds) >= equivTableSize){
	    equivTableSize += cc_size_incr;
	    rgEquivTable = (D_uint32*)realloc(rgEquivTable,
					      sizeof(D_uint32)*
					      equivTableSize);
	    for(int i =equivTableSize-cc_size_incr; i <equivTableSize; ++i)
	      rgEquivTable[i] = i;
	  }
	  if(pSrc8[0] == bgVal8){ // first pixel in row y
	    pDst[0] = 0;
	  }
	  else{
	    if(0 != pPrevDst[0])
	      pDst[0] = pPrevDst[0];
	    else if(0 != pPrevDst[1])
	      pDst[0] = pPrevDst[1];
	    else{
	      ++numBlobIds;
	      pDst[0] = numBlobIds;
	    }
	  }
	  // the rest of the pixels in row y
	  for(int x = 1; x < wm1; ++x){
	    if(pSrc8[x] == bgVal8){
	      pDst[x] = 0;
	      continue;
	    }
	    lblNW = rgEquivTable[pPrevDst[x-1]];
	    lblN = rgEquivTable[pPrevDst[x]];
	    lblNE = rgEquivTable[pPrevDst[x+1]];
	    lblW = rgEquivTable[pDst[x-1]];
	    if(0 == (lblNW | lblN | lblNE | lblW)){ // new blobID
	      ++numBlobIds;
	      pDst[x] = numBlobIds;
	      continue;
	    }
	    // figure out how many distinct blobs there are (and their ids)
	    numCurBlobs = 0;
	    if(0 != lblNW){
	      rgCurBlobs[numCurBlobs] = lblNW;
	      ++numCurBlobs;
	    }
	    if(0 != lblN){
	      if(lblN != lblNW){
		rgCurBlobs[numCurBlobs] = lblN;
		++numCurBlobs;
	      }
	    }
	    if(0 != lblNE){
	      if((lblNE != lblNW) && (lblNE != lblN)){
		rgCurBlobs[numCurBlobs] = lblNE;
		++numCurBlobs;
	      }
	    }
	    if(0 != lblW){
	      if((lblW != lblNW) && (lblW != lblN) && (lblW != lblNE)){
		rgCurBlobs[numCurBlobs] = lblW;
		++numCurBlobs;
	      }
	    }
#ifdef DEBUG
	    if(0 == numCurBlobs)
	      fprintf(stderr, "logic error! (%s:%d)\n", __FILE__,__LINE__);
#endif
	    if(1 == numCurBlobs){ // a single blobID is touching, so use it
	      pDst[x] = rgCurBlobs[0];
	      continue;
	    }
	    // multiple blob IDs.  Need to update equivalence table
	    minLabel = rgCurBlobs[0];
	    if(rgCurBlobs[1] < minLabel)
	      minLabel = rgCurBlobs[1];
	    if(numCurBlobs > 2){
	      if(rgCurBlobs[2] < minLabel)
		minLabel = rgCurBlobs[2];
	      if(numCurBlobs > 3){ // but I dont think numCurBlobs CAN be > 3
		if(rgCurBlobs[3] < minLabel)
		  minLabel = rgCurBlobs[3];
	      }
	    }
#ifdef DEBUG
	    if(minLabel == 0){
	      fprintf(stderr,"minLabel==0!\n");
	      fprintf(stderr,"numCurBlobs=%d rgCurBlobs=%d %d %d %d\n",
		      numCurBlobs, rgCurBlobs[0], rgCurBlobs[1],
		      rgCurBlobs[2], rgCurBlobs[3]);
	    }
#endif
	    oldLbl = pPrevDst[x-1];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
	      oldLbl = newLbl;
	    }
	    oldLbl = pPrevDst[x];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
		oldLbl = newLbl;
	    }
	    oldLbl = pPrevDst[x+1];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
		oldLbl = newLbl;
	    }
	    oldLbl = pDst[x-1];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
	      oldLbl = newLbl;
	    }
	    pDst[x] = minLabel;
	  }//end for(x...
	  // do the last pixel in row y
	  for(int x = wm1; x==wm1; ++x){ // last pixel in row y
	    if(pSrc8[x] == bgVal8){
	      pDst[x] = 0;
	      continue;
	    }
	    lblNW = rgEquivTable[pPrevDst[x-1]];
	    lblN = rgEquivTable[pPrevDst[x]];
	    lblW = rgEquivTable[pDst[x-1]];
	    if(0 == (lblNW | lblN | lblW)){ // new blobID
	      ++numBlobIds;
	      pDst[x] = numBlobIds;
	      continue;
	    }
	    // figure out how many distinct blobs there are (and their ids)
	    numCurBlobs = 0;
	    if(0 != lblNW){
	      rgCurBlobs[numCurBlobs] = lblNW;
	      ++numCurBlobs;
	    }
	    if(0 != lblN){
	      if(lblN != lblNW){
		rgCurBlobs[numCurBlobs] = lblN;
		++numCurBlobs;
	      }
	    }
	    if(0 != lblW){
	      if((lblW != lblNW) && (lblW != lblN)){
		rgCurBlobs[numCurBlobs] = lblW;
		++numCurBlobs;
	      }
	    }
#ifdef DEBUG
	    if(0 == numCurBlobs)
	      fprintf(stderr, "logic error! (%s:%d)\n", __FILE__,__LINE__);
#endif
	    if(1 == numCurBlobs){ // a single blobID is touching, so use it
	      pDst[x] = rgCurBlobs[0];
	      continue;
	    }
	    // multiple blob IDs.  Need to update equivalence table
	    minLabel = rgCurBlobs[0];
	    if(rgCurBlobs[1] < minLabel)
	      minLabel = rgCurBlobs[1];
	    if(numCurBlobs > 2){// but I dont think numCurBlobs CAN be > 2
	      if(rgCurBlobs[2] < minLabel)
		minLabel = rgCurBlobs[2];
	    }
#ifdef DEBUG
	    if(minLabel == 0){
	      fprintf(stderr,"minLabel==0!\n");
	      fprintf(stderr,"numCurBlobs=%d rgCurBlobs=%d %d %d %d\n",
		      numCurBlobs, rgCurBlobs[0], rgCurBlobs[1],
		      rgCurBlobs[2], rgCurBlobs[3]);
	    }
#endif
	    oldLbl = pPrevDst[x-1];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
	      oldLbl = newLbl;
	    }
	    oldLbl = pPrevDst[x];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
		oldLbl = newLbl;
	    }
	    oldLbl = pPrevDst[x+1];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
		oldLbl = newLbl;
	    }
	    oldLbl = pDst[x-1];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
	      oldLbl = newLbl;
	    }
	    pDst[x] = minLabel;
	  }//end last pixel in row y
	  pPrevDst = pDst;
	  pSrc8 += w;
	  pDst += w;
	}// end for(y...


	// resolve equivalences
	rgBlobIDs = (D_uint32*)calloc(numBlobIds+1, sizeof(D_uint32));
	D_CHECKPTR(rgBlobIDs);
	for(unsigned int i = 0; i <= (unsigned int)numBlobIds; ++i){
	  if(rgEquivTable[i] != i){
	    unsigned int j;
	    j = rgEquivTable[i];
	    while(rgEquivTable[j] != j)
	      j = rgEquivTable[j];
	    if(rgEquivTable[i] != j){
	      rgEquivTable[i] = j;
	    }
	  }
	  rgBlobIDs[rgEquivTable[i]] = 1;
	}

	// count number of CCs and prepare mapping for second pass through img
	{
	  int ccCount = 0;

	  rgBlobIDs[0] = 0;
	  ccCount = 0;
	  for(int i = 1; i <= numBlobIds; ++i){
	    if(rgBlobIDs[i] > 0){
	      ++ccCount;
	      rgBlobIDs[i] = ccCount;
	    }
	  }
	  for(int i = 1; i <= numBlobIds; ++i){
	    rgEquivTable[i] = rgBlobIDs[rgEquivTable[i]];
	  }
	  (*numCCs) = ccCount;
	}
	
	// do second pass to map the CC labels to be sequential
	pDst = imgDst.dataPointer_u32();
	for(int y = 0, idx = 0; y < h; ++y){
	  for(int x = 0; x < w; ++x, ++idx){
	    pDst[idx] = rgEquivTable[pDst[idx]];
	  }
	}

	free(rgBlobIDs);

      }// end block within case statement
      break;
    default:
      fprintf(stderr, "DConnectedComponentLabeler::getCCimage_() not "
	      "implemented for image type %d\n", imgSrc.getImageType());
      exit(1);
  }

  free(rgEquivTable);
}












///sets imgDst to a CC map of pixels with value val
/**all pixels that did not have value val will be labeled with
 * component ID 0 (background), whether they were originally
 * background or not.  Only adjacent pixels with value val are
 * considered to be connected.  This may not be what you want.  Only
 * 8-connected neighborhoods are currently implemented, and only for
 * images of type DImage_u8. The number of components with value val
 * is stored in *numCCs. (non-val pixels are labeled as 0, and
 * component 0 does not count as a component in the number of
 * components)
 */
void DConnectedComponentLabeler::getCCimageForVal_(DImage &imgDst,
						   const DImage &imgSrc,
						   int *numCCs,
						   D_uint32 val,
						   bool f8connected,
						   bool fSequentialIDs){
  int w, h, wm1;
  D_uint32 *pDst;
  D_uint16 *pSrc16;
  D_uint32 *pSrc32;
  int numBlobIds = 0; // blob IDs used so far (before merging equivalent IDs)
  D_uint32 lblW, lblNW, lblN, lblNE;
  D_uint32 rgCurBlobs[4]; //distinct blob IDs for current top/left neighborhood
  int numCurBlobs;// num of distinct blob IDs in current top/left neighborhood
  D_uint32 minLabel;
  D_uint32 newLbl;
  D_uint32 oldLbl;
  D_uint32 *rgEquivTable;
  D_uint32 *rgBlobIDs;
  int equivTableSize = 0;
  int cc_size_incr = CC_SIZE_INCR;

  if(!f8connected){
    fprintf(stderr, "DConnectedComponentLabeler::getCCimageForVal_() "
	    "only 8-connected is implemented\n");
    exit(1);
  }

  w = imgSrc.width();
  wm1 = w-1;
  h = imgSrc.height();

  if(cc_size_incr < w){
    cc_size_incr = w;
  }

  if((w != imgDst.width()) || (h != imgDst.height()) ||
     (imgDst.getImageType() != DImage::DImage_u32))
    imgDst.create(w, h, DImage::DImage_u32);
  pDst = imgDst.dataPointer_u32();
  
  rgEquivTable = (D_uint32*)malloc(sizeof(D_uint32)*cc_size_incr);
  D_CHECKPTR(rgEquivTable);
  equivTableSize = cc_size_incr;
  for(int i = 0; i < cc_size_incr; ++i)
    rgEquivTable[i] = i;

  switch(imgSrc.getImageType()){
    case DImage::DImage_u8:
      {
	D_uint8 val8;
	D_uint8 *pSrc8;
	D_uint32 *pPrevDst;
	pSrc8 = imgSrc.dataPointer_u8();
	val8 = (D_uint8)val;
	// do first row (any time the current pixel isn't equal to the previous
	// and equals val, then new blobID
	if(pSrc8[0] != val8){ // first pixel (on first row)
	  pDst[0] = 0;
	}
	else{
	  ++numBlobIds;
	  pDst[0] = numBlobIds;
	}
	// the rest of the pixels on first row
	for(int x = 1; x < w; ++x){
	  if(pSrc8[x] != val8){
	    pDst[x] = 0;
	    continue;
	  }
	  if(pSrc8[x-1] != val8){
	    ++numBlobIds;
	  }
	  pDst[x] = numBlobIds;
	}

	// now do the rest of the rows
	pPrevDst = pDst;
	pSrc8 += w;
	pDst += w;
	for(int y = 1; y < h; ++y){
	  if((w+numBlobIds) >= equivTableSize){
	    equivTableSize += cc_size_incr;
	    rgEquivTable = (D_uint32*)realloc(rgEquivTable,
					      sizeof(D_uint32)*
					      equivTableSize);
	    for(int i =equivTableSize-cc_size_incr; i <equivTableSize; ++i)
	      rgEquivTable[i] = i;
	  }
	  if(pSrc8[0] != val8){ // first pixel in row y
	    pDst[0] = 0;
	  }
	  else{
	    if(0 != pPrevDst[0])
	      pDst[0] = pPrevDst[0];
	    else if(0 != pPrevDst[1])
	      pDst[0] = pPrevDst[1];
	    else{
	      ++numBlobIds;
	      pDst[0] = numBlobIds;
	    }
	  }
	  // the rest of the pixels in row y
	  for(int x = 1; x < wm1; ++x){
	    if(pSrc8[x] != val8){
	      pDst[x] = 0;
	      continue;
	    }
	    lblNW = rgEquivTable[pPrevDst[x-1]];
	    lblN = rgEquivTable[pPrevDst[x]];
	    lblNE = rgEquivTable[pPrevDst[x+1]];
	    lblW = rgEquivTable[pDst[x-1]];
	    if(0 == (lblNW | lblN | lblNE | lblW)){ // new blobID
	      ++numBlobIds;
	      pDst[x] = numBlobIds;
	      continue;
	    }
	    // figure out how many distinct blobs there are (and their ids)
	    numCurBlobs = 0;
	    if(0 != lblNW){
	      rgCurBlobs[numCurBlobs] = lblNW;
	      ++numCurBlobs;
	    }
	    if(0 != lblN){
	      if(lblN != lblNW){
		rgCurBlobs[numCurBlobs] = lblN;
		++numCurBlobs;
	      }
	    }
	    if(0 != lblNE){
	      if((lblNE != lblNW) && (lblNE != lblN)){
		rgCurBlobs[numCurBlobs] = lblNE;
		++numCurBlobs;
	      }
	    }
	    if(0 != lblW){
	      if((lblW != lblNW) && (lblW != lblN) && (lblW != lblNE)){
		rgCurBlobs[numCurBlobs] = lblW;
		++numCurBlobs;
	      }
	    }
#ifdef DEBUG
	    if(0 == numCurBlobs)
	      fprintf(stderr, "logic error! (%s:%d)\n", __FILE__,__LINE__);
#endif
	    if(1 == numCurBlobs){ // a single blobID is touching, so use it
	      pDst[x] = rgCurBlobs[0];
	      continue;
	    }
	    // multiple blob IDs.  Need to update equivalence table
	    minLabel = rgCurBlobs[0];
	    if(rgCurBlobs[1] < minLabel)
	      minLabel = rgCurBlobs[1];
	    if(numCurBlobs > 2){
	      if(rgCurBlobs[2] < minLabel)
		minLabel = rgCurBlobs[2];
	      if(numCurBlobs > 3){ // but I dont think numCurBlobs CAN be > 3
		if(rgCurBlobs[3] < minLabel)
		  minLabel = rgCurBlobs[3];
	      }
	    }
#ifdef DEBUG
	    if(minLabel == 0){
	      fprintf(stderr,"minLabel==0!\n");
	      fprintf(stderr,"numCurBlobs=%d rgCurBlobs=%d %d %d %d\n",
		      numCurBlobs, rgCurBlobs[0], rgCurBlobs[1],
		      rgCurBlobs[2], rgCurBlobs[3]);
	    }
#endif
	    oldLbl = pPrevDst[x-1];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
	      oldLbl = newLbl;
	    }
	    oldLbl = pPrevDst[x];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
		oldLbl = newLbl;
	    }
	    oldLbl = pPrevDst[x+1];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
		oldLbl = newLbl;
	    }
	    oldLbl = pDst[x-1];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
	      oldLbl = newLbl;
	    }
	    pDst[x] = minLabel;
	  }//end for(x...
	  // do the last pixel in row y
	  for(int x = wm1; x==wm1; ++x){ // last pixel in row y
	    if(pSrc8[x] != val8){
	      pDst[x] = 0;
	      continue;
	    }
	    lblNW = rgEquivTable[pPrevDst[x-1]];
	    lblN = rgEquivTable[pPrevDst[x]];
	    lblW = rgEquivTable[pDst[x-1]];
	    if(0 == (lblNW | lblN | lblW)){ // new blobID
	      ++numBlobIds;
	      pDst[x] = numBlobIds;
	      continue;
	    }
	    // figure out how many distinct blobs there are (and their ids)
	    numCurBlobs = 0;
	    if(0 != lblNW){
	      rgCurBlobs[numCurBlobs] = lblNW;
	      ++numCurBlobs;
	    }
	    if(0 != lblN){
	      if(lblN != lblNW){
		rgCurBlobs[numCurBlobs] = lblN;
		++numCurBlobs;
	      }
	    }
	    if(0 != lblW){
	      if((lblW != lblNW) && (lblW != lblN)){
		rgCurBlobs[numCurBlobs] = lblW;
		++numCurBlobs;
	      }
	    }
#ifdef DEBUG
	    if(0 == numCurBlobs)
	      fprintf(stderr, "logic error! (%s:%d)\n", __FILE__,__LINE__);
#endif
	    if(1 == numCurBlobs){ // a single blobID is touching, so use it
	      pDst[x] = rgCurBlobs[0];
	      continue;
	    }
	    // multiple blob IDs.  Need to update equivalence table
	    minLabel = rgCurBlobs[0];
	    if(rgCurBlobs[1] < minLabel)
	      minLabel = rgCurBlobs[1];
	    if(numCurBlobs > 2){// but I dont think numCurBlobs CAN be > 2
	      if(rgCurBlobs[2] < minLabel)
		minLabel = rgCurBlobs[2];
	    }
#ifdef DEBUG
	    if(minLabel == 0){
	      fprintf(stderr,"minLabel==0!\n");
	      fprintf(stderr,"numCurBlobs=%d rgCurBlobs=%d %d %d %d\n",
		      numCurBlobs, rgCurBlobs[0], rgCurBlobs[1],
		      rgCurBlobs[2], rgCurBlobs[3]);
	    }
#endif
	    oldLbl = pPrevDst[x-1];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
	      oldLbl = newLbl;
	    }
	    oldLbl = pPrevDst[x];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
		oldLbl = newLbl;
	    }
	    oldLbl = pPrevDst[x+1];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
		oldLbl = newLbl;
	    }
	    oldLbl = pDst[x-1];
	    while(oldLbl > minLabel){
	      newLbl = rgEquivTable[oldLbl];
	      rgEquivTable[oldLbl] = minLabel;
	      oldLbl = newLbl;
	    }
	    pDst[x] = minLabel;
	  }//end last pixel in row y
	  pPrevDst = pDst;
	  pSrc8 += w;
	  pDst += w;
	}// end for(y...


	// resolve equivalences
	rgBlobIDs = (D_uint32*)calloc(numBlobIds+1, sizeof(D_uint32));
	D_CHECKPTR(rgBlobIDs);
	for(unsigned int i = 0; i <= (unsigned int)numBlobIds; ++i){
	  if(rgEquivTable[i] != i){
	    unsigned int j;
	    j = rgEquivTable[i];
	    while(rgEquivTable[j] != j)
	      j = rgEquivTable[j];
	    if(rgEquivTable[i] != j){
	      rgEquivTable[i] = j;
	    }
	  }
	  rgBlobIDs[rgEquivTable[i]] = 1;
	}

	// count number of CCs and prepare mapping for second pass through img
	{
	  int ccCount = 0;

	  rgBlobIDs[0] = 0;
	  ccCount = 0;
	  for(int i = 1; i <= numBlobIds; ++i){
	    if(rgBlobIDs[i] > 0){
	      ++ccCount;
	      rgBlobIDs[i] = ccCount;
	    }
	  }
	  for(int i = 1; i <= numBlobIds; ++i){
	    rgEquivTable[i] = rgBlobIDs[rgEquivTable[i]];
	  }
	  (*numCCs) = ccCount;
	}
	
	// do second pass to map the CC labels to be sequential
	pDst = imgDst.dataPointer_u32();
	for(int y = 0, idx = 0; y < h; ++y){
	  for(int x = 0; x < w; ++x, ++idx){
	    pDst[idx] = rgEquivTable[pDst[idx]];
	  }
	}

	free(rgBlobIDs);

      }// end block within case statement
      break;
    default:
      fprintf(stderr, "DConnectedComponentLabeler::getCCimage_() not "
	      "implemented for image type %d\n", imgSrc.getImageType());
      exit(1);
  }

  free(rgEquivTable);
}

///Fill in DConnectedComponentInfo for each CC (including background)
/**The caller provides an already-allocated array, rgCCInfo, of
 * DConnectedComponentInfo objects for this function to fill in with
 * the information about the connected components in imgCCMap.
 * rgCCInfo must be large enough to hold at least numComponentsPlusBG
 * DConnectedComponentInfo objects, where numComponentsPlusBG is the
 * number of components in imgCCMap (including the background
 * component, which has label 0), therefore the number of
 * non-background components plus 1 slot for the bg component.
 */
void DConnectedComponentLabeler::getCCInfoFromCCimage(DImage &imgCCMap,
						      DConnectedComponentInfo
						      *rgCCInfo,
						      int numComponentsPlusBG){
  D_uint32 lbl;
  int w, h;
  D_uint32 *p32;
  
  w = imgCCMap.width();
  h = imgCCMap.height();
  p32 = imgCCMap.dataPointer_u32();
  for(int i = 0; i < numComponentsPlusBG; ++i){
    rgCCInfo[i].clear();
  }
  for(int y = 0, idx = 0; y < h; ++y){
    for(int x = 0; x < w; ++x, ++idx){
      lbl = p32[idx];
      if(lbl >= numComponentsPlusBG){
	fprintf(stderr, "DConnectedComponentLabeler::getCCInfoFromCCimage() "
		"numComponentsPlusBG=%d, but found label %d in imgCCMap\n",
		numComponentsPlusBG, lbl);
	abort();
      }
      if(0 == rgCCInfo[lbl].pixels){ // initialize the first time CC is seen
	rgCCInfo[lbl].label = lbl;
	rgCCInfo[lbl].bbLeft = rgCCInfo[lbl].bbRight = x;
	rgCCInfo[lbl].bbTop = rgCCInfo[lbl].bbBottom = y;
	rgCCInfo[lbl].startX = x;
	rgCCInfo[lbl].startY = y;
	rgCCInfo[lbl].pixels = 1;
	rgCCInfo[lbl].centroidX = x;
	rgCCInfo[lbl].centroidY = y;
      }	
      else{
	if(x < rgCCInfo[lbl].bbLeft)
	  rgCCInfo[lbl].bbLeft = x;
	if(x > rgCCInfo[lbl].bbRight)
	  rgCCInfo[lbl].bbRight = x;
	if(y < rgCCInfo[lbl].bbTop)
	  rgCCInfo[lbl].bbTop = y;
	if(y > rgCCInfo[lbl].bbBottom)
	  rgCCInfo[lbl].bbBottom = y;
	++(rgCCInfo[lbl].pixels);
	rgCCInfo[lbl].centroidX += x;
	rgCCInfo[lbl].centroidY += y;
      } /* end else */
    }//end for(x...
  }//end for(y...

  // finish calculating centroid of each CC
  for(int i=0; i < numComponentsPlusBG; ++i){
    if(rgCCInfo[i].pixels > 0){
      rgCCInfo[i].centroidX /= rgCCInfo[i].pixels;
      rgCCInfo[i].centroidY /= rgCCInfo[i].pixels;
    }
  }
}

DImage DConnectedComponentLabeler::getRGBImageFromCCImage(DImage &imgSrc,
							  bool fKeepNumbers){
  DImage imgDst;
  int w, h;
  D_uint32 *p32;
  D_uint8 *p8;
  const int NUMCOLORS = 33;
  D_uint8 rgRGBvals[NUMCOLORS*3] = {0,255,0,//greys,red intentionally skipped
				    0,0,255,
				    255,255,0,
				    255,0,255,
				    0,255,255,

				    255,192,255,
				    192,255,255,
				    255,255,192,
				    
				    128,192,192,
				    192,192,128,
				    192,128,192,
				    
				    64, 0, 128,
				    64,128,0,
				    128, 64, 0,
				    128,0,64,
				    0, 128, 64,
				    0,64,128,
				    
				    0, 128, 0,
				    0, 0, 128,
				    128,128,0,
				    128,0,128,
				    
				    255,128,0,
				    255,0,128,
				    128,255,0,
				    128,0,255,
				    0,128,255,
				    0,255,128,
				    
				    255,219,168,
				    168,255,219,
				    219,168,255,
				    219,255,168,
				    168,219,255,
				    255,168,255};
  
  w = imgSrc.width();
  h = imgSrc.height();
  imgDst.create(w, h, DImage::DImage_RGB);
  p8 = imgDst.dataPointer_u8();
  p32 = imgSrc.dataPointer_u32();
  if(fKeepNumbers){
    for(int y=0, idx=0; y < h; ++y){
      for(int x = 0; x < w; ++x, ++idx){
	if(0 == p32[idx]){
	  p8[idx*3] = p8[idx*3+1] = p8[idx*3+2] = 0;
	}
	else{
	  p8[idx*3] = (D_uint8)(p32[idx] >> 16);
	  p8[idx*3+1] = 0x80 | (D_uint8)(p32[idx] >> 8);
	  p8[idx*3+2] = (D_uint8)(p32[idx]);
	}
      }
    }
  }
  else{
    for(int y=0, idx=0; y < h; ++y){
      for(int x = 0; x < w; ++x, ++idx){
	D_uint8 colorIdx;
	if(0 == p32[idx]){
	  p8[idx*3] = p8[idx*3+1] = p8[idx*3+2] = 0;
	}
	else{
	  colorIdx = p32[idx] % NUMCOLORS;
	  p8[idx*3] = rgRGBvals[colorIdx*3];
	  p8[idx*3+1] = rgRGBvals[colorIdx*3+1];
	  p8[idx*3+2] = rgRGBvals[colorIdx*3+2];
	}
      }
    }
  }
  return imgDst;
}
