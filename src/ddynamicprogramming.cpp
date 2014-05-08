#include "ddynamicprogramming.h"
#include "dmath.h"
#include "dcolorspace.h"

//whichOne says which is the min.  Ties go to c so diagonal is preferred.
static inline double min3double(double a, double b, double c, int *whichOne){
  if(a < b){
    if(a < c){
      (*whichOne) = 1;
      return a;
    }
    (*whichOne) = 3;
    return c;
  }
  if(b < c){
    (*whichOne) = 2;
    return b;
  }
  (*whichOne) = 3;
  return c;
}






///Find the DP alignment cost (and alignment if desired) between 2 feature vects
/** This finds optimal alignment of fv1 to fv2 according to the
    Dynamic Time Warping algorithm used by Rath and Manmatha.  Rath
    and Manmatha used 15 for bandRadius, so we calculate that a good
    value may be avgWordWidth / 11.  If pathLen is set to not NULL, it
    will be filled with the length of the path.  If rgPath is not
    NULL, it must already be allocated with enough space to hold any
    possible DP path (can just set it to
    (1+fv1.vectLen)*(1+fv2.vectLen) for convenience.  The first
    *pathLen slots will be filled with the left-to-right mapping as
    follows: 0=match (pixel columns aligned); 1=move pixel over a
    column in A but not B; 2=move pixel over a column in B but not A;
    Both rgPrev and rgTable are used locally and will be
    allocated/freed locally if NULL.  If not NULL (for example if this
    function is being called a lot and you want to just reuse some
    memory), they must already be allocated to at least
    (1+fv1.vectLen)*(1+fv2.vectLen).

    TODO: nonDiagonalCost is empirically chosen right now.  The higher
    this cost, the more linear the DP warped will be because it makes
    it more costly to stretch or compress (go east or south) than to
    progress linearly through the DP table.
**/


double DDynamicProgramming::findDPAlignment(const DFeatureVector &fv1,
					    const DFeatureVector &fv2,
					    int bandRadius, double bandCost,
					    double nonDiagonalCost,
					    int *pathLen, int *rgPath,
					    int *rgPrev, double *rgTable){
  // parameter variables:
  //   int *rgPrev; /* table of which direction was best (previous) */
  //   double *rgTable; /* Dynamic Time Warping table */
  int i, j, k; /* i=table row, j=table col, k=feature number */
  int Wa, Wb; /* Wa=width of word a (same as feature len), likewise for b */
  int WaPlus1, WbPlus1, WbPlus2;
  double d; /* current cost of aligning features i of word a with features j of
	      word b.  This is added to minimum previoud cost (min of W,NW,N)*/
  int idx;
  int dims;
  int whichOne;
  bool fAllocTable = false;
  bool fAllocPrev = false;
  int pLen;// path length used locally, copied back to *pathLen if not NULL
  Wa = fv1.vectLen;
  Wb = fv2.vectLen;
  WaPlus1 = Wa+1;
  WbPlus1 = Wb+1;
  WbPlus2 = Wb+2;
  dims = fv1.dimensions;

  if(bandRadius < 2){
    fprintf(stderr, "Warning: DDynamicProgramming::findDPAlignment() changing "
	    "bandRadius to 2 (was set to %d)\n",bandRadius);
    bandRadius = 2;
  }
  if((fv1.dimensions) != (fv2.dimensions)){
    fprintf(stderr, "DDynamicProgramming::findDPAlignment() dimensions of fv1 "
	    "(%d) and fv2 (%d) must be equal!\n", fv1.dimensions,
	    fv2.dimensions);
    exit(1);
  }
  if((fv1.fDataIsDouble)!=(fv2.fDataIsDouble)){
    fprintf(stderr, "DDynamicProgramming::findDPAlignment() data types must "
	    "match!\n");
    exit(1);
  }
  if((fv1.fDataGroupedByDim)!=(fv2.fDataGroupedByDim)){
    fprintf(stderr, "DDynamicProgramming::findDPAlignment() data for fv1 and "
	    "fv2 must be both grouped by dimension or both interleaved.  They "
	    "currently don't match!\n");
    exit(1);
  }

  if(NULL == rgTable){
    rgTable = (double*)malloc(sizeof(double) * WaPlus1 * WbPlus1);
    D_CHECKPTR(rgTable);
    fAllocTable = true;
  }
  if(NULL == rgPrev){
    rgPrev = (int*)malloc(sizeof(int) * WaPlus1 * WbPlus1);
    D_CHECKPTR(rgPrev);
    fAllocPrev = true;
  }
    /* initialize the first row and column of table */

  for(j = 0; j < WbPlus1; ++j){
    rgPrev[j] = 1; /* prev is left (West) for top row */
    rgTable[j] = j + bandCost;
  }
  for(i = 0; i < WaPlus1; ++i){
    rgPrev[i*WbPlus1] = 2; /* prev is up (North) for left col */
    rgTable[i*WbPlus1] = i + bandCost;
  }
  rgTable[0] = 0;

  /* we are doing this one-based, even though the features are 0-based,
     that is why we have i-1 and j-1 in the feature index, because we
     have the initial row and col 0 initialized above so we don't have to
     check for row/col 0 every time through the loop */
  if(fv1.fDataGroupedByDim){
    if(fv1.fDataIsDouble){
      for(i = 1; i < WaPlus1; ++i){
	for(j = 1; j < WbPlus1; ++j){
	  d= 0.;
	  idx = i*WbPlus1+j;
	  for(k = 0; k < dims; ++k){
	    d += ((fv1.pDbl[Wa*k+i-1] - fv2.pDbl[Wb*k+j-1]) *
		  (fv1.pDbl[Wa*k+i-1] - fv2.pDbl[Wb*k+j-1]));
	  }
	  // rgTable[idx] = d + min3double(nonDiagonalCost + rgTable[idx-1]+
	  // 				checkBandCost(i-1,j,Wa,Wb,bandRadius,
	  // 					      bandCost),
	  // 				nonDiagonalCost + 
	  // 				rgTable[idx-WbPlus1] +
	  // 				checkBandCost(i,j-1,Wa,Wb, bandRadius,
	  // 					      bandCost),
	  // 				rgTable[idx-WbPlus2]+
	  // 				checkBandCost(i-1,j-1,Wa,Wb,bandRadius,
	  // 					      bandCost),
	  // 				&whichOne);
	  double bc;//bandcost
	  bc = checkBandCost(i,j,Wa,Wb,bandRadius, bandCost);
	  rgTable[idx] = d + min3double(nonDiagonalCost + rgTable[idx-1]+bc,
	  				nonDiagonalCost + rgTable[idx-WbPlus1] + bc,
	  				rgTable[idx-WbPlus2]+ bc,
	  				&whichOne);
	  rgPrev[idx] = whichOne;
	}//for j
      }//for i
    }//if(fv1.fDataIsDouble)
    else{// data is float, so use float arrays (pFlt) in DFeatureVectors
      fprintf(stderr,"NYI!!! (%s:%d)\n",__FILE__,__LINE__);
      exit(1);
      for(i = 1; i < WaPlus1; ++i){
	for(j = 1; j < WbPlus1; ++j){
	  d= 0.;
	  idx = i*WbPlus1+j;
	  for(k = 0; k < dims; ++k){
	    d += ((double)((fv1.pFlt[Wa*k+i-1] - fv2.pFlt[Wb*k+j-1]) *
			      (fv1.pFlt[Wa*k+i-1] - fv2.pFlt[Wb*k+j-1])));
	  }
	  rgTable[idx] = d + min3double(nonDiagonalCost +rgTable[idx-1]+
					checkBandCost(i-1,j,Wa,Wb,bandRadius,
						      bandCost),
					nonDiagonalCost +
					rgTable[idx-WbPlus1] +
					checkBandCost(i,j-1,Wa,Wb, bandRadius,
						      bandCost),
					rgTable[idx-WbPlus2]+
					checkBandCost(i-1,j-1,Wa,Wb,bandRadius,
						      bandCost),
					&whichOne);
	  rgPrev[idx] = whichOne;
	}//for j
      }//for i
    }//else
  }//if(fv1.fGroupedByDim)
  else{// data is not grouped by dimension.  interleaved instead
      fprintf(stderr,"NYI!!! (%s:%d)\n",__FILE__,__LINE__);
      exit(1);
    if(fv1.fDataIsDouble){
      for(i = 1; i < WaPlus1; ++i){
	for(j = 1; j < WbPlus1; ++j){
	  d= 0.;
	  idx = i*WbPlus1+j;
	  for(k = 0; k < dims; ++k){
	    d += ((fv1.pDbl[(i-1)*dims+k] - fv2.pDbl[(j-1)*dims+k]) *
		  (fv1.pDbl[(i-1)*dims+k] - fv2.pDbl[(j-1)*dims+k]));
	  }
	  rgTable[idx] = d + min3double(nonDiagonalCost +rgTable[idx-1]+
					checkBandCost(i-1,j,Wa,Wb,bandRadius,
						      bandCost),
					nonDiagonalCost +
					rgTable[idx-WbPlus1] +
					checkBandCost(i,j-1,Wa,Wb, bandRadius,
						      bandCost),
					rgTable[idx-WbPlus2]+
					checkBandCost(i-1,j-1,Wa,Wb,bandRadius,
						      bandCost),
					&whichOne);
	  rgPrev[idx] = whichOne;
	}//for j
      }//for i
    }//if(fv1.fDataIsDouble)
    else{// data is float, so use float arrays (pFlt) in DFeatureVectors
      for(i = 1; i < WaPlus1; ++i){
	for(j = 1; j < WbPlus1; ++j){
	  d= 0.;
	  idx = i*WbPlus1+j;
	  for(k = 0; k < dims; ++k){
	    d += ((double)((fv1.pFlt[(i-1)*dims+k] - fv2.pFlt[(j-1)*dims+k]) *
			       (fv1.pFlt[(i-1)*dims+k] - fv2.pFlt[(j-1)*dims+k])));
	  }
	  rgTable[idx] = d + min3double(nonDiagonalCost +rgTable[idx-1]+
					checkBandCost(i-1,j,Wa,Wb,bandRadius,
						      bandCost),
					nonDiagonalCost +
					rgTable[idx-WbPlus1] +
					checkBandCost(i,j-1,Wa,Wb, bandRadius,
						      bandCost),
					rgTable[idx-WbPlus2]+
					checkBandCost(i-1,j-1,Wa,Wb,bandRadius,
						      bandCost),
					&whichOne);
	  rgPrev[idx] = whichOne;
	}//for j
      }//for i
    }//else
  }
  if(NULL != pathLen)
    (*pathLen) = 0;
  pLen = 0;
  i = Wa;
  j = Wb;
  while((i != 0) || (j != 0)){
    if( (i<0) || (j<0)){
      fprintf(stderr,"ERROR! i=%d j=%d (%s:%d)\n",i,j,__FILE__,__LINE__);
      abort();
    }
    idx = i*WbPlus1+j;
    // ++(*pathLen);
    ++pLen;
    if(1 == rgPrev[idx]){ /* west (left) */
      --j;
    }
    else if(2 == rgPrev[idx]){ /* north (up) */
      --i;
    }
    else{ /* diagonal (NorthWest) */
      --j;
      --i;
    }
  }

  if(NULL != rgPath){// copy the backward path, but in reverse so it is forward
    int pathidx;
    i = Wa;
    j = Wb;
    // pathidx = (*pathLen)-1;
    pathidx = pLen-1;
    while((i != 0) || (j != 0)){
      if( (i<0) || (j<0)){
	fprintf(stderr,"ERROR! i=%d j=%d (%s:%d)\n",i,j,__FILE__,__LINE__);
	abort();
      }
      idx = i*WbPlus1+j;
      if(1 == rgPrev[idx]){ /* east (right) */
	if(pathidx < 0){
	  fprintf(stderr,"pathidx=%d\n",pathidx);
	  abort();
	}
	rgPath[pathidx] = 1;
	--j;
      }
      else if(2 == rgPrev[idx]){ /* south (down) */
	if(pathidx < 0){
	  fprintf(stderr,"pathidx=%d\n",pathidx);
	  abort();
	}
	rgPath[pathidx] = 2;
	--i;
      }
      else{ /* diagonal (SouthEast) */
	if(pathidx < 0){
	  fprintf(stderr,"pathidx=%d\n",pathidx);
	  abort();
	}
	rgPath[pathidx] = 0;
	--j;
	--i;
      }
      --pathidx;
    }
  }

  // if (0 == (*pathLen)){
  if (0 == pLen){
    fprintf(stderr, "ERROR! (%s:%d)\n", __FILE__, __LINE__);
    // (*pathLen) = 1;
    pLen = 1;
  }

  d = rgTable[WaPlus1*WbPlus1-1] / pLen;

  if(NULL != pathLen)
    (*pathLen) = pLen;

  if(fAllocTable)
    free(rgTable);
  if(fAllocPrev)
    free(rgPrev);
  return d;
}

//Do a piecewise linear warp of img1 to img2 using previously calculated rgPath
/**This function assumes that rgPath has been calculated on DFeatureVectors extracted from img1 and img2, respectively, and having the same length as the width (or height if fVertical is true) as img1 and img2. If a column is squished, then the min value of any column mapping to it is used.  This function is mainly for debug.**/
DImage DDynamicProgramming::piecewiseLinearWarpDImage(DImage &img1,
						      int warpToLength,
						      int pathLen, int *rgPath,
						      bool fVertical){
  DImage imgWarp;
  DImage::DImageType imgType;
  int w, h;
  D_uint8 *p8dst, *p8src;

  w = img1.width();
  h = img1.height();
  imgType = img1.getImageType();
  p8src = img1.dataPointer_u8();
  

  if((DImage::DImage_RGB != imgType) && (DImage::DImage_u8 != imgType)){
    fprintf(stderr, "DDynamicProgramming::piecewiseLinearWarpDImage() only "
	    "implemented for 8-bit GS and RGB images\n");
    exit(1);
  }

  if(!fVertical){
    //TODO: For just computing morph costs with distance maps, we may want to 
    //just create a sideways stretched image so memory writes will be
    //sequential?
    int x1, x2;
    x1 = x2 = 0;
    imgWarp.create(warpToLength, h, imgType);
    p8dst = imgWarp.dataPointer_u8();
    if(DImage::DImage_RGB == imgType){
      imgWarp.fill(255,255,255);
      for(int i=0; i < pathLen; ++i){
	if(0 == rgPath[i]){
	  for(int y = 0; y < h; ++y){
	    p8dst[3*(y*warpToLength+x2)] = p8src[3*(y*w+x1)];
	    p8dst[3*(y*warpToLength+x2)+1] = p8src[3*(y*w+x1)+1];
	    p8dst[3*(y*warpToLength+x2)+2] = p8src[3*(y*w+x1)+2];
	  }
	  ++x1;
	  ++x2;
	}
	else if(1 == rgPath[i]){//east (stretch)
	  for(int y = 0; y < h; ++y){
	    p8dst[3*(y*warpToLength+x2)] = p8src[3*(y*w+x1)];
	    p8dst[3*(y*warpToLength+x2)+1] = p8src[3*(y*w+x1)+1];
	    p8dst[3*(y*warpToLength+x2)+2] = p8src[3*(y*w+x1)+2];
	  }
	  ++x2;
	}
	else if(2 == rgPath[i]){//south (squish!)
	  for(int y = 0; y < h; ++y){
	    float H1, H2, S1, S2, V1, V2;
	    DColorSpace::getHSVFromRGB(p8src[3*(y*w+x1)],p8src[3*(y*w+x1)+1],
				       p8src[3*(y*w+x1)+2], &H1, &S1, &V1);
	    DColorSpace::getHSVFromRGB(p8dst[3*(y*warpToLength+x2)],
				       p8dst[3*(y*warpToLength+x2)+1],
				       p8dst[3*(y*warpToLength+x2)+2],
				       &H2, &S2, &V2);
	    if(V1 < V2){//use the darker of the color already there vs. new one
	      p8dst[3*(y*warpToLength+x2)] = p8src[3*(y*w+x1)];
	      p8dst[3*(y*warpToLength+x2)+1] = p8src[3*(y*w+x1)+1];
	      p8dst[3*(y*warpToLength+x2)+2] = p8src[3*(y*w+x1)+2];
	    }
	  }
	  ++x1;
	}
	else{
	    fprintf(stderr, "logic error! (%s:%d)\n", __FILE__, __LINE__);
	    exit(1);
	}
      }
    }
    else{
      imgWarp.fill(255.);
      for(int i=0; i < pathLen; ++i){
		if(0 == rgPath[i]){
		  for(int y = 0; y < h; ++y){
		    p8dst[y*warpToLength+x2] = p8src[y*w+x1];
		  }
		  ++x1;
		  ++x2;
		}
		else if(1 == rgPath[i]){//east (stretch)
		  for(int y = 0; y < h; ++y){
		    p8dst[y*warpToLength+x2] = p8src[y*w+x1];
		  }
		  ++x2;
		}
		else if(2 == rgPath[i]){//south (squish!)
		  for(int y = 0; y < h; ++y){
		    //use the darker of the color already there vs. new one
		    if(p8src[y*w+x1] < p8dst[y*warpToLength+x2]){
			p8dst[y*warpToLength+x2] = p8src[y*w+x1];
		    }
		  }
		  ++x1;
		}
		else{
		  fprintf(stderr, "logic error! (%s:%d)\n", __FILE__, __LINE__);
		  exit(1);
		}
		if (x1>=w)
		{
			//fprintf(stderr, "DP warping out of bounds x1=%d, w=%d (%s:%d)\n",x1,w, __FILE__, __LINE__);
			x1 = w-1;
		}
		if (x2>=warpToLength)
		{
			//fprintf(stderr, "DP warping out of bounds x2=%d, warpToLength=%d (%s:%d)\n",x2,warpToLength, __FILE__, __LINE__);
			x2 = warpToLength-1;
		}
      }//for i
    }//else
  }
  else{//verticle
    //TODO test
    //Underconstruction
    
	int y1, y2;
	y1 = y2 = 0;
	imgWarp.create(w, warpToLength, imgType);
	p8dst = imgWarp.dataPointer_u8();
	if(DImage::DImage_RGB == imgType){
		imgWarp.fill(255,255,255);
		for(int i=0; i < pathLen; ++i){
			if(0 == rgPath[i]){
			  for(int x = 0; x < w; ++x){
			    p8dst[3*(y2*w+x)] = p8src[3*(y1*w+x)];
			    p8dst[3*(y2*w+x)+1] = p8src[3*(y1*w+x)+1];
			    p8dst[3*(y2*w+x)+2] = p8src[3*(y1*w+x)+2];
			  }
			  ++y1;
			  ++y2;
			}
			else if(1 == rgPath[i]){//east (stretch)
			  for(int x = 0; x < w; ++x){
			    p8dst[3*(y2*w+x)] = p8src[3*(y1*w+x)];
			    p8dst[3*(y2*w+x)+1] = p8src[3*(y1*w+x)+1];
			    p8dst[3*(y2*w+x)+2] = p8src[3*(y1*w+x)+2];
			  }
			  ++y2;
			}
			else if(2 == rgPath[i]){//south (squish!)
			  for(int x = 0; x < w; ++x){
			    float H1, H2, S1, S2, V1, V2;
			    DColorSpace::getHSVFromRGB(p8src[3*(y1*w+x)],p8src[3*(y1*w+x)+1],
							 p8src[3*(y1*w+x)+2], &H1, &S1, &V1);
			    DColorSpace::getHSVFromRGB(p8dst[3*(y2*w+x)],
							 p8dst[3*(y2*w+x)+1],
							 p8dst[3*(y2*w+x)+2],
							 &H2, &S2, &V2);
			    if(V1 < V2){//use the darker of the color already there vs. new one
				p8dst[3*(y2*w+x)] = p8src[3*(y1*w+x)];
				p8dst[3*(y2*w+x)+1] = p8src[3*(y1*w+x)+1];
				p8dst[3*(y2*w+x)+2] = p8src[3*(y1*w+x)+2];
			    }
			  }
			  ++y1;
			}
			else{
			    fprintf(stderr, "logic error! (%s:%d)\n", __FILE__, __LINE__);
			    exit(1);
			}
		}
	}
	else{
		imgWarp.fill(255.);
		for(int i=0; i < pathLen; ++i){
			if(0 == rgPath[i]){
			  for(int x = 0; x < w; ++x){
			    p8dst[y2*w+x] = p8src[y1*w+x];
			  }
			  ++y1;
			  ++y2;
			}
			else if(1 == rgPath[i]){//east (stretch)
			  for(int x = 0; x < w; ++x){
			    p8dst[y2*w+x] = p8src[y1*w+x];
			  }
			  ++y2;
			}
			else if(2 == rgPath[i]){//south (squish!)
			  for(int x = 0; x < w; ++x){
			    //use the darker of the color already there vs. new one
			    if(p8src[y1*w+x] < p8dst[y2*w+x]){
				p8dst[y2*w+x] = p8src[y1*w+x];
			    }
			  }
			  ++y1;
			}
			else{
			  fprintf(stderr, "logic error! (%s:%d)\n", __FILE__, __LINE__);
			  exit(1);
			}
			if (y1>=h)
			{
				//fprintf(stderr, "DP warping out of bounds y1=%d, h=%d (%s:%d)\n",y1,h, __FILE__, __LINE__);
				y1 = h-1;
			}
			if (y2>=warpToLength)
			{
				//fprintf(stderr, "DP warping out of bounds y2=%d, warpToLength=%d (%s:%d)\n",y2,warpToLength, __FILE__, __LINE__);
				y2 = warpToLength-1;
			}
		}//for i
	}//else
  }
  return imgWarp;
}

/**Fills rgCoords[x] with the mapping from x in fvect0 to x' in fvect1
   Expects that rgCoords is already allocated to (w0+1)*sizeof(double) **/
void DDynamicProgramming::getCoord0MappingsToCoord1(int w0, int w1,
						    double *rgCoords,
						    int pathLen, int *rgPath){
  int x0, x1;
  x0 = x1 = 0;
  int i;
  int rgOffsX0[3] = {1,0,1};
  int rgOffsX1[3] = {1,1,0};

  //  printf("pathLen=%d\n",pathLen);
  // fprintf(stderr, "TODO:getCoord0MappingsToCoord1() use averages and interpolation for E,S.  \n");
  i = 0;
  while((0 == x0)&&(i<pathLen)){// first column is y=-1 since we had extra row/col in table
    x0 += rgOffsX0[rgPath[i]];
    x1 += rgOffsX1[rgPath[i]];
    ++i;
  }
  for(; i < pathLen; ++i){
    rgCoords[x0-1] = (x1>0)?(x1-1):0;
    //    printf("rgCoords[%d]=%.0f\n",x0-1,rgCoords[x0-1]);
    if((rgPath[i]<0)||(rgPath[i]>2)){
      fprintf(stderr, "logic error! (%s:%d)\n", __FILE__, __LINE__);
      exit(1);
    }
    x0 += rgOffsX0[rgPath[i]];
    x1 += rgOffsX1[rgPath[i]];
    if((x0<0)||(x0>w0)){
      fprintf(stderr,"ERROR! DDynamicProgramming::getCoord0MappingsToCoord1(): x0=%d but w0=%d\n", x0, w0);
      abort();
    }
    if((x1<0)||(x1>w1)){
      fprintf(stderr,"ERROR! DDynamicProgramming::getCoord0MappingsToCoord1(): x1=%d but w1=%d\n", x1, w1);
      abort();
    }
    // if(0 == rgPath[i]){
    //   ++x0;
    //   ++x1;
    // }
    // else if(1 == rgPath[i]){//east (stretch)
    //   ++x1;
    // }
    // else if(2 == rgPath[i]){//south (squish!)
    //   ++x0;
    // }
    // else{
    //   fprintf(stderr, "logic error! (%s:%d)\n", __FILE__, __LINE__);
    //   exit(1);
    // }
    // if(x0 < w0){
    //   if(x1 < w1){
    // 	rgCoords[x0-1] = x1-1;
    //   }
    //   else{
    // 	fprintf(stderr,"ERROR! DDynamicProgramming::getCoord0MappingsToCoord1(): x1=%d but w1=%d\n", x1, w1);
    //   }
    // }
    // else
    //   fprintf(stderr,"ERROR! DDynamicProgramming::getCoord0MappingsToCoord1(): x0=%d but w0=%d\n", x0, w0);
  }//for i
  // printf("rgCoords[0]=%lf changing to 0\n",rgCoords[0]);
  //  rgCoords[0] = 0;// anchor the beginning


  while(x0 <= w0){ //anchor the end to last point found
    rgCoords[x0-1] = (x1>0)?(x1-1):0;
    ++x0;
  }
}


/**Fills rgCoords[x] with the mapping from x' in fvect1 to x in fvect0
   Expects that rgCoords is already allocated to w1*sizeof(double) **/
void DDynamicProgramming::getCoord1MappingsToCoord0(int w1, double *rgCoords,
						    int pathLen, int *rgPath){
  int x0, x1;
  x0 = x1 = 0;

  // fprintf(stderr, "TODO:getCoord1MappingsToCoord0() use averages and interpolation for E,S\n");
  
  for(int i=0; i < pathLen; ++i){
    if(x1 <= w1)
      rgCoords[x1] = x0;
    else
      fprintf(stderr,"ERROR! DDynamicProgramming::getCoord1MappingsToCoord0(): x1=%d but w1=%d\n", x1, w1);
    if(0 == rgPath[i]){
      ++x0;
      ++x1;
    }
    else if(1 == rgPath[i]){//east (stretch)
      ++x1;
    }
    else if(2 == rgPath[i]){//south (squish!)
      ++x0;
    }
    else{
      fprintf(stderr, "logic error! (%s:%d)\n", __FILE__, __LINE__);
      exit(1);
    }
  }//for i
}


DImage DDynamicProgramming::visualizeDPTable(int *rgPrev, int Wa, int Wb,
					     int bandRadius){
  fprintf(stderr, "DDynamicProgramming::visualizeDPTable() NYI!\n");
  exit(1);
}


void DDynamicProgramming::debugImages(int Wa, int Wb, int pathLen,
				      int *rgPath, int *rgPrev,
				      double *rgTable){
  DImage imgPath;
  DImage imgPrev;
  DImage imgTable;
  DImage imgMapping0to1;
  DImage imgMapping1to0;
  DImage imgPathForFig;//shows path and sakoe-chiba band of radius 15
  int i, j, k; /* i=table row, j=table col, k=feature number */
  //int Wa, Wb; /* Wa=width of word a (same as feature len), likewise for b */
  int WaPlus1, WbPlus1;
  D_uint8 *p8;

  WaPlus1 = Wa+1;
  WbPlus1 = Wb+1;

  imgPrev.create(WbPlus1, WaPlus1, DImage::DImage_RGB);
  imgTable.create(WbPlus1, WaPlus1, DImage::DImage_RGB);
  imgTable.clear();

#if 1
  //prev
  imgPrev.clear();
  p8 = imgPrev.dataPointer_u8();
  for(int y=0, idx=0; y < WaPlus1; ++y){
    for(int x=0; x < WbPlus1; ++x, ++idx){
      if(rgPrev[idx] == 1){//west
	imgPrev.drawPixel(x,y,255,0,0);
      }
      else if(rgPrev[idx] == 2){//up
	imgPrev.drawPixel(x,y,0,0,255);
      }
      else{//diag
      }
    }
  }
  imgPrev.save("/tmp/prev.ppm");
#endif


#if 1
  //path

  printf("pathLen=%d path:",pathLen);
  for(int i=0; i < pathLen; ++i)
    printf("%d", rgPath[i]);
  printf("\n");

  imgPath = imgPrev;
  int x0, x1;
  x0 = x1 = 0;
  for(int i=0; i < pathLen; ++i){
    if(0 == rgPath[i]){
      ++x0;
      ++x1;
    }
    else if(1 == rgPath[i]){//east (stretch)
      ++x1;
    }
    else if(2 == rgPath[i]){//south (squish!)
      ++x0;
    }
    else{
      fprintf(stderr, "logic error! (%s:%d)\n", __FILE__, __LINE__);
      exit(1);
    }
    if((x0<WaPlus1)&&(x1<WbPlus1))
      imgPath.drawPixel(x1,x0,255,255,0,0.5);
    else
      printf("DDynamicProgramming::debugImages() pixel %d,%d out of bounds "
	     "(Wa=%d,Wb=%d)\n",x0,x1,Wa,Wb);
  }//for i
  imgPath.save("/tmp/path.ppm");
#endif

#if 1
  //table
  for(int y=0, idx=0; y < WaPlus1; ++y){
    for(int x=0; x < WbPlus1; ++x, ++idx){
      int val;
      val = (int)(30*rgTable[idx]);
      if(val < 0)
	imgTable.drawPixel(x,y,255,0,0);
      else if(val <= 255)
	imgTable.drawPixel(x,y,val,0,0);
      else if(val <= 255*2)
	imgTable.drawPixel(x,y,0,val,0);
      else if(val <= 255*3)
	imgTable.drawPixel(x,y,0,0,val);
      else if(val <= 255*4)
	imgTable.drawPixel(x,y,val,0,val);
      else if(val <= 255*5)
	imgTable.drawPixel(x,y,val,val,0);
      else if(val <= 255*6)
	imgTable.drawPixel(x,y,0,val,val);
      else if(val <= 255*7)
	imgTable.drawPixel(x,y,val,val,val);
      else
	imgTable.drawPixel(x,y,255,255,255);
    }
  }
  imgTable.save("/tmp/table.ppm");
#endif


#if 1
  //path for figure

  printf("pathLen=%d path:",pathLen);
  for(int i=0; i < pathLen; ++i)
    printf("%d", rgPath[i]);
  printf("\n");

  imgPathForFig = imgPrev;
  imgPathForFig.fill(0,0,0);

  for(int y=0, idx=0; y < WaPlus1; ++y){
    for(int x=0; x < WbPlus1; ++x, ++idx){
      int val;
      val = (int)(rgTable[idx]);
      if((val >= 999.) || (val < 0.))
	imgPathForFig.drawPixel(x,y,255,255,255,0.0);
      else{
	imgPathForFig.drawPixel(x,y,164,164,164,0.0);
      }
    }
  }

  //int x0, x1;
  x0 = x1 = 0;
  for(int i=0; i < pathLen; ++i){
    if(0 == rgPath[i]){
      ++x0;
      ++x1;
    }
    else if(1 == rgPath[i]){//east (stretch)
      ++x1;
    }
    else if(2 == rgPath[i]){//south (squish!)
      ++x0;
    }
    else{
      fprintf(stderr, "logic error! (%s:%d)\n", __FILE__, __LINE__);
      exit(1);
    }
    if((x0<WaPlus1)&&(x1<WbPlus1))
      imgPathForFig.drawPixel(x1,x0,0,0,0,0.0);
    else
      printf("DDynamicProgramming::debugImages() pixel %d,%d out of bounds "
	     "(Wa=%d,Wb=%d)\n",x0,x1,Wa,Wb);
  }//for i

  imgPathForFig.save("/tmp/pathForFig.ppm");
#endif


#if 1
  //DP mapping from 0 to 1
  double *rgCoords;
  
  rgCoords = (double*)malloc(sizeof(double)*WaPlus1);
  D_CHECKPTR(rgCoords);
  getCoord0MappingsToCoord1(Wa,Wb,rgCoords,pathLen, rgPath);
  DImage imgMap0to1;

  imgMap0to1.create((Wa>Wb)?(Wa+1):(Wb+1), 100, DImage::DImage_RGB);
  imgMap0to1.clear();
  int rgColors[3*14]={255,0,0,  0,255,0,  0,0,255,  255,255,0,  255,0,255,  
		      0,255,255,  255,255,255, 128,0,0,  0,128,0,  0,0,128,
		      128,128,0,  128,0,128,  0,128,128,  128,128,128};
  // printf("makeing map0to1.ppm with Wa=%d\n",Wa);
  for(int x=0; x < Wa; ++x){
    int rval,gval,bval;
    // printf("  x=%d",x); fflush(stdout);
    // printf("  rgCoords[%d]=%f\n",x,rgCoords[x]);
    rval = rgColors[(x%14)*3];
    gval = rgColors[(x%14)*3+1];
    bval = rgColors[(x%14)*3+2];
    imgMap0to1.drawLine(x,0,rgCoords[x],99, rval,gval,bval, 0.2);
  }
  imgMap0to1.save("/tmp/map0to1.ppm");
#endif
}

