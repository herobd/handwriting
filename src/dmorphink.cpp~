#include "dmorphink.h"
#include "dprofile.h"
#include "dfeaturevector.h"
#include "ddynamicprogramming.h"
#include "ddistancemap.h"
#include "dmedialaxis.h"
#include "dtimer.h"
#include "dwordfeatures.h"
#include <string.h>

#include <signal.h>

#define SUBSAMPLE 1
#define SAVE_IMAGES 0
#define NORMALIZE_VERT_PROF 0
#define USE_DISTANCE_FROM_DP_COST 0
#define OUTPUT_DEBUG_FINAL_COST 0

#define MAXDISTCOLORS 500 /*for heatmaps*/

int DEBUG_stopAtStep = 1;

//#define COMPENSATE_FOR_OOB_DISTMAP 1 //we now do this all the time

DMorphInk::DMorphInk(){
  avgMovement0to1 = 0.;
  avgMovement1to0 = 0.;
  numPointRows = 0;
  numPointCols = 0;
  rgPoints0X = NULL;
  rgPoints0Y = NULL;
  rgPlacePoints0 = NULL;
  rgPoints1X = NULL;
  rgPoints1Y = NULL;
  rgPointsDPX = NULL;
  rgPointsDPY = NULL;
  rgPointsPrevX = NULL;
  rgPointsPrevY = NULL;
  ma0prevPosX = NULL;
  ma0prevPosY = NULL;
  rgMA0X = NULL;
  rgMA0Y = NULL;
  rgMA1X = NULL;
  rgMA1Y = NULL;
#if NEW_WARP
  rgMA0r = NULL;
  rgMA0c = NULL;
  rgMA0s = NULL;
  rgMA0t = NULL;
  rgMA0quadW = NULL;
  rgMA0quadH = NULL;
  rgTempMA0idx = NULL;
#endif
  rgTempMA0X = NULL;
  rgTempMA0Y = NULL;
  lenMA0=0;
  lenMA1=0;

  pimg0 = NULL;
  pimg1 = NULL;
  w0 = h0 = w1 = h1 = 0;
  meshLevel = 0;
  initialColSpacing = 0;
  initialRowSpacing = 0;
  curColSpacing = 0.;
  curRowSpacing = 0.;
  warpCostDP = 0.;
  warpCostDPhoriz = 0.;
  fOnlyDoOneDirection = false;
  fOnlyDoCoarseAlignment = false;
  fFaintMesh = false;
}

DMorphInk::~DMorphInk(){
  if(numPointRows>0){
    free(rgPoints0X);
    free(rgPoints0Y);
    free(rgPlacePoints0);
    free(rgPoints1X);
    free(rgPoints1Y);
    free(rgPointsDPX);
    free(rgPointsDPY);
    free(rgPointsPrevX);
    free(rgPointsPrevY);
  }
  if(rgMA0X != NULL){
    free(rgMA0X);
    free(rgMA0Y);
    free(rgMA1X);
    free(rgMA1Y);
#if NEW_WARP
    free(rgMA0r);
    free(rgMA0c);
    free(rgMA0s);
    free(rgMA0t);
    free(rgMA0quadW);
    free(rgMA0quadH);
    free(rgTempMA0idx);
#endif
    free(rgDPMA0X);
    free(rgDPMA0Y);
    free(rgTempMA0X);
    free(rgTempMA0Y);
  }
  rgPoints0X = rgPoints0Y = rgPoints1X = rgPoints1Y =
    rgPointsDPX = rgPointsDPY = rgPointsPrevX = rgPointsPrevY = NULL;
  rgPlacePoints0 = NULL;
  if(NULL != ma0prevPosX)
    free(ma0prevPosX);
  if(NULL != ma0prevPosY)
    free(ma0prevPosY);
  ma0prevPosX = ma0prevPosY = NULL;
  numPointRows = numPointCols = 0;
  pimg0 = pimg1 = NULL;
  w0 = h0 = w1 = h1 = 0;
  meshLevel = 0;
  initialColSpacing = 0;
  initialRowSpacing = 0;
  curColSpacing = 0.;
  curRowSpacing = 0.;
}

void DMorphInk::init(const DImage &src0, const DImage &src1, bool fMakeCopies,
		     int initialMeshSpacing, int bandRadius,
		     double nonDiagonalDPcost){
  w0 = src0.width();
  h0 = src0.height();
  w1 = src1.width();
  h1 = src1.height();


  if((w0 < 1) || (h0 < 1)){
    fprintf(stderr, "DMorphInk::init() called with image0 width %d or height %d "
	    "less than 1\n", w0, h0);
    abort();
  }
  if((w1 < 1) || (h1 < 1)){
    fprintf(stderr, "DMorphInk::init() called with image0 width %d or height %d "
	    "less than 1\n", w0, h0);
    abort();
  }

  if(fMakeCopies){
    img0tmp = src0;
    img1tmp = src1;
    pimg0 = &img0tmp;
    pimg1 = &img1tmp;
  }
  else{
    pimg0 = (DImage*)&src0;
    pimg1 = (DImage*)&src1;
  }

  //create distance maps for each image
  // imgDist0 = DDistanceMap::getDistFromInkBitonal(*pimg0,1000,-1000);
  // imgDist1 = DDistanceMap::getDistFromInkBitonal(*pimg1,1000,-1000);
  DDistanceMap::getDistFromInkBitonal_(imgDist0, *pimg0,1000,-1000);
  DDistanceMap::getDistFromInkBitonal_(imgDist1, *pimg1,1000,-1000);
#if SAVE_IMAGES
  {
    DImage imgD;
    imgD = imgDist0.convertedImgType(DImage::DImage_u8);
    printf("saving /tmp/dist0.pgm\n");
    imgD.save("/tmp/dist0.pgm");
    printf("saving /tmp/dist1.pgm\n");
    imgD = imgDist1.convertedImgType(DImage::DImage_u8);
    imgD.save("/tmp/dist1.pgm");
  }
#endif
  //don't leave memory hanging if we call init() more than once:
  if(rgMA0X != NULL){
    free(rgMA0X);
    free(rgMA0Y);
    free(rgMA1X);
    free(rgMA1Y);
#if NEW_WARP
    free(rgMA0r);
    free(rgMA0c);
    free(rgMA0s);
    free(rgMA0t);
    free(rgMA0quadW);
    free(rgMA0quadH);
    free(rgTempMA0idx);
#endif
    free(rgDPMA0X);
    free(rgDPMA0Y);
    free(rgTempMA0X);
    free(rgTempMA0Y);
  }

  imgMA0 = imgInk0 = DMedialAxis::getMedialAxisImageFromDistMap(imgDist0,true,
								&lenMA0);
  imgMA1 = DMedialAxis::getMedialAxisImageFromDistMap(imgDist1,true,&lenMA1);
#if SAVE_IMAGES
  printf("saving /tmp/medial_axis1.pgm\n");
  imgMA1.save("/tmp/medial_axis1.pgm");
  printf("saving /tmp/medial_axis0.pgm\n");
  imgMA0.save("/tmp/medial_axis0.pgm");
#endif


  //make lists of MA0 pixels
  rgMA0X = (double*)malloc(sizeof(double)*(1+lenMA0));
  D_CHECKPTR(rgMA0X);
  rgMA0Y = (double*)malloc(sizeof(double)*(1+lenMA0));
  D_CHECKPTR(rgMA0Y);
#if NEW_WARP
  rgMA0r = (int*)malloc(sizeof(double)*(1+lenMA0));
  D_CHECKPTR(rgMA0r);
  rgMA0c = (int*)malloc(sizeof(double)*(1+lenMA0));
  D_CHECKPTR(rgMA0c);
  rgMA0s = (double*)malloc(sizeof(double)*(1+lenMA0));
  D_CHECKPTR(rgMA0s);
  rgMA0t = (double*)malloc(sizeof(double)*(1+lenMA0));
  D_CHECKPTR(rgMA0t);
  rgMA0quadW = (double*)malloc(sizeof(double)*(1+lenMA0));
  D_CHECKPTR(rgMA0quadW);
  rgMA0quadH = (double*)malloc(sizeof(double)*(1+lenMA0));
  D_CHECKPTR(rgMA0quadH);
  rgTempMA0idx = (int*)malloc(sizeof(double)*(1+lenMA0));
  D_CHECKPTR(rgTempMA0idx);
#endif

  rgTempMA0X = (double*)malloc(sizeof(double)*lenMA0);
  D_CHECKPTR(rgTempMA0X);
  rgTempMA0Y = (double*)malloc(sizeof(double)*lenMA0);
  D_CHECKPTR(rgTempMA0Y);

  rgDPMA0X = (double*)malloc(sizeof(double)*(1+lenMA0));
  D_CHECKPTR(rgDPMA0X);
  rgDPMA0Y = (double*)malloc(sizeof(double)*(1+lenMA0));
  D_CHECKPTR(rgDPMA0Y);
  D_uint8 *p0;
  p0 = imgMA0.dataPointer_u8();
  for(int y=0, idx=0, ma_pt=0; y < h0; ++y){
    for(int x=0; x<w0; ++x, ++idx){
      if(p0[idx] > 0){
	rgMA0X[ma_pt] = x;
	rgMA0Y[ma_pt] = y;
	++ma_pt;//MedialAxis point
      }
    }
  }
  //make lists of MA1 pixels
  rgMA1X = (double*)malloc(sizeof(double)*(1+lenMA1));
  D_CHECKPTR(rgMA1X);
  rgMA1Y = (double*)malloc(sizeof(double)*(1+lenMA1));
  D_CHECKPTR(rgMA1Y);
  D_uint8 *p1;
  p1 = imgMA1.dataPointer_u8();
  for(int y=0, idx=0, ma_pt=0; y < h1; ++y){
    for(int x=0; x<w1; ++x, ++idx){
      if(p1[idx] > 0){
	rgMA1X[ma_pt] = x;
	rgMA1Y[ma_pt] = y;
	++ma_pt;//MedialAxis point
      }
    }
  }



  // //TODO::get medial axis image of img0
  // fprintf(stderr,"TODO: try just using the medial axis instead of all ink\n");
  // imgInk0.create(w0, h0, DImage::DImage_u8);
  // D_uint8 *psrc, *pdst;
  // psrc = pimg0->dataPointer_u8();
  // pdst = imgInk0.dataPointer_u8();
  // for(int y=0, idx=0; y < h0; ++y){
  //   for(int x=0; x < w0; ++x, ++idx){
  //     if(psrc[idx] == 0x00)
  // 	pdst[idx] = 0x01;
  //     else
  // 	pdst[idx] = 0x00;
  //   }
  // }

  imgMA1.invertGrayscale();
  // imgDist1 = DDistanceMap::getDistFromInkBitonal(imgMA1,1000,-1000);
  DDistanceMap::getDistFromInkBitonal_(imgDist1, imgMA1,1000,-1000);
  imgMA1.invertGrayscale();
#if SAVE_IMAGES
  imgMA1.save("/tmp/imgMA1.ppm");
  imgDist1.save("/tmp/imgDist1.ppm");
#endif

  resetMeshes(initialMeshSpacing,initialMeshSpacing,
	      bandRadius,nonDiagonalDPcost);
  // printf("init() has run to just after resetMeshes().  calling exit()\n");
  // exit(1);
}


//This is essentially the initail alignment of meshes (?)
///reset meshes. uses col/row spacing on mesh0, uses DP to align mesh1
/**The last column and row may be different size than the rest so we
   can stay on pixel boundaries in mesh0. If -1 is specified, the
   initialColSpacing or initialRowSpacing (respectively) will be used
   - this assumes that they have been previously set, of course.*/
void DMorphInk::resetMeshes(int columnSpacing, int rowSpacing,
			    int bandRadius, double nonDiagonalDPcost){
  int numPoints;
  //int lastColW;
  int numMeshPointCols;
  int numMeshPointRows;

  if(-1==columnSpacing)
    columnSpacing = initialColSpacing;
  if(-1==rowSpacing)
    rowSpacing = initialRowSpacing;

  numMeshPointCols = 1 + (w0+columnSpacing-1)/columnSpacing;
  numMeshPointRows = 1 + (h0+rowSpacing-1)/rowSpacing;

  // release memory if already being used
  if(numPointRows>0){
    free(rgPoints0X);
    free(rgPoints0Y);
    free(rgPlacePoints0);
    free(rgPoints1X);
    free(rgPoints1Y);
    free(rgPointsDPX);
    free(rgPointsDPY);
    free(rgPointsPrevX);
    free(rgPointsPrevY);
  }
  numPoints = numMeshPointCols * numMeshPointRows;
  // printf("A: numPoints=%d\n",numPoints);
  //allocate memory
  rgPoints0X = (double*)malloc(sizeof(double)*numPoints);
  D_CHECKPTR(rgPoints0X);
  rgPoints0Y = (double*)malloc(sizeof(double)*numPoints);
  D_CHECKPTR(rgPoints0Y);
  rgPlacePoints0 = (bool*)malloc(sizeof(bool)*numPoints);
  D_CHECKPTR(rgPlacePoints0);
  rgPoints1X = (double*)malloc(sizeof(double)*numPoints);
  D_CHECKPTR(rgPoints1X);
  rgPoints1Y = (double*)malloc(sizeof(double)*numPoints);
  D_CHECKPTR(rgPoints1Y);
  rgPointsDPX = (double*)malloc(sizeof(double)*numPoints);
  D_CHECKPTR(rgPointsDPX);
  rgPointsDPY = (double*)malloc(sizeof(double)*numPoints);
  D_CHECKPTR(rgPointsDPY);
  rgPointsPrevX = (double*)malloc(sizeof(double)*numPoints);
  D_CHECKPTR(rgPointsPrevX);
  rgPointsPrevY = (double*)malloc(sizeof(double)*numPoints);
  D_CHECKPTR(rgPointsPrevY);
  // ma0prevPosX = (double*)malloc(sizeof(double)*lenMA0);
  // D_CHECKPTR(ma0prevPosX);
  // ma0prevPosY = (double*)malloc(sizeof(double)*lenMA0);
  // D_CHECKPTR(ma0prevPosY);
  // for(int i=0; i < lenMA0; ++i){
  //   ma0prevPosX[i] = ma0prevPosY[i] = -999.;
  // }

  numPointRows=numMeshPointRows;
  numPointCols=numMeshPointCols;
  //set the x,y position of each control point in mesh0
  for(int py=0, idx=0; py < numPointRows; ++py){
    for(int px = 0; px < numPointCols; ++px, ++idx){
      rgPoints0X[idx] = px * columnSpacing;
      if(rgPoints0X[idx] >= w0)//last column may not be as wide
	rgPoints0X[idx] = w0-1;
      rgPoints0Y[idx] = py * rowSpacing;
      if(rgPoints0Y[idx] >= h0)//last row may not be as tall
	rgPoints0Y[idx] = h0-1;
      // rgPointsDPX[idx] =  rgPoints0X[idx];
      // rgPointsDPY[idx] =  rgPoints0Y[idx];
       rgPointsDPX[idx] =  10*px;
       rgPointsDPY[idx] =  10*py;
    }
  }
  //do DP x-alignment of img0 to img1 to decide how to set warp1 x-coords
  //  fprintf(stderr, "TODO: clean up this code!  variable names confusing,etc!\n");
  DProfile prof1, prof2;

  DFeatureVector fv1, fv2;
  fv1 = DWordFeatures::extractWordFeatures(*pimg0,true,false,true,true,127);
  fv2 = DWordFeatures::extractWordFeatures(*pimg1,true,false,true,true,127);
  double dblDPcost;
  int pathLen;
  int *rgPath;
  int *rgPrev;
  double *rgTable;
  rgPath = (int*)malloc(sizeof(int)*(w0+2)*(w1+2));
  D_CHECKPTR(rgPath);
  rgPrev = (int*)malloc(sizeof(int)*(w0+2)*(w1+2));
  D_CHECKPTR(rgPrev);
  rgTable = (double*)malloc(sizeof(double)*(w0+2)*(w1+2));
  D_CHECKPTR(rgTable);

  // printf("fv1.vectLen=%d fv2.vectLen=%d\n", fv1.vectLen, fv2.vectLen);

  dblDPcost = DDynamicProgramming::findDPAlignment(fv1, fv2, bandRadius,1000.,
						   nonDiagonalDPcost, &pathLen,
						   rgPath, rgPrev, rgTable);
  // {
  //   DDynamicProgramming::debugImages(w0,w1,pathLen,rgPath,rgPrev,rgTable);
  //   // features of w0:
  //   for(int k=0; k < 4; ++k){
  //     FILE *ffv;
  //     char stFFVfile[1024];
  //     sprintf(stFFVfile,"/tmp/fv1_%d.dat",k);
  //     ffv = fopen(stFFVfile,"wb");
  //     if(!ffv){fprintf(stderr,"badd file open!\n");exit(1);}
  //     for(int x=0; x<w0; ++x){
  // 	fprintf(ffv,"%d %lf\n",x,fv1.pDbl[x+k*fv1.vectLen]);
  //     }
  //     fclose(ffv);
  //   }
  // }
  // fprintf(stderr,"done calling debugImages()\n");
  //  printf("horizontal DP cost was %lf\n", dblDPcost);
  warpCostDP = dblDPcost;
  double *rgXMappings0to1;
  rgXMappings0to1 = (double*)malloc(sizeof(double)*(w0+1));
  D_CHECKPTR(rgXMappings0to1);
  DDynamicProgramming::getCoord0MappingsToCoord1(w0,w1,rgXMappings0to1,pathLen,
						 rgPath);
#if SAVE_IMAGES
  {
    DImage imgDPwarp;
    imgDPwarp = DDynamicProgramming::piecewiseLinearWarpDImage(*pimg0, w1,
							       pathLen, rgPath,
							       false);
    imgDPwarp.save("/tmp/dpwarp.pgm");
  }
#endif


  //do DP y-alignment of img0 to img1 to decide how to set warp1 y-coords
  //The features are built JUST ON THE VERTICAL PROFILE
  DProfile Vprof1, Vprof2;
  Vprof1.getImageVerticalProfile(*pimg0);
  Vprof2.getImageVerticalProfile(*pimg1);
  double *pVProf1, *pVProf2;
  pVProf1 = Vprof1.dataPointer();
  pVProf2 = Vprof2.dataPointer();

#if NORMALIZE_VERT_PROF
  {
    double max1=0., max2=0.;
    for(int i=0, len=Vprof1.dataLen(); i < len; ++i){
      if(pVProf1[i] > max1)
	max1 = pVProf1[i];
    }
    if(max1 != 0.){
      for(int i=0, len=Vprof1.dataLen(); i < len; ++i){
	pVProf1[i] /= max1;
      }
    }

    for(int i=0, len=Vprof2.dataLen(); i < len; ++i){
      if(pVProf2[i] > max1)
	max2 = pVProf2[i];
    }
    if(max2 != 0.){
      for(int i=0, len=Vprof2.dataLen(); i < len; ++i){
	pVProf2[i] /= max2;
      }
    }

  }

#endif

  DFeatureVector Vfv1, Vfv2;
  Vfv1.setData_dbl(pVProf1, h0, 1, true, true, true);
  Vfv2.setData_dbl(pVProf2, h1, 1, true, true, true);
  double dblDPcostV;
  int pathLenV;
  int *rgPathV;
  int *rgPrevV;
  double *rgTableV;
  rgPathV = (int*)malloc(sizeof(int)*(h0+2)*(h1+2));
  rgPrevV = (int*)malloc(sizeof(int)*(h0+2)*(h1+2));
  rgTableV = (double*)malloc(sizeof(double)*(h0+2)*(h1+2));
  dblDPcostV = DDynamicProgramming::findDPAlignment(Vfv1, Vfv2,
						    bandRadius*2/3+1, 1000.,
						    nonDiagonalDPcost,&pathLenV,
						    rgPathV, rgPrevV, rgTableV);
  // {
  //   DDynamicProgramming::debugImages(w0,w1,pathLen,rgPath,rgPrev,rgTable);
  //   // features of w0:
  //   for(int k=0; k < 1; ++k){
  //     FILE *ffv;
  //     char stFFVfile[1024];
  //     sprintf(stFFVfile,"/tmp/Vfv1_%d.dat",k);
  //     ffv = fopen(stFFVfile,"wb");
  //     if(!ffv){fprintf(stderr,"badd file open!\n");exit(1);}
  //     for(int x=0; x<h0; ++x){
  // 	fprintf(ffv,"%d %lf\n",x,Vfv1.pDbl[x+k*fv1.vectLen]);
  //     }
  //     fclose(ffv);
  //   }
  // }
  warpCostDPv = dblDPcostV;
  //  printf("vertical DP cost was %lf\n", dblDPcostV);
  
  double *rgYMappings0to1;
  rgYMappings0to1 = (double*)malloc(sizeof(double)*(h0+1));
  D_CHECKPTR(rgYMappings0to1);
  DDynamicProgramming::getCoord0MappingsToCoord1(h0,h1,rgYMappings0to1,pathLenV,
						 rgPathV);

  //TODO: fix this so I don't get bad mappings.  I commented out warnings below
  // for now
  //now use DP mappings to set control points of warp1
  for(int py=0, idx=0; py < numPointRows; ++py){
    for(int px = 0; px < numPointCols; ++px, ++idx){
      if( (((int)rgPoints0X[idx]) < 0) ||(((int)rgPoints0X[idx])>=w0)){
	fprintf(stderr,"error! (%s:%d) value=%d w0=%d\n",__FILE__,__LINE__,
		(int)rgPoints0X[idx],w0);
	abort();
      }
      rgPoints1X[idx] = rgXMappings0to1[(int)rgPoints0X[idx]];
           // rgPointsDPX[idx] = rgPoints1X[idx];
            // rgPointsDPX[idx] = rgPoints0X[idx];
      if(rgPoints1X[idx] >= w1){//in case of bad mapping
	fprintf(stderr, "DMorphInk::resetMeshes() bad x mapping! changing x from %f to %d(%s:%d) w0=%d w1=%d px=%d py=%d idx=%d\n",
	 	rgPoints1X[idx],w1-1,__FILE__,__LINE__,w0,w1,px,py,idx);
	rgPoints1X[idx] = w1-1;
	exit(1);
      }
      if( (((int)rgPoints0Y[idx]) < 0) ||(((int)rgPoints0Y[idx])>=h0)){
	fprintf(stderr,"error! (%s:%d) value=%d h0=%d\n",__FILE__,__LINE__,
		(int)rgPoints0Y[idx],h0);
	abort();
      }
      rgPoints1Y[idx] = rgYMappings0to1[(int)rgPoints0Y[idx]];
           // rgPointsDPY[idx] = rgPoints1Y[idx];
            // rgPointsDPY[idx] = rgPoints0Y[idx];
      if(rgPoints1Y[idx] >= h1){//in case of bad mapping
	rgPoints1Y[idx] = h1-1;
	fprintf(stderr, "DMorphInk::resetMeshes() bad y mapping! changing y from %f to %d(%s:%d)\n",
	 	rgPoints1Y[idx],h1-1,__FILE__,__LINE__);
	exit(1);
      }
    }
  }

  //assign the Medial Axis points their corresponding quad indexes
  //assignMA0qidxs();
  // free memory used for DP warping, etc.
  free(rgPath);
  free(rgPrev);
  free(rgTable);
  free(rgPathV);
  free(rgPrevV);
  free(rgTableV);
  free(rgXMappings0to1);
  free(rgYMappings0to1);
  this->initialColSpacing = columnSpacing;
  this->initialRowSpacing = rowSpacing;
  curColSpacing = initialColSpacing;
  curRowSpacing = initialRowSpacing;
  meshLevel = 0;

#if SPEED_TEST
  for(int i=0; i < numPointCols*numPointRows; ++i)
    rgPlacePoints0[i]=false;
  for(int i=0; i < lenMA0; ++i){
    double s,t;
    int r,c;
    s = rgMA0X[i] / curColSpacing;
    t = rgMA0Y[i] / curRowSpacing;
    c = (int)s;
    if(c == (numPointCols-1))//if last pt, make it s=1.0 for c=numPointsCols-2
      --c;
    r = (int)t;
    if(r == (numPointRows-1))//if last pt, make it t=1.0 for r=numPointsRows-2
      --r;
    rgPlacePoints0[r*numPointCols+c]=true;
    if((r+1) < numPointRows){
      rgPlacePoints0[(r+1)*numPointCols+c]=true;
      if((c+1) < numPointCols){
	rgPlacePoints0[(r+1)*numPointCols+(c+1)]=true;
      }
    }
    if((c+1) < numPointCols){
      rgPlacePoints0[r*numPointCols+(c+1)]=true;
    }
  }  
#endif //SPEED_TEST

#if NEW_WARP
  for(int i=0; i < lenMA0; ++i){
    //this code is what used to happen in warpPoint.  now we
    //pre-calculate it and store it in arrays so it doesn't have to be
    //done every time warpPoint is called.  We do this both when we reset the
    //meshes and when we refine the meshes.  Those are the only times they need
    //to be computed because they stay the same for the entire refinement level
    double s,t;
    int r,c;
    int idx;
    s = rgMA0X[i] / curColSpacing;
    t = rgMA0Y[i] / curRowSpacing;
    c = (int)s;
    if(c == (numPointCols-1))//if last pt, make it s=1.0 for c=numPointsCols-2
      --c;
    r = (int)t;
    if(r == (numPointRows-1))//if last pt, make it t=1.0 for r=numPointsRows-2
      --r;
    idx = r*numPointCols+c;
    s -= c;
    t -= r;
    double quadW, quadH;
    quadW = curColSpacing;
    quadH = curRowSpacing;
    if(c==(numPointCols-2)){
      quadW = rgPoints0X[idx+1] - rgPoints0X[idx];
      if(quadW>0.)//changed from >=1. to > 0.(changed after dissertation)
	s = (rgMA0X[i]/*was x*/ - rgPoints0X[idx]) / quadW;
      else
	s = 0.;
    }
    if(r==(numPointRows-2)){
      quadH = rgPoints0Y[idx+numPointCols] - rgPoints0Y[idx];
      if(quadH>0.)//changed from >=1. to >0. (after dissertation)
	t = (rgMA0Y[i]/*was y*/ - rgPoints0Y[idx]) / quadH;
      else
	t = 0.;
    }
    if(s<0.)
      s=0.;
    if(t<0.)
      t=0.;
    if(s>1.)
      s=1.;
    if(t>1.)
      t=1.;
    //save it for later
    rgMA0r[i] = r;
    rgMA0c[i] = c;
    rgMA0s[i] = s;
    rgMA0t[i] = t;
    rgMA0quadW[i] = quadW;
    rgMA0quadH[i] = quadH;
  }
#endif

  //now save the DP-warped MedialAxis0 locations as starting point for
  //cost calculation of "distance moved"
#if SAVE_IMAGES
  DImage imgStartPoints;
  imgStartPoints.create(w1,h1,DImage::DImage_u8);
  imgStartPoints.fill(255);
#endif
  for(int i=0; i < lenMA0; ++i){
    double xp,yp;
     if(warpPoint(rgMA0X[i], rgMA0Y[i], &xp, &yp)){
       // if(warpPoint(rgMA0X[i], rgMA0Y[i], &xp, &yp)||true){
      rgDPMA0X[i] = xp;
      rgDPMA0Y[i] = yp;
    }
    else{
      fprintf(stderr, "DMorphInk::resetMeshes() DP warp point%d out of bounds! (%.lf,%.lf) w0=%d h0=%d w1=%d h1=%d\n",i,xp,yp,w0,h0,w1,h1);
	  exit(1);
    //   rgDPMA0X[i] = 0.;
    //   rgDPMA0Y[i] = 0.;
     }
#if SAVE_IMAGES
    if((xp>=0.) && (yp>=0.) && (xp<w1) &&(yp<h1)){
      imgStartPoints.drawPixel((int)xp,(int)yp,0);
    }
#endif
  }
#if SAVE_IMAGES
  imgStartPoints.save("/tmp/startpoints.pgm");
#endif

}


//TODO test
void DMorphInk::morphOneWay(	const DImage &srcFrom, 
							const DImage &srcTo, 
							int bandWidthDP, 
							double nonDiagonalCostDP, 
							int meshSpacingStatic,
				  			int numRefinementsStatic, 
				  			double meshDiv) {
	int numImprovesPerRefinement = 3;
	int meshSpacing = (int)(srcFrom.height() / meshDiv);
	if(meshSpacing < 4)
	  	meshSpacing = 4;
	if(-1 != meshSpacingStatic)
	 	meshSpacing = meshSpacingStatic;
	init(srcFrom, srcTo, false, meshSpacing, bandWidthDP,nonDiagonalCostDP);
	if(fOnlyDoCoarseAlignment){
	}
	else{
	  	int numRefinements = 0;
	  	while(meshSpacing > 16){
	    		meshSpacing /=2;
	    		++numRefinements;
	  	}
	  	
	  	if(-1 != numRefinementsStatic)
	    		numRefinements = numRefinementsStatic;    
	   
	   	//TODO Brian: Make this refined until no improvment found?
	  	for(int ref=0; ref <= numRefinements; ++ref){
	   		for(int imp=0; imp < numImprovesPerRefinement; ++imp){
				improveMorph();
	    		}
	    		saveCurrentMorph();
	    		if(ref < numRefinements){
				refineMeshes();
	    		}
	  	}
	}
}

void DMorphInk::saveCurrentMorph()
{
	//do stuff
	if (DEBUG_stopAtStep)
		raise(SIGINT);
}





///Returns the "difference" metric between two words.
/**Use this function to just get the cost of comparing two words to
   each other.  Don't need to call init() or anything because this
   function will do that.  However, it assumes that the images are
   cropped, size-normalized, deslanted, and baselines aligned, etc.
   This will do an automorph from 0 to 1 and then from 1 to 0 and add
   the two costs together to get the final cost. bandWidth is the
   width of the Dynamic Programming band (Rath and Manmatha used 15
   for the George Washington letters. This is the final version of the
   function used in the dissertation, which uses both a dist map from
   red and blue in each calculation.*/
double DMorphInk::getWordMorphCost(const DImage &src0,
				   const DImage &src1,
				   int bandWidthDP,
				   double nonDiagonalCostDP,
				   int meshSpacingStatic,
				   int numRefinementsStatic,
				   double meshDiv,
				   double lengthMismatchPenalty){
  double cost = 0.;
  double cost2 = 0.;

  double lenPen = 0.;
  double wLong = 0., wShort=0.;
  if(src0.width() > src1.width()){
    wLong = src0.width();
    wShort = src1.width();
  }
  else{
    wLong = src1.width();
    wShort = src0.width();
  }
  lenPen = lengthMismatchPenalty*(wLong-wShort)/wLong;
  

  //get cost to morph from src0 to src1
  morphOneWay(	src0, 
			src1, 
			bandWidthDP, 
			nonDiagonalCostDP, 
			meshSpacingStatic,
			numRefinementsStatic, 
			meshDiv);
  cost = getCost() + lenPen;
  warpCostDPfull = warpCostDP + warpCostDPv;
  warpCostDPhoriz = warpCostDP;

  if(fOnlyDoOneDirection)
    return cost;


  //get cost to morph from src1 to src0
  morphOneWay(	src1, 
			src0, 
			bandWidthDP, 
			nonDiagonalCostDP, 
			meshSpacingStatic,
			numRefinementsStatic, 
			meshDiv);
  cost2 = getCost() + lenPen;
  warpCostDPfull += warpCostDP + warpCostDPv;
  warpCostDPhoriz += warpCostDP;


  return cost + cost2;
  // return cost;
}


///same as getWordMorphCost except that it forces numImprovesPerRefinement to 1 and numRefinements to 0 for a fast pass at the data.
double DMorphInk::getWordMorphCostFast(const DImage &src0,
				       const DImage &src1,
				       int bandWidthDP,
				       double nonDiagonalCostDP,
				       int meshSpacingStatic,
				       int numRefinementsStatic,
				       double meshDiv,
				       double lengthMismatchPenalty){
  int meshSpacing;
  int numRefinements;
  int numImprovesPerRefinement = 3;
  double cost = 0.;
  double cost2 = 0.;

  double lenPen = 0.;
  double wLong = 0., wShort=0.;
  if(src0.width() > src1.width()){
    wLong = src0.width();
    wShort = src1.width();
  }
  else{
    wLong = src1.width();
    wShort = src0.width();
  }
  lenPen = lengthMismatchPenalty*(wLong-wShort)/wLong;
  

  //get cost to morph from src0 to src1
  
  meshSpacing = (int)(src0.height() / meshDiv);
  if(meshSpacing < 4)
    meshSpacing = 4;
  if(-1 != meshSpacingStatic)
    meshSpacing = meshSpacingStatic;
  init(src0, src1, false, meshSpacing, bandWidthDP,nonDiagonalCostDP);

  if(fOnlyDoCoarseAlignment){
  }
  else{
    numRefinements = 0;
    while(meshSpacing > 16){
      meshSpacing /=2;
      ++numRefinements;
    }
    if(-1 != numRefinementsStatic)
      numRefinements = numRefinementsStatic;
#if SPEED_TEST3
    numImprovesPerRefinement = 1;
    numRefinements = 0;
#endif
    for(int ref=0; ref <= numRefinements; ++ref){
      for(int imp=0; imp < numImprovesPerRefinement; ++imp){
#if SPEED_TEST2
	improveMorphFast();
#else
	improveMorph();
#endif
      }
      if(ref < numRefinements){
	refineMeshes();
      }
    }
  }
  cost = getCost() + lenPen;
  warpCostDPfull = warpCostDP + warpCostDPv;
  warpCostDPhoriz = warpCostDP;

  if(fOnlyDoOneDirection)
    return cost;


  //get cost to morph from src0 to src1
  meshSpacing = (int)(src1.height() / meshDiv);
  if(meshSpacing < 4)
    meshSpacing = 4;
  if(-1 != meshSpacingStatic)
    meshSpacing = meshSpacingStatic;

  init(src1, src0, false, meshSpacing, bandWidthDP,nonDiagonalCostDP);

  if(fOnlyDoCoarseAlignment){
  }
  else{
    numRefinements = 0;
    while(meshSpacing > 16){
      meshSpacing /=2;
      ++numRefinements;
    }
    if(-1 != numRefinementsStatic)
      numRefinements = numRefinementsStatic;
#if SPEED_TEST3
    numImprovesPerRefinement = 1;
    numRefinements = 0;
#endif
    for(int ref=0; ref <= numRefinements; ++ref){
      for(int imp=0; imp < numImprovesPerRefinement; ++imp){
#if SPEED_TEST2
	improveMorphFast();
#else
	improveMorph();
#endif
      }
      if(ref < numRefinements){
	refineMeshes();
      }
    }
  }
  cost2 = getCost() + lenPen;
  warpCostDPfull += warpCostDP + warpCostDPv;
  warpCostDPhoriz += warpCostDP;


  return cost + cost2;
  // return cost;
}



/**Calculates total cost as the avg distance MedialAxis0 pixels travel
(from their original DP-warped positions) plus the avg cost for MA1
pixels in distance map created from warped MA0.  */
double DMorphInk::getCost(){
  double costForward = 0.;
  double costBackwards = 0.;
  //int w, h;
  double costTotal = 0.;
  DImage imgWarped;
  DImage imgWarpedDist;
  D_uint8 *p8;
  signed int *ps32;
  int xpMin, xpMax, ypMin, ypMax;
  int wDM, hDM;//widht, height of the distance map created in this function
  
  xpMin = w1;
  xpMax = 0;
  ypMin = h1;
  ypMax = 0;
  ps32 = (signed int*)imgDist1.dataPointer_u32();
  for(int i=0; i < lenMA0; ++i){
    double xp,yp;
#if NEW_WARP
    warpPointNew(rgMA0s[i], rgMA0t[i], rgMA0r[i]*numPointCols+rgMA0c[i],
		 &xp, &yp/*, rgMA0X[i], rgMA0Y[i]*/);
    {//no longer using if. just put curly brace to start the block
#else
    if(warpPoint(rgMA0X[i], rgMA0Y[i], &xp, &yp)){
#endif
      int ixp, iyp;
      int addDistX, addDistY; // if the position is off the distmap, compensate
      ixp=(int)xp;
      iyp=(int)yp;
      if(ixp < xpMin)
	xpMin = ixp;
      if(iyp < ypMin)
	ypMin = iyp;
      if(ixp > xpMax)
	xpMax = ixp;
      if(iyp > ypMax)
	ypMax = iyp;
      addDistX = addDistY = 0;
      //added this to fix the position and cost if ixp,iyp out of distmap
#if 1 /* COMPENSATE_FOR_OOB_DISTMAP*/
      if(ixp < 0){
	addDistX = 0-ixp;
	ixp = 0;
      }
      else if(ixp >= w1){
	addDistX = ixp-w1+1;
	ixp = w1-1;
      }
      if(iyp < 0){
	addDistY = 0-iyp;
	iyp = 0;
      }
      else if(iyp >= h1){
	addDistY = iyp-h1+1;
	iyp = h1-1;
      }
#endif /*COMPENSATE_FOR_OOB_DISTMAP*/
      if((ixp>=0)&&(ixp<w1)&&(iyp>=0)&&(iyp<h1))
	costForward += (double)(ps32[w1*iyp+ixp]);
      else{
	fprintf(stderr,"THIS SHOULDN'T HAPPEN! (%s:%d)\n",__FILE__,__LINE__);
	abort();
      }
      costForward += addDistX + addDistY;
    }
    // else{
    //   printf("warpPoint() returned false! 2\n");
    // }
  }
  if(lenMA0 > 0)
    costForward /=(lenMA0);

  if(xpMin>0)
    xpMin = 0;
  if(ypMin > 0)
    ypMin = 0;
  if(xpMax < w1)
    xpMax = w1;
  if(ypMax < h1)
    ypMax = h1;
  wDM = xpMax - xpMin+1;
  hDM = ypMax - ypMin+1;
  // imgWarped.create(w1+40,h1+40,DImage::DImage_u8);
  // imgWarped.create(w1,h1,DImage::DImage_u8);
  imgWarped.create(wDM,hDM,DImage::DImage_u8);
  imgWarped.fill(255);
  p8 = imgWarped.dataPointer_u8();
  for(int i=0; i < lenMA0; ++i){
    double xp,yp;
#if NEW_WARP
    warpPointNew(rgMA0s[i], rgMA0t[i], rgMA0r[i]*numPointCols+rgMA0c[i],
		 &xp, &yp/*, rgMA0X[i], rgMA0Y[i]*/);
    {//no longer using if.  just start block with curly brace
#else
    if(warpPoint(rgMA0X[i], rgMA0Y[i], &xp, &yp)){
#endif
      int ixp, iyp;
      // ixp=(int)xp+20;
      // iyp=(int)yp+20;
      // if((ixp>=0)&&(ixp<(w1+40))&&(iyp>=0)&&(iyp<(h1+40)))
      // 	p8[iyp*(w1+40)+ixp]=0x00;
      // ixp=(int)xp;
      // iyp=(int)yp;
      // if((ixp>=0)&&(ixp<(w1))&&(iyp>=0)&&(iyp<(h1)))
	// p8[iyp*(w1)+ixp]=0x00;
      ixp=(int)xp;
      iyp=(int)yp;
      // if((ixp>=0)&&(ixp<wDM)&&(iyp>=0)&&(iyp<hDM))
      if((ixp>=xpMin)&&(ixp<=xpMax)&&(iyp>=ypMin)&&(iyp<=ypMax))
	p8[(iyp-ypMin)*wDM+ixp-xpMin]=0x00;
      else{
	fprintf(stderr,"getCost() rgMA0[%d] mapped to %d,%d wDM=%d hDM=%d\n",
		i,ixp,iyp,wDM,hDM);
	// fprintf(stderr,"getCost() rgMA0[%d] mapped to %d,%d w1=%d h1=%d\n",
	// 	i,ixp,iyp,w1,h1);
	//	abort();
      }

    }
    // else{
    //   printf("warpPoint() returned false! 2\n");
    // }
  }

  DDistanceMap::getDistFromInkBitonal_(imgWarpedDist, imgWarped,1000,-1000);

  ps32 = (signed int*)imgWarpedDist.dataPointer_u32();
  int numMA1pixels;
  numMA1pixels =0;
  for(int i=0; i < lenMA1; ++i){
    int ixp, iyp;
    // ixp = 20+(int)rgMA1X[i];
    // iyp = 20+(int)rgMA1Y[i];
    // ixp = (int)rgMA1X[i];
    // iyp = (int)rgMA1Y[i];
    // if((ixp>=0)&&(iyp>=0)&&(ixp<(w1+40))&&(iyp<(h1+40))){
    //   costBackwards += (double)ps32[iyp*(w1+40)+ixp];
    ixp=(int)rgMA1X[i];
    iyp=(int)rgMA1Y[i];
    // if((ixp>=0)&&(iyp>=0)&&(ixp<wDM)&&(iyp<hDM)){
    if((ixp>=xpMin)&&(ixp<=xpMax)&&(iyp>=ypMin)&&(iyp<=ypMax)){
      costBackwards += (double)ps32[(iyp-ypMin)*wDM+ixp-xpMin];
      ++numMA1pixels;
    }
    else{
      printf("getCost() out of bounds ixp,iyp=%d,%d wDM=%d hDM=%d\n",ixp,iyp,wDM,hDM);
    }
  }
  if(numMA1pixels > 0)
    costBackwards /= numMA1pixels;
  // if(lenMA1 > 0)
  //   costBackwards /=(lenMA1);//we also want this to be with respect to word 0

  costTotal = costForward + costBackwards;

  return costTotal;
}



/**Iterates through each point within a rectangular region around the
   current control point looking for lowest cost position to place the
   control point. */
void DMorphInk::improveMorph(){
  int radiusX, radiusY;

  //try vertex in each valid position within a square about one-third
  //the mesh spacing (positions are only valid if they don't cause
  //crossing mesh edges or messed up vertices and they are within
  //the image boundaries)


  radiusX = curColSpacing*.4;
  radiusY = curRowSpacing*.4;
  if(radiusX<1)
    radiusX = 1;
  if(radiusY<1)
    radiusY = 1;
  // printf("improving morph with radiusX=%d radiusY=%d\n",radiusX, radiusY);
  for(int r=0, idx=0; r < numPointRows; ++r){
    	for(int c=0; c < numPointCols; ++c, ++idx){
      	double curControlX, curControlY;//current x,y (in img1) of control point
      	double Vcost;//cost at curControlX,curControlY
      	double bestX, bestY;//best (lowest cost) x,y found so far
      	double bestCost;//best (lowest) cost
      	int minY, maxY, minX, maxX;//bounds of search area
      	curControlX = rgPoints1X[r*numPointCols+c];
      	curControlY = rgPoints1Y[r*numPointCols+c];
#if SPEED_TEST
      	if(!rgPlacePoints0[r*numPointCols+c]){
			//printf("skipping r=%d c=%d\n",r, c);
			continue;
      	}
#endif 	//SPEED_TEST

      	//decide which MedialAxis0 pixels are within the four quads of this ctl pt
      	tempMA0len = 0;
      	// for(int i=0; i < lenMA0; ++i){
      	for(int i=0; i < lenMA0; i+=SUBSAMPLE){
			int maRow;
			int maCol;
			maCol = (int)(rgMA0X[i] / curColSpacing);
			maRow = (int)(rgMA0Y[i] / curRowSpacing);
			//is medial axis point within the four quads this control pt affects?
			//if((maRow>=(r-1)) && (maRow<=r) && (maCol>=(c-1)) && (maCol<=c)){
			if(((maRow==(r-1)) || (maRow==r)) && ((maCol==(c-1)) || (maCol==c))){
#if NEW_WARP
	  			rgTempMA0idx[tempMA0len] = i;
#else//also need these if debugging
	  			rgTempMA0X[tempMA0len] = rgMA0X[i];
	  			rgTempMA0Y[tempMA0len] = rgMA0Y[i];
#endif
	  			++tempMA0len;
			}
      	}


      	Vcost = bestCost =
#if NEW_WARP
		getVertexPositionCostNew(r, c, curControlX, curControlY);
#else
		getVertexPositionCost(r, c, curControlX, curControlY);
#endif

      	bestX = curControlX;
      	bestY = curControlY;

      	minY = curControlY-radiusY;
      	// if(minY<0)
      	// 	minY=0;
      	maxY = curControlY+radiusY;
      	// if(maxY>=h1)
      	// 	maxY=h1-1;
      	minX = curControlX-radiusX;
      	// if(minX<0)
      	// 	minX=0;
      	maxX = curControlX+radiusX;
      	// if(maxX>=w1)
      	// 	maxX=w1-1;
      		//Constrain the search area with volleyball position rules. (top
      	//corners must stay above bottom corners, left corners left of
      	//right corners).  TODO:I might want to see about relaxing this
      	//so they just can't cross the line between the two constraining
      	//corners instead of either corner, but then I'd have to do the
      	//check within the tight loop.  Either way, I can still get
      	//concave quads, but I don't think quads that cross edges
      	if(r>0){
			if(minY < rgPoints1Y[(r-1)*numPointCols+c])
	  			minY = rgPoints1Y[(r-1)*numPointCols+c];
	  			
			if((c>0) && (minY < rgPoints1Y[(r-1)*numPointCols+c-1]))
	  			minY = rgPoints1Y[(r-1)*numPointCols+c-1];
	  			
			if((c<(numPointCols-1))&&(minY < rgPoints1Y[(r-1)*numPointCols+c+1]))
	  			minY = rgPoints1Y[(r-1)*numPointCols+c+1];
      	}
      	
      	if(r<(numPointRows-1)){
			if(maxY > rgPoints1Y[(r+1)*numPointCols+c])
	 	 		maxY = rgPoints1Y[(r+1)*numPointCols+c];
	 	 		
			if((c>0) && (maxY > rgPoints1Y[(r+1)*numPointCols+c-1]))
	  			maxY = rgPoints1Y[(r+1)*numPointCols+c-1];
	  			
			if((c<(numPointCols-1))&&(maxY > rgPoints1Y[(r+1)*numPointCols+c+1]))
	 	 		maxY = rgPoints1Y[(r+1)*numPointCols+c+1];
      	}
      	
      	if(c>0){
			if(minX < rgPoints1X[r*numPointCols+c-1])
	  			minX = rgPoints1X[r*numPointCols+c-1];
	  			
			if((r>0)&&(minX < rgPoints1X[(r-1)*numPointCols+c-1]))
	  			minX = rgPoints1X[(r-1)*numPointCols+c-1];
	  			
			if((r<numPointRows-1)&&(minX < rgPoints1X[(r+1)*numPointCols+c-1]))
	  			minX = rgPoints1X[(r+1)*numPointCols+c-1];
      	}
      	
      	if(c<(numPointCols-1)){
			if(maxX > rgPoints1X[r*numPointCols+c+1])
	  			maxX = rgPoints1X[r*numPointCols+c+1];
	  			
			if((r>0)&&(maxX > rgPoints1X[(r-1)*numPointCols+c+1]))
	    			maxX = rgPoints1X[(r-1)*numPointCols+c+1];
	    			
			if((r<(numPointRows-1))&&(maxX > rgPoints1X[(r+1)*numPointCols+c+1]))
	  			maxX = rgPoints1X[(r+1)*numPointCols+c+1];
      	}
      	// for(int ty=minY; ty <= maxY; ++ty){
      	// 	for(int tx=minX; tx <= maxX; ++tx){
      	for(int ty=minY; ty <= maxY; ++ty){
			for(int tx=minX; tx <= maxX; ++tx){
	 		 	double curCost;
#if NEW_WARP
	  			curCost = getVertexPositionCostNew(r, c, tx, ty);
#else
	  			curCost = getVertexPositionCost(r, c, tx, ty);
#endif

	  			curCost += sqrt((tx-curControlX)*(tx-curControlX)+(ty-curControlY)*(ty-curControlY))/100.;

	  			if(curCost < bestCost){
	    				bestCost = curCost;
	    				bestX = tx;
	    				bestY = ty;
	  			}
			}//for tx
      	}//for ty
      	rgPoints1X[r*numPointCols+c] = bestX;
      	rgPoints1Y[r*numPointCols+c] = bestY;
    	}
  }
  // free(rgTempMA0X);
  // free(rgTempMA0Y);
}

/**similar to improveMorph except that it doesn't check every
position's placement cost. it does coarse-to-fine. */
void DMorphInk::improveMorphFast(){
  int radiusX, radiusY;

  int searchSkip = 1;

  //try vertex in each valid position within a square about one-third
  //the mesh spacing (positions are only valid if they don't cause
  //crossing mesh edges or messed up vertices and they are within
  //the image boundaries)


  radiusX = curColSpacing*.4;
  radiusY = curRowSpacing*.4;
  if(radiusX<1)
    radiusX = 1;
  if(radiusY<1)
    radiusY = 1;

#if 0
  if(radiusX < radiusY)
    searchSkip = radiusX / 4;
  else
    searchSkip = radiusY / 4;
  if(searchSkip > 4)
    searchSkip = 4;
  if(searchSkip < 1)
    searchSkip = 1;
#endif
  // printf("%d/%d ", radiusX,radiusY);
  searchSkip = 4;
  if(searchSkip > radiusX)
    searchSkip = radiusX;
  if(searchSkip > radiusY)
    searchSkip = radiusY;

  // printf("improving morph with radiusX=%d radiusY=%d\n",radiusX, radiusY);
  for(int r=0, idx=0; r < numPointRows; ++r){
    for(int c=0; c < numPointCols; ++c, ++idx){
      double Vx, Vy;//current x,y (in img1) of control point
      double Vcost;//cost at Vx,Vy
      double bestX, bestY;//best (lowest cost) x,y found so far
      double bestCost;//best (lowest) cost
      int minY, maxY, minX, maxX;//bounds of search area
      Vx = rgPoints1X[r*numPointCols+c];
      Vy = rgPoints1Y[r*numPointCols+c];
#if SPEED_TEST
      if(!rgPlacePoints0[r*numPointCols+c]){
	//printf("skipping r=%d c=%d\n",r, c);
	continue;
      }
#endif //SPEED_TEST

      //decide which MedialAxis0 pixels are within the four quads of this ctl pt
      tempMA0len = 0;
      // for(int i=0; i < lenMA0; ++i){
      for(int i=0; i < lenMA0; i+=SUBSAMPLE){
	int maRow;
	int maCol;
	maCol = (int)(rgMA0X[i] / curColSpacing);
	maRow = (int)(rgMA0Y[i] / curRowSpacing);
	//is medial axis point within the four quads this control pt affects?
	//	if((maRow>=(r-1)) && (maRow<=r) && (maCol>=(c-1)) && (maCol<=c)){
	if(((maRow==(r-1)) || (maRow==r)) && ((maCol==(c-1)) || (maCol==c))){
#if NEW_WARP
	  rgTempMA0idx[tempMA0len] = i;
#else
	  rgTempMA0X[tempMA0len] = rgMA0X[i];
	  rgTempMA0Y[tempMA0len] = rgMA0Y[i];
#endif
	  ++tempMA0len;
	}
      }


      Vcost = bestCost =
#if NEW_WARP
	getVertexPositionCostNew(r, c, Vx, Vy);
#else
	getVertexPositionCost(r, c, Vx, Vy);
#endif

      bestX = Vx;
      bestY = Vy;

      minY = Vy-radiusY;
      // if(minY<0)
      // 	minY=0;
      maxY = Vy+radiusY;
      // if(maxY>=h1)
      // 	maxY=h1-1;
      minX = Vx-radiusX;
      // if(minX<0)
      // 	minX=0;
      maxX = Vx+radiusX;
      // if(maxX>=w1)
      // 	maxX=w1-1;
      //Constrain the search area with volleyball position rules. (top
      //corners must stay above bottom corners, left corners left of
      //right corners).  TODO:I might want to see about relaxing this
      //so they just can't cross the line between the two constraining
      //corners instead of either corner, but then I'd have to do the
      //check within the tight loop.  Either way, I can still get
      //concave quads, but I don't think quads that cross edges
      if(r>0){
	if(minY < rgPoints1Y[(r-1)*numPointCols+c])
	  minY = rgPoints1Y[(r-1)*numPointCols+c];
	if((c>0) && (minY < rgPoints1Y[(r-1)*numPointCols+c-1]))
	  minY = rgPoints1Y[(r-1)*numPointCols+c-1];
	if((c<(numPointCols-1))&&(minY < rgPoints1Y[(r-1)*numPointCols+c+1]))
	  minY = rgPoints1Y[(r-1)*numPointCols+c+1];
      }
      if(r<(numPointRows-1)){
	if(maxY > rgPoints1Y[(r+1)*numPointCols+c])
	  maxY = rgPoints1Y[(r+1)*numPointCols+c];
	if((c>0) && (maxY > rgPoints1Y[(r+1)*numPointCols+c-1]))
	  maxY = rgPoints1Y[(r+1)*numPointCols+c-1];
	if((c<(numPointCols-1))&&(maxY > rgPoints1Y[(r+1)*numPointCols+c+1]))
	  maxY = rgPoints1Y[(r+1)*numPointCols+c+1];
      }
      if(c>0){
	if(minX < rgPoints1X[r*numPointCols+c-1])
	  minX = rgPoints1X[r*numPointCols+c-1];
	if((r>0)&&(minX < rgPoints1X[(r-1)*numPointCols+c-1]))
	  minX = rgPoints1X[(r-1)*numPointCols+c-1];
	if((r<numPointRows-1)&&(minX < rgPoints1X[(r+1)*numPointCols+c-1]))
	  minX = rgPoints1X[(r+1)*numPointCols+c-1];
      }
      if(c<(numPointCols-1)){
	if(maxX > rgPoints1X[r*numPointCols+c+1])
	  maxX = rgPoints1X[r*numPointCols+c+1];
	if((r>0)&&(maxX > rgPoints1X[(r-1)*numPointCols+c+1]))
	    maxX = rgPoints1X[(r-1)*numPointCols+c+1];
	if((r<(numPointRows-1))&&(maxX > rgPoints1X[(r+1)*numPointCols+c+1]))
	  maxX = rgPoints1X[(r+1)*numPointCols+c+1];
      }
      // for(int ty=minY; ty <= maxY; ++ty){
      // 	for(int tx=minX; tx <= maxX; ++tx){
      int searchSkipTmp;
      searchSkipTmp = searchSkip;
      while(searchSkipTmp >= 1){
	for(int ty=minY; ty <= maxY; ty+=searchSkipTmp){
	  for(int tx=minX; tx <= maxX; tx+=searchSkipTmp){
	    double curCost;
#if NEW_WARP
	    curCost = getVertexPositionCostNew(r, c, tx, ty);
#else
	    curCost = getVertexPositionCost(r, c, tx, ty);
#endif
	    
	    curCost += sqrt((tx-Vx)*(tx-Vx)+(ty-Vy)*(ty-Vy))/100.;
	    
	    if(curCost < bestCost){
	      bestCost = curCost;
	      bestX = tx;
	      bestY = ty;
	    }
	  }//for tx
	}//for ty
	if((bestX-searchSkipTmp+1)>minX)
	  minX = bestX-searchSkipTmp+1;
	if((bestX+searchSkipTmp-1)<maxX)
	  maxX = bestX+searchSkipTmp-1;
	if(maxX<minX)
	  maxX=minX;
	if((bestY-searchSkipTmp+1)>minY)
	  minY = bestY-searchSkipTmp+1;
	if((bestY+searchSkipTmp-1)<maxY)
	  maxY = bestY+searchSkipTmp-1;
	if(maxY<minY)
	  maxY=minY;
	// if(searchSkipTmp==1)
	//   searchSkipTmp = 0;
	// else
	//   searchSkipTmp = 1;
	searchSkipTmp /= 2;
      }//while(searchSkipTmp >= 1)
      rgPoints1X[r*numPointCols+c] = bestX;
      rgPoints1Y[r*numPointCols+c] = bestY;
    }
  }
  // free(rgTempMA0X);
  // free(rgTempMA0Y);
}










//What does this do?

// this version has been modified to only split the last row/col when they are
// larger than curRowSpacing/curColSpacing and not split in half, but so all
// but the very last row/col remain the same size in mesh0
void DMorphInk::refineMeshes(){
  int numPointColsNew;
  int numPointRowsNew;
  double *rgPoints0Xnew;//x-coords of mesh0 quad vertices (row-major order)
  double *rgPoints0Ynew;//y-coords
  double *rgPoints1Xnew;//same, but for mesh1
  double *rgPoints1Ynew;//y-coords for mesh1
  // everything should stay the same except for numPointRows, numPointCols, and

  //  fprintf(stderr,"refineMeshes()\n");

  if((numPointCols < 2)||(numPointRows < 2)){
    fprintf(stderr, "DMorphInk::refineMeshes() called without a valid mesh!\n");
    exit(1);
  }

  bool fSplitLastCol = false;
  bool fSplitLastRow = false;
  if((rgPoints0X[numPointCols-1]-rgPoints0X[numPointCols-2])>(curColSpacing/2.))
    fSplitLastCol = true;
  if((rgPoints0Y[(numPointRows-1)*numPointCols]-
      rgPoints0Y[(numPointRows-2)*numPointCols])>(curRowSpacing/2.))
    fSplitLastRow = true;

  numPointColsNew = numPointCols + (numPointCols-2);
  if(fSplitLastCol)
    ++numPointColsNew;
  numPointRowsNew = numPointRows + (numPointRows-2);
  if(fSplitLastRow)
    ++numPointRowsNew;

  rgPoints0Xnew =
    (double*)malloc(sizeof(double)*numPointColsNew*numPointRowsNew);
  D_CHECKPTR(rgPoints0Xnew);
  rgPoints0Ynew =
    (double*)malloc(sizeof(double)*numPointColsNew*numPointRowsNew);
  D_CHECKPTR(rgPoints0Ynew);
  rgPoints1Xnew =
    (double*)malloc(sizeof(double)*numPointColsNew*numPointRowsNew);
  D_CHECKPTR(rgPoints1Xnew);
  rgPoints1Ynew =
    (double*)malloc(sizeof(double)*numPointColsNew*numPointRowsNew);
  D_CHECKPTR(rgPoints1Ynew);

  for(int r=0; r < numPointRowsNew; ++r){
    for(int c=0; c < numPointColsNew; ++c){
      rgPoints0Xnew[r*numPointColsNew+c] = rgPoints0Ynew[r*numPointColsNew+c] =
      rgPoints1Xnew[r*numPointColsNew+c] = rgPoints1Ynew[r*numPointColsNew+c] =
	-999.;
    }
  }
  for(int r=0; r < numPointRows; ++r){
    for(int c=0; c < numPointCols; ++c){
      bool fSplitRight, fSplitBelow;
      fSplitRight = ((c < (numPointCols-2))||
			  (fSplitLastCol&&(c==(numPointCols-2))));
      fSplitBelow = ((r < (numPointRows-2))||
			  (fSplitLastRow&&(r==(numPointRows-2))));
      int rNew,cNew;//row,column in the rgPointsNew arrays
      rNew = r*2;
      cNew = c*2;
      if((r==(numPointRows-1)) && !fSplitLastRow)
	--rNew;
      if((c==(numPointCols-1)) && !fSplitLastCol)
	--cNew;

      //Mesh0 --------------------------
      //this point
      rgPoints0Xnew[rNew*numPointColsNew+cNew] = rgPoints0X[r*numPointCols+c];
      rgPoints0Ynew[rNew*numPointColsNew+cNew] = rgPoints0Y[r*numPointCols+c];
      //midpoint to right
      // if(c < (numPointCols-1)){
      if(fSplitRight){
	rgPoints0Xnew[rNew*numPointColsNew+cNew+1] =
	  (rgPoints0X[r*numPointCols+c] + curColSpacing / 2.);
	//(rgPoints0X[r*numPointCols+c] + rgPoints0X[r*numPointCols+c+1]) / 2.;
	rgPoints0Ynew[rNew*numPointColsNew+cNew+1] =
	  rgPoints0Y[r*numPointCols+c];//all y's on row are same
	//(rgPoints0Y[r*numPointCols+c] + rgPoints0Y[r*numPointCols+c+1]) / 2.;
      }
      //midpoint below
      if(fSplitBelow){
	rgPoints0Xnew[(rNew+1)*numPointColsNew+cNew] =
	  rgPoints0X[r*numPointCols+c]; //all x's in this column are the same
	//(rgPoints0X[r*numPointCols+c]+rgPoints0X[(r+1)*numPointCols+c]) / 2.;
	rgPoints0Ynew[(rNew+1)*numPointColsNew+cNew] =
	  (rgPoints0Y[r*numPointCols+c] + curRowSpacing / 2.);
	//(rgPoints0Y[r*numPointCols+c]+ rgPoints0Y[(r+1)*numPointCols+c]) / 2.;
      }
      //midpoint below and right (don't need bilinear interpolation for mesh 0
      //because it is regular and rectangular)
      if(fSplitRight && fSplitBelow){
	rgPoints0Xnew[(rNew+1)*numPointColsNew+cNew+1] = 
	  rgPoints0Xnew[(rNew)*numPointColsNew+cNew+1];//same as x above it
	rgPoints0Ynew[(rNew+1)*numPointColsNew+cNew+1] =
	  rgPoints0Ynew[(rNew+1)*numPointColsNew+cNew];//same as y left of it
      }
    }
  }
#if 1
  for(int r=0; r < numPointRows; ++r){
    for(int c=0; c < numPointCols; ++c){
      bool fSplitRight, fSplitBelow;
      double sSplit, tSplit; // the s,t coordinate of the point added
      fSplitRight = ((c < (numPointCols-2))||
			  (fSplitLastCol&&(c==(numPointCols-2))));
      fSplitBelow = ((r < (numPointRows-2))||
			  (fSplitLastRow&&(r==(numPointRows-2))));
      int rNew,cNew;//row,column in the rgPointsNew arrays
      rNew = r*2;
      cNew = c*2;
      if((r==(numPointRows-1)) && !fSplitLastRow)
	--rNew;
      if((c==(numPointCols-1)) && !fSplitLastCol)
	--cNew;
      sSplit = 0.5;
      if(c==(numPointCols-2))
	sSplit = (rgPoints0Xnew[rNew*numPointColsNew+cNew+1] -
		  rgPoints0X[r*numPointCols+c]) /
	  (rgPoints0X[r*numPointCols+c+1] - rgPoints0X[r*numPointCols+c]);
      tSplit = 0.5;
      if(r==(numPointRows-2))
	tSplit = (rgPoints0Ynew[(rNew+1)*numPointColsNew+cNew] - 
		  rgPoints0Y[r*numPointCols+c]) / 
	  (rgPoints0Y[(r+1)*numPointCols+c] - rgPoints0Y[r*numPointCols+c]);

      //Mesh1 --------------------------
      //this point
      rgPoints1Xnew[rNew*numPointColsNew+cNew] = rgPoints1X[r*numPointCols+c];
      rgPoints1Ynew[rNew*numPointColsNew+cNew] = rgPoints1Y[r*numPointCols+c];
      //midpoint to right
      if(fSplitRight){
	rgPoints1Xnew[rNew*numPointColsNew+cNew+1] =
	  (1.-sSplit)*rgPoints1X[r*numPointCols+c] +
	  sSplit*rgPoints1X[r*numPointCols+c+1];
	rgPoints1Ynew[rNew*numPointColsNew+cNew+1] =
	  (1.-sSplit)*rgPoints1Y[r*numPointCols+c] +
	  sSplit*rgPoints1Y[r*numPointCols+c+1];
      }
      //midpoint below
      if(fSplitBelow){
	rgPoints1Xnew[(rNew+1)*numPointColsNew+cNew] =
	  (1.-tSplit)*rgPoints1X[r*numPointCols+c] +
	  tSplit*rgPoints1X[(r+1)*numPointCols+c];
	rgPoints1Ynew[(rNew+1)*numPointColsNew+cNew] =
	  (1.-tSplit)*rgPoints1Y[r*numPointCols+c] +
	  tSplit*rgPoints1Y[(r+1)*numPointCols+c];
      }
      //midpoint below and right (use bilinear interpolation)
      if((c <(numPointCols-1)) &&(r < (numPointRows-1))){
	double Xtop, Ytop, Xbot, Ybot;
	Xtop = rgPoints1Xnew[rNew*numPointColsNew+cNew+1];
	Ytop = rgPoints1Ynew[rNew*numPointColsNew+cNew+1];
	Xbot = (1.-sSplit)*rgPoints1X[(r+1)*numPointCols+c] +
	  sSplit*rgPoints1X[(r+1)*numPointCols+c+1];
	Ybot = (1.-sSplit)*rgPoints1Y[(r+1)*numPointCols+c] +
	  sSplit*rgPoints1Y[(r+1)*numPointCols+c+1];
	rgPoints1Xnew[(rNew+1)*numPointColsNew+cNew+1] =
	  (1.-tSplit)*Xtop + tSplit*Xbot;
	rgPoints1Ynew[(rNew+1)*numPointColsNew+cNew+1] =
	  (1.-tSplit)*Ytop + tSplit*Ybot;
      }
    }
  }
#endif
  bool fBad = false;
  for(int r=0; r < numPointRowsNew; ++r){
    for(int c=0; c < numPointColsNew; ++c){
      if((-999. == rgPoints0Xnew[r*numPointColsNew+c]) ||
	 (-999. == rgPoints0Ynew[r*numPointColsNew+c]) ||
	 (-999. == rgPoints1Xnew[r*numPointColsNew+c]) ||
	 (-999. == rgPoints1Ynew[r*numPointColsNew+c])){
	printf("refined point r=%d,c=%d was not assigned a value! x0=%lf,y0=%lf x1=%lf,y1=%lf\n",
	       r,c,rgPoints0Xnew[r*numPointColsNew+c],rgPoints0Ynew[r*numPointColsNew+c],rgPoints1Xnew[r*numPointColsNew+c],rgPoints1Ynew[r*numPointColsNew+c]);
	printf("numPointCols=%d numPointColsNew=%d\nnumPointRows=%d numPointRowsNew=%d\n",numPointCols,numPointColsNew,numPointRows,numPointRowsNew);
	fBad = true;
      }
    }
  }
  if(fBad){
    FILE *fout;
    fout = fopen("/tmp/meshes.txt","wb");
    fprintf(fout,"numPointCols=%d numPointColsNew=%d\nnumPointRows=%d numPointRowsNew=%d\n",numPointCols,numPointColsNew,numPointRows,numPointRowsNew);
    fprintf(fout,"curRowSpacing=%lf curColSpacing=%lf\n",curRowSpacing,curColSpacing);
    fprintf(fout,"last row height=%lf-%lf=%lf\n",rgPoints0Y[(numPointRows-1)*numPointCols],rgPoints0Y[(numPointCols-2)*numPointCols],(rgPoints0Y[(numPointRows-1)*numPointCols]- rgPoints0Y[(numPointCols-2)*numPointCols]));
    fprintf(fout,"fSplitLastRow=%d fSplitLastCol=%d\n",(int)fSplitLastRow, (int)fSplitLastCol);
    fprintf(fout,"BEFORE REFINE: mesh0\n");
    for(int r=0; r < numPointRows; ++r){
      for(int c=0; c < numPointCols; ++c){
	bool fSplitRight, fSplitBelow;
	fSplitRight = ((c < (numPointCols-2))||
		       (fSplitLastCol&&(c==(numPointCols-2))));
	fSplitBelow = ((r < (numPointRows-2))||
		       (fSplitLastRow&&(r==(numPointRows-2))));
	fprintf(fout,"  %c%c(%6.1lf,%6.1lf)",fSplitRight?'R':'r',fSplitBelow?'B':'b',rgPoints0X[r*numPointCols+c],
	       rgPoints0Y[r*numPointCols+c]);
      }
      fprintf(fout,"\n");
    }

    fprintf(fout,"BEFORE REFINE: mesh1\n");
    for(int r=0; r < numPointRows; ++r){
      for(int c=0; c < numPointCols; ++c){
	bool fSplitRight, fSplitBelow;
	fSplitRight = ((c < (numPointCols-2))||
		       (fSplitLastCol&&(c==(numPointCols-2))));
	fSplitBelow = ((r < (numPointRows-2))||
		       (fSplitLastRow&&(r==(numPointRows-2))));
	fprintf(fout,"  %c%c(%6.1lf,%6.1lf)",fSplitRight?'R':'r',fSplitBelow?'B':'b',rgPoints1X[r*numPointCols+c],
	       rgPoints1Y[r*numPointCols+c]);
      }
      fprintf(fout,"\n");
    }

    fprintf(fout,"mesh0:\n");
    for(int r=0; r < numPointRowsNew; ++r){
      for(int c=0; c < numPointColsNew; ++c){
	fprintf(fout,"  (%6.1lf,%6.1lf)",rgPoints0Xnew[r*numPointColsNew+c],
		rgPoints0Ynew[r*numPointColsNew+c]);
      }
      fprintf(fout,"\n");
    }
    fprintf(fout,"\nmesh1:\n");
    for(int r=0; r < numPointRowsNew; ++r){
      for(int c=0; c < numPointColsNew; ++c){
	fprintf(fout,"  (%6.1lf,%6.1lf)",rgPoints1Xnew[r*numPointColsNew+c],
		rgPoints1Ynew[r*numPointColsNew+c]);
      }
      fprintf(fout,"\n");
    }
    fclose(fout);
    exit(1);
  }

  numPointCols = numPointColsNew;
  numPointRows = numPointRowsNew;

  free(rgPoints0X);
  free(rgPoints0Y);
  free(rgPoints1X);
  free(rgPoints1Y);

  rgPoints0X = rgPoints0Xnew;
  rgPoints0Y = rgPoints0Ynew;
  rgPoints1X = rgPoints1Xnew;
  rgPoints1Y = rgPoints1Ynew;
  curColSpacing /= 2.;
  curRowSpacing /= 2.;
  ++meshLevel;


  free(rgPlacePoints0);
  rgPlacePoints0 = (bool*)malloc(sizeof(bool)*numPointCols*numPointRows);
  D_CHECKPTR(rgPlacePoints0);
#if SPEED_TEST
  for(int i=0; i < numPointCols*numPointRows; ++i)
    rgPlacePoints0[i]=false;
  for(int i=0; i < lenMA0; ++i){
    double s,t;
    int r,c;
    s = rgMA0X[i] / curColSpacing;
    t = rgMA0Y[i] / curRowSpacing;
    c = (int)s;
    if(c == (numPointCols-1))//if last pt, make it s=1.0 for c=numPointsCols-2
      --c;
    r = (int)t;
    if(r == (numPointRows-1))//if last pt, make it t=1.0 for r=numPointsRows-2
      --r;
    rgPlacePoints0[r*numPointCols+c]=true;
    if((r+1) < numPointRows){
      rgPlacePoints0[(r+1)*numPointCols+c]=true;
      if((c+1) < numPointCols){
	rgPlacePoints0[(r+1)*numPointCols+(c+1)]=true;
      }
    }
    if((c+1) < numPointCols){
      rgPlacePoints0[r*numPointCols+(c+1)]=true;
    }
  }
#endif//SPEED_TEST



#if NEW_WARP
  for(int i=0; i < lenMA0; ++i){
    //this code is what used to happen in warpPoint.  now we
    //pre-calculate it and store it in arrays so it doesn't have to be
    //done every time warpPoint is called.  We do this both when we reset the
    //meshes and when we refine the meshes.  Those are the only times they need
    //to be computed because they stay the same for the entire refinement level
    double s,t;
    int r,c;
    int idx;
    s = rgMA0X[i] / curColSpacing;
    t = rgMA0Y[i] / curRowSpacing;
    c = (int)s;
    if(c == (numPointCols-1))//if last pt, make it s=1.0 for c=numPointsCols-2
      --c;
    r = (int)t;
    if(r == (numPointRows-1))//if last pt, make it t=1.0 for r=numPointsRows-2
      --r;
    idx = r*numPointCols+c;
    s -= c;
    t -= r;
    double quadW, quadH;
    quadW = curColSpacing;
    quadH = curRowSpacing;
    if(c==(numPointCols-2)){
      quadW = rgPoints0X[idx+1] - rgPoints0X[idx];
      if(quadW>0.)//changed from >=1. to > 0.(changed after dissertation)
	s = (rgMA0X[i]/*was x*/ - rgPoints0X[idx]) / quadW;
      else
	s = 0.;
    }
    if(r==(numPointRows-2)){
      quadH = rgPoints0Y[idx+numPointCols] - rgPoints0Y[idx];
      if(quadH>0.)//changed from >=1. to >0. (after dissertation)
	t = (rgMA0Y[i]/*was y*/ - rgPoints0Y[idx]) / quadH;
      else
	t = 0.;
    }
    if(s<0.)
      s=0.;
    if(t<0.)
      t=0.;
    if(s>1.)
      s=1.;
    if(t>1.)
      t=1.;
    //save it for later
    rgMA0r[i] = r;
    rgMA0c[i] = c;
    rgMA0s[i] = s;
    rgMA0t[i] = t;
    rgMA0quadW[i] = quadW;
    rgMA0quadH[i] = quadH;
  }
#endif

  //assign the Medial Axis points their corresponding quad indexes
  //assignMA0qidxs();
}






//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//   The functions below are for visualization / debug only
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------





#if 0
//this estimate may no longer be valid (I didn't double-check, so I just
//commented it out.
//this is a cost estimate (entire word) using D1 instead of D0' (for video)
/**It is not exactly the same calculation used in the morphing because
 * here I divide by the number of pixels (lenMA0) instead of smoothing
 * with (lenMA0+1) and because in the actual morphing algorithm we do
 * only the MA pixels within the four quads adjacent to each single
 * control point instead of using the entire thing.
 */
// double DMorphInk::getCost_fast(){
//   double costDist = 0.;
//   signed int *ps32;


//   ps32 = (signed int*)imgDist1.dataPointer_u32();
//   for(int i=0; i < lenMA0; ++i){
//     double xp,yp;
//     if(warpPoint(rgMA0X[i], rgMA0Y[i], &xp, &yp)){
//       int ixp, iyp;
//       ixp=(int)xp;
//       iyp=(int)yp;
//       if((ixp>=0)&&(ixp<w1)&&(iyp>=0)&&(iyp<h1))
// 	costDist += (double)(ps32[w1*iyp+ixp]);
//     }
//     else{
//       printf("getCost_fast() warpPoint==false\n");
//     }
//   }
//   costDist /=(lenMA0+1);//we want this to be with respect to word 0
//   return costDist;
// }
#endif





/// Used for DEBUG/visualization only.
/**  Returns the morph at some point between 0. and 1. (inclusive).
    If fWarpOnly is true, there will be no fade/blend, just warping of
    the pixels from src0 to the mesh of src1.  ***This is a simplified
    function*** It assumes that the mesh at t=0 is rectangular and
    does not do backward mapping for pixel values and instead does
    super-sampled forward mapping.*/
DImage DMorphInk::getMorphAtTime(double dblTime, bool fWarpOnly,
				 bool fDrawMesh, double pad_val,
				 bool fOverlayImg1, bool fShowMedialAxis0,
				 bool fShowMedialAxis1, bool fNoMorph){
  int w, h;
  double *rgXs, *rgYs;//the warped mesh at time=dblTime
  DImage imgMorphed;
  D_uint8 *pSrc0, *pSrc1;
  D_uint8 *pDst;

  DImage imgQidx;

  if((dblTime < 0.) || (dblTime > 1.)){
    fprintf(stderr,
  	    "DMorphInk::getMorphAtTime() dblTime must be >= 0. and <= 1.\n");
    exit(1);
  }
  w = (w0 >= w1) ? w0 : w1;
  h = (h0 >= h1) ? h0 : h1;
  if(pimg0->getImageType() != pimg1->getImageType()){
    fprintf(stderr,
	    "DMorphInk::getMorphAtTime() types of the two images being morphed must be the same\n");
    return imgMorphed;
  }
  if((pimg0->getImageType() != DImage::DImage_u8) && 
     (pimg0->getImageType() != DImage::DImage_RGB)){
    fprintf(stderr,
	    "DMorphInk::getMorphAtTime() only supports 8-bit RGB or GS images\n");
    return imgMorphed;
  }
  if(pimg0->getImageType() != DImage::DImage_u8){
    fprintf(stderr,
	    "DMorphInk::getMorphAtTime() only supports 8-bit GS images so far\n");
    exit(1);
  }


  imgMorphed.create(w,h,pimg0->getImageType());
  pDst = imgMorphed.dataPointer_u8();
  pSrc0 = pimg0->dataPointer_u8();
  pSrc1 = pimg1->dataPointer_u8();

  imgQidx.create(w,h,DImage::DImage_u8);
  imgQidx.fill(255);

  // define rgXs and rgYs as the intermediate mesh point locations at t=dblTime
  rgXs = (double*)malloc(sizeof(double)*numPointRows*numPointCols);
  D_CHECKPTR(rgXs);
  rgYs = (double*)malloc(sizeof(double)*numPointRows*numPointCols);
  D_CHECKPTR(rgYs);
  for(int r=0, qidx=0; r < numPointRows; ++r){
    for(int c=0; c < numPointCols; ++c,++qidx){
      rgXs[qidx] = (1.0-dblTime)*rgPoints0X[qidx]+dblTime*rgPoints1X[qidx];
      rgYs[qidx] = (1.0-dblTime)*rgPoints0Y[qidx]+dblTime*rgPoints1Y[qidx];
    }
  }

  // do forward mapping (This is the main warp)
  if(!fNoMorph){
    for(int r=0, qidx=0; r < (numPointRows-1); ++r){
      for(int c=0; c < (numPointCols-1); ++c, ++qidx){
	int idxTL; //point index of top-left vertex of quad with index qidx
	idxTL = qidx + (qidx / (numPointCols-1));
	
	for(double t=0.; t < 1.; t+=0.001){
	  for(double s=0.; s < 1.; s+=0.001){
	    double Xtop, Ytop, Xbot, Ybot, Px, Py;
	    D_uint8 val;
	    Xtop =(1.-s)*rgXs[idxTL] + s*rgXs[idxTL+1];
	    Ytop =(1.-s)*rgYs[idxTL] + s*rgYs[idxTL+1];
	    Xbot =(1.-s)*rgXs[idxTL+numPointCols]+s*rgXs[idxTL+numPointCols+1];
	    Ybot =(1.-s)*rgYs[idxTL+numPointCols]+s*rgYs[idxTL+numPointCols+1];
	    Px = (1.-t)*Xtop + t*Xbot;
	    Py = (1.-t)*Ytop + t*Ybot;
	    int Sx, Sy; //x,y in source image (assume rectangular quad in mesh)
	    Sx = (1.-s)*rgPoints0X[idxTL]+s*rgPoints0X[idxTL+1];
	    Sy = (1.-t)*rgPoints0Y[idxTL]+t*rgPoints0Y[idxTL+numPointCols];
	    val = pSrc0[Sy*w0+Sx];
	    if((Py < h)&&(Px < w))
	      pDst[((int)Py)*w+(int)Px] = val;
	  }
	}
      }
    }
  }
  else{
    imgMorphed.fill(255);
  }
  if(fOverlayImg1){
    if(imgMorphed.getImageType() == DImage::DImage_u8)
      imgMorphed = imgMorphed.convertedImgType(DImage::DImage_RGB);
    D_uint8 *p8src;
    p8src = pimg1->dataPointer_u8();
    for(int y=0, idx=0; y < h1; ++y){
      for(int x=0; x < w1; ++x, ++idx){
	if(p8src[idx] < 255){
	  imgMorphed.drawPixel(x,y,(int)(p8src[idx]),0,0,0.5);
	}
      }
    }
  }
  if(fDrawMesh){
    double meshTrans;
    meshTrans = 0.2;
    if(fFaintMesh)
      meshTrans = 0.7;
    if(imgMorphed.getImageType() == DImage::DImage_u8)
      imgMorphed = imgMorphed.convertedImgType(DImage::DImage_RGB);
    for(int r=0, meshIdx=0; r < numPointRows; ++r){
      for(int c=0; c < numPointCols; ++c, ++meshIdx){
	if(c>0)
	  imgMorphed.drawLine(rgXs[r*numPointCols+c],rgYs[r*numPointCols+c],
			      rgXs[r*numPointCols+c-1],rgYs[r*numPointCols+c-1],
			      // 255,0,0,0.2);
			      255,255,0,meshTrans);
	if(c<(numPointCols-1))
	  imgMorphed.drawLine(rgXs[r*numPointCols+c],rgYs[r*numPointCols+c],
			      rgXs[r*numPointCols+c+1],rgYs[r*numPointCols+c+1],
			      // 0,255,0,0.2);
			      255,255,0,meshTrans);
	if(r>0)
	  imgMorphed.drawLine(rgXs[r*numPointCols+c],rgYs[r*numPointCols+c],
			      rgXs[(r-1)*numPointCols+c],
			      rgYs[(r-1)*numPointCols+c],
			      // 0,0,255,0.2);
			      255,255,0,meshTrans);
	if(r<(numPointRows-1))
	  imgMorphed.drawLine(rgXs[r*numPointCols+c],rgYs[r*numPointCols+c],
			      rgXs[(r+1)*numPointCols+c],
			      rgYs[(r+1)*numPointCols+c],
			      // 255,0,255,0.2);
			      255,255,0,meshTrans);
      }
    }
  }
  if(fShowMedialAxis0){
    if(imgMorphed.getImageType() == DImage::DImage_u8)
      imgMorphed = imgMorphed.convertedImgType(DImage::DImage_RGB);
    // for(int ma_pt=0; ma_pt < lenMA0; ++ma_pt){
    for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
      double xp, yp;
      if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,dblTime))
	imgMorphed.drawPixel((int)xp,(int)yp,0,255,0,0.);
    }
  }
  if(fShowMedialAxis1){
    if(imgMorphed.getImageType() == DImage::DImage_u8)
      imgMorphed = imgMorphed.convertedImgType(DImage::DImage_RGB);
    // for(int ma_pt=0; ma_pt < lenMA1; ++ma_pt){
    for(int ma_pt=0; ma_pt < lenMA1; ma_pt+=SUBSAMPLE){
      imgMorphed.drawPixel((int)rgMA1X[ma_pt],(int)rgMA1Y[ma_pt],255,255,0,0.);
    }
  }



  free(rgXs);
  free(rgYs);

  imgQidx.save("/tmp/qidx.pgm");

  return imgMorphed;
}

DImage DMorphInk::getImageWithMesh(int meshNum){
  DImage imgRet;
  double *rgXs, *rgYs;

  if(0 == meshNum){
    imgRet = (*pimg0);
    rgXs = rgPoints0X;
    rgYs = rgPoints0Y;
  }
  else{
    imgRet = (*pimg1);
    rgXs = rgPoints1X;
    rgYs = rgPoints1Y;
  }


  if(imgRet.getImageType() == DImage::DImage_u8)
    imgRet = imgRet.convertedImgType(DImage::DImage_RGB);
  for(int r=0, meshIdx=0; r < numPointRows; ++r){
    for(int c=0; c < numPointCols; ++c, ++meshIdx){
      if(c>0)
	imgRet.drawLine(rgXs[r*numPointCols+c],rgYs[r*numPointCols+c],
			rgXs[r*numPointCols+c-1],rgYs[r*numPointCols+c-1],
			255,0,0,0.2);
      if(c<(numPointCols-1))
	imgRet.drawLine(rgXs[r*numPointCols+c],rgYs[r*numPointCols+c],
			rgXs[r*numPointCols+c+1],rgYs[r*numPointCols+c+1],
			0,255,0,0.2);
      if(r>0)
	imgRet.drawLine(rgXs[r*numPointCols+c],rgYs[r*numPointCols+c],
			rgXs[(r-1)*numPointCols+c],
			rgYs[(r-1)*numPointCols+c],
			0,0,255,0.2);
      if(r<(numPointRows-1))
	imgRet.drawLine(rgXs[r*numPointCols+c],rgYs[r*numPointCols+c],
			rgXs[(r+1)*numPointCols+c],
			rgYs[(r+1)*numPointCols+c],
			255,0,255,0.2);
    }
  }

  return imgRet;
}






DImage DMorphInk::getMedialAxisWithMesh(int meshNum){
  DImage imgRet;
  double *rgXs, *rgYs;

  if(0 == meshNum){
    imgRet.create(w0,h0,DImage::DImage_RGB);
    imgRet.fill(255,255,255);
    rgXs = rgPoints0X;
    rgYs = rgPoints0Y;
  }
  else{
    imgRet.create(w1,h1,DImage::DImage_RGB);
    imgRet.fill(255,255,255);
    rgXs = rgPoints1X;
    rgYs = rgPoints1Y;
  }


  for(int r=0, meshIdx=0; r < numPointRows; ++r){
    for(int c=0; c < numPointCols; ++c, ++meshIdx){
      if(c>0)
	imgRet.drawLine(rgXs[r*numPointCols+c],rgYs[r*numPointCols+c],
			rgXs[r*numPointCols+c-1],rgYs[r*numPointCols+c-1],
			244,184,0,0.);
      if(c<(numPointCols-1))
	imgRet.drawLine(rgXs[r*numPointCols+c],rgYs[r*numPointCols+c],
			rgXs[r*numPointCols+c+1],rgYs[r*numPointCols+c+1],
			244,184,0,0.);
      if(r>0)
	imgRet.drawLine(rgXs[r*numPointCols+c],rgYs[r*numPointCols+c],
			rgXs[(r-1)*numPointCols+c],
			rgYs[(r-1)*numPointCols+c],
			244,184,0,0.);
      if(r<(numPointRows-1))
	imgRet.drawLine(rgXs[r*numPointCols+c],rgYs[r*numPointCols+c],
			rgXs[(r+1)*numPointCols+c],
			rgYs[(r+1)*numPointCols+c],
			244,184,0,0.2);
    }
  }

  if(0 == meshNum){
    for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
      //      double xp, yp;
      // if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,1.))
      // 	imgRet.drawPixel((int)xp,(int)yp,0,0,0,0.);
      imgRet.drawPixel((int)rgMA0X[ma_pt],(int)rgMA0Y[ma_pt],255,0,0,0.);
    }
  }
  else{
    for(int ma_pt=0; ma_pt < lenMA1; ma_pt+=SUBSAMPLE){
      imgRet.drawPixel((int)rgMA1X[ma_pt],(int)rgMA1Y[ma_pt],0,0,255,0.);
    }
  }

  return imgRet;
}


DImage DMorphInk::getWarpedMedialAxis0_A0prime(){
  DImage imgRet;
  double *rgXs, *rgYs;
  int w, h;
  w = w0;
  h = h0;
  //  if(pimg1->width() > w)
  if(w < w1)
    w = w1;
  if(h < h1)
    h = h1;
  
  imgRet.create(w,h,DImage::DImage_RGB);
  imgRet.fill(255,255,255);


  rgXs = rgPoints0X;
  rgYs = rgPoints0Y;
  for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
    double xp, yp;
    if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,1.))
      imgRet.drawPixel((int)xp,(int)yp,255,0,0,0.);
  }
  return imgRet;
}

DImage DMorphInk::getBothMedialAxes(bool fWarpPoints0, bool fShowGreenDistance){
  DImage imgRet;
  double *rgXs, *rgYs;
  int w, h;
  w = w0;
  h = h0;
  //  if(pimg1->width() > w)
  if(w < w1)
    w = w1;
  if(h < h1)
    h = h1;
  
  imgRet.create(w,h,DImage::DImage_RGB);
  imgRet.fill(255,255,255);

  if(fShowGreenDistance){
    int qlen;
    int *rgQx;
    int *rgQy;
    DImage imgAxis0, imgAxis1;
    DImage imgD0, imgD1;
    DImage imgGreen;
    D_uint8 *p8Green;
    imgAxis0.create(w,h,DImage::DImage_u8);
    imgAxis0.fill(255.);
    imgAxis1.create(w,h,DImage::DImage_u8);
    imgAxis1.fill(255.);
    imgGreen.create(w,h,DImage::DImage_RGB);
    imgGreen.fill(255.,255.,255.);
    p8Green = imgGreen.dataPointer_u8();
    qlen = 0;
    rgQx = new int[w*h];
    D_CHECKPTR(rgQx);
    rgQy = new int[w*h];
    D_CHECKPTR(rgQy);
    for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
      double xp, yp;
      if(fWarpPoints0){
	if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,1.))
	  imgAxis0.drawPixel((int)xp,(int)yp,0);
      }
      else{
	imgAxis0.drawPixel((int)rgMA0X[ma_pt],(int)rgMA0Y[ma_pt],0);
      }
    }
    DDistanceMap::getDistFromInkBitonal_(imgD0, imgAxis0);
    for(int ma_pt=0; ma_pt < lenMA1; ma_pt+=SUBSAMPLE){
      int px, py;
      px = (int)rgMA1X[ma_pt];
      py = (int)rgMA1Y[ma_pt];
      imgAxis1.drawPixel((int)rgMA1X[ma_pt],(int)rgMA1Y[ma_pt],0);
    }
    DDistanceMap::getDistFromInkBitonal_(imgD1, imgAxis1);

    for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
      double xp, yp;
      int ixp, iyp;
      ixp = -1; iyp = -1;
      if(fWarpPoints0){
	if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,1.)){
	  ixp = (int)xp;
	  iyp =(int)yp;
	}
      }
      else{
	ixp = (int)rgMA0X[ma_pt];
	iyp = (int)rgMA0Y[ma_pt];
      }
      if((ixp>=0) && (iyp>=0)){
	D_sint32 curDist;
	D_sint32 *ps32;
	ps32 = (D_sint32*)imgD1.dataPointer_u32();
	curDist=ps32[iyp*w+ixp];
	if((curDist >= 0) && (ixp>=0) && (ixp < w) && (iyp>=0) && (iyp < h) &&
	   (p8Green[3*(iyp*w+ixp)]==255)&&(p8Green[3*(iyp*w+ixp)+1]==255)&&
	   (p8Green[3*(iyp*w+ixp)+2]==255)){//white
	  imgGreen.drawPixel(ixp,iyp,255,255,0,0.);
	  rgQx[qlen] = ixp;
	  rgQy[qlen] = iyp;
	  ++qlen;
	}
      }
    }
    while(qlen > 0){
      int ixp, iyp;
      D_sint32 curDist;
      D_sint32 *ps32;
      ps32 = (D_sint32*)imgD1.dataPointer_u32();
      ixp = rgQx[qlen-1];
      iyp = rgQy[qlen-1];
      curDist=ps32[iyp*w+ixp];
      imgGreen.drawPixel(ixp,iyp,255,255,0,0.);
      --qlen;
      if((iyp>0)&&(ps32[(iyp-1)*w+ixp] == curDist-1) &&
	 (p8Green[3*((iyp-1)*w+ixp)]==255)&&
	 (p8Green[3*((iyp-1)*w+ixp)+1]==255)&&
	 (p8Green[3*((iyp-1)*w+ixp)+2]==255)){
	imgGreen.drawPixel(ixp,iyp-1,255,255,0,0.);
	rgQx[qlen] = ixp;
	rgQy[qlen] = iyp-1;
	++qlen;
      }
      if((iyp<(h-1))&&(ps32[(iyp+1)*w+ixp] == curDist-1) &&
	 (p8Green[3*((iyp+1)*w+ixp)]==255)&&
	 (p8Green[3*((iyp+1)*w+ixp)+1]==255)&&
	 (p8Green[3*((iyp+1)*w+ixp)+2]==255)){
	imgGreen.drawPixel(ixp,iyp+1,255,255,0,0.);
	rgQx[qlen] = ixp;
	rgQy[qlen] = iyp+1;
	++qlen;
      }
      if((ixp>0)&&(ps32[iyp*w+ixp-1] == curDist-1) &&
	 (p8Green[3*(iyp*w+ixp-1)]==255)&&
	 (p8Green[3*(iyp*w+ixp-1)+1]==255)&&
	 (p8Green[3*(iyp*w+ixp-1)+2]==255)){
	imgGreen.drawPixel(ixp-1,iyp,255,255,0,0.);
	rgQx[qlen] = ixp-1;
	rgQy[qlen] = iyp;
	++qlen;
      }
      if((ixp<(w-1))&&(ps32[iyp*w+ixp+1] == curDist-1) &&
	 (p8Green[3*(iyp*w+ixp+1)]==255)&&
	 (p8Green[3*(iyp*w+ixp+1)+1]==255)&&
	 (p8Green[3*(iyp*w+ixp+1)+2]==255)){
	imgGreen.drawPixel(ixp+1,iyp,255,255,0,0.);
	rgQx[qlen] = ixp+1;
	rgQy[qlen] = iyp;
	++qlen;
      }
    }

    for(int ma_pt=0; ma_pt < lenMA1; ma_pt+=SUBSAMPLE){
      int ixp, iyp;
      ixp = (int)rgMA1X[ma_pt];
      iyp = (int)rgMA1Y[ma_pt];
      if((ixp>=0) && (iyp>=0)){
	D_sint32 curDist;
	D_sint32 *ps32;
	ps32 = (D_sint32*)imgD0.dataPointer_u32();
	curDist=ps32[iyp*w+ixp];
	if((curDist >= 0) && (ixp>=0) && (ixp < w) && (iyp>=0) && (iyp < h) &&
	   (p8Green[3*(iyp*w+ixp)]==255)&&(p8Green[3*(iyp*w+ixp)+1]==255)&&
	   (p8Green[3*(iyp*w+ixp)+2]==255)){//white
	  imgGreen.drawPixel(ixp,iyp,255,255,0,0.);
	  rgQx[qlen] = ixp;
	  rgQy[qlen] = iyp;
	  ++qlen;
	}
      }
    }
    while(qlen > 0){
      int ixp, iyp;
      D_sint32 curDist;
      D_sint32 *ps32;
      ps32 = (D_sint32*)imgD0.dataPointer_u32();
      ixp = rgQx[qlen-1];
      iyp = rgQy[qlen-1];
      curDist=ps32[iyp*w+ixp];
      imgGreen.drawPixel(ixp,iyp,255,255,0,0.);
      --qlen;
      if((iyp>0)&&(ps32[(iyp-1)*w+ixp] == curDist-1) &&
	 (p8Green[3*((iyp-1)*w+ixp)]==255)&&
	 (p8Green[3*((iyp-1)*w+ixp)+1]==255)&&
	 (p8Green[3*((iyp-1)*w+ixp)+2]==255)){
	imgGreen.drawPixel(ixp,iyp-1,255,255,0,0.);
	rgQx[qlen] = ixp;
	rgQy[qlen] = iyp-1;
	++qlen;
      }
      if((iyp<(h-1))&&(ps32[(iyp+1)*w+ixp] == curDist-1) &&
	 (p8Green[3*((iyp+1)*w+ixp)]==255)&&
	 (p8Green[3*((iyp+1)*w+ixp)+1]==255)&&
	 (p8Green[3*((iyp+1)*w+ixp)+2]==255)){
	imgGreen.drawPixel(ixp,iyp+1,255,255,0,0.);
	rgQx[qlen] = ixp;
	rgQy[qlen] = iyp+1;
	++qlen;
      }
      if((ixp>0)&&(ps32[iyp*w+ixp-1] == curDist-1) &&
	 (p8Green[3*(iyp*w+ixp-1)]==255)&&
	 (p8Green[3*(iyp*w+ixp-1)+1]==255)&&
	 (p8Green[3*(iyp*w+ixp-1)+2]==255)){
	imgGreen.drawPixel(ixp-1,iyp,255,255,0,0.);
	rgQx[qlen] = ixp-1;
	rgQy[qlen] = iyp;
	++qlen;
      }
      if((ixp<(w-1))&&(ps32[iyp*w+ixp+1] == curDist-1) &&
	 (p8Green[3*(iyp*w+ixp+1)]==255)&&
	 (p8Green[3*(iyp*w+ixp+1)+1]==255)&&
	 (p8Green[3*(iyp*w+ixp+1)+2]==255)){
	imgGreen.drawPixel(ixp+1,iyp,255,255,0,0.);
	rgQx[qlen] = ixp+1;
	rgQy[qlen] = iyp;
	++qlen;
      }
    }
    // 	while(curDist > 0){
    // 	  if((iyp>0)&&(ps32[(iyp-1)*w+ixp] == curDist-1)){
    // 	    --iyp;
    // 	  }
    // 	  else if((iyp<(h-1))&&(ps32[(iyp+1)*w+ixp] == curDist-1)){
    // 	    ++iyp;
    // 	  }
    // 	  else if((ixp>0)&&(ps32[iyp*w+ixp-1] == curDist-1)){
    // 	    --ixp;
    // 	  }
    // 	  else if((ixp<(w-1))&&(ps32[iyp*w+ixp+1] == curDist-1)){
    // 	    ++ixp;
    // 	  }
    // 	  if((curDist > 0) && (ixp>=0) && (ixp < w) && (iyp>=0) && (iyp < h))
    // 	    imgGreen.drawPixel(ixp,iyp,0,255,0,0.);
    // 	  --curDist;
    // 	}
    //   }
    // }
    imgRet = imgGreen;
    imgGreen.save("/tmp/green.ppm");
    delete [] rgQx;
    delete [] rgQy;
  }




  rgXs = rgPoints0X;
  rgYs = rgPoints0Y;
  for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
    double xp, yp;
    if(fWarpPoints0){
      if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,1.))
	imgRet.drawPixel((int)xp,(int)yp,255,0,0,0.);
    }
    else{
      imgRet.drawPixel((int)rgMA0X[ma_pt],(int)rgMA0Y[ma_pt],255,0,0,0.);
    }
  }
  rgXs = rgPoints1X;
  rgYs = rgPoints1Y;
  D_uint8 *p8;
  p8 = imgRet.dataPointer_u8();
  for(int ma_pt=0; ma_pt < lenMA1; ma_pt+=SUBSAMPLE){
    int px, py;
    px = (int)rgMA1X[ma_pt];
    py = (int)rgMA1Y[ma_pt];
    //int pidx;
    //    pidx = py*w+px;
    // if((p8[3*pidx]==255)&&(p8[3*pidx+1]==0)&&(p8[3*pidx+2]==0))
    //   imgRet.drawPixel((int)rgMA1X[ma_pt],(int)rgMA1Y[ma_pt],0,0,255,0.5);
    // else
      imgRet.drawPixel((int)rgMA1X[ma_pt],(int)rgMA1Y[ma_pt],0,0,255,0.3);
  }

  return imgRet;
}





void DMorphInk::setPreviousPoints(){
  int numPoints;
  numPoints = numPointCols*numPointRows;
  if(NULL != rgPointsPrevX)
    free(rgPointsPrevX);
  if(NULL != rgPointsPrevY)
    free(rgPointsPrevY);
  if(NULL != ma0prevPosX)
    free(ma0prevPosX);
  if(NULL != ma0prevPosY)
    free(ma0prevPosY);

  rgPointsPrevX = (double*)malloc(sizeof(double)*numPoints);
  D_CHECKPTR(rgPointsPrevX);
  rgPointsPrevY = (double*)malloc(sizeof(double)*numPoints);
  D_CHECKPTR(rgPointsPrevY);
  if(lenMA0 < 1)
    fprintf(stderr,"setPreviousPoints() warning: lenMA0 < 1! (%d)\n",lenMA0);
  ma0prevPosX = (double*)malloc(sizeof(double)* lenMA0);
  D_CHECKPTR(ma0prevPosX);
  ma0prevPosY = (double*)malloc(sizeof(double)* lenMA0);
  D_CHECKPTR(ma0prevPosY);
  for(int i=0; i < numPoints; ++i){
    rgPointsPrevX[i] = rgPoints1X[i];
    rgPointsPrevY[i] = rgPoints1Y[i];
  }
  for(int ma_pt=0; ma_pt < lenMA0; ++ma_pt){
    double xp, yp;
    if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,1.0)){
      ma0prevPosX[ma_pt] = xp;
      ma0prevPosY[ma_pt] = yp;
    }
    else{
      ma0prevPosX[ma_pt] = -999.;
      ma0prevPosY[ma_pt] = -999.;
    }
  }
}

/// Returns the morph at some point between 0. and 1. (inclusive)
/** If fWarpOnly is true, there will be no fade/blend, just warping of the
    pixels from src0 to the mesh of src1.  ***This is a simplified function***  It assumes that the mesh at t=0 is rectangular and does not do backward mapping for pixel values and instead does super-sampled forward mapping.*/
DImage DMorphInk::getMorphOverlayAtTime(double dblTime, 
					bool fDrawMesh, double pad_val,
					bool fShowMedialAxis0,
					bool fShowMedialAxis1, int scale,
					bool fStartAtDP,
					bool fStartAtPreviousEnd){
  int w, h;
  double *rgXs, *rgYs;//the warped mesh at time=dblTime
  DImage imgMorphed;
  D_uint8 *pSrc0, *pSrc1;
  D_uint8 *pDst;
  //  double *rgPoints0Xcopy, *rgPoints0Ycopy;
  int numPoints;

  numPoints = numPointCols*numPointRows;
  // printf("B: numPoints=%d\n",numPoints);

  if((dblTime < 0.) || (dblTime > 1.)){
    fprintf(stderr,
  	    "DMorphInk::getMorphOverlayAtTime() dblTime must be >= 0. and <= 1.\n");
    exit(1);
  }
  w = (w0 >= w1) ? w0 : w1;
  h = (h0 >= h1) ? h0 : h1;
  if(pimg0->getImageType() != pimg1->getImageType()){
    fprintf(stderr,
	    "DMorphInk::getMorphOverlayAtTime() types of the two images being morphed must be the same\n");
    return imgMorphed;
  }
  if((pimg0->getImageType() != DImage::DImage_u8) && 
     (pimg0->getImageType() != DImage::DImage_RGB)){
    fprintf(stderr,
	    "DMorphInk::getMorphOverlayAtTime() only supports 8-bit RGB or GS images\n");
    return imgMorphed;
  }
  if(pimg0->getImageType() != DImage::DImage_u8){
    fprintf(stderr,
	    "DMorphInk::getMorphOverlayAtTime() only supports 8-bit GS images so far\n");
    exit(1);
  }


  imgMorphed.create(w*scale,h*scale,pimg0->getImageType());
  pDst = imgMorphed.dataPointer_u8();
  pSrc0 = pimg0->dataPointer_u8();
  pSrc1 = pimg1->dataPointer_u8();


  // define rgXs and rgYs as the intermediate mesh point locations at t=dblTime
  rgXs = (double*)malloc(sizeof(double)*numPointRows*numPointCols);
  D_CHECKPTR(rgXs);
  rgYs = (double*)malloc(sizeof(double)*numPointRows*numPointCols);
  D_CHECKPTR(rgYs);
  for(int r=0, qidx=0; r < numPointRows; ++r){
    for(int c=0; c < numPointCols; ++c,++qidx){
      if(fStartAtDP){
	rgXs[qidx] = (1.0-dblTime)*rgPointsDPX[qidx]+dblTime*rgPoints1X[qidx];
	rgYs[qidx] = (1.0-dblTime)*rgPointsDPY[qidx]+dblTime*rgPoints1Y[qidx];
      }
      else if(fStartAtPreviousEnd){
	rgXs[qidx] = (1.0-dblTime)*rgPointsPrevX[qidx]+dblTime*rgPoints1X[qidx];
	rgYs[qidx] = (1.0-dblTime)*rgPointsPrevY[qidx]+dblTime*rgPoints1Y[qidx];
      }
      else{
	rgXs[qidx] = (1.0-dblTime)*rgPoints0X[qidx]+dblTime*rgPoints1X[qidx];
	rgYs[qidx] = (1.0-dblTime)*rgPoints0Y[qidx]+dblTime*rgPoints1Y[qidx];
      }
    }
  }

  imgMorphed.fill(255);
  if(fDrawMesh){
    double meshTrans;
    meshTrans = 0.2;
    if(fFaintMesh)
      meshTrans = 0.7;
    if(imgMorphed.getImageType() == DImage::DImage_u8)
      imgMorphed = imgMorphed.convertedImgType(DImage::DImage_RGB);
    for(int r=0, meshIdx=0; r < numPointRows; ++r){
      for(int c=0; c < numPointCols; ++c, ++meshIdx){
	if(c>0){
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    imgMorphed.drawLine(rgXs[r*numPointCols+c]*scale,
				rgYs[r*numPointCols+c]*scale+linew,
				rgXs[r*numPointCols+c-1]*scale,
				rgYs[r*numPointCols+c-1]*scale+linew,
				// 255,0,0,0.2);
				//255,255,0,meshTrans, true);
				198,198,0,meshTrans, true);
	  }
	}
	if(c<(numPointCols-1)){
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    imgMorphed.drawLine(rgXs[r*numPointCols+c]*scale,
				rgYs[r*numPointCols+c]*scale+linew,
				rgXs[r*numPointCols+c+1]*scale,
				rgYs[r*numPointCols+c+1]*scale+linew,
				// 0,255,0,0.2);
				//255,255,0,meshTrans, true);
				198,198,0,meshTrans, true);
	  }
	}
	if(r>0){
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    imgMorphed.drawLine(rgXs[r*numPointCols+c]*scale+linew,
				rgYs[r*numPointCols+c]*scale,
				rgXs[(r-1)*numPointCols+c]*scale+linew,
				rgYs[(r-1)*numPointCols+c]*scale,
				// 0,0,255,0.2);
				//255,255,0,meshTrans, true);
				198,198,0,meshTrans, true);
	  }
	}
	if(r<(numPointRows-1)){
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    imgMorphed.drawLine(rgXs[r*numPointCols+c]*scale+linew,
				rgYs[r*numPointCols+c]*scale,
				rgXs[(r+1)*numPointCols+c]*scale+linew,
				rgYs[(r+1)*numPointCols+c]*scale,
				// 255,0,255,0.2);
				//255,255,0,meshTrans, true);
				198,198,0,meshTrans, true);
	  }
	}
      }
    }
  }


  if(fShowMedialAxis0){
    if(imgMorphed.getImageType() == DImage::DImage_u8)
      imgMorphed = imgMorphed.convertedImgType(DImage::DImage_RGB);

    if(fStartAtDP){
      for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
	double xp, yp;
	if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,1.0)){
	  xp = (1.-dblTime)*rgPointsDPX[ma_pt] + dblTime * xp;
	  yp = (1.-dblTime)*rgPointsDPY[ma_pt] + dblTime * yp;
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    imgMorphed.drawLine(xp*scale-scale/2,yp*scale+linew,
			       xp*scale+scale/2,yp*scale+linew,255,0,0,0.,true);
	  }
	}
      }
    }
    else if(fStartAtPreviousEnd){
      for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
	double xp, yp;
	if(NULL == ma0prevPosX){
	  fprintf(stderr,"ma0prevPosX is NULL!\n");
	  exit(1);
	}
	if((warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,1.0)) &&
	   (ma0prevPosX[ma_pt]!=-999.)){
	  xp = (1.-dblTime)*ma0prevPosX[ma_pt] + dblTime * xp;
	  yp = (1.-dblTime)*ma0prevPosY[ma_pt] + dblTime * yp;
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    imgMorphed.drawLine(xp*scale-scale/2,yp*scale+linew,
			       xp*scale+scale/2,yp*scale+linew,255,0,0,0.,true);
	  }
	}
      }
    }
    else{
      for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
	double xp, yp;
	if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,dblTime)){
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    imgMorphed.drawLine(xp*scale-scale/2,yp*scale+linew,
			       xp*scale+scale/2,yp*scale+linew,255,0,0,0.,true);
	  }
	}
      }
    }
  }


  // if(fShowMedialAxis0){
  //   if(imgMorphed.getImageType() == DImage::DImage_u8)
  //     imgMorphed = imgMorphed.convertedImgType(DImage::DImage_RGB);

  //   if(fStartAtDP){
  //     rgPoints0Xcopy = new double[numPoints];
  //     rgPoints0Ycopy = new double[numPoints];
  //     for(int i=0; i < numPoints; ++i){
  // 	rgPoints0Xcopy[i] = rgPoints0X[i];
  // 	rgPoints0Ycopy[i] = rgPoints0Y[i];
  //     }
  //     for(int i=0; i < numPoints; ++i){
  // 	rgPoints0X[i] = rgPointsDPX[i];
  // 	rgPoints0Y[i] = rgPointsDPY[i];
  //     }
  //   }
  //   else if(fStartAtPreviousEnd){
  //     rgPoints0Xcopy = new double[numPoints];
  //     rgPoints0Ycopy = new double[numPoints];
  //     for(int i=0; i < numPoints; ++i){
  // 	rgPoints0Xcopy[i] = rgPoints0X[i];
  // 	rgPoints0Ycopy[i] = rgPoints0Y[i];
  //     }
  //     for(int i=0; i < numPoints; ++i){
  // 	rgPoints0X[i] = rgPointsPrevX[i];
  // 	rgPoints0Y[i] = rgPointsPrevY[i];
  //     }
  //   }

  //   // for(int ma_pt=0; ma_pt < lenMA0; ++ma_pt){
  //   for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
  //     double xp, yp;
  //     //rgDPMA0X
  //     //if(warpPointAtTime(rgDPMA0X[ma_pt],rgDPMA0Y[ma_pt],&xp,&yp,dblTime)){
  //     if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,dblTime)){
  // 	for(int linew=-1*scale/2; linew<=scale/2; ++linew){
  // 	  //	  imgMorphed.drawPixel((int)xp,(int)yp,255,0,0,0.);
  // 	  imgMorphed.drawLine(xp*scale-scale/2,yp*scale+linew,
  // 			      xp*scale+scale/2,yp*scale+linew, 255,0,0,0.,true);
  // 	}
  //     }
  //   }
  //   if(fStartAtDP|fStartAtPreviousEnd){
  //     for(int i=0; i < numPoints; ++i){
  // 	rgPoints0X[i] = rgPoints0Xcopy[i];
  // 	rgPoints0Y[i] = rgPoints0Ycopy[i];
  //     }
  //     delete [] rgPoints0Xcopy;
  //     delete [] rgPoints0Ycopy;
  //     delete [] ma0prevPosX;
  //     delete [] ma0prevPosY;
  //   }
  // }


  if(fShowMedialAxis1){
    if(imgMorphed.getImageType() == DImage::DImage_u8)
      imgMorphed = imgMorphed.convertedImgType(DImage::DImage_RGB);
    // for(int ma_pt=0; ma_pt < lenMA1; ++ma_pt){
    for(int ma_pt=0; ma_pt < lenMA1; ma_pt+=SUBSAMPLE){
      for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	//	imgMorphed.drawPixel((int)rgMA1X[ma_pt],(int)rgMA1Y[ma_pt],0,0,255,0.);
	imgMorphed.drawLine(rgMA1X[ma_pt]*scale-scale/2,
			    rgMA1Y[ma_pt]*scale+linew,
			    rgMA1X[ma_pt]*scale+scale/2,
			    rgMA1Y[ma_pt]*scale+linew,
			    0,0,255,0.,true);
      }
    }
  }


  free(rgXs);
  free(rgYs);

  return imgMorphed;
}



DImage DMorphInk::getLinearHeatmapAtTimeFromDP(double dblTime,
					       bool fShowMedialAxisDP,
					       bool fShowMA0prime,
					       bool fShowMedialAxisTarget,
					       int scale, double *avgDistMoved){
  int w, h;
  DImage imgHeatmap;
  D_uint8 *pDst;
  double sumDistMoved = 0.;
  if((dblTime < 0.) || (dblTime > 1.)){
    fprintf(stderr,
  	    "DMorphInk::getLinearHeatmapAtTimeFromDP() dblTime must be >= 0. and <= 1.\n");
    exit(1);
  }
  w = (w0 >= w1) ? w0 : w1;
  h = (h0 >= h1) ? h0 : h1;
  if(pimg0->getImageType() != pimg1->getImageType()){
    fprintf(stderr,
	    "DMorphInk:::getLinearHeatmapAtTimeFromDP() types of the two images being morphed must be the same\n");
    return imgHeatmap;
  }
  if((pimg0->getImageType() != DImage::DImage_u8) && 
     (pimg0->getImageType() != DImage::DImage_RGB)){
    fprintf(stderr,
	    "DMorphInk:::getLinearHeatmapAtTimeFromDP() only supports 8-bit RGB or GS images\n");
    return imgHeatmap;
  }
  if(pimg0->getImageType() != DImage::DImage_u8){
    fprintf(stderr,
	    "DMorphInk:::getLinearHeatmapAtTimeFromDP() only supports 8-bit GS images so far\n");
    exit(1);
  }
  int rgDistColors[MAXDISTCOLORS*3];
  int stepsize = 255/((MAXDISTCOLORS+3)/4);
  for(int cc=0; cc <= MAXDISTCOLORS/4; ++cc){
    rgDistColors[cc*3+0] = 0;
    rgDistColors[cc*3+1] = cc*stepsize;
    rgDistColors[cc*3+2] = 255;
  }
  for(int cc=0; cc <= MAXDISTCOLORS/4; ++cc){
    rgDistColors[(MAXDISTCOLORS/4+cc)*3+0] = 0;
    rgDistColors[(MAXDISTCOLORS/4+cc)*3+1] = 255;
    rgDistColors[(MAXDISTCOLORS/4+cc)*3+2] = 255-cc*stepsize;
    }
  for(int cc=0; cc <= MAXDISTCOLORS/4; ++cc){
    rgDistColors[(MAXDISTCOLORS/2+cc)*3+0] = cc*stepsize;
    rgDistColors[(MAXDISTCOLORS/2+cc)*3+1] = 255;
    rgDistColors[(MAXDISTCOLORS/2+cc)*3+2] = 0;
  }
  for(int cc=0; cc < MAXDISTCOLORS/4; ++cc){
    rgDistColors[(MAXDISTCOLORS*3/4+cc)*3+0] = 255;
    rgDistColors[(MAXDISTCOLORS*3/4+cc)*3+1] = 255-cc*stepsize;
    rgDistColors[(MAXDISTCOLORS*3/4+cc)*3+2] = 0;
  }

  imgHeatmap.create(w*scale,h*scale,pimg0->getImageType());
  pDst = imgHeatmap.dataPointer_u8();
  imgHeatmap.fill(255);

  if(imgHeatmap.getImageType() == DImage::DImage_u8)
    imgHeatmap = imgHeatmap.convertedImgType(DImage::DImage_RGB);

  if(fShowMedialAxisTarget){
    // for(int ma_pt=0; ma_pt < lenMA1; ++ma_pt){
    for(int ma_pt=0; ma_pt < lenMA1; ma_pt+=SUBSAMPLE){
      for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	//	imgHeatmap.drawPixel((int)rgMA1X[ma_pt],(int)rgMA1Y[ma_pt],0,0,255,0.);
	imgHeatmap.drawLine(rgMA1X[ma_pt]*scale-scale/2,
			    rgMA1Y[ma_pt]*scale+linew,
			    rgMA1X[ma_pt]*scale+scale/2,
			    rgMA1Y[ma_pt]*scale+linew,
			    150,150,150,0.,true);
      }
    }
  }
  
  if(fShowMedialAxisDP){
    // for(int ma_pt=0; ma_pt < lenMA1; ++ma_pt){
    for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
      for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	//	imgHeatmap.drawPixel((int)rgMA1X[ma_pt],(int)rgMA1Y[ma_pt],0,0,255,0.);
	imgHeatmap.drawLine(rgDPMA0X[ma_pt]*scale-scale/2,
			    rgDPMA0Y[ma_pt]*scale+linew,
			    rgDPMA0X[ma_pt]*scale+scale/2,
			    rgDPMA0Y[ma_pt]*scale+linew,
			    200,200,200,0.,true);
      }
    }
  }

  if(fShowMA0prime){
    for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
      double xp, yp;
      if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,1.0)){
	for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	  //	imgHeatmap.drawPixel((int)rgMA1X[ma_pt],(int)rgMA1Y[ma_pt],0,0,255,0.);
	  
	  imgHeatmap.drawLine(xp*scale-scale/2,
			      yp*scale+linew,
			      xp*scale+scale/2,
			      yp*scale+linew,
			      32,32,32,0.,true);
	}
      }
    }
  }
  double distanceMoved;
  int distMoved;
  int numMA0pts;
  numMA0pts = 0;
  for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE, ++numMA0pts){
    double xp, yp, xpp, ypp;
    if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,1.0)){
      xpp = (1.-dblTime)*rgDPMA0X[ma_pt] + dblTime*xp;
      ypp = (1.-dblTime)*rgDPMA0Y[ma_pt] + dblTime*yp;

      distanceMoved = sqrt((xpp - rgDPMA0X[ma_pt])*(xpp - rgDPMA0X[ma_pt]) +
			   (ypp - rgDPMA0Y[ma_pt])*(ypp - rgDPMA0Y[ma_pt]));
      sumDistMoved += distanceMoved;
      distMoved = (int)(distanceMoved*60);
      if(distMoved > (MAXDISTCOLORS-1))
	distMoved = MAXDISTCOLORS-1;
      if(distMoved < 0)
	distMoved = 0;
      for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	imgHeatmap.drawLine(xpp*scale-scale/2,ypp*scale+linew,
			    xpp*scale+scale/2,ypp*scale+linew,
			    rgDistColors[distMoved*3],
			    rgDistColors[distMoved*3+1],
			    rgDistColors[distMoved*3+2],
			    0.,true);
      }	
	
	
      // for(int linew=-1*scale/2; linew<=scale/2; ++linew){
      // 	imgHeatmap.drawLine(xpp*scale-scale/2,
      // 			    ypp*scale+linew,
      // 			      xpp*scale+scale/2,
      // 			      ypp*scale+linew,
      // 			      255,0,0,0.,true);
      // }
    }
  }
  if((NULL != avgDistMoved) && (numMA0pts > 0))
    (*avgDistMoved) = sumDistMoved / numMA0pts;

  return imgHeatmap;
}




/// linearly interpolated the positions of MA0 starting at DP-warped positions
/** dblTime should be on the range [0.,1.] The positions of MA0 are
    linearly interpolated from the DP position to the final morphed
    position. If fShowMedialAxisDP is true, the original positions will be
    shown in gray.  Dark blue is 0 movement.  Red is most movement.
 */
DImage DMorphInk::getMovementHeatmapAtTime(double dblTime, 
					   bool fDrawMesh,
					   bool fShowMedialAxisDP,
					   bool fShowMedialAxisTarget,
					   int scale,
					   bool fStartAtDP,
					   bool fStartAtPreviousEnd){
  int w, h;
  double *rgXs, *rgYs;//the warped mesh at time=dblTime
  DImage imgHeatmap;
  D_uint8 *pDst;
  double *rgPoints0Xcopy = NULL;
  double *rgPoints0Ycopy = NULL;
  int numPoints;

  numPoints = numPointCols*numPointRows;
  // printf("B: numPoints=%d\n",numPoints);

  if((dblTime < 0.) || (dblTime > 1.)){
    fprintf(stderr,
  	    "DMorphInk::getMovementHeatmapAtTime() dblTime must be >= 0. and <= 1.\n");
    exit(1);
  }
  w = (w0 >= w1) ? w0 : w1;
  h = (h0 >= h1) ? h0 : h1;
  if(pimg0->getImageType() != pimg1->getImageType()){
    fprintf(stderr,
	    "DMorphInk:::getMovementHeatmapAtTime() types of the two images being morphed must be the same\n");
    return imgHeatmap;
  }
  if((pimg0->getImageType() != DImage::DImage_u8) && 
     (pimg0->getImageType() != DImage::DImage_RGB)){
    fprintf(stderr,
	    "DMorphInk:::getMovementHeatmapAtTime() only supports 8-bit RGB or GS images\n");
    return imgHeatmap;
  }
  if(pimg0->getImageType() != DImage::DImage_u8){
    fprintf(stderr,
	    "DMorphInk:::getMovementHeatmapAtTime() only supports 8-bit GS images so far\n");
    exit(1);
  }


  imgHeatmap.create(w*scale,h*scale,pimg0->getImageType());
  pDst = imgHeatmap.dataPointer_u8();


  // define rgXs and rgYs as the intermediate mesh point locations at t=dblTime
  rgXs = (double*)malloc(sizeof(double)*numPointRows*numPointCols);
  D_CHECKPTR(rgXs);
  rgYs = (double*)malloc(sizeof(double)*numPointRows*numPointCols);
  D_CHECKPTR(rgYs);
  for(int r=0, qidx=0; r < numPointRows; ++r){
    for(int c=0; c < numPointCols; ++c,++qidx){
      if(fStartAtDP){
	if(fDrawMesh){
	  fprintf(stderr,"getMovementHeatmapAtTime() TURNING OFF fDrawMesh because fStartAtDP is true\n");
	  fDrawMesh = false;
	}
      }
      else if(fStartAtPreviousEnd){
	rgXs[qidx] = (1.0-dblTime)*rgPointsPrevX[qidx]+dblTime*rgPoints1X[qidx];
	rgYs[qidx] = (1.0-dblTime)*rgPointsPrevY[qidx]+dblTime*rgPoints1Y[qidx];
      }
      else{
	rgXs[qidx] = (1.0-dblTime)*rgPoints0X[qidx]+dblTime*rgPoints1X[qidx];
	rgYs[qidx] = (1.0-dblTime)*rgPoints0Y[qidx]+dblTime*rgPoints1Y[qidx];
      }
    }
  }

  imgHeatmap.fill(255);
  if(fDrawMesh){
    double meshTrans;
    meshTrans = 0.2;
    if(fFaintMesh)
      meshTrans = 0.7;
    if(imgHeatmap.getImageType() == DImage::DImage_u8)
      imgHeatmap = imgHeatmap.convertedImgType(DImage::DImage_RGB);
    for(int r=0, meshIdx=0; r < numPointRows; ++r){
      for(int c=0; c < numPointCols; ++c, ++meshIdx){
	if(c>0){
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    imgHeatmap.drawLine(rgXs[r*numPointCols+c]*scale,
				rgYs[r*numPointCols+c]*scale+linew,
				rgXs[r*numPointCols+c-1]*scale,
				rgYs[r*numPointCols+c-1]*scale+linew,
				// 255,0,0,0.2);
				255,255,0,meshTrans, true);
	  }
	}
	if(c<(numPointCols-1)){
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    imgHeatmap.drawLine(rgXs[r*numPointCols+c]*scale,
				rgYs[r*numPointCols+c]*scale+linew,
				rgXs[r*numPointCols+c+1]*scale,
				rgYs[r*numPointCols+c+1]*scale+linew,
				// 0,255,0,0.2);
				255,255,0,meshTrans, true);
	  }
	}
	if(r>0){
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    imgHeatmap.drawLine(rgXs[r*numPointCols+c]*scale+linew,
				rgYs[r*numPointCols+c]*scale,
				rgXs[(r-1)*numPointCols+c]*scale+linew,
				rgYs[(r-1)*numPointCols+c]*scale,
				// 0,0,255,0.2);
				255,255,0,meshTrans, true);
	  }
	}
	if(r<(numPointRows-1)){
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    imgHeatmap.drawLine(rgXs[r*numPointCols+c]*scale+linew,
				rgYs[r*numPointCols+c]*scale,
				rgXs[(r+1)*numPointCols+c]*scale+linew,
				rgYs[(r+1)*numPointCols+c]*scale,
				// 255,0,255,0.2);
				255,255,0,meshTrans, true);
	  }
	}
      }
    }
  }
  if(fShowMedialAxisTarget){
    if(imgHeatmap.getImageType() == DImage::DImage_u8)
      imgHeatmap = imgHeatmap.convertedImgType(DImage::DImage_RGB);
    // for(int ma_pt=0; ma_pt < lenMA1; ++ma_pt){
    for(int ma_pt=0; ma_pt < lenMA1; ma_pt+=SUBSAMPLE){
      for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	//	imgHeatmap.drawPixel((int)rgMA1X[ma_pt],(int)rgMA1Y[ma_pt],0,0,255,0.);
	imgHeatmap.drawLine(rgMA1X[ma_pt]*scale-scale/2,
			    rgMA1Y[ma_pt]*scale+linew,
			    rgMA1X[ma_pt]*scale+scale/2,
			    rgMA1Y[ma_pt]*scale+linew,
			    150,150,150,0.,true);
      }
    }
  }

  if(fShowMedialAxisDP){
    if(imgHeatmap.getImageType() == DImage::DImage_u8)
      imgHeatmap = imgHeatmap.convertedImgType(DImage::DImage_RGB);
    // for(int ma_pt=0; ma_pt < lenMA1; ++ma_pt){
    for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
      for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	//	imgHeatmap.drawPixel((int)rgMA1X[ma_pt],(int)rgMA1Y[ma_pt],0,0,255,0.);
	imgHeatmap.drawLine(rgDPMA0X[ma_pt]*scale-scale/2,
			    rgDPMA0Y[ma_pt]*scale+linew,
			    rgDPMA0X[ma_pt]*scale+scale/2,
			    rgDPMA0Y[ma_pt]*scale+linew,
			    200,200,200,0.,true);
      }
    }
  }



  // if(fShowMedialAxis0){
  if(true){
    if(fStartAtPreviousEnd){
      rgPoints0Xcopy = new double[numPoints];
      rgPoints0Ycopy = new double[numPoints];
      for(int i=0; i < numPoints; ++i){
	rgPoints0Xcopy[i] = rgPoints0X[i];
	rgPoints0Ycopy[i] = rgPoints0Y[i];
      }
      for(int i=0; i < numPoints; ++i){
	rgPoints0X[i] = rgPointsPrevX[i];
	rgPoints0Y[i] = rgPointsPrevY[i];
      }
    }

    if(imgHeatmap.getImageType() == DImage::DImage_u8)
      imgHeatmap = imgHeatmap.convertedImgType(DImage::DImage_RGB);
    // for(int ma_pt=0; ma_pt < lenMA0; ++ma_pt){
    int rgDistColors[MAXDISTCOLORS*3];
    int stepsize = 255/((MAXDISTCOLORS+3)/4);
    for(int cc=0; cc <= MAXDISTCOLORS/4; ++cc){
      rgDistColors[cc*3+0] = 0;
      rgDistColors[cc*3+1] = cc*stepsize;
      rgDistColors[cc*3+2] = 255;
    }
    for(int cc=0; cc <= MAXDISTCOLORS/4; ++cc){
      rgDistColors[(MAXDISTCOLORS/4+cc)*3+0] = 0;
      rgDistColors[(MAXDISTCOLORS/4+cc)*3+1] = 255;
      rgDistColors[(MAXDISTCOLORS/4+cc)*3+2] = 255-cc*stepsize;
    }
    for(int cc=0; cc <= MAXDISTCOLORS/4; ++cc){
      rgDistColors[(MAXDISTCOLORS/2+cc)*3+0] = cc*stepsize;
      rgDistColors[(MAXDISTCOLORS/2+cc)*3+1] = 255;
      rgDistColors[(MAXDISTCOLORS/2+cc)*3+2] = 0;
    }
    for(int cc=0; cc < MAXDISTCOLORS/4; ++cc){
      rgDistColors[(MAXDISTCOLORS*3/4+cc)*3+0] = 255;
      rgDistColors[(MAXDISTCOLORS*3/4+cc)*3+1] = 255-cc*stepsize;
      rgDistColors[(MAXDISTCOLORS*3/4+cc)*3+2] = 0;
    }
    // DImage imgHMPalette;
    // imgHMPalette.create(1,MAXDISTCOLORS,DImage::DImage_RGB);
    // for(int y=0; y < MAXDISTCOLORS; ++y){
    //   imgHMPalette.drawPixel(0,y,rgDistColors[y*3],rgDistColors[y*3+1],
    // 			     rgDistColors[y*3+2]);
    // }
    //    imgHMPalette.save("/tmp/hmpal.ppm");


    for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
      double xp, yp;
      double distanceMoved;
      int distMoved;
      //rgDPMA0X
      //if(warpPointAtTime(rgDPMA0X[ma_pt],rgDPMA0Y[ma_pt],&xp,&yp,dblTime)){
      if(fStartAtDP){
	warpPoint(rgDPMA0X[ma_pt], rgDPMA0Y[ma_pt], &xp, &yp);
	printf("  xp=%.0lf yp=%.0lf   ",xp,yp);
	// xp = (1.0-dblTime)*rgDPMA0X[ma_pt] + dblTime*xp;
	// yp = (1.0-dblTime)*rgDPMA0Y[ma_pt] + dblTime*yp;
	printf("  xpt=%.0lf ypt=%.0lf  dblTime=%lf\n",xp,yp, dblTime);
	distanceMoved = sqrt((xp - rgDPMA0X[ma_pt])*(xp - rgDPMA0X[ma_pt]) +
			     (yp - rgDPMA0Y[ma_pt])*(yp - rgDPMA0Y[ma_pt]));
	distMoved = (int)(distanceMoved*60);
	if(distMoved > (MAXDISTCOLORS-1))
	  distMoved = MAXDISTCOLORS-1;
	if(distMoved < 0)
	  distMoved = 0;
	for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	  imgHeatmap.drawLine(xp*scale-scale/2,yp*scale+linew,
			      xp*scale+scale/2,yp*scale+linew,
			      rgDistColors[distMoved*3],
			      rgDistColors[distMoved*3+1],
			      rgDistColors[distMoved*3+2],
			      0.,true);
	}							       
      }
      else if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,dblTime)){
	distanceMoved = sqrt((xp - rgDPMA0X[ma_pt])*(xp - rgDPMA0X[ma_pt]) +
			     (yp - rgDPMA0Y[ma_pt])*(yp - rgDPMA0Y[ma_pt]));
	distMoved = (int)(distanceMoved*60);
	if(distMoved > (MAXDISTCOLORS-1))
	  distMoved = MAXDISTCOLORS-1;
	if(distMoved < 0)
	  distMoved = 0;
	for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	  //	  imgHeatmap.drawPixel((int)xp,(int)yp,255,0,0,0.);
	  imgHeatmap.drawLine(xp*scale-scale/2,yp*scale+linew,
			      xp*scale+scale/2,yp*scale+linew,
			      rgDistColors[distMoved*3],
			      rgDistColors[distMoved*3+1],
			      rgDistColors[distMoved*3+2],
			      0.,true);
	}
      }
    }
    if(fStartAtPreviousEnd){
      for(int i=0; i < numPoints; ++i){
	rgPoints0X[i] = rgPoints0Xcopy[i];
	rgPoints0Y[i] = rgPoints0Ycopy[i];
      }
      delete [] rgPoints0Xcopy;
      delete [] rgPoints0Ycopy;
      rgPoints0Xcopy = NULL;
      rgPoints0Ycopy = NULL;
    }
  }




  free(rgXs);
  free(rgYs);

  return imgHeatmap;
}




DImage DMorphInk::getBothMedialAxesWithMesh(bool fWarpPoints0, int mesh0or1){
  DImage imgRet;
  double *rgXs, *rgYs;
  int w, h;
  int scale = 8;

  w = w0;
  h = h0;
  //  if(pimg1->width() > w)
  if(w < w1)
    w = w1;
  if(h < h1)
    h = h1;
  
  imgRet.create(w*8,h*8,DImage::DImage_RGB);
  imgRet.fill(255,255,255);

  int numMeshPointCols;
  int numMeshPointRows;

  numMeshPointRows = numPointRows;
  numMeshPointCols = numPointCols;
  if(mesh0or1 == 0){
    int numPoints;

    //        initialColSpacing /=2;
    //        initialRowSpacing /=2;
#if 0
    numMeshPointCols = 1 + (w0-1+initialColSpacing-1)/initialColSpacing;
    numMeshPointRows = 1 + (h0-1+initialRowSpacing-1)/initialRowSpacing;

    numPoints = numMeshPointCols * numMeshPointRows;

    rgXs = (double*)malloc(sizeof(double)*numPoints);
    D_CHECKPTR(rgXs);
    rgYs = (double*)malloc(sizeof(double)*numPoints);
    D_CHECKPTR(rgYs);

    //set the x,y position of each control point in mesh0
    for(int py=0, idx=0; py < numMeshPointRows; ++py){
      for(int px = 0; px < numMeshPointCols; ++px, ++idx){
	rgXs[idx] = px * initialColSpacing;
	if(rgXs[idx] >= w0)//last column may not be as wide
	  rgXs[idx] = w0-1;
	rgYs[idx] = py * initialRowSpacing;
	if(rgYs[idx] >= h0)//last row may not be as tall
	  rgYs[idx] = h0-1;
	printf("row %d col %d: x,y= %lf,%lf\n",py,px,rgXs[idx],rgYs[idx]);
      }
    }
#else
    numMeshPointCols = 1 + (w0-1+initialColSpacing-1)/initialColSpacing;
    numMeshPointRows = 1 + (h0-1+initialRowSpacing-1)/initialRowSpacing;

    numMeshPointCols += numMeshPointCols-1;
    numMeshPointRows += numMeshPointRows-1;

    numPoints = numMeshPointCols * numMeshPointRows;

    rgXs = (double*)malloc(sizeof(double)*numPoints);
    D_CHECKPTR(rgXs);
    rgYs = (double*)malloc(sizeof(double)*numPoints);
    D_CHECKPTR(rgYs);

    //set the x,y position of each control point in mesh0
    for(int py=0, idx=0; py < numMeshPointRows; ++py){
      for(int px = 0; px < numMeshPointCols; ++px, ++idx){
	rgXs[idx] = px * initialColSpacing/2.;
	if(rgXs[idx] >= w0)//last column may not be as wide
	  rgXs[idx] = w0-1;
	rgYs[idx] = py * initialRowSpacing/2.;
	if(rgYs[idx] >= h0)//last row may not be as tall
	  rgYs[idx] = h0-1;
	printf("row %d col %d: x,y= %lf,%lf\n",py,px,rgXs[idx],rgYs[idx]);
      }
    }

#endif
    printf("numMeshPointRows=%d numMeshPointCols=%d xlast=%lf ylast=%lf\n",
	   numMeshPointRows,numMeshPointCols,rgXs[numMeshPointCols-1],
	   rgYs[numMeshPointCols*numMeshPointRows-1]);
  }
  else if(mesh0or1 == 1){
    rgXs = rgPoints1X;
    rgYs = rgPoints1Y;
    printf("h0=%d initialRowSpacing=%d\n",h0, initialRowSpacing);
    printf("numPointRows=%d numPointCols=%d xlast=%lf ylast=%lf\n",
	   numPointRows,numPointCols,rgXs[numPointRows-1],
	   rgYs[numPointCols*numPointRows-1]);
	   
  }
  else{
    fprintf(stderr,"DMorphInk::getBothMedialAxesWithMesh() mesh0or1 must be 0 or 1 (was %d)\n",mesh0or1);
    exit(1);
  }

  // DRAW THE MESH
  if(true){
    double meshTrans;
    meshTrans = 0.0;
    if(fFaintMesh)
      meshTrans = 0.0;
    int r1,g1,b1,r2,g2,b2, rr,gg,bb;
    r1=199; g1=199; b1=0;
    r2=199; g2=199; b2=0;
    for(int r=0, meshIdx=0; r < numMeshPointRows; ++r){
      for(int c=0; c < numMeshPointCols; ++c, ++meshIdx){
	if(c>0){
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    if(1==r%2){
	      rr=r2; gg=g2; bb=b2;
	    }
	    else{
	      rr=r1; gg=g1; bb=b1;
	    }
	    imgRet.drawLine(rgXs[r*numMeshPointCols+c]*scale,
				rgYs[r*numMeshPointCols+c]*scale+linew,
				rgXs[r*numMeshPointCols+c-1]*scale,
				rgYs[r*numMeshPointCols+c-1]*scale+linew,
			    rr,gg,bb,meshTrans, true);
				// 128,128,128,meshTrans, true);
	  }
	}
	if(c<(numMeshPointCols-1)){
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    if(1==r%2){
	      rr=r2; gg=g2; bb=b2;
	    }
	    else{
	      rr=r1; gg=g1; bb=b1;
	    }
	    imgRet.drawLine(rgXs[r*numMeshPointCols+c]*scale,
				rgYs[r*numMeshPointCols+c]*scale+linew,
				rgXs[r*numMeshPointCols+c+1]*scale,
				rgYs[r*numMeshPointCols+c+1]*scale+linew,
			    rr,gg,bb,meshTrans, true);
				// 128,128,128,meshTrans, true);
	  }
	}
	if(r>0){
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    if(1==c%2){
	      rr=r2; gg=g2; bb=b2;
	    }
	    else{
	      rr=r1; gg=g1; bb=b1;
	    }
	    imgRet.drawLine(rgXs[r*numMeshPointCols+c]*scale+linew,
				rgYs[r*numMeshPointCols+c]*scale,
				rgXs[(r-1)*numMeshPointCols+c]*scale+linew,
				rgYs[(r-1)*numMeshPointCols+c]*scale,
			    rr,gg,bb,meshTrans, true);
				// 128,128,128,meshTrans, true);
	  }
	}
	if(r<(numMeshPointRows-1)){
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    if(1==c%2){
	      rr=r2; gg=g2; bb=b2;
	    }
	    else{
	      rr=r1; gg=g1; bb=b1;
	    }
	    imgRet.drawLine(rgXs[r*numMeshPointCols+c]*scale+linew,
				rgYs[r*numMeshPointCols+c]*scale,
				rgXs[(r+1)*numMeshPointCols+c]*scale+linew,
				rgYs[(r+1)*numMeshPointCols+c]*scale,
			    rr,gg,bb,meshTrans, true);
				// 128,128,128,meshTrans, true);
	  }
	}
      }
    }
  }
  // scale the mesh image down

  if(mesh0or1 == 0){
    free(rgXs);
    free(rgYs);
  }
  DImage imgRetTmp;
  imgRetTmp = imgRet;
  imgRetTmp.scaledDownPow2_(imgRet, 2);


  for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
    double xp, yp;
    if(fWarpPoints0){
      if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,1.)){
	imgRet.drawPixel(2*(int)xp,2*(int)yp,255,0,0,0.);
	imgRet.drawPixel(2*(int)xp+1,2*(int)yp,255,0,0,0.);
	imgRet.drawPixel(2*(int)xp,2*(int)yp+1,255,0,0,0.);
	imgRet.drawPixel(2*(int)xp+1,2*(int)yp+1,255,0,0,0.);
      }
    }
    else{
      imgRet.drawPixel(2*(int)rgMA0X[ma_pt],2*(int)rgMA0Y[ma_pt],255,0,0,0.);
      imgRet.drawPixel(2*(int)rgMA0X[ma_pt]+1,2*(int)rgMA0Y[ma_pt],255,0,0,0.);
      imgRet.drawPixel(2*(int)rgMA0X[ma_pt],2*(int)rgMA0Y[ma_pt]+1,255,0,0,0.);
      imgRet.drawPixel(2*(int)rgMA0X[ma_pt]+1,2*(int)rgMA0Y[ma_pt]+1,255,0,0,0.);
    }
  }

  D_uint8 *p8;
  p8 = imgRet.dataPointer_u8();
  for(int ma_pt=0; ma_pt < lenMA1; ma_pt+=SUBSAMPLE){
    int px, py, pidx;
    px = (int)rgMA1X[ma_pt];
    py = (int)rgMA1Y[ma_pt];
    pidx = 2*py*2*w+2*px;
    //clear mesh pixels since we're using transparency, but not red medial axis
    if((p8[3*pidx]!=255)||(p8[3*pidx+1]!=0)||(p8[3*pidx+2]!=0)){
      imgRet.drawPixel(2*(int)rgMA1X[ma_pt],2*(int)rgMA1Y[ma_pt],255,255,255,0.);
      imgRet.drawPixel(2*(int)rgMA1X[ma_pt]+1,2*(int)rgMA1Y[ma_pt],255,255,255,0.);
      imgRet.drawPixel(2*(int)rgMA1X[ma_pt],2*(int)rgMA1Y[ma_pt]+1,255,255,255,0.);
      imgRet.drawPixel(2*(int)rgMA1X[ma_pt]+1,2*(int)rgMA1Y[ma_pt]+1,255,255,255,0.);
    }
    //draw blue with transparency
    imgRet.drawPixel(2*(int)rgMA1X[ma_pt],2*(int)rgMA1Y[ma_pt],0,0,255,0.3);
    imgRet.drawPixel(2*(int)rgMA1X[ma_pt]+1,2*(int)rgMA1Y[ma_pt],0,0,255,0.3);
    imgRet.drawPixel(2*(int)rgMA1X[ma_pt],2*(int)rgMA1Y[ma_pt]+1,0,0,255,0.3);
    imgRet.drawPixel(2*(int)rgMA1X[ma_pt]+1,2*(int)rgMA1Y[ma_pt]+1,0,0,255,0.3);
  }

  return imgRet;
}



DImage DMorphInk::getWarpedMedialAxesWithSmoothMesh(bool fShowMA0,
						    bool fShowMA1){
  DImage imgRet;
  double *rgXs, *rgYs;
  int w3, h3;
  int scale = 20;
  double invScale = 1./5.;
  //  int scale = 1;
  //  double invScale = 1.;

  int h1,w1,h2,w2;
  w1 = pimg0->width();
  w2 = pimg1->width();
  h1 = pimg0->height();
  h2 = pimg1->height();

  w3 = w1;
  if(w3 < w2)
    w3 = w2;
  h3 = h1;
  if(h3 < h2)
    h3 = h2;

  imgRet.create(w3*scale,h3*scale,DImage::DImage_RGB);
  imgRet.fill(255,255,255);
  rgXs = rgPoints1X;
  rgYs = rgPoints1Y;

  double meshTrans = 0.;

  for(int r=0, meshIdx=0; r < numPointRows; ++r){
    for(int c=0; c < numPointCols; ++c, ++meshIdx){
      if(c>0){
	for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	  imgRet.drawLine(rgXs[r*numPointCols+c]*scale,
			      rgYs[r*numPointCols+c]*scale+linew,
			      rgXs[r*numPointCols+c-1]*scale,
			      rgYs[r*numPointCols+c-1]*scale+linew,
			      // 198,198,0,meshTrans, true);
			      244,184,0,meshTrans, true);
	}
      }
      if(c<(numPointCols-1)){
	for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	  imgRet.drawLine(rgXs[r*numPointCols+c]*scale,
			      rgYs[r*numPointCols+c]*scale+linew,
			      rgXs[r*numPointCols+c+1]*scale,
			      rgYs[r*numPointCols+c+1]*scale+linew,
			      // 198,198,0,meshTrans, true);
			      244,184,0,meshTrans, true);
	}
      }
      if(r>0){
	for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	  imgRet.drawLine(rgXs[r*numPointCols+c]*scale+linew,
			      rgYs[r*numPointCols+c]*scale,
			      rgXs[(r-1)*numPointCols+c]*scale+linew,
			      rgYs[(r-1)*numPointCols+c]*scale,
			      // 198,198,0,meshTrans, true);
			      244,184,0,meshTrans, true);
	}
      }
      if(r<(numPointRows-1)){
	for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	  imgRet.drawLine(rgXs[r*numPointCols+c]*scale,
			  rgYs[r*numPointCols+c]*scale,
			  rgXs[(r+1)*numPointCols+c]*scale+linew,
			  rgYs[(r+1)*numPointCols+c]*scale,
			  // 198,198,0,meshTrans, true);
			  244,184,0,meshTrans, true);
	}
      }
    }
  }
  //  imgRet = imgRet.scaled(invScale,invScale,DImage::DImageTransSmooth);
  // imgRet.save("/tmp/_imgRet.ppm");
  // char stCmd[1024];
  // sprintf(stCmd,"pnmscale %lf /tmp/_imgRet.ppm > /tmp/imgRet.ppm",invScale);
  // system(stCmd);
  // imgRet.load("/tmp/imgRet.ppm");

  if(fShowMA0){
    for(int ma_pt=0; ma_pt < lenMA0; ma_pt+=SUBSAMPLE){
      double xp, yp;
      if(warpPoint(rgMA0X[ma_pt],rgMA0Y[ma_pt], &xp, &yp)){
	int x,y;
	x = (int)xp;
	y = (int)yp;
	// imgRet.drawPixel(x,y,255,0,0,0.);
	for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	  imgRet.drawLine(x*scale-scale/2,y*scale+linew,
			  x*scale+scale/2,y*scale+linew,255,0,0,0.,true);
	}
      }
    }
  }
  if(fShowMA1){
    D_uint8 *p8;
    p8 = imgRet.dataPointer_u8();
    for(int ma_pt=0; ma_pt < lenMA1; ma_pt+=SUBSAMPLE){
      //if the image pixel is red, draw purple, otherwise draw blue
      int x, y;
      x = (int)(rgMA1X[ma_pt]);
      y = (int)(rgMA1Y[ma_pt]);
      if(fShowMA0){
	if(( (p8[3*(y*scale*w3*scale+x*scale)+0]==(D_uint8)255) &&
	     (p8[3*(y*scale*w3*scale+x*scale)+1]==(D_uint8)0) &&
	     (p8[3*(y*scale*w3*scale+x*scale)+2]==(D_uint8)0))){//draw purple (76,0,17
	  //if(( (p8[3*(y*scale*w3+x*scale)+0]==(D_uint8)244) &&
	  //(p8[3*(y*scale*w3+x*scale)+1]==(D_uint8)184) &&
	  //(p8[3*(y*scale*w3+x*scale)+2]==(D_uint8)0))){//draw purple (76,0,178)
	  //imgRet.drawPixel(x,y,76,0,178,0.);
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    imgRet.drawLine(x*scale-scale/2,y*scale+linew,
			    x*scale+scale/2,y*scale+linew,76,0,178,0.,true);
	  }
	}
	else{//draw blue
	  //	  imgRet.drawPixel(x,y,0,0,255,0.);
	  for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	    imgRet.drawLine(x*scale-scale/2,y*scale+linew,
			    x*scale+scale/2,y*scale+linew,0,0,255,0.,true);
	  }
	  // p8[3*(y*scale*w3+x*scale)+0] = 0;
	  // p8[3*(y*scale*w3+x*scale)+1] = 255;
	  // p8[3*(y*scale*w3+x*scale)+2] = 0;
	}
      }
      else{//draw blue
	//	imgRet.drawPixel(x,y,0,0,255,0.);
	for(int linew=-1*scale/2; linew<=scale/2; ++linew){
	  imgRet.drawLine(x*scale-scale/2,y*scale+linew,
			  x*scale+scale/2,y*scale+linew,0,0,255,0.,true);
	}
      }
    }
  }

  imgRet.save("/tmp/_imgRet.ppm");
  char stCmd[1024];
  int itmp;
  sprintf(stCmd,"pnmscale %lf /tmp/_imgRet.ppm > /tmp/imgRet.ppm",invScale);
  itmp = system(stCmd);
  imgRet.load("/tmp/imgRet.ppm");

    // if(warpPointAtTime(rgMA0X[ma_pt],rgMA0Y[ma_pt],&xp,&yp,1.))
    //   imgRet.drawPixel((int)xp,(int)yp,255,0,0,0.);
  return imgRet;
}






