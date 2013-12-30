#ifndef DMORPHINK_H
#define DMORPHINK_H

#include "dimage.h"
#include "dmath.h"
#include <stdio.h>

#define SPEED_TEST 1
#define SPEED_TEST2 1
#define SPEED_TEST3 1



#define DEBUG_WARP_POINT 0

#define NEW_WARP 1 /*attempt at faster warpPoint function & re-tooling for it*/

//#define COMPENSATE_FOR_OOB_DISTMAP 1 //we now do this all the time

//#define CHECK_OLD_WARP 1 /*I fixed the difference between warpPoint and warpPointAtTime after the dissertation.  They both use >0. now instead of warpPoint using >=1.0 like it used to.*/


///This class provides ink-morphing functionality to morph one word image to another.
/** This class provides ink-morphing functionality to morph one word
    image (src0) to another (src1).  Each image has a mesh associated
    with it.  mesh0 is assumed to be rectangular and mesh1 will get
    warped.  Mesh0 is evenly spaced every N pixels and Mesh1 is given
    initial alignment using Dynamic Programming (DP) mapping of the
    control points.

    compute medial axis of img0.  compute distance map of img1.

    The "improve morph" algorithm then consists of the
    following:  for each control point find the best new
    location of any within the perturbation area.  The best location
    is that with lowest morph cost.  Morph cost calculated as sum of
    distance map values at warped position of each medial axis pixel.
    
    Refine mesh algorithm is just to split each existing quad into 4 quads
    in both meshes without moving any control points. (Just insert a control
    point at the midpoint of each mesh edge)
    
    Overall process for auto-morph:
    improve morph a few times, then refine once, then improve morph a few.
    Repeat until total cost isn't changing much or for a set number of times.

**/


class DMorphInk{
public:
  DMorphInk();
  ~DMorphInk();
  void morphOneWay(	const DImage &srcFrom, 
				const DImage &srcTo, 
				int bandWidthDP, 
				double nonDiagonalCostDP, 
				int meshSpacingStatic,
				int numRefinementsStatic, 
				double meshDiv);
  void saveCurrentMorph(int iteration);
  void saveCurrentMAImagesAndControlPoints(int id);
  double getWordMorphCost(const DImage &src0, const DImage &src1,
			  int bandWidthDP = 15,
			  double nonDiagonalCostDP=0.,
			  int meshSpacingStatic=-1,
			  int numRefinementsStatic=-1,
			  double meshDiv=4.0,
			  double lengthMismatchPenalty=0.0);
  double getWordMorphCostFast(const DImage &src0, const DImage &src1,
					   int bandWidthDP = 15,
					   double nonDiagonalCostDP=0.,
					   int meshSpacingStatic=-1,
					   int numRefinementsStatic=-1,
					   double meshDiv=4.0,
					   double lengthMismatchPenalty=0.0);

  void init(const DImage &src0, const DImage &src1, bool fMakeCopies = true,
	    int initialMeshSpacing=50, int bandRadius=15,
	    double nonDiagonalDPcost = 0.);
  void resetMeshes(int columnSpacing, int rowSpacing,int bandRadius=15,
		   double nonDiagonalDPcost = 0.);
  void refineMeshes();
  void improveMorph();
  void improveMorphFast();
  void improveMorphSimAnn();
  		double acceptProbability(double newCost, double oldCost, double temp);
  		double get_new_temp(int time);
  double getVertexPositionCost(int r, int c, double x1, double y1);
  bool warpPoint(double x, double y, double *xp, double *yp, bool fDebugPrint=false);
#if NEW_WARP
  double getVertexPositionCostNew(int r, int c, double x1, double y1);
  void warpPointNew(double s, double t, int meshPointIdx,double *xp, double *yp
		    /*debug: ,double x,double y, int ii=-1, int l=-1,
		    int rgTempMA0idx_i=-1, int r=-1, int c=-1*/);
#endif//NEW_WARP
  double getCost();

  // bool isCoordInQuad(double x, double y, int quadIdx,
  // 		     double *rgMeshXs, double *rgMeshYs);
  // int getQuadIdxFromCoord(double x, double y, double *rgXs, double *rgYs,
  // 			  int lastQuadIdx = -1);
  // double bilinearInterpDistMap(double xp, double yp, signed int *psDist,
  // 			       int w, int h);
  //-------------------------------------------------------------------------
  // functions for visualization / debug
  //-------------------------------------------------------------------------
  DImage getImageWithMesh(int meshNum);
  DImage getMedialAxisWithMesh(int meshNum);
  DImage getBothMedialAxes(bool fWarpPoints0, bool fShowGreenDistance = false);
  DImage getBothMedialAxesWithMesh(bool fWarpPoints0, int mesh0or1);
  DImage getWarpedMedialAxis0_A0prime();
  DImage getMorphAtTime(double dblTime, bool fWarpOnly = false,
			bool fDrawMesh = false, double pad_val = 128.,
			bool fOverlayImg1=false, bool fShowMedialAxis0 = false,
			bool fShowMedialAxis1 = false, bool fNoMorph=false);
  DImage getMorphOverlayAtTime(double dblTime, 
			       bool fDrawMesh = false, double pad_val = 128.,
			       bool fShowMedialAxis0 = false,
			       bool fShowMedialAxis1 = false, int scale=1,
			       bool fStartAtDP = false,
			       bool fStartAtPreviousEnd = false);
  DImage getMovementHeatmapAtTime(double dblTime, 
				  bool fDrawMesh = false,
				  bool fShowMedialAxisDP = false,
				  bool fShowMedialAxisTarget = false,
				  int scale=1,
				  bool fStartAtDP = false,
				  bool fStartAtPreviousEnd = false);
  DImage getLinearHeatmapAtTimeFromDP(double dblTime,
				      bool fShowMedialAxisDP = false,
				      bool fShowMA0prime = false,
				      bool fShowMedialAxisTarget = false,
				      int scale = 1, double *avgDistMoved=NULL);
  DImage getWarpedMedialAxesWithSmoothMesh(bool fShowMA0=true,
					   bool fShowMA1=true);
  void setPreviousPoints();
  DImage getMorphAtTimeForOtherImages(DImage &img0, DImage &img1);

  bool warpPointAtTime(double x,double y, double *xp, double *yp, double t=1.0);

  //  double getCost_fast();//for videos

  // void flowInkMA0to1(

  double warpCostDPfull;//cost of horizontal+vertical both directions
  //the next two will only have the most recent costs (second direction unless
  //fOnlyDoOneDirection is true in which case only the first direction).
  double warpCostDP;//cost of horizontal DP warping calculated by resetMeshes()
  double warpCostDPv;//cost of vertical DP warping calculated by resetMeshes()
  double warpCostDPhoriz;//cost of horizontal DP warping (forward+backward)

  double avgMovement0to1;
  double avgMovement1to0;
  double distmapCost0to1;
  double distmapCost1to0;
  double cost0to1;
  double cost1to0;

  bool fOnlyDoOneDirection;
  bool fOnlyDoCoarseAlignment;

  //private:
  int numPointRows;//numPointRows and numPointCols may differ, but must be kept
  int numPointCols;//consistent between the two meshes
  
  double *rgPoints0X;//x-coords of mesh0 quad vertices (row-major order)
  double *rgPoints0Y;//y-coords
  double *rgPoints1X;//same, but for mesh1
  double *rgPoints1Y;//y-coords for mesh1
  double *rgPointsDPX;//rgPoints1X copied right after coarse alignment with DP
  double *rgPointsDPY;//rgPoints1Y copied right after coarse alignment with DP
  double *rgPointsPrevX;
  double *rgPointsPrevY;
  double *ma0prevPosX;//prev positions of each medial axis pixel from ma0
  double *ma0prevPosY;//set to -999. if warpPointAtTime returns false
  double *rgMA0X;//x-coords for medial axis 0
  double *rgMA0Y;//y-coords for medial axis 0
  double *rgMA1X;//x-coords for medial axis 1
  double *rgMA1Y;//y-coords for medial axis 1
#if NEW_WARP
  //these arrays hold values that used to be calculated for every warpPoint call
  int *rgMA0r;
  int *rgMA0c;
  double *rgMA0s;
  double *rgMA0t;
  double *rgMA0quadW;
  double *rgMA0quadH;
  //this array holds the indexes for rgMA0 points that are in 4 adjacent quads
  //instead of copying everything into rgTempMA0X/Y, etc, we now just
  //copy the indexes and look them up in the original full arrays rgMA0X, etc.
  int *rgTempMA0idx;
#endif

  double *rgDPMA0X;//DP-warped coords of MedialAxis0
  double *rgDPMA0Y;//DP-warped coords of MedialAxis0

  double *rgTempMA0X;//used by improveMorph,getVertexPositionCost2 
  double *rgTempMA0Y;//used by improveMorph,getVertexPositionCost2 
  int tempMA0len;//used by improveMorph,getVertexPositionCost2 
  int lenMA0;//number of MA0 (medial axis 0) points
  int lenMA1;//number of MA1 (medial axis 1) points
  DImage *pimg0,*pimg1;//pointers to the source images (or copies if applicable)
  DImage img0tmp, img1tmp;//if images are copied by init(), point to these
  DImage imgDist0, imgDist1;//distance maps (DImage_u32) for each image
  DImage imgMA0, imgMA1;//medialAxis images
  DImage imgInk0;//medial axis (or all ink) in image0
  int w0, h0, w1, h1; //width,height of the two images (don't have to be equal)
  int initialColSpacing, initialRowSpacing;
  double curColSpacing, curRowSpacing;
  int meshLevel;
  double *rgBackmapX1toX0;
  double *rgBackmapY1toY0;
  bool fFaintMesh;


  bool *rgPlacePoints0;//if false, no points in 4 adjacent quads, so don't place
  
  
  	struct neighbor_point {
		double cost;
		int x;
		int y;
		bool operator<(const neighbor_point& other) const
			{ return cost > other.cost; }
	};
};



///This function is used for visualization / debug only.
/**given point x, y in mesh0, find the warped coordinate at time=t [0. to 1.].*/
inline bool DMorphInk::warpPointAtTime(double x, double y, double *xp,
				       double *yp, double dblTime){
  int r, c, idx;
  double Xtop, Ytop, Xbot, Ybot;
  double s,t;
  double Ax, Ay, Bx, By, Cx, Cy, Dx, Dy;//quad vertices at dblTime (NW,NE,SE,SW)

  if((x<0)||(y<0)||(x>(w0-1))||(y>(h0-1))){
    // fprintf(stderr,
    // 	    "DMorphInk::warpPoint() point %lf,%lf out of img0 bounds "
    // 	    "(0,0-%d,%d)\n",x,y,w0-1,h0-1);
    (*xp)=0.;
    (*yp)=0.;
    return false;
  }
  if((dblTime < 0.)||(dblTime > 1.0)){
    fprintf(stderr,"DMorphInk::warpPointAtTime() dblTime must be 0. to 1.\n");
    exit(1);
  }

  //r,c are easy to find because mesh0 is evenly spaced
  s = x / curColSpacing;
  t = y / curRowSpacing;
  c = (int)s;
  if(c == (numPointCols-1))//if last pixel, make it s=1.0 for c=numPointsCols-2
    --c;
  r = (int)t;
  if(r == (numPointRows-1))//if last pixel, make it t=1.0 for r=numPointsRows-2
    --r;
  idx = r*numPointCols+c;
  s -= c;
  t -= r;
  
  double quadW;
  double quadH;
  quadW = curColSpacing;
  quadH = curRowSpacing;
  if(c==(numPointCols-2)){
    //quadW = rgPoints0X[idx+1] - rgPoints0X[idx];
    //if(quadW>=1.)
    //  s = (x - rgPoints0X[idx]) / quadW;
    quadW = (w0-1)-(c*curColSpacing);
    if(quadW>0.)
      s = (x - c*curColSpacing) / quadW;
    else
      s = 0.;
  }
  if(r==(numPointRows-2)){
    // quadH = rgPoints0Y[idx+numPointCols] - rgPoints0Y[idx];
    // if(quadH>=1.)
    //   t = (y - rgPoints0Y[idx]) / quadH;
    // else
    //   t = 0.;
    quadH = (h0-1) - (r*curRowSpacing);
    if(quadH>0.)
      t = (y - r*curRowSpacing) / quadH;
    else
      t = 0.;
  }
  if((s > 1.)||(t>1.)){
    fprintf(stderr,"WPAT: bad s or t!   x,y=%lf,%lf  r=%d c=%d  idx=%d numPointsCols=%d  s=%lf t=%lf  quadW=%lf quadH=%lf\n",x,y,r,c,idx,numPointCols,s,t,quadW,quadH);
    fprintf(stderr,"      vertices0: %lf,%lf  %lf,%lf  %lf,%lf  %lf,%lf\n",
	    rgPoints0X[idx],rgPoints0Y[idx],
	    rgPoints0X[idx+1],rgPoints0Y[idx+1],
	    rgPoints0X[idx+numPointCols+1],rgPoints0Y[idx+numPointCols+1],
	    rgPoints0X[idx+numPointCols],rgPoints0Y[idx+numPointCols]);
    fprintf(stderr,"      vertices1: %lf,%lf  %lf,%lf  %lf,%lf  %lf,%lf\n",
	    rgPoints1X[idx],rgPoints1Y[idx],
	    rgPoints1X[idx+1],rgPoints1Y[idx+1],
	    rgPoints1X[idx+numPointCols+1],rgPoints1Y[idx+numPointCols+1],
	    rgPoints1X[idx+numPointCols],rgPoints1Y[idx+numPointCols]);
    fprintf(stderr,"      CLAMPING VALUE to range 0. to 1.\n");
  }
  if(s<0.)
    s=0.;
  if(t<0.)
    t=0.;
  if(s>1.)
    s=1.;
  if(t>1.)
    t=1.;

  //vertices are in clockwise order (upper-left, up-rt, lower-rt, lower-left)
  //calculate the quad vertices at time=dblTime:
  Ax = (1.0-dblTime)*rgPoints0X[idx]+dblTime*rgPoints1X[idx];
  Ay = (1.0-dblTime)*rgPoints0Y[idx]+dblTime*rgPoints1Y[idx];
  Bx = (1.0-dblTime)*rgPoints0X[idx+1]+dblTime*rgPoints1X[idx+1];
  By = (1.0-dblTime)*rgPoints0Y[idx+1]+dblTime*rgPoints1Y[idx+1];
  Cx = (1.0-dblTime)*rgPoints0X[idx+numPointCols+1]+
    dblTime*rgPoints1X[idx+numPointCols+1];
  Cy = (1.0-dblTime)*rgPoints0Y[idx+numPointCols+1]+
    dblTime*rgPoints1Y[idx+numPointCols+1];
  Dx = (1.0-dblTime)*rgPoints0X[idx+numPointCols]+
    dblTime*rgPoints1X[idx+numPointCols];
  Dy = (1.0-dblTime)*rgPoints0Y[idx+numPointCols]+
    dblTime*rgPoints1Y[idx+numPointCols];
  //bilinearly interpolate the point x,y to get it's warped coordinates xp,yp:
  Xtop = (1.-s)*Ax + s*Bx;
  Ytop = (1.-s)*Ay + s*By;
  Xbot = (1.-s)*Dx + s*Cx;
  Ybot = (1.-s)*Dy + s*Cy;
  (*xp) = (1.-t)*Xtop + t*Xbot;
  (*yp) = (1.-t)*Ytop + t*Ybot;
  return true;
}

/**given point x, y in mesh0, find the warped coordinate at time=t [0. to 1.].*/
inline bool DMorphInk::warpPoint(double x, double y, double *xp, double *yp,bool fDebugPrint){
  int r, c, idx;
  double s,t;

#if DEBUG_WARP_POINT
  if(fDebugPrint){
    printf("warpPoint() inputs: x=%lf y=%lf\n",x,y);
  }

  if((x<0)||(y<0)||(x>(w0-1))||(y>(h0-1))){
    (*xp)=x;
    (*yp)=y;
    fprintf(stderr,
	    "warpPoint(x,y) out of image0. x=%lf,y=%lf w0=%d,h0=%d\n",
	    x,y,w0,h0);
    return false;
  }
#endif

  //r,c are easy to find because mesh0 is evenly spaced
  s = x / curColSpacing;
  t = y / curRowSpacing;
  // c = (int)s;
  // r = (int)t;
  c = (int)s;
  if(c == (numPointCols-1))//if last pixel, make it s=1.0 for c=numPointsCols-2
    --c;
  r = (int)t;
  if(r == (numPointRows-1))//if last pixel, make it t=1.0 for r=numPointsRows-2
    --r;
  idx = r*numPointCols+c;
  s -= c;
  t -= r;
#if DEBUG_WARP_POINT
  if(fDebugPrint){
    printf("  s=%lf t=%lf  c=%d r=%d\n",s,t, c, r);
  }
#endif

  double quadW;
  double quadH;
  quadW = curColSpacing;
  quadH = curRowSpacing;
  if(c==(numPointCols-2)){
    quadW = rgPoints0X[idx+1] - rgPoints0X[idx];
    if(quadW>0.)//changed from >=1. to > 0.(changed after dissertation)
      s = (x - rgPoints0X[idx]) / quadW;
    else //Brian: prevents control point cross-over?
      s = 0.;
// #if CHECK_OLD_WARP
//     {
//       double otherQuadW;
//       double otherS;
//       otherQuadW = (w0-1)-(c*curColSpacing);
//       if(abs(otherQuadW-quadW) > 0.0001){
// 	fprintf(stderr,"CHECK_OLD_WARP found discrepancy: quadW=%lf otherQuadW=%lf\n",quadW,otherQuadW);
//       }
//       if(quadW>0.)
// 	otherS = (x - c*curColSpacing) / quadW;
//       else
// 	otherS = 0.;
//       if(abs(otherS-s) > 0.0001){
// 	fprintf(stderr,"CHECK_OLD_WARP found discrepancy: s=%lf otherS=%lf\n",s,
// 		otherS);
//       }
//     }
// #endif
  }
  if(r==(numPointRows-2)){
    quadH = rgPoints0Y[idx+numPointCols] - rgPoints0Y[idx];
    if(quadH>0.)//changed from >=1. to >0. (after dissertation)
      t = (y - rgPoints0Y[idx]) / quadH;
    else
      t = 0.;
  }
#if DEBUG_WARP_POINT
  if((s > 1.)||(t>1.)){
    fprintf(stderr,"warpPoint(): bad s or t!!!   x,y=%lf,%lf  r=%d c=%d  idx=%d numPointsCols=%d  s=%lf t=%lf  quadW=%lf quadH=%lf (CLAMPING)\n",
	   x,y,r,c,idx,numPointCols,s,t,quadW,quadH);
  }
#endif
  if(s<0.)
    s=0.;
  if(t<0.)
    t=0.;
  if(s>1.)
    s=1.;
  if(t>1.)
    t=1.;
     
#if DEBUG_WARP_POINT
  if(fDebugPrint){
    printf("  quadW=%lf quadH=%lf\n",quadW, quadH);
    printf("  corrected s=%lf t=%lf\n",s,t);
  }
#endif

  //bilinearly interpolate the point x,y to get its warped coordinates xp,yp:
  double Xtop, Ytop, Xbot, Ybot;
  Xtop = (1.-s)*rgPoints1X[idx] + s*rgPoints1X[idx+1];
  Ytop = (1.-s)*rgPoints1Y[idx] + s*rgPoints1Y[idx+1];
  Xbot = (1.-s)*rgPoints1X[idx+numPointCols] + s*rgPoints1X[idx+numPointCols+1];
  Ybot = (1.-s)*rgPoints1Y[idx+numPointCols] + s*rgPoints1Y[idx+numPointCols+1];
  (*xp) = (1.-t)*Xtop + t*Xbot;
  (*yp) = (1.-t)*Ytop + t*Ybot;
#if DEBUG_WARP_POINT
  if(fDebugPrint){
    printf("  Xtop=%lf Ytop=%lf  Xbot=%lf Ybot=%lf  xp=%lf yp=%lf\n",
	   Xtop, Ytop, Xbot, Ybot, *xp, *yp);
  }
#endif
  return true;
}


//cost of a vertex position (using distance map of warped ink pixels in 4 quads)
//inlining this function sped us up about 10%
inline double DMorphInk::getVertexPositionCost(int r, int c,
					       double x1, double y1){
  double sumCost = 0.;
  int inkPxls = 0;//count of ink pixels in quads
  signed int *psDist;
  int idxVertex;

  idxVertex = r*numPointCols+c;//index in control point arrays of this vertex
  psDist = (signed int*)imgDist1.dataPointer_u32();

  for(int i=0; i < tempMA0len; ++i){
    double xp, yp;
    int ixp, iyp;
    double vxTmp, vyTmp; // store the vertex before changing it to set back
    
    vxTmp = rgPoints1X[idxVertex];
    vyTmp = rgPoints1Y[idxVertex];
    rgPoints1X[idxVertex] = x1;
    rgPoints1Y[idxVertex] = y1;
    warpPoint(rgTempMA0X[i], rgTempMA0Y[i], &xp, &yp);
    rgPoints1X[idxVertex] = vxTmp;
    rgPoints1Y[idxVertex] = vyTmp;
    ixp=(int)xp;
    iyp=(int)yp;

    int addDistX, addDistY; // if the position is off the distmap, compensate
    addDistX = addDistY = 0;
#if 1 /*COMPENSATE_FOR_OOB_DISTMAP*/
    if(ixp < 0){
      addDistX = 0-ixp;
      ixp = 0;
    }
    else if(ixp >  (w1-1)){
      addDistX = ixp-w1+1;
      ixp = w1-1;
    }
    if(iyp < 0){
      addDistY = 0-iyp;
      iyp = 0;
    }
    else if(iyp > (h1-1)){
      addDistY = iyp-h1+1;
      iyp = h1-1;
    }
#endif /*COMPENSATE_FOR_OOB_DISTMAP*/

    sumCost += psDist[w1*iyp+ixp] + addDistX + addDistY;
    ++inkPxls;
    
  }
  if(sumCost > 0)
    sumCost = sumCost/(inkPxls);
  //sumCost = sumCost/(1+inkPxls);

  return sumCost;
}







#if NEW_WARP

#include <execinfo.h>
inline void print_stack(){
  void **rgStackFrames;
  int numFrames;

  rgStackFrames= new void*[100];
  numFrames = backtrace(rgStackFrames, 100);
  char **rgStackFrameSymbols;
  
  rgStackFrameSymbols = backtrace_symbols(rgStackFrames, numFrames);
  fprintf(stderr,"stack trace:\n");
  for(int i=0; i < numFrames; ++i){
    fprintf(stderr,"  %s\n",rgStackFrameSymbols[i]);
  }
  free(rgStackFrameSymbols);
  delete [] rgStackFrames;
}


/**given point x, y in mesh0, find the warped coordinate at time=t [0. to 1.].*/
//instead of x,y use s,t, and idx
//this new version speeds us up by 33% from where we were when using the old warpPoint function
inline void DMorphInk::warpPointNew(double s, double t, int meshPointIdx,
				    double *xp, double *yp//,
				    /*debug:double x,double y,
				    int ii,int l,int rgTempMA0idx_i,
				    int r, int c*/){
  //FYI: meshPointIdx = r*numPointCols+c;

  //bilinearly interpolate the point x,y to get its warped coordinates xp,yp:
  double Xtop, Ytop, Xbot, Ybot;
  double oneMs, oneMt;//1-s and 1-t
  oneMs = 1.-s;
  oneMt = 1.-t;
  Xtop = (oneMs)*rgPoints1X[meshPointIdx] + s*rgPoints1X[meshPointIdx+1];
  Ytop = (oneMs)*rgPoints1Y[meshPointIdx] + s*rgPoints1Y[meshPointIdx+1];
  Xbot = (oneMs)*rgPoints1X[meshPointIdx+numPointCols] +
    s*rgPoints1X[meshPointIdx+numPointCols+1];
  Ybot = (oneMs)*rgPoints1Y[meshPointIdx+numPointCols] +
    s*rgPoints1Y[meshPointIdx+numPointCols+1];
  (*xp) = (oneMt)*Xtop + t*Xbot;
  (*yp) = (oneMt)*Ytop + t*Ybot;

#if 0
  {//debug
    double xp2, yp2;
    if(!warpPoint(x,y,&xp2,&yp2)){
      fprintf(stderr,"warpPointNew() debug check with warpPoint: warpPoint returned false! xp=%lf,yp=%lf  xp2=%lf,yp2=%lf\n",*xp,*yp,xp2,yp2);
      print_stack();
      exit(1);
    }
    if((xp2!=(*xp))||(yp2!=(*yp))){
      fprintf(stderr,"warpPointNew() debug check with warpPoint: xp=%lf,yp=%lf  xp2=%lf,yp2=%lf  s=%lf,t=%lf x=%lf,y=%lf ii=%d l=%d rgTempMA0idx_i=%d r=%d c=%d\n",
	      *xp,*yp,xp2,yp2,s,t,x,y,ii,l,rgTempMA0idx_i,r,c);
      fprintf(stderr,"numPointRows=%d numPointCols=%d\n",numPointRows,
	      numPointCols);
      print_stack();
      exit(1);
    }
  }
#endif

  //  return true;
}



//cost of a vertex position (using distance map of warped ink pixels in 4 quads)
//inlining this function speeds us up about 10%
/*
r: row in rgPoints1X,Y
c: columm in rgPoints1X,Y
x1, y1: the location of the current control point (in img1)
*/
inline double DMorphInk::getVertexPositionCostNew(int r, int c,
						  double x1, double y1){
  double sumCost = 0.;
  int inkPxls = 0;//count of ink pixels in quads
  signed int *psDist;
  int idxVertex;

  idxVertex = r*numPointCols+c;//index in control point arrays of this vertex
  psDist = (signed int*)imgDist1.dataPointer_u32();

  for(int i=0; i < tempMA0len; ++i){
    double xp, yp;
    int ixp, iyp;
    double vxTmp, vyTmp; // store the vertex before changing it to set back
    
    vxTmp = rgPoints1X[idxVertex];
    vyTmp = rgPoints1Y[idxVertex];
    rgPoints1X[idxVertex] = x1;
    rgPoints1Y[idxVertex] = y1;
    // warpPoint(rgTempMA0X[i], rgTempMA0Y[i], &xp, &yp);
    warpPointNew(rgMA0s[rgTempMA0idx[i]], rgMA0t[rgTempMA0idx[i]],
		 rgMA0r[rgTempMA0idx[i]]*numPointCols+rgMA0c[rgTempMA0idx[i]],
		 &xp, &yp/*debug: ,rgTempMA0X[i], rgTempMA0Y[i],i,lenMA0,
			   rgTempMA0idx[i],r,c*/);
    rgPoints1X[idxVertex] = vxTmp;
    rgPoints1Y[idxVertex] = vyTmp;
    ixp=(int)xp;
    iyp=(int)yp;

    int addDistX, addDistY; // if the position is off the distmap, compensate
    addDistX = addDistY = 0;
#if 1 /*COMPENSATE_FOR_OOB_DISTMAP*/
    if(ixp < 0){
      addDistX = 0-ixp;
      ixp = 0;
    }
    else if(ixp >  (w1-1)){
      addDistX = ixp-w1+1;
      ixp = w1-1;
    }
    if(iyp < 0){
      addDistY = 0-iyp;
      iyp = 0;
    }
    else if(iyp > (h1-1)){
      addDistY = iyp-h1+1;
      iyp = h1-1;
    }
#endif /*COMPENSATE_FOR_OOB_DISTMAP*/

    sumCost += psDist[w1*iyp+ixp] + addDistX + addDistY;
    ++inkPxls;
    
  }
  if(sumCost > 0)//sumCost==0 if inkPxls is zero, and some other times too
    sumCost = sumCost/(inkPxls);
  //sumCost = sumCost/(1+inkPxls);

#if 0
  {//debug
    if(sumCost != getVertexPositionCost(r,c,x1,y1)){
      fprintf(stderr,"getVertexPositionCostNew() returns different value than "
	      "getVertexPositionCost()\n");
      exit(1);
    }
  }
#endif

  return sumCost;
}
#endif //NEW_WARP



#endif
