#ifndef DDYNAMIC_PROGRAMMING_H
#define DDYNAMIC_PROGRAMMING_H

#include "dfeaturevector.h"
#include "dimage.h"

class DDynamicProgramming{
public:
  static double findDPAlignment(const DFeatureVector &fv1,
				const DFeatureVector &fv2,
				int bandRadius, double bandCost=1000.,
				double nonDiagonalCost = 0.,
				int *pathLen = NULL,
				int *rgPath = NULL,
				int *rgPrev = NULL,
				double *rgTable = NULL);
  static DImage piecewiseLinearWarpDImage(DImage &img1, int warpToLength,
					  int pathLen, int *rgPath,
					  bool fVertical);
  static DImage visualizeDPTable(int *rgPrev, int Wa, int Wb, int bandRadius=0);
  static void debugImages(int Wa, int Wb, int pathLen,
			  int *rgPath, int *rgPrev, double *rgTable);
  static void getCoord0MappingsToCoord1(int w0, int w1, double *rgCoords,
					int pathLen, int *rgPath);
  static void getCoord1MappingsToCoord0(int w1, double *rgCoords,
					int pathLen, int *rgPath);
  static double checkBandCost(int i,int j,int numRows,int numCols,
			      int bandRadius, double bandCost = 1000.);
};

///adds high cost to anything outside of Sakoe-Chiba band
inline double DDynamicProgramming::checkBandCost(int i,int j,int numRows,
						 int numCols, int bandRadius,
						 double bandCost){
  double iExpect, jExpect;

  iExpect = numRows * j / (double)numCols;
  if(i > (iExpect + bandRadius))
    return bandCost;
  jExpect = ((numRows-1)+numCols * i) /  (double)numRows;
  if(j > (jExpect + bandRadius))
    return bandCost;
  return 0.;
}



#endif
