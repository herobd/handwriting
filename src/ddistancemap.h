#ifndef DDISTANCEMAP_H
#define DDISTANCEMAP_H

#include "dimage.h"

class DDistanceMap{
public:

  static DImage getDistFromInkBitonal(DImage &src, int maxDist=240,
				      int maxNegDist=-12);
  static void getDistFromInkBitonal_(DImage &dst, DImage &src, int maxDist=240,
				     int maxNegDist=-12);

};
inline DImage DDistanceMap::getDistFromInkBitonal(DImage &src, int maxDist,
					   int maxNegDist){
  DImage imgResult;
  getDistFromInkBitonal_(imgResult, src, maxDist, maxNegDist);
  return imgResult;
}

#endif
