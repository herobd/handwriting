#ifndef DCONNECTEDCOMPDISTANCE_H
#define DCONNECTEDCOMPDISTANCE_H

#include "dimage.h"

class DConnectedComponentDistance{
public:
  static void getCCDist_(DImage &imgCCDist, const DImage &imgCCs);
  static void getNearestCC_(DImage &imgNearestCC, const DImage &imgCCs);
  static void getCCDistAndNearestCC_(DImage &imgCCDist, DImage &imgNearestCC,
				     const DImage &imgCCs);
};



#endif
