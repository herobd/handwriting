#ifndef DMEDIALAXIS_H
#define DMEDIALAXIS_H

#include "dimage.h"

class DMedialAxis{
 public:
  static DImage getMedialAxisImageFromDistMap(DImage &imgDistMap,
					      bool fThin=true,
					      int *numPoints=NULL);
  static DImage getMedialAxisImage(DImage &img, bool fZeroIsInk=true,
				   int thresholdVal=127);
  static DImage colorizeMedialAxisImage(DImage &imgMA, int r, int g, int b,
					int bgR=255, int bgG=255, int bgB=255);
};


#endif
