#ifndef DMORPHOLOGY_H
#define DMORPHOLOGY_H

#include "dimage.h"


class DMorphology{
public:
  // static DStructuringElement SEL_RECT3x3; // default structuring element
  static void dilate3x3_(DImage &imgDst, DImage &imgSrc, bool fBlackIsFG=true);
  static void erode3x3_(DImage &imgDst, DImage &imgSrc, bool fBlackIsFG=true);
};

#endif
