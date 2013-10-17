#ifndef DCONNECTEDCOMPINFO_H
#define DCONNECTEDCOMPINFO_H

#include "ddefs.h"

class DConnectedComponentInfo{
public:
  DConnectedComponentInfo();
  ~DConnectedComponentInfo();
  void clear();

  D_uint32 label;///<CC label (background is label 0)
  int bbLeft;///<bounding box left
  int bbTop;///<bounding box Top
  int bbRight;///<bounding box right
  int bbBottom;///<bounding box bottom
  int startX;///<x-coord of first pixel (traversing image in row-major order)
  int startY;///<x-coord of first pixel (traversing image in row-major order)
  int pixels;///<Number of pixels belonging to the CC
  int centroidX;
  int centroidY;
};

inline void DConnectedComponentInfo::clear(){
  pixels = 0;
  bbLeft = bbTop = bbRight = bbBottom = startX = startY = -1;
  centroidX = centroidY = -1;
  label = 0xffffffff;
}

#endif
