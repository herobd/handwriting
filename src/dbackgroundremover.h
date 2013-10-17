#ifndef DBACKGROUNDREMOVER_H
#define DBACKGROUNDREMOVER_H

#include "dimage.h"

class DProgress; // forward declaration

class DBackgroundRemover{
public:
  static void removeBackground(DImage &imgDst, const DImage &imgSrc,
			       DImage *pImgBG = NULL, int radiusX = -1,
			       int radiusY = -1,DProgress *pProg = NULL,
			       int numThreads = 8);
  static void removeLowResBackground(DImage &imgDst, const DImage &imgSrc,
				     DImage *pImgBG = NULL, int radiusX = -1,
				     int radiusY = -1,DProgress *pProg = NULL,
				     int numThreads = 1,
				     int bgMinW = 100, int bgMinH = 100,
				     int bgMaxW = 800, int bgMaxH = 800);

};

#endif
