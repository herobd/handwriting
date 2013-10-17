#ifndef DSLANTANGLE_H
#define DSLANTANGLE_H

#include "dimage.h"
#include "drect.h"

class DSlantAngle{
public:

  static double getTextlineSlantAngleDeg(DImage &imgBW, int rlThresh,
					 int x0=0,int y0=0,int x1=-1,int y1=-1,
					 double *weight=NULL,
					 unsigned int *rgSlantHist=NULL,
					 DImage *imgAngleHist=NULL);
  static double getAllTextlinesSlantAngleDeg(DImage &imgBW, int rlThresh,
   					     DRect *rgTextlineRects,
					     int numRects, int numThreads=12,
					     double *rgAngles=NULL,
					     double *rgWeights=NULL);
private:
  static void* getSlant_thread_func(void *params);
};



#endif
