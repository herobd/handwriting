#ifndef DGLOBALSKEW_H
#define DGLOBALSKEW_H

#include "dimage.h"

/**This class provides functionality for detecting the global skew of
   a document image. */

class DGlobalSkew{
public:
  static double getSkewAng_var(const DImage &img, double angFirstDeg = -45.,
			       double angLastDeg = 45.,
			       double angDeltaDeg = .1);
  static double getSkewAng_fast(const DImage &img, double angFirstDeg = -45.,
				double angLastDeg = 45.,
 				double angDeltaDeg = .1,
				double angDeltaDegLowRes = .5,
				bool fUseEdges = true);
//   static double getSkewAng_fast_3(DImage &img, double angFirstDeg = -45.,
// 				double angLastDeg = 45.,
//  				double angDeltaDeg = .1,
// 				double angDeltaDegLowRes = .5,
// 				bool fUseEdges = true);

  static double getSkewAng_postl(const DImage &img, double angFirstDeg = -45.,
				double angLastDeg = 45.,
 				double angDeltaDeg = .1,
				double angDeltaDegLowRes = .5,
				bool fUseEdges = true);

};


#endif
