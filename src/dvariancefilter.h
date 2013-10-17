#ifndef DVARIANCEFILTER_H
#define DVARIANCEFILTER_H

#include "dimage.h"

/**This class provides a convenient mechanism for calculating the
 * local variance (or standard deviation) of an image.  For each
 * pixel, the variance (or standard deviation) is calculated within
 * the surrounding local window of size (2*radiusX+1) by (2*radiusY+1)
 * pixels.  If radiusX or radiusY is -1, then the size will default to
 * max(3*imageWidth/200, 3*imageHeight/200).  If fAlreadyPadded is
 * false, then the image is padded by replicating the edge pixels
 * before running the filter.  The destination image, imgDst, will be
 * of type DImage_dbl_multi, and will have the same number of channels
 * as imgSrc.  If imgSrc was an RGB image, the channels will be split
 * in imgDst (all R values in channel 0, all green in 1, all blue in
 * 2, so RRRRRR...GGGGGGG...BBBBBBB... instead of interleaved
 * RGBRGBRGBRGBRGB...).
 *
 * This implementation uses  integral images for speed, as described in
 * "Efficient Implementation of Local Adaptive Thresholding Techniques
 * Using Integral Images" by Shafait, Keysers, and Breuel that was to
 * appear in Proc. Document Recognition and Retrieval XV, IST/SPIE
 * Annual Symposium, San Jose, CA, January 2008. It uses more memory
 * than a convolution approach since the integral image must be
 * created and uses long longs or long doubles to prevent overflow.*/


class DVarianceFilter{
public:

  
  static void varianceImage_(DImage &imgDst, const DImage &imgSrc,
			     bool fAlreadyPadded = false,
			     int radiusX = -1, int radiusY = -1);
  static void standardDeviationImage_(DImage &imgDst, const DImage &imgSrc,
				      bool fAlreadyPadded = false,
				      int radiusX = -1, int radiusY = -1);

};


#endif
