#ifndef DEDGEDETECTOR_H
#define DEDGEDETECTOR_H

#include "dimage.h"

class DEdgeDetector{
  public:
  
  static void sobel_(DImage &imgGradMag, const DImage &imgSrc,
		    DImage *imgDir=NULL);
  static void sobelHorizEdges_(DImage &imgGradMag, const DImage &imgSrc);
  static void sobelVertEdges_(DImage &imgGradMag, const DImage &imgSrc);
  static void sobelHorizEdgesKeepSign_(DImage &imgGrad, const DImage &imgSrc);
  static void sobelVertEdgesKeepSign_(DImage &imgGrad, const DImage &imgSrc);
  static void sobel_RGB_(DImage &imgGradMag, const DImage &imgSrc,
			 DImage *imgDir=NULL);
  static void sobel_RGB_combined_(DImage &imgGradMag, const DImage &imgSrc,
				  DImage *imgDir=NULL);

  static void prewittHorizEdges_(DImage &imgGradMag, const DImage &imgSrc);
  static void prewittVertEdges_(DImage &imgGradMag, const DImage &imgSrc);
  static void prewittHorizEdgesKeepSign_(DImage &imgGrad,const DImage &imgSrc);
  static void prewittVertEdgesKeepSign_(DImage &imgGrad, const DImage &imgSrc);

  static void longPrewittHorizEdgesKeepSign_(DImage &imgGrad,
					     const DImage &imgSrc,
					     int radius = 2);
  static void longPrewittVertEdgesKeepSign_(DImage &imgGrad,
					    const DImage &imgSrc,
					    int radius = 2);


//   static void laplace3x3_(DImage &imgGradMag, DImage &imgSrc);
  static void laplace3x3Zeros_(DImage &imgGradMag, DImage &imgSrc,
			       bool fNearestOnly = true);

  static void convertGradDirsToEdgeDirs(DImage &imgDirs,
					bool fQuads1and2 = true);

  static void nonMaximalSuppression_(DImage &imgDst, DImage &imgSrc,
				     DImage &imgDirs);

  static void showDirsInRGBImage_(DImage &imgDst, DImage &imgGradMag,
				  DImage &imgDirs, int thresh = 20);
//   static void convertRadianDirsToDegs(DImage &imgDirs);
//   static void limitDirsToQuadrants1and2(DImage &imgDirs);
//   static void limitDirsToQuadrants1and4(DImage &imgDirs);

//   static void cannyEdgeDetect();///<do canny edge detection
//   static void edgeDetectRGB();///<detect edges in each channel, result is RGB
//   static void edgeDetectMagRGB(bool fUseCIELab = true);///< use the color distance in computing the edge strength, result is grayscale and gradient direction
};

#endif
