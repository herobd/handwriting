#ifndef DTHRESHOLDER_H
#define DTHRESHOLDER_H

#include "dimage.h"

class DProgress; // forward declaration

///This class provides functionality for thresholding images
/**Unless stated otherwise, image pixels EQUAL to a threshold value
   are set to zero when performing global thresholding, just
   like the pixels less than the threshold value.  Only pixels greater
   than the threshold are set to 255.*/

class DThresholder{
public:

  static void threshImage_(DImage &imgDst, const DImage &imgSrc, double tval);
  static void threshImageSpecial_(DImage &imgDst, const DImage &imgSrc,
				  double t1, double t2,
				  double lowVal = 0., double midVal=-1.,
				  double highVal = 255.);
  static void otsuThreshImage_(DImage &imgDst, const DImage &imgSrc);
  static void niblackThreshImage_(DImage &imgDst, const DImage &imgSrc,
				  int radius = 7, double K=-0.2,
				  DProgress *pProg = NULL, int numThreads = 1);
  static void sauvolaNiblackThreshImage_(DImage &imgDst, const DImage &imgSrc,
					 int window = 15, double R=128.,
					 double K=0.5,
					 DProgress *pProg = NULL,
					 int numThreads = 1);
  static void sauvolaThreshImage_(DImage &imgDst, const DImage &imgSrc,
				  int window = 15, double R=128.,
				  double K=0.5,
				  DProgress *pProg = NULL, int numThreads = 1);
  static void milewskiThreshImage_precollaboration(DImage &imgDst, const DImage &imgSrc,
				   int strokeThickness=5, int P = -1,
				   int kappa=10, double sineConst = 0.5,
				   DProgress *pProg=NULL, int numThreads = 1);
  static void milewskiThreshImage_(DImage &imgDst, const DImage &imgSrc,
				   int strokeThickness=5, int P = -1,
				   int kappa=10, double sineConst = 0.5,
				   int voteThresh = 4,
				   DProgress *pProg=NULL, int numThreads = 1);
  static void milewskiThreshImage3_(DImage &imgDst, const DImage &imgSrc,
				   int strokeThickness=5, int P = -1,
				   int kappa=10, double sineConst = 0.5,
				   int voteThresh = 4,
				   DProgress *pProg=NULL, int numThreads = 1);
  static void ccThreshImage_(DImage &imgDst, const DImage &imgSrc,
			     DProgress *pProg = NULL, int numThreads = 12,
			     int *tval=NULL);
			     
  
  static double getOtsuThreshVal(const DImage &img);
  static double getOtsuThreshVal(unsigned int *pHist, int histLen);
  static double getOtsuThreshVal2(unsigned int *pHist, int histLen);
  static double getOtsuThreshVal3(unsigned int *pHist, int histLen);
  static double getOtsuThreshVal(unsigned char *pImgData, int w, int h);
  static double getCCThreshVal(const DImage &img, DProgress *pProg = NULL,
			       int numThreads = 1);
  static double getCCThreshVal2(const DImage &img, DProgress *pProg = NULL,
			       int numThreads = 1);

private:
  static void* getCCThreshVal_thread_func(void *params);
};

#endif
