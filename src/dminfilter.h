#ifndef DMINFILTER_H
#define DMINFILTER_H

#include "dimage.h"


class DProgress; // forward declaration

///This class provides min filter functionality for DImage objects
class DMinFilter {
public:

  enum DMinFiltType{
    DMinFilt_slow,///<Slow method
    DMinFilt_circle,///<Histogram (based on Huang) method, circular kernel
    DMinFilt_square///<Histogram (based on Huang) method, square kernel
  };

  DMinFilter();
  DMinFilter(int radiusX, int radiusY, DMinFiltType filtType);
  ~DMinFilter();

  void setRadius(int radiusX, int radiusY);
  int getRadiusX();
  int getRadiusY();
  void setType(DMinFiltType filtType);
  DMinFiltType getType();
  void setNumThreads(int numThreads);
  int getNumThreads();
  
  DImage filterImage(const DImage &imgSrc, bool fAlreadyPadded = false,
		     DProgress *pProg = NULL);
  void filterImage_(DImage &imgDst,
		    const DImage &imgSrc, bool fAlreadyPadded = false,
		    DProgress *pProg = NULL);

  static void minFilterImage(DImage &imgDst, const DImage &imgSrc,
				bool fAlreadyPadded = false,
				int radiusX = 1, int radiusY = 1,
				DMinFiltType filtType = DMinFilt_circle,
				DProgress *pProg = NULL, int numThreads = 1);

private:
  DMinFiltType _minFiltType;
  int _radiusX;
  int _radiusY;
  D_uint8 *rgKern;
  int *rgRightEdge;
  int _numThreads;

  static void fill_square_kern_offsets(int radiusX, int radiusY,
				       unsigned char *rgKernFull,
				       int *rgRightEdge, int *numKernPxls);
  static void fill_circle_kern_offsets(int radiusX, int radiusY,
				       unsigned char *rgKernFull,
				       int *rgRightEdge, int *numKernPxls);
  static void minFiltHuang_u8(DImage &imgDst, const DImage &imgSrc,
			      int radiusX, int radiusY,
			      int wKern, int hKern, D_uint8 *rgKern,
			      int numKernPxls, int *rgRightEdge,
			      DProgress *pProg = NULL,
			      int progStart = 0, int progMax = 1,
			      int threadNumber = 0, int numThreads = 1);
  static void minFiltHuang_u8_square(DImage &imgDst, const DImage &imgSrc,
				     int radiusX, int radiusY,
				     int wKern, int hKern, D_uint8 *rgKern,
				     int numKernPxls,
				     DProgress *pProg = NULL,
				     int progStart = 0, int progMax = 1,
				     int threadNumber = 0, int numThreads = 1);
  static void* DMinFilter_Huang8threadWrap(void* params);

  /// copy constructor is private so nobody can use it
  DMinFilter(const DMinFilter &src);
};



#endif
