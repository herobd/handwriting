#ifndef DMAXFILTER_H
#define DMAXFILTER_H

#include "dimage.h"


class DProgress; // forward declaration

///This class provides max filter functionality for DImage objects
class DMaxFilter {
public:

  enum DMaxFiltType{
    DMaxFilt_slow,///<Slow method
    DMaxFilt_circle,///<Histogram (based on Huang) method, circular kernel
    DMaxFilt_square///<Histogram (based on Huang) method, square kernel
  };

  DMaxFilter();
  DMaxFilter(int radiusX, int radiusY, DMaxFiltType filtType);
  ~DMaxFilter();

  void setRadius(int radiusX, int radiusY);
  int getRadiusX();
  int getRadiusY();
  void setType(DMaxFiltType filtType);
  DMaxFiltType getType();
  void setNumThreads(int numThreads);
  int getNumThreads();
  
  DImage filterImage(const DImage &imgSrc, bool fAlreadyPadded = false,
		     DProgress *pProg = NULL);
  void filterImage_(DImage &imgDst,
		    const DImage &imgSrc, bool fAlreadyPadded = false,
		    DProgress *pProg = NULL);

  static void maxFilterImage(DImage &imgDst, const DImage &imgSrc,
				bool fAlreadyPadded = false,
				int radiusX = 1, int radiusY = 1,
				DMaxFiltType filtType = DMaxFilt_circle,
				DProgress *pProg = NULL, int numThreads = 1);

private:
  DMaxFiltType _maxFiltType;
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
  static void maxFiltHuang_u8(DImage &imgDst, const DImage &imgSrc,
			      int radiusX, int radiusY,
			      int wKern, int hKern, D_uint8 *rgKern,
			      int numKernPxls, int *rgRightEdge,
			      DProgress *pProg = NULL,
			      int progStart = 0, int progMax = 1,
			      int threadNumber = 0, int numThreads = 1);
  static void maxFiltHuang_u8_square(DImage &imgDst, const DImage &imgSrc,
				     int radiusX, int radiusY,
				     int wKern, int hKern, D_uint8 *rgKern,
				     int numKernPxls,
				     DProgress *pProg = NULL,
				     int progStart = 0, int progMax = 1,
				     int threadNumber = 0, int numThreads = 1);
  static void* DMaxFilter_Huang8threadWrap(void* params);

  /// copy constructor is private so nobody can use it
  DMaxFilter(const DMaxFilter &src);
};



#endif
