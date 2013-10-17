#ifndef DMEDIANFILTER_H
#define DMEDIANFILTER_H

#include "dimage.h"


class DProgress; // forward declaration

///This class provides median filter functionality for DImage objects
class DMedianFilter {
public:

  enum DMedFiltType{
    DMedFilt_default, ///<Choose whatever the default is for the image type
    DMedFilt_slow,///<Slow method
    DMedFilt_Huang_circle,///<A variation on Huang's method, circular kernel
    DMedFilt_Huang_square,///<A variation on Huang's method, square kernel
    DMedFilt_separable///<Separable approximation (not the true median)
  };

  DMedianFilter();
  DMedianFilter(int radiusX, int radiusY, DMedFiltType filtType);
  ~DMedianFilter();

  void setRadius(int radiusX, int radiusY);
  int getRadiusX();
  int getRadiusY();
  void setType(DMedFiltType filtType);
  DMedFiltType getType();
  void setNumThreads(int numThreads);
  int getNumThreads();
  
  DImage filterImage(const DImage &imgSrc, bool fAlreadyPadded = false,
		     DProgress *pProg = NULL);
  void filterImage_(DImage &imgDst,
		    const DImage &imgSrc, bool fAlreadyPadded = false,
		    DProgress *pProg = NULL);

  static void medianFilter_separable_u8(DImage &imgDst, const DImage &imgSrc,
					 bool fAlreadyPadded, int radius);


  static void medianFilterImage(DImage &imgDst, const DImage &imgSrc,
				bool fAlreadyPadded = false,
				int radiusX = 1, int radiusY = 1,
				DMedFiltType filtType = DMedFilt_default,
				DProgress *pProg = NULL, int numThreads = 1);

private:
  DMedFiltType _medFiltType;
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
  static void medFiltHuang_u8(DImage &imgDst, const DImage &imgSrc,
			      int radiusX, int radiusY,
			      int wKern, int hKern, D_uint8 *rgKern,
			      int numKernPxls, int *rgRightEdge,
			      DProgress *pProg = NULL,
			      int progStart = 0, int progMax = 1,
			      int threadNumber = 0, int numThreads = 1);
  static void medFiltHuang_u8_square(DImage &imgDst, const DImage &imgSrc,
				     int radiusX, int radiusY,
				     int wKern, int hKern, D_uint8 *rgKern,
				     int numKernPxls,
				     DProgress *pProg = NULL,
				     int progStart = 0, int progMax = 1,
				     int threadNumber = 0, int numThreads = 1);
  static void* DMedianFilter_Huang8threadWrap(void* params);


  /// copy constructor is private so nobody can use it
  DMedianFilter(const DMedianFilter &src);
};



#endif
