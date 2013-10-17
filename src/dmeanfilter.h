#ifndef DMEANFILTER_H
#define DMEANFILTER_H

#include "dimage.h"


class DProgress; // forward declaration

///This class provides mean filter functionality for DImage objects
class DMeanFilter {
public:

  enum DMeanFiltType{
    DMeanFilt_circle,///Circular/elliptical kernel
    DMeanFilt_square///square/rectangle kernel
  };

  DMeanFilter();
  ~DMeanFilter();

  void setRadii(int radiusX, int radiusY);
  void setType(DMeanFiltType filtType);
  DMeanFiltType getType();
  
  DImage filterImage(const DImage &imgSrc, bool fAlreadyPadded = false,
		     DProgress *pProg = NULL);
  void filterImage_(DImage &imgDst,
		    const DImage &imgSrc, bool fAlreadyPadded = false,
		    DProgress *pProg = NULL);

  static void meanFilterImage(DImage &imgDst, const DImage &imgSrc,
			      bool fAlreadyPadded = false,
			      int radiusX = 1, int radiusY = 1,
			      DMeanFiltType filtType = DMeanFilt_circle,
			      DProgress *pProg = NULL);
  static void meanFilterUseIntegralImage_(DImage &imgDst, const DImage &imgSrc,
					  int radiusX, int radiusY,
					  bool fAlreadyPadded = false);

private:
  DMeanFiltType _meanFiltType;
  int _radiusX, _radiusY;
  D_uint8 *rgKern;
  int *rgRightEdge;

  static void fill_square_kern_offsets(int radiusX, int radiusY,
				       unsigned char *rgKernFull,
				       int *rgRightEdge, int *numKernPxls);
  static void fill_circle_kern_offsets(int radiusX, int radiusY,
				       unsigned char *rgKernFull,
				       int *rgRightEdge, int *numKernPxls);
  static void meanFilt_u8(DImage &imgDst, const DImage &imgSrc,
			  int radiusX, int radiusY, D_uint8 *rgKern,
			  int numKernPxls, int *rgRightEdge,
			  DProgress *pProg = NULL,
			  int progStart = 0, int progMax = 1);
  /// copy constructor is private so nobody can use it
  DMeanFilter(const DMeanFilter &src);
};



#endif
