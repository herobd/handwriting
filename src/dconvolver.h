#ifndef DCONVOLVER_H
#define DCONVOLVER_H

#include "dimage.h"
#include "dkernel2d.h"

class DProgress; // forward declaration

class DConvolver{
  public:


  static void convolve_(DImage &imgDst, const DImage &imgSrc, DKernel2D &kern, 
			bool fSrcAlreadyPadded=false,
			bool fResize=false, bool fConvertBack=false,
			bool fUseDoublePrec=false, DProgress *pProg=NULL);

private:
  static void convSepFlt(DImage &imgDst, const DImage &imgSrc, DKernel2D &kern,
			 bool fResize, DProgress *pProg);
  static void convSepFltOld(DImage &imgDst, const DImage &imgSrc, DKernel2D &kern,
			 bool fResize, DProgress *pProg);
  static void convSepDbl(DImage &imgDst, const DImage &imgSrc, DKernel2D &kern,
			 bool fResize, DProgress *pProg);
  static void convRegFlt(DImage &imgDst, const DImage &imgSrc, DKernel2D &kern,
			 bool fResize, DProgress *pProg);
  static void convRegDbl(DImage &imgDst, const DImage &imgSrc, DKernel2D &kern,
			 bool fResize, DProgress *pProg);
};


#endif

