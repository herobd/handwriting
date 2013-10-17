#ifndef DKERNEL2D_H
#define DKERNEL2D_H

#include "dimage.h"
#include <stdio.h>

class DKernel2D{
public:
  DKernel2D();
  DKernel2D(const DKernel2D &src);
  ~DKernel2D();
  const DKernel2D& operator=(const DKernel2D &src);
  const DKernel2D& operator=(const DImage &src);

  DImage toDImage();

  bool isInitialized();
  bool isSeparable();
  void setRect(int radiusX, int radiusY, bool fSeparable = true);
  void setCirc(int radiusX, int radiusY);
  void setGauss(int radiusX, int radiusY, bool fSeparable = true);
  void setLaplace();
  void setLoG(int radiusX, int radiusY, bool fSeparable = true);
  void setFakeLoG(int radiusX, int radiusY);
  void setData_flt(float *rgData, int w, int h, bool fIsSeparable);
  void setData_dbl(double *rgData, int w, int h, bool fIsSeparable);

  
  void scaleValues(double scaleBy);
  void setNumSigmas(double sigmasX, double sigmasY);
  float *getData_flt();
  double *getData_dbl();
  int getRadiusX();
  int getRadiusY();
  int getWidth();
  int getHeight();

  void print(FILE *fout = stdout);

private:
  void deleteBuffers();
  void allocBuffers(int w, int h);
  bool fSep; ///true iff this is a separable kernel
  int radiusX; /// # of pixels on each side of the central pixel
  int radiusY;
  float *rgData_flt;/// regular kernel values: (2*radiusX+1)*(2*radiusY+1)
  double *rgData_dbl;/// regular kernel values: (2*radiusX+1)*(2*radiusY+1)
  float *rgSep_flt; /// separable kernel values:(2*radiusX+1) + (2*radiusY+1)
  double *rgSep_dbl; /// separable kernel values:(2*radiusX+1) + (2*radiusY+1)
  double numSigsX;
  double numSigsY;
};


#endif

