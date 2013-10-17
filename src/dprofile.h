#ifndef DPROFILE_H
#define DPROFILE_H

#include "dimage.h"

class DProfile{
public:
  DProfile();
  DProfile(const DProfile &src);
  ~DProfile();
//   const DProfile& operator=(const DProfile &src); // assignment operator

  void getImageProfile(const DImage &img, double axisAngDeg = 0.,
		       bool fNormalize = true, bool fInterp = false);
  void getImageVerticalProfile(const DImage &img, bool fNormalize = true);
  void getImageHorizontalProfile(const DImage &img, bool fNormalize = true);
  void getAngledVertProfile(const DImage &img, double ang,
			    int fNormalize = true);
  void getVertMaxRunlengthProfile(const DImage &img,
				  D_uint32 rgbVal/*0xffffff*/,
				  bool fNormalize = true);
  void getHorizMaxRunlengthProfile(const DImage &img,
				   D_uint32 rgbVal/*=0xffffff*/,
				   bool fNormalize = true);
  void getVertAvgRunlengthProfile(const DImage &img,
				  D_uint32 rgbVal/*=0xffffff*/,
				  bool fNormalize = true);
  void getHorizAvgRunlengthProfile(const DImage &img,
				   D_uint32 rgbVal/*=0xffffff*/,
				   bool fNormalize = true);

  void copyFromData(double *rgData, int dataLen);
  void copyFromData(float *rgData, int dataLen);
  void copyFromData(int *rgData, int dataLen);
  void smoothAvg(int windowRadius);
//   void smoothGauss(double sigma, double numSigmas=2.);
//   void scaleMaxToVal(double maxVal=100.);

  void saveGnuplot(const char *stPath);
  DImage toDImage(int numPixels=100, bool fVertical=false, D_uint8 fg=0xff,
		  D_uint8 bg=0x80);
  
  double* dataPointer(); // return rgData
  double max();
  double min();
  int dataLen(); // return len
protected:
  double *rgProf;
  int len;
};

inline double* DProfile::dataPointer(){return rgProf;}
inline int DProfile::dataLen(){return len;}

#endif
