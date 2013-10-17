#include "dkernel2d.h"
#include "dmemalign.h"
#include "dimage.h"
#include "dinstancecounter.h"
#include <string.h>

#define DEFAULT_NUM_SIGMAS 3.0

///Default Constructor
DKernel2D::DKernel2D(){
  DInstanceCounter::addInstance("DKernel2D");
  fSep = false;
  radiusX = 0;
  radiusY = 0;
  rgData_flt = NULL;
  rgData_dbl = NULL;
  rgSep_flt = NULL;
  rgSep_dbl = NULL;
  numSigsX = DEFAULT_NUM_SIGMAS;
  numSigsY = DEFAULT_NUM_SIGMAS;
}

///Destructor
DKernel2D::~DKernel2D(){
  DInstanceCounter::removeInstance("DKernel2D");
  deleteBuffers();
}

///Copy Construstor (makes a deep copy of the data)
DKernel2D::DKernel2D(const DKernel2D &src){
  int w, h;
  DInstanceCounter::addInstance("DKernel2D");
  fSep = src.fSep;
  radiusX = src.radiusX;
  radiusY = src.radiusY;
  numSigsX = src.numSigsX;
  numSigsY = src.numSigsY;
  w = 2*radiusX+1;
  h = 2*radiusY+1;
  if((w < 1) || (h < 1)){
    rgData_flt = NULL;
    rgData_dbl = NULL;
    rgSep_flt = NULL;
    rgSep_dbl = NULL;
  }
  else{
    allocBuffers(w, h);
    memcpy(rgData_flt, src.rgData_flt, sizeof(float)*w*h);
    memcpy(rgData_dbl, src.rgData_dbl, sizeof(double)*w*h);
    memcpy(rgSep_flt, src.rgSep_flt, sizeof(float)*(w+h));
    memcpy(rgSep_dbl, src.rgSep_dbl, sizeof(double)*(w+h));
  }
}

///Assignment operator (deep-copies the src kernel)
const DKernel2D& DKernel2D::operator=(const DKernel2D &src){
  int w, h;
  if(this != &src){
    fSep = src.fSep;
    radiusX = src.radiusX;
    radiusY = src.radiusY;
    numSigsX = src.numSigsX;
    numSigsY = src.numSigsY;
    w = 2*radiusX+1;
    h = 2*radiusY+1;
    deleteBuffers();
    if((w < 1) || (h < 1)){
      rgData_flt = NULL;
      rgData_dbl = NULL;
      rgSep_flt = NULL;
      rgSep_dbl = NULL;
    }
    else{
      allocBuffers(w, h);
      memcpy(rgData_flt, src.rgData_flt, sizeof(float)*w*h);
      memcpy(rgData_dbl, src.rgData_dbl, sizeof(double)*w*h);
      memcpy(rgSep_flt, src.rgSep_flt, sizeof(float)*(w+h));
      memcpy(rgSep_dbl, src.rgSep_dbl, sizeof(double)*(w+h));
    }
  }
  return *this;
}

///Assignment from a DImage object (creates a kernel based on srcImg contents)
const DKernel2D& DKernel2D::operator=(const DImage &srcImg){
  int w, h;
  DImage imgDbl;
  double *pDataDbl;
  double *pDataDst;
  float *pDataFlt;

  if((0 == (srcImg.width() & 1)) || (0 == (srcImg.height() & 1))){
    fprintf(stderr, "DKernel2D::operator=() source image w,h must be odd\n");
    exit(1);
  }
  if(1 != srcImg.numChannels()){
    fprintf(stderr, "DKernel2D::operator=() source image numChannels != 1\n");
    exit(1);
  }
  if(srcImg.getImageType() == DImage::DImage_cmplx){
    fprintf(stderr, "DKernel2D::operator=() DImage_cmplx not supported\n");
    exit(1);
  }
  srcImg.convertedImgType_(imgDbl, DImage::DImage_dbl_multi, 1);
  pDataDbl = imgDbl.dataPointer_dbl();
  fSep = false;
  radiusX = srcImg.width() / 2;
  radiusY = srcImg.height() / 2;
  numSigsX = DEFAULT_NUM_SIGMAS;
  numSigsY = DEFAULT_NUM_SIGMAS;
  w = 2*radiusX+1;
  h = 2*radiusY+1;
  deleteBuffers();
  allocBuffers(w, h);
  pDataDst = rgData_dbl;
  pDataFlt = rgData_flt;
  
  for(int y = 0; y < h; ++y){
    memcpy(pDataDst, pDataDbl, sizeof(double) * w);
    for(int x = 0; x < w; ++x)
      pDataFlt[x] = (float)(pDataDst[x]);
    pDataDst += w;
    pDataDbl += w;
    pDataFlt += w;
  }

  return *this;
}

///Returns a DImage object (of type DImage_dbl_multi) with the kernel contents
/**Non-separable kernels are formatted as one might expect.  Separable kernels
 * show the two separable functions in a "+" pattern in the image.
 */
DImage DKernel2D::toDImage(){
  int w, h;
  DImage dst;
  double *pDst;

  w = 2*radiusX+1;
  h = 2*radiusY+1;

  dst.create(w, h, DImage::DImage_dbl_multi, 1);
  pDst = dst.dataPointer_dbl(0);

  if(!fSep){ // non-separable: just copy the data to the image
    memcpy(pDst, rgData_dbl, sizeof(double)*w*h);
  }
  else{ // separable: make a "+" with the two separable 1-d kernels
    memset(pDst, 0, sizeof(double)*w*h);
    for(int y = 0; y < h; ++y){
      pDst[y*w+w/2] = rgSep_dbl[w+y];
    }
    pDst += h/2*w;
    for(int x = 0; x < w; ++x){
      pDst[x] = rgSep_dbl[x];
    }
  }

  return dst;
}

///used internally by DKernel2D
void DKernel2D::deleteBuffers(){
  if(rgData_flt != NULL){
    daligned_free((void*)rgData_flt);
    rgData_flt = NULL;
  }
  if(rgData_dbl != NULL){
    daligned_free((void*)rgData_dbl);
    rgData_dbl = NULL;
  }
  if(rgSep_flt != NULL){
    daligned_free((void*)rgSep_flt);
    rgSep_flt = NULL;
  }
  if(rgSep_dbl != NULL){
    daligned_free((void*)rgSep_dbl);
    rgSep_dbl = NULL;
  }
}

///used internally by DKernel2D
void DKernel2D::allocBuffers(int w, int h){
  rgData_flt = (float*)daligned_malloc(sizeof(float)*w*h, 16);
  if(!rgData_flt){
    fprintf(stderr, "DKernel2D::allocBuffers() no memory\n");
    exit(1);
  }
  rgData_dbl = (double*)daligned_malloc(sizeof(double)*w*h, 16);
  if(!rgData_flt){
    fprintf(stderr, "DKernel2D::allocBuffers() no memory\n");
    exit(1);
  }
  rgSep_flt = (float*)daligned_malloc(sizeof(float)*(w+h), 16);
  if(!rgSep_flt){
    fprintf(stderr, "DKernel2D::allocBuffers() no memory\n");
    exit(1);
  }
  rgSep_dbl = (double*)daligned_malloc(sizeof(double)*(w+h), 16);
  if(!rgSep_dbl){
    fprintf(stderr, "DKernel2D::allocBuffers() no memory\n");
    exit(1);
  }
}

///Returns true iff the kernel is a separable kernel, false otherwise
/**The result tells you whether the kernel IS currently a separable
   kernel, not whether the kernel COULD BE a separable kernel.  (i.e.,
   while a separable gaussian kernel can be created, that doesn't mean
   that any particular gausiian kernel has been constructed as a
   separable kernel.  Only those that are created to be separable can
   be convolved separably.) */
bool DKernel2D::isSeparable(){
  return fSep;
}

///Returns true iff the kernel is initialized with data, false otherwise
bool DKernel2D::isInitialized(){
  return ((rgData_flt != NULL) || (rgData_dbl != NULL) ||
	  (rgSep_flt != NULL) || (rgSep_dbl != NULL));
}

///Scales the values of each position in the kernel by scaleBy
/**Multiplies each kernel value by scaleBy */
void DKernel2D::scaleValues(double scaleBy){
  int len;
  int w, h;
  float fltScaleBy;
  
  w = radiusX*2+1;
  h = radiusY*2+1;
  len = w*h;
  fltScaleBy = (float)scaleBy;

  if(rgData_flt){
    for(int i = 0; i < len; ++i)
      rgData_flt[i] *= fltScaleBy;
  }
  if(rgData_dbl){
    for(int i = 0; i < len; ++i)
      rgData_dbl[i] *= scaleBy;
  }
  if(rgSep_flt){
    for(int i = 0; i < (w+h); ++i)
      rgSep_flt[i] *= fltScaleBy;
  }
  if(rgSep_dbl){
    for(int i = 0; i < (w+h); ++i)
      rgSep_dbl[i] *= scaleBy;
  }
}

///Set how many sigmas to include within kernel area for gaussian-based kernels
/**By default, when creation methods such as setGuass() and setLoG()
 * are called, the kernel is created such that 3.0 standard deviations
 * fit within the kernel created with a given radiusX and radiusY.  If
 * you want more or less standard deviations to fit within the kernel,
 * you can adjust how many here.  Note that sigmasX and sigmasY are
 * NOT the actual sigma (standard deviation) of the gaussian function,
 * but rather the NUMBER of sigmas that will fit, so sigmaX = (radiusX
 * / sigmasX).
 */
void DKernel2D::setNumSigmas(double sigmasX, double sigmasY){
  numSigsX = sigmasX;
  numSigsY = sigmasY;
}

///get a pointer to float data for this kernel
/** The kernel data is stored redundantly so that either float or
 *  doubles can be used from the same kernel.  For non-separable
 *  kernels, the pointer points to a buffer that is big enough for
 *  getWidth()*getHeight() floats.  For separable kernels, the pointer
 *  points to a buffer that is big enough for getWidth()+getHeight()
 *  floats (add instead of multiply).
 */
float* DKernel2D::getData_flt(){
  if(fSep)
    return rgSep_flt;
  return rgData_flt;
}

///get a pointer to double data for this kernel
/** The kernel data is stored redundantly so that either float or
 *  doubles can be used from the same kernel.  For non-separable
 *  kernels, the pointer points to a buffer that is big enough for
 *  getWidth()*getHeight() doubles.  For separable kernels, the pointer
 *  points to a buffer that is big enough for getWidth()+getHeight()
 *  doubles (add instead of multiply).
 */
double* DKernel2D::getData_dbl(){
  if(fSep)
    return rgSep_dbl;
  return rgData_dbl;
}

///Returns the kernel radius along the x-axis (kernel width=radiusX*2+1)
int DKernel2D::getRadiusX(){
  return radiusX;
}

///Returns the kernel radius along the y-axis (kernel height=radiusY*2+1)
int DKernel2D::getRadiusY(){
  return radiusY;
}
///Returns the width of the kernel, which is radiusX*2+1
/**Note that the kernel width and height are always odd because the
 * radius does not include the central pixel.
 */
int DKernel2D::getWidth(){
  return radiusX*2+1;
}
///Returns the width of the kernel, which is radiusX*2+1
/**Note that the kernel width and height are always odd because the
 * radius does not include the central pixel.
 */
int DKernel2D::getHeight(){
  return radiusY*2+1;
}

///Set the kernel to a rectangular (averaging) filter kernel
/**radiusX and radiusY specify how far the kernel extends from the
   central pixel.  fSeparable specifies whether this should be a
   separable kernel, which allows dconvolver to perform convolution
   faster than standard non-separable convolution.  The value of each
   pixel of the kernel is scaled equally so that the entire filter
   sums to 1.
 */
void DKernel2D::setRect(int radiusX, int radiusY, bool fSeparable){
  int w, h;
  float fltVal;
  double dblVal;
  w = radiusX*2+1;
  h = radiusY*2+1;
  this->radiusX = radiusX;
  this->radiusY = radiusY;
  fSep = fSeparable;
  deleteBuffers();
  allocBuffers(w,h);


  if(fSeparable){
    dblVal = 1./w;
    fltVal = (float)dblVal;
    for(int x = 0; x < w; ++x){
      rgSep_flt[x] = fltVal;
      rgSep_dbl[x] = dblVal;
    }
    dblVal = 1./h;
    fltVal = (float)dblVal;
    for(int y = 0; y < h; ++y){
      rgSep_flt[w+y] = fltVal;
      rgSep_dbl[w+y] = dblVal;
    }
  }
  else{
    dblVal = 1./(w*h);
    fltVal = (float)dblVal;
    for(int y = 0, idx = 0; y < h; ++y){
      for(int x = 0; x < w; ++x, ++idx){
	rgData_flt[idx] = fltVal;
	rgData_dbl[idx] = dblVal;
      }
    }
  }
}
///Set the kernel to a circular (actually elliptical) averaging filter kernel
/**radiusX and radiusY specify how far the kernel extends from the
   central pixel.  Circular kernels cannot be separable.  The value of
   each pixel of the kernel within the elliptical area is scaled
   equally so that the entire filter sums to 1.  Areas outside the
   ellipse are set to zero.
 */
void DKernel2D::setCirc(int radiusX, int radiusY){
  int w, h;
  double dblSum;
  double rx2,ry2;
  w = radiusX*2+1;
  h = radiusY*2+1;
  this->radiusX = radiusX;
  this->radiusY = radiusY;
  fSep = false;
  deleteBuffers();
  allocBuffers(w,h);

  dblSum = 0.;
  rx2=(0.5+radiusX)*(0.5+radiusX);
  ry2=(0.5+radiusY)*(0.5+radiusY);
  for(int y = -radiusY, idx = 0; y <= radiusY; ++y){
    for(int x = -radiusX; x <= radiusX; ++x, ++idx){
      if(((x*x/rx2)+(y*y/ry2)) <= 1.){
	rgData_dbl[idx] = 1.;
	dblSum += 1.;
      }
      else
	rgData_dbl[idx] = 0.;
    }
  }
  for(int y = -radiusY, idx = 0; y <= radiusY; ++y){
    for(int x = -radiusX; x <= radiusX; ++x, ++idx){
      rgData_dbl[idx] /= dblSum;
      rgData_flt[idx] = rgData_dbl[idx];
    }
  }
}
///Set the kernel to an axis-aligned Gaussian filter kernel
/**fSeparable should be set to true if you want a separable kernel,
 * which is faster than a non-separable kernel when convolving.  (for
 * an NxN separable kernel, convolution is a 2*N operation instead of
 * N*N operation for each image pixel).  By default, 3 standard
 * deviations (sigmas) fit under the gaussian curve on each side of
 * the central pixel.  If you need more (or fewer) standard deviations
 * to fit within the kernel, use setNumSigmas().
*/
void DKernel2D::setGauss(int radiusX, int radiusY, bool fSeparable){
  int w, h;
  double xp, yp;
  double one_2piSigXsigY;
  double coef;
  double coefX, coefY;
  double sigmaX, sigmaY;
  double dblSum;

  w = radiusX*2+1;
  h = radiusY*2+1;
  this->radiusX = radiusX;
  this->radiusY = radiusY;
  fSep = fSeparable;
  deleteBuffers();
  allocBuffers(w,h);
  sigmaX = radiusX / numSigsX;
  sigmaY = radiusY / numSigsY;

  coefX = -0.5 / sigmaX / sigmaX;
  coefY = -0.5 / sigmaY / sigmaY;
  if(fSeparable){
    one_2piSigXsigY = 1. / (2* M_PI * sigmaX * sigmaY);
    dblSum = 0.;
    for(int x = 0; x < w; ++x){
      xp = x-radiusX;
      rgSep_dbl[x] = one_2piSigXsigY * exp(xp*xp * coefX);
      dblSum += rgSep_dbl[x];
    }
    for(int x = 0; x < w; ++x){
      rgSep_dbl[x] /= dblSum;
      rgSep_flt[x] = (float)(rgSep_dbl[x]);
    }
    dblSum = 0.;
    for(int y = 0; y < h; ++y){
      yp = y-radiusY;
      rgSep_dbl[w+y] = one_2piSigXsigY * exp(yp*yp * coefY);
      dblSum += rgSep_dbl[w+y];
    }
    for(int y = 0; y < h; ++y){
      rgSep_dbl[w+y] /= dblSum;
      rgSep_flt[w+y] = (float)(rgSep_dbl[w+y]);
    }
  }
  else{
    one_2piSigXsigY = 1. / (2* M_PI * sigmaX * sigmaY);
    dblSum = 0.;
    for(int y = 0, idx = 0; y < h; ++y){
      yp = y-radiusY;
      for(int x = 0; x < w; ++x, ++idx){
	xp = x-radiusX;
	rgData_dbl[idx] =
	  one_2piSigXsigY * exp(coefX*xp*xp + coefY*yp*yp);
	dblSum += rgData_dbl[idx];
      }
    }
    for(int idx = 0; idx < w*h; ++idx){
      rgData_dbl[idx] /= dblSum;
      rgData_flt[idx] = (float)(rgData_dbl[idx]);
    }
  }
}

///Set the kernel to a 3x3 Laplace (-8 in the center, 1 everywhere else)
void DKernel2D::setLaplace(){
  int w, h;
  w = 3;
  h = 3;
  this->radiusX = 1;
  this->radiusY = 1;
  fSep = false;
  deleteBuffers();
  allocBuffers(w,h);
  for(int idx = 0; idx < 9; ++idx){
    rgData_dbl[idx] = 1.;
    rgData_flt[idx] = 1.;
  }
  rgData_dbl[4] = -8.;
  rgData_flt[4] = -8.;
}

///Set the kernel to a Laplacian of Gaussian (2nd derivative of the Gaussian)
void DKernel2D::setLoG(int radiusX, int radiusY, bool fSeparable){
  int w, h;
  double xp, yp;
  double one_2piSigXsigY;
  double sigmaX, sigmaY;
  double sigmaX2, sigmaY2;
  double sigmaX4, sigmaY4;
  double ksum = 0.;

  sigmaX = radiusX / numSigsX;
  sigmaY = radiusY / numSigsY;
  sigmaX2 = sigmaX * sigmaX;
  sigmaY2 = sigmaY * sigmaY;
  sigmaX4 = sigmaX2 * sigmaX2;
  sigmaY4 = sigmaY2 * sigmaY2;

  w = radiusX*2+1;
  h = radiusY*2+1;
  this->radiusX = radiusX;
  this->radiusY = radiusY;
  fSep = fSeparable;
  deleteBuffers();
  allocBuffers(w,h);
  one_2piSigXsigY = 1. / (2* M_PI * sigmaX * sigmaY);

  if(fSeparable){
    fprintf(stderr, "DKernel2D::setLoG() separable kernel not implemented\n");
    exit(1);
  }
  else{
    one_2piSigXsigY = 1. / (2* M_PI * sigmaX * sigmaY);
    for(int y = 0, idx = 0; y < h; ++y){
      yp = y-radiusY;
      for(int x = 0; x < w; ++x, ++idx){
	xp = x-radiusX;
	rgData_dbl[idx] =
	  one_2piSigXsigY * exp(-0.5*(xp*xp/sigmaX2 + yp*yp/sigmaY2)) *
	  ((xp*xp/sigmaX4 - 1./sigmaX2) + (yp*yp/sigmaY4 - 1./sigmaY2));
	rgData_flt[idx] = (float)(rgData_dbl[idx]);
	ksum += rgData_dbl[idx];
      }
    }
  }
  printf("LoG ksum=%f\n", ksum);
//   // get rid of the bias (ksum)
//   ksum /= (w*h);
//   for(int i = 0; i < w*h; ++i){
//     rgData_dbl[i] -= ksum;
//     rgData_flt[i] -= ksum;
//   }
//   ksum = 0.;
//   for(int i = 0; i < w*h; ++i){
//     ksum += rgData_dbl[i];
//   }
//   printf("after adjusting, ksum=%f\n", ksum);

}

///Set the kernel to a fake LaplacianOfGaussian
void DKernel2D::setFakeLoG(int radiusX, int radiusY){
  int w, h;
  double xp, yp;
  double one_2piSigXsigY;
  double sigmaX, sigmaY;
  double sigmaX2, sigmaY2;
  double sigmaX4, sigmaY4;
  double XoverY;
  double ksum = 0.;

  XoverY = (double)radiusX / (double)radiusY;

  sigmaX = radiusX / numSigsX;
  sigmaY = sigmaX;
//   sigmaY = radiusY / numSigsY;
  sigmaX2 = sigmaX * sigmaX;
  sigmaY2 = sigmaY * sigmaY;
  sigmaX4 = sigmaX2 * sigmaX2;
  sigmaY4 = sigmaY2 * sigmaY2;

  w = radiusX*2+1;
  h = radiusY*2+1;
  this->radiusX = radiusX;
  this->radiusY = radiusY;
  fSep = false;
  deleteBuffers();
  allocBuffers(w,h);
  one_2piSigXsigY = 1. / (2* M_PI * sigmaX * sigmaY);

  one_2piSigXsigY = 1. / (2* M_PI * sigmaX * sigmaY);
  for(int y = 0, idx = 0; y < h; ++y){
    yp = (y-radiusY)*XoverY;
    for(int x = 0; x < w; ++x, ++idx){
      xp = x-radiusX;
      rgData_dbl[idx] =
	one_2piSigXsigY * exp(-0.5*(xp*xp/sigmaX2 + yp*yp/sigmaY2)) *
	((xp*xp/sigmaX4 - 1./sigmaX2) + (yp*yp/sigmaY4 - 1./sigmaY2));
      rgData_flt[idx] = (float)(rgData_dbl[idx]);
      ksum += rgData_dbl[idx];
    }
  }
  printf("Fake LoG ksum=%f\n", ksum);
//   // get rid of the bias (ksum)
//   ksum /= (w*h);
//   for(int i = 0; i < w*h; ++i){
//     rgData_dbl[i] -= ksum;
//     rgData_flt[i] -= ksum;
//   }
//   ksum = 0.;
//   for(int i = 0; i < w*h; ++i){
//     ksum += rgData_dbl[i];
//   }
//   printf("after adjusting, ksum=%f\n", ksum);
}

///set the kernel data to the values specified in rgData
/**The internal data arrays are allocated and the values from rgData
 * are copied into the internal data arrays.  The size of the buffers
 * depends on whether fIsSeparable is true or false.  If fIsSeparable
 * is false, then the buffers hold w*h floats/doubles (multiply).  If
 * fIsSeparable is true, then the buffers only hold w+h floats/doubles
 * (add).  Not all kernels can be made separable, but if you use separable
 * kernels when possible, convolution is typically much faster.  When
 * calling either setData function, the internal arrays for both float and
 * double values are allocated and initialized to the proper values.  Both
 * functions are provided as a convenience to the caller so the data being
 * copied can be either float or double.
 *
 * For separable kernels, the horizontal component of the kernel is
 * first (w values), followed by the vertical component (h values).
 */
void DKernel2D::setData_flt(float *rgData, int w, int h, bool fIsSeparable){
  if((0 == (w & 1)) || (0 == (h & 1))){
    fprintf(stderr, "DKernel2D::setData_flt() w,h must be odd\n");
    exit(1);
  }
  this->radiusX = w/2;
  this->radiusY = h/2;
  this->fSep = fIsSeparable;
  deleteBuffers();
  allocBuffers(w,h);
  
  for(int y =0, idx=0; y < h; ++y){
    for(int x = 0; x < w; ++x, ++idx){
      rgData_flt[idx] = rgData[idx];
      rgData_dbl[idx] = (double)(rgData[idx]);
    }
  }
}
///set the kernel data to the values specified in rgData
/**The internal data arrays are allocated and the values from rgData
 * are copied into the internal data arrays.  The size of the buffers
 * depends on whether fIsSeparable is true or false.  If fIsSeparable
 * is false, then the buffers hold w*h floats/doubles (multiply).  If
 * fIsSeparable is true, then the buffers only hold w+h floats/doubles
 * (add).  Not all kernels can be made separable, but if you use separable
 * kernels when possible, convolution is typically much faster.  When
 * calling either setData function, the internal arrays for both float and
 * double values are allocated and initialized to the proper values.  Both
 * functions are provided as a convenience to the caller so the data being
 * copied can be either float or double.
 *
 * For separable kernels, the horizontal component of the kernel is
 * first (w values), followed by the vertical component (h values).
 */
void DKernel2D::setData_dbl(double *rgData, int w, int h, bool fIsSeparable){
  if((0 == (w & 1)) || (0 == (h & 1))){
    fprintf(stderr, "DKernel2D::setData_dbl() w,h must be odd\n");
    exit(1);
  }
  this->radiusX = w/2;
  this->radiusY = h/2;
  this->fSep = fIsSeparable;
  deleteBuffers();
  allocBuffers(w,h);
  for(int y =0, idx=0; y < h; ++y){
    for(int x = 0; x < w; ++x, ++idx){
      rgData_flt[idx] = (float)(rgData[idx]);
      rgData_dbl[idx] = rgData[idx];
    }
  }
}


///Print the contents of the kernel (to 3 decimal places) to stdout (or a file)
void DKernel2D::print(FILE *fout){
  int w, h;
  double sum = 0.;

  if(NULL == rgData_dbl){
    fprintf(fout, "kernel not initialized\n");
    return;
  }


  w = radiusX*2+1;
  h = radiusY*2+1;
  if(fSep){
    for(int x = 0; x < w; ++x){
      fprintf(fout,"%0.03f ", rgSep_dbl[x]);
      sum += rgSep_dbl[x];
    }
    printf("\n");
    for(int y = 0; y < h; ++y){
      fprintf(fout,"   %0.03f\n", rgSep_dbl[w+y]);
      sum += rgSep_dbl[y];
    }
  }
  else{
    for(int y = 0; y < h; ++y){
      for(int x = 0; x < w; ++x){
	fprintf(fout,"%0.03f ", rgData_dbl[y*w+x]);
	sum += rgData_dbl[y*w+x];
      }
      fprintf(fout,"\n");
    }
  }
  fprintf(fout, "sum of kernel values=%.20f\n", sum);
}
