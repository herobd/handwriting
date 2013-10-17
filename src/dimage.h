#ifndef DIMAGE_H
#define DIMAGE_H

#include <stdio.h>
#include <stdlib.h>
#include "ddefs.h"
#include <complex>
#include <map>
#include <vector>
#include <string>
#include "dsize.h"

/// This class is the internal storage structure for images of various data formats
/** The DImage class encapsulates storage of image data within memory.
 * Various image types are supported (unsigned char, float, complex,
 * etc.)
 *
 * Conversion between data types is supported
 * when it seems to make sense.  Basic manipulation (such as scaling,
 * rotation, drawing lines into the image, etc.) is also supported.
 *
 * Image data can be directly accessed and manipulated via the pointers
 * returned by the dataPointer_<type> functions.
 *
 * Most manipulation functions that return a DImage object also have a
 * version of the function with the same name, followed by an
 * underscore ("_") character.  These functions use the first
 * parameter as the destination image instead of creating a new image
 * object.  This eliminates an extra copy that otherwise must be done
 * to return a DImage object when the destination image is not passed
 * (by reference) as a parameter.
 */



class DImage {
  friend class DImageIO;
public:
  
  enum DImageType{
    DImage_u8 = 1,///< unsigned 8-bit (usually grayscale)
    DImage_u16,///< unsigned 16-bit (usually high dynamic range GS)
    DImage_u32,///< unsigned 32-bit (usually a component map, etc.)
    DImage_RGB,///< unsigned 8-bit RGB
    DImage_RGB_16,///< unsigned 16-bit per channel RGB)
    DImage_flt_multi, ///< float data, may be multi-channel (non-interleaved)
    DImage_dbl_multi, ///< double data, may be multi-channel (non-interleaved)
    DImage_cmplx ///< std::complex<double> (as defined in <complex>)
  };

  enum DFileFormat{
    DFileFormat_pnm = 0, ///< PNM format corresponding to the image
    DFileFormat_pbm_plain, ///< PBM plain-text format (bitonal images only)
    DFileFormat_pgm_plain, ///< PGM plain-text (color is converted)
    DFileFormat_ppm_plain, ///< PPM plain-text (BW/grayscale goes to RGB)
    DFileFormat_pbm, ///< PBM format forced (for bitonal images only)
    DFileFormat_pgm, ///< PGM format forced (color is converted)
    DFileFormat_ppm, ///< PPM format forced (BW/grayscale goes to RGB)
    DFileFormat_jpg, ///< JPEG
    DFileFormat_png, ///< PNG
    DFileFormat_tiff, ///< TIFF file
    DFileFormat_gif, ///< GIF format
    DFileFormat_bmp, ///< Windows Bitmap (.bmp) format
    DFileFormat_gnuplot, ///< gnuplot splot format
    DFileFormat_unknown///<a format that is not recognized
  };

  enum DImageTransformMode{
    DImageTransSmooth = 1,///<interpolate pixel values (doesn't preserve vals)
    DImageTransSample///<sample pixel values (jagged, but fast; preserves vals)
  };

  enum CompositeMode{ // TODO: add other composite modes as needed
    CompositeMode_Source = 1,
    CompositeMode_SourceOver
  };

  enum DImagePadMethod{
    DImagePadValue = 0,///<Fill pad area with specified value
    DImagePadReplicate,///<Replicate the pixels on the image edge into the pad
    DImagePadMirror,///<Mirror the pixels across the edge
    DImagePadWrap,///<Wrap around to the values on opposite side of image
    DImagePadNone///<Leave the pad area uninitialized
  };

  DImage();
  DImage(const DImage &src, bool fCopyDataBuffer = true);
  DImage(const char *stPath);
  ~DImage();
  const DImage& operator=(const DImage &src); // assignment operator

  bool load(const char *stPath, DFileFormat fmt=DFileFormat_unknown);
  bool save(const char *stPath, DFileFormat fmt = DFileFormat_pnm,
	    bool fSaveProps=true, bool fSaveComments = true);
  bool create(int w, int h, DImageType imgType,int numChannels = -1,
	      D_AllocationMethod allocMeth = AllocationMethod_malloc);
  void fill(double val);
  void fill(int R, int G, int B);
  void fill(std::complex<double> complexVal);
  void clear();
  void setDataBuffer(void *pBuf, D_AllocationMethod allocMeth);
  void* releaseDataBuffer();
  void deallocateBuffer(); // use appropriate deallocation for pData
  

  DImage copy(int x, int y, int w, int h) const;
  void copy_(DImage &imgDst, int x, int y, int w, int h,
	     D_AllocationMethod allocMeth = AllocationMethod_src) const;
  void pasteFromImage(int dstX, int dstY,
		      const DImage &imgSrc, int srcX, int srcY,int w,int h,
		      DImage::CompositeMode cmode = CompositeMode_Source);
  
  void setDataRange(double min, double max,
		    float perctMin=0., float perctMax=0.);
  void capDataRange(double min, double max);
  void setDataRangeLog(double base=10.);
  void getDataRange(double *min, double *max);
  void addValueToPixels(int val);
  void addValueToPixels(double val);
  void dividePixelsByValue(double val);
  void multiplyPixelsByValue(double val);

  DImage multiplyByImage(const DImage imgOther);
  void multiplyByImage_(DImage &imgDst, const DImage imgOther);

  DImage addToImage(const DImage imgOther);
  void addToImage_(DImage &imgDst, const DImage imgOther);

  void drawPixel(int x, int y, int GS, float transparency=0.);
  void drawPixel(int x, int y, int R, int G, int B, float transparency=0.);



  void drawLine(int x0, int y0, int x1, int y1, int GS, float transparency=0.,
		bool fIgnoreClipWarning=false);
  void drawLine(int x0, int y0, int x1, int y1, int R, int G, int B,
		float transparency = 0., bool fIgnoreClipWarning=false);
  void drawText(const char *st, int x, int y, int pointSize = 32,
		int gs_val=0, double transparency=0.,
		const char *stFont = NULL, int *w=NULL, int *h=NULL,
		double offsPctX=0., double offsPctY=0.);
  void drawTextRGB(const char *st, int x, int y, int pointSize = 32,
		   int R=0, int G=0, int B=0, double transparency=0.,
		   const char *stFont = NULL, int *w=NULL, int *h=NULL,
		   double offsPctX=0., double offsPctY=0.);
  void getTextSize(const char *st, int pointSize = 32,
		   const char *stFont = NULL, int *w=NULL, int *h=NULL);
  void drawText_unicode(std::string str, int x, int y, int pointSize = 32,
			int gs_val=0, double transparency=0.,
			const char *stFont = NULL, int *w=NULL, int *h=NULL,
			double offsPctX=0., double offsPctY=0.);
  void drawText_unicode(std::string str, int x, int y, int pointSize = 32,
			int R=0, int G=0, int B=0, double transparency=0.,
			const char *stFont = NULL, int *w=NULL, int *h=NULL,
			double offsPctX=0., double offsPctY=0.);
  
  

#if 0
  void drawPixelAntiAliased(double x, double y, int GS, float transparency=0.);
  void drawPixelAntiAliased(double x, double y, int R, int G, int B,
			    float transparency=0.);
  void drawLineAntiAliased(double x0, double y0, double x1, double y1, int GS,
			   float transparency=0.,
			   bool fIgnoreClipWarning=false);
  void drawLineAntiAliased(double x0, double y0, double x1, double y1,
			   int R, int G, int B,
			   float transparency = 0.,
			   bool fIgnoreClipWarning=false);
#endif
  void drawRect(int x0, int y0, int x1, int y1, int GS,
		float transparency = 0., bool fill=false);
  void drawRect(int x0, int y0, int x1, int y1, int R, int G, int B,
		float transparency = 0., bool fill=false);

  DImage resized(int newW, int newH, double pad_val) const;
  DImage resized(int newW, int newH,
		 int pad_R, int pad_G, int pad_B) const;
  void resized_(DImage &imgDst, int newW, int newH, double pad_val) const;
  void resized_(DImage &imgDst, int newW, int newH,
	       int pad_R, int pad_G, int pad_B) const;

  DImage padEdges(int left, int right, int top, int bottom,
		  DImagePadMethod padMethod = DImagePadValue,
		  double val=0.) const;
  void padEdges_(DImage &imgDst, int left, int right, int top, int bottom,
		 DImagePadMethod padMethod = DImagePadValue,
		 double val=0.) const;

  void invertGrayscale();

  void splitRGB(DImage &imgR, DImage &imgG, DImage &imgB) const;
  void combineRGB(DImage &imgR, DImage &imgG, DImage &imgB);


  void setLogicalSize(int w, int h);
  
  DImage scaled(int dstW, int dstH,
		DImageTransformMode mode = DImageTransSmooth) const;
  DImage scaled(double scaleX, double scaleY,
		DImageTransformMode mode = DImageTransSmooth) const;

  void scaled_(DImage &imgDst, int dstW, int dstH,
	       DImageTransformMode mode = DImageTransSmooth) const;
  void scaled_(DImage &imgDst, double scaleX, double scaleY,
	       DImageTransformMode mode = DImageTransSmooth) const;
  void scaledDownPow2_(DImage &imgDst, int numHalves,
		       DImageTransformMode mode = DImageTransSmooth) const;
  void scaledUpPow2_(DImage &imgDst, int numPows, int dstW=-1, int dstH=-1,
		     DImageTransformMode mode = DImageTransSample) const;
  int scaledDownNumHalves(int minW=100, int minH=100,
			  int maxW=800, int maxH=800) const;

  DImage translated(double tx, double ty, double pad_val,
		    DImageTransformMode mode = DImageTransSmooth) const;
  DImage translated(double tx, double ty, int pad_R, int pad_G, int pad_B,
		    DImageTransformMode mode = DImageTransSmooth) const;
  void translated_(DImage &imgDst, double tx, double ty, double pad_val,
		   DImageTransformMode mode = DImageTransSmooth) const;
  void translated_(DImage &imgDst, double tx, double ty,
		   int pad_R, int pad_G, int pad_B,
		   DImageTransformMode mode = DImageTransSmooth) const;

  
  DImage rotated(double ang, double pad_val, bool fResize,
		 DImageTransformMode mode = DImageTransSmooth) const;
  DImage rotated(double ang, int pad_R, int pad_G, int pad_B,
		 bool fResize,
		 DImageTransformMode mode = DImageTransSmooth) const;
  void rotated_(DImage &imgDst, double ang, double pad_val, bool fResize,
		DImageTransformMode mode = DImageTransSmooth) const;
  void rotated_(DImage &imgDst, double ang, int pad_R, int pad_G, int pad_B,
		bool fResize,
		DImageTransformMode mode = DImageTransSmooth) const;
  DSize getRotatedSize(double ang) const;
  void rotate90_(DImage &imgDst, int ang) const;
  
  DImage shearedH(double ang, double pad_val, bool fResize,
		  DImageTransformMode mode = DImageTransSmooth) const;
  DImage shearedH(double ang, int pad_R, int pad_G, int pad_B,
		  bool fResize,
		  DImageTransformMode mode = DImageTransSmooth) const;
  void shearedH_(DImage &imgDst, double ang, double pad_val, bool fResize,
		 DImageTransformMode mode = DImageTransSmooth) const;
  void shearedH_(DImage &imgDst, double ang, int pad_R, int pad_G, int pad_B,
		 bool fResize,
		 DImageTransformMode mode = DImageTransSmooth) const;
  int getShearedHWidth(double ang) const;
  

  DImage shearedV(double ang, double pad_val,
		  bool fResize,
		  DImageTransformMode mode = DImageTransSmooth) const;
  DImage shearedV(double ang, int pad_R, int pad_G, int pad_B,
		  bool fResize,
		  DImageTransformMode mode = DImageTransSmooth) const;
  void shearedV_(DImage &imgDst, double ang, double pad_val,
		 bool fResize,
		 DImageTransformMode mode = DImageTransSmooth) const;
  void shearedV_(DImage &imgDst, double ang, int pad_R, int pad_G, int pad_B,
		 bool fResize,
		 DImageTransformMode mode = DImageTransSmooth) const;
  int getShearedVHeight(double ang) const;

  void setPixel(int x, int y, double val, int chan=0);
  void setPixel(int x, int y, int Rval, int Gval, int Bval);

  DImage convertedImgType(DImageType imgType, int numDstChannels = -1,
			  D_uint32 srcChannelMask = 0xffffffff,
			  D_AllocationMethod allocMeth =
			  AllocationMethod_src) const;
  void convertedImgType_(DImage &imgDst, DImageType imgType,
			 int numDstChannels =-1,
			 D_uint32 srcChannelMask = 0xffffffff,
			 D_AllocationMethod allocMeth =
			 AllocationMethod_src) const;

  void setProperty(const std::string propName, std::string propVal);
  std::string getPropertyVal(std::string propName);
  std::string getPropertyValByIndex(unsigned int propNum);
  int getNumProperties() const;
  void clearProperties();

  void addComment(std::string stComment);
  std::string getCommentByIndex(unsigned int idx) const;
  int getNumComments() const;
  void clearComments();

  void copyProperties(DImage &src);
  void copyComments(DImage &src);

  // accessor methods
  int width() const {return (NULL==pData)?0:_w;}///<return width of the image
  int height() const {return (NULL==pData)?0:_h;}///<return height of the image
  ///returns the number of channels in the image
  int numChannels() const {return (NULL==pData)?0:_numChan;}
  ///returns true if data is interleaved (ex: RGBRGBRGB instead of RRRGGGBBB)
  bool channelsAreInterleaved() const {return _fInterleaved;}
  ///returns the type of image this is. (DImage_u8 if image is still undefined)
  DImageType getImageType() const {return _imgType;}
  ///returns the memory allocation method used for the data of this image
  D_AllocationMethod getAllocMethod() const {return _allocMethod;}

  static void setAlignment(int alignment);
  static int getAlignment();
  static DFileFormat getImageFileFormat(const char *stPath);
  static DFileFormat getImageFileFormat(FILE *fin);
  static DImage blendImages(DImage img1, DImage img2, double t=0.5);

  // for the following functions, using channel is an error unless the image
  // is in DImageColor_MULTI color mode (separate, noninterleaved channels)
  D_uint8* dataPointer_u8() const;
  D_uint16* dataPointer_u16() const;
  D_uint32* dataPointer_u32() const;
  float* dataPointer_flt(int channel = 0) const;
  double* dataPointer_dbl(int channel = 0) const;
  std::complex<double>* dataPointer_cmplx() const;
    
private:
  int _w; // width 
  int _h; // height
  int _actualW; // the real allocated width (not the logical width) of image
  int _actualH; // the real allocated height (not the logical height) of image
  int _numChan; // number of channels
  DImageType _imgType; // type of image
//   DImageDataType _dtype; // base type of data
//   DImageColorType _colorType; // color format
  bool _fInterleaved; // true iff channels are interleaved
  D_AllocationMethod _allocMethod; // see enum in ddefs.h
  D_uint8 *pData; // pointer to data buffer (needs cast to proper type for use)
  size_t _dataSize; // size of data buffer (may be incorrect if pData==NULL)
  size_t _sampleSize; // sizeof(data type)
  static int _data_alignment; // alignment boundary for data allocation
  static int _drawTextFileNum; // used to create temporary files in drawText
  // properties (saved as comments)
  std::map<const std::string, std::string> mapProps;
  std::vector<std::string> vectComments;

  // _internal_DImage_assign is called by the copy constructor,
  // assignment operator, and a few other member functions
  static void _internal_DImage_assign(DImage &dst, const DImage &src,
				      bool fCopyData, bool fCopyProps,
				      bool fCopyComments);
};


inline void DImage::setPixel(int x, int y, double val, int chan){
  int idx;
  if((x >= _w) || (y >= _h) || (x < 0) || (y < 0)){
    fprintf(stderr, "DImage::setPixel(%d,%d) out of bounds w=%d h=%d\n",
	    x, y, _w, _h);
    return;
    // abort();
  }
  if((chan >= _numChan) && (DImage_cmplx != _imgType)){
    fprintf(stderr, "DImage::setPixel() chan=%d, image only has %d channels\n",
	    chan, _numChan);
    abort();
  }
  idx = chan * _w * _h + _w * y + x;
  switch(_imgType){
    case DImage_u8:
      {
	pData[idx] = (D_uint8)val;
      }
      break;
    case DImage_u16:
      {
	((D_uint16*)pData)[idx] = (D_uint16)val;
      }
      break;
    case DImage_u32:
	((D_uint32*)pData)[idx] = (D_uint32)val;
      break;
    case DImage_RGB:
    case DImage_RGB_16:
      fprintf(stderr,
	      "DImage::setPixel() wrong function called for RGB image\n");
      abort();
      break;
    case DImage_flt_multi:
	((float*)pData)[idx] = (float)val;
      break;
    case DImage_dbl_multi:
	((double*)pData)[idx] = (double)val;
      break;
    case DImage_cmplx:
      if((chan < 0) || (chan > 1)){
	fprintf(stderr, "DImage::setPixel() chan=%d. Must be 0(real part) or "
		"1(imaginary part)\n", chan);
	abort();
      }
      {
	std::complex<double> *pcmplx;
	
	idx = _w * y + x;
	pcmplx = (std::complex<double> *)pData;
	if(0 == chan){
	  pcmplx[idx] = val;
	}
	else{
	  pcmplx[idx] = std::complex<double>(0.,val);
	}
      }
	
      break;
    default:
      fprintf(stderr, "DImage::setPixel() bas case value\n");
      abort();
  }

}

inline void DImage::setPixel(int x, int y, int Rval, int Gval, int Bval){
  int idx;
  if((x >= _w) || (y >= _h) || (x < 0) || (y < 0)){
    fprintf(stderr, "DImage::setPixel(%d,%d) out of bounds w=%d h=%d\n",
	    x, y, _w, _h);
    return;
    // abort();
  }
  idx = 3 * (_w * y + x);
  if(_imgType == DImage_RGB){
    pData[idx] = (D_uint8)Rval;
    pData[idx+1] = (D_uint8)Gval;
    pData[idx+2] = (D_uint8)Bval;
  }
  else if(_imgType == DImage_RGB_16){
    D_uint16 *p16;
    p16 = (D_uint16*)pData;
    p16[idx] = (D_uint16)Rval;
    p16[idx+1] = (D_uint16)Gval;
    p16[idx+2] = (D_uint16)Bval;
  }
  else{
    fprintf(stderr, "DImage::setPixel() RGB function called for"
	    " non-RGB image\n");
    abort();
  }
}


inline D_uint8* DImage::dataPointer_u8() const{
#ifdef DEBUG
  if((DImage_u8 != _imgType) && (DImage_RGB != _imgType))
    fprintf(stderr, "DImage::dataPointer_u8() on wrong type of image %d\n",
	    (int)_imgType);
#endif
  return (D_uint8*)pData;
}
inline D_uint16* DImage::dataPointer_u16() const{
#ifdef DEBUG
  if((DImage_u16 != _imgType) && (DImage_RGB_16 != _imgType))
    fprintf(stderr, "DImage::dataPointer_u16() on wrong type of image\n");
#endif
  return (D_uint16*)pData;
}
inline D_uint32* DImage::dataPointer_u32() const{
#ifdef DEBUG
  if(DImage_u32 != _imgType)
    fprintf(stderr, "DImage::dataPointer_u32() on wrong type of image\n");
#endif
  return (D_uint32*)pData;
}
inline float* DImage::dataPointer_flt(int channel) const{
#ifdef DEBUG
  if(DImage_flt_multi != _imgType)
    fprintf(stderr, "DImage::dataPointer_flt() on wrong type of image\n");
  if((channel < 0) || (channel >= _numChan))
    fprintf(stderr, "DImage::dataPointer_flt() on invalid channel\n");
#endif
  return (float*)&(pData[_w*_h*channel*sizeof(float)]);
}
inline double* DImage::dataPointer_dbl(int channel) const{
#ifdef DEBUG
  if(DImage_dbl_multi != _imgType)
    fprintf(stderr, "DImage::dataPointer_dbl() on wrong type of image\n");
  if((channel < 0) || (channel >= _numChan))
    fprintf(stderr, "DImage::dataPointer_dbl() on invalid channel\n");
#endif
  return (double*)&(pData[_w*_h*channel*sizeof(double)]);
}
inline std::complex<double>* DImage::dataPointer_cmplx() const{
  return (std::complex<double>*)pData;
}


#endif
