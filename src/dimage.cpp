#include "dmemalign.h"
#include "dimage.h"
#include "dinstancecounter.h"
#include "dimageio.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <complex>
#include <math.h>
#include "drect.h"
int DImage::_data_alignment = 16; // by default, 8-byte aligned
int DImage::_drawTextFileNum = 0;
//TODO:The methods do not necessarily protect against using the source
// and destination being the same object, (therefore problems with create()ing
// the source data to use as dest data could arise without warning if the
// person using this library isn't careful).  This should be protected against
// by checking to see if src==dst, and if so, detaching the data before 
// overwriting it, and then cleaning up the detached data.

//TODO: instead of abort() on errors, I need to have an error,warning,and
// debug reporting mechanism in the library so that applications can handle
// the errors however they want, instead of the library just ending the
// application whenever it feels like it (due to errors). (Although, the
// app shouldn't be calling in a manner that generates errors to begin with).



inline double DImage_radFromDeg(double degrees){
  return (M_PI * degrees) / 180.;
}
inline double DImage_degFromRad(double radians){
  return 180.*radians/M_PI;
}
inline double DImage_minDbl(double a, double b){
  return (a < b) ? a : b;
}
inline double DImage_maxDbl(double a, double b){
  return (a > b) ? a : b;
}

/// default constructor
/** Creates an empty DImage object with no data buffer (data buffer is NULL)*/
DImage::DImage(){
  DInstanceCounter::addInstance("DImage");
//   fprintf(stderr, ">>>DImage default constructor this=%p\n", this);
  pData = NULL;
  // default values for undefined image:
  _w = 0;
  _h = 0;
  _actualW = 0;
  _actualH = 0;
  _dataSize = 0;
  _sampleSize = sizeof(D_uint8);
  _imgType = DImage_u8;
  _fInterleaved = true;
  _allocMethod = AllocationMethod_malloc;
  _numChan = 1;
}

///copy constructor
/**Makes a deep copy of src in the new DImage object*/
// DImage::DImage(const DImage &src, bool fCopyDataBuffer):mapProps(src.mapProps),vectComments(src.vectComments){
//   DInstanceCounter::addInstance("DImage");
// //   fprintf(stderr, ">>>DImage copy constructor this=%p\n", this);
//   pData = NULL;
//   _internal_DImage_assign(*this, src, fCopyDataBuffer, false, false);
//   //  *this = src; // call the assignment operator member function to copy src
// }
DImage::DImage(const DImage &src, bool fCopyDataBuffer){
  DInstanceCounter::addInstance("DImage");
//   fprintf(stderr, ">>>DImage copy constructor this=%p\n", this);
  pData = NULL;
  _internal_DImage_assign(*this, src, fCopyDataBuffer, false, false);
  mapProps = src.mapProps;
  vectComments = src.vectComments;
  //  *this = src; // call the assignment operator member function to copy src
}


///Constructor that reads image data from an existing file named stPath
/** If stPath references a valid (and readable) image file of a
 *  supported format, it will be loaded into the new DImage object*/
DImage::DImage(const char *stPath){
  DInstanceCounter::addInstance("DImage");
//   fprintf(stderr, ">>>DImage filename constructor this=%p stPath=%s\n",
// 	  this, stPath);
  pData = NULL;
  _w = 0;
  _h = 0;
  _actualW = 0;
  _actualH = 0;
  _dataSize = 0;
  _sampleSize = sizeof(D_uint8);
  _imgType = DImage_u8;
  _fInterleaved = true;
  _allocMethod = AllocationMethod_malloc;
  this->load(stPath);
}

///Destructor
DImage::~DImage(){
  DInstanceCounter::removeInstance("DImage");
//   fprintf(stderr, ">>>DImage destructor this=%p\n", this);
  deallocateBuffer();
  _w = 0;
  _h = 0;
  _actualW = 0;
  _actualH = 0;
  _dataSize = 0;
  _sampleSize = sizeof(D_uint8);
  _imgType = DImage_u8;
  _fInterleaved = true;
  _allocMethod = AllocationMethod_malloc;
}

///Assignment operator - makes a deep copy
const DImage& DImage::operator=(const DImage &src){
//   fprintf(stderr,">>>DImage::operator=() this=%p\n", this);
  if(this != &src){ // don't assign self (ex: 'imgA = imgA')
    //copy data,props,comments
    _internal_DImage_assign(*this,src,true,true,true);
  }
  return *this;
}

// This function is called by the copy constructor, the assignment operator,
// and a few other functions
// to copy the member variables, etc. from the src object to dst object.
// fCopyProps should be false from copy constructor to avoid doing a
// second copy of the properties container (std::map).
// fCopyProps should be true from assignment operator since it hasn't 
// been copied.
void DImage::_internal_DImage_assign(DImage &dst, const DImage &src,
				     bool fCopyData, bool fCopyProps,
				     bool fCopyComments){
//   fprintf(stderr, ">>>_internal_DImage_assign() fCopyData=%d fCopyProps=%d "
// 	  "fCopyComments=%d\n", fCopyData, fCopyProps, fCopyComments);
  // TODO:copy stuff here! (make sure I haven't forgotten anything)
#ifdef DEBUG    
  if((&src) == (&dst)){
    fprintf(stderr,"DImage::_internal_DImage_assign() called with src==dst\n");
    return;
  }
#endif
  dst.deallocateBuffer();// if this has data, clear it first
  dst._w = src._w;
  dst._h = src._h;
  dst._actualW = src._actualW;
  dst._actualH = src._actualH;
  dst._numChan = src._numChan;
  dst._imgType = src._imgType;
  dst._fInterleaved = src._fInterleaved;
  dst._allocMethod = src._allocMethod;
  dst._dataSize = src._dataSize;
  dst._sampleSize = src._sampleSize;
  
  
  if(fCopyProps)
    dst.mapProps = src.mapProps;

  if(fCopyComments)
    dst.vectComments = src.vectComments;

  if((fCopyData) && (NULL != src.pData)){
    if(AllocationMethod_daligned == dst._allocMethod)
      dst.pData =
	(D_uint8*)daligned_malloc(dst._dataSize, DImage::_data_alignment);
    else if(AllocationMethod_malloc == dst._allocMethod)
      dst.pData = (D_uint8*)malloc(dst._dataSize);
    else if(AllocationMethod_new == dst._allocMethod)
      dst.pData = new D_uint8[dst._dataSize];
    else{
      fprintf(stderr,"DImage::_internal_DImage_assign() "
	      "invalid allocation method used\n");
    }
    
    // TODO: instead of just checking if it worked or not, maybe we should
    // (in the case of failure) also tell the memory manager to save some of
    // the managed images or buffers out to disk and then try allocating again.
    // if I change this, I should change the code in create(), also.
    D_CHECKPTR(dst.pData);
    memcpy(dst.pData, src.pData, dst._dataSize);
  }

}

/// copy the properties (but not comments) of src to this DImage
/** Any previously existing properties (but not comments) for this
    DImage will be discarded */
void DImage::copyProperties(DImage &src){
  this->mapProps = src.mapProps;
}

/// copy the comments (but not properties) of src to this DImage
/** Any previously existing comments (but not properties) for this
    DImage will be discarded. */
void DImage::copyComments(DImage &src){
  this->vectComments = src.vectComments;
}

/// create a DImage object of a specified type with w by h buffer
/** This function does not initialize the pixel data.  Call the
 * appropriate fill() routine if you want the pixel data initialized.
 */
bool DImage::create(int w, int h, DImageType imgType, int numChannels,
		    D_AllocationMethod allocMeth){
  size_t bufSize;
  bool retVal = true;

  this->_imgType = imgType;

  _fInterleaved = true;//interleave all types except for DImage_<flt/dbl>_multi
  _numChan = 1;// all types are 1 channel except the RGB and multi-chan types
  switch(imgType){
    case DImage_u8:
      _sampleSize=sizeof(D_uint8);
      break;
    case DImage_u16:
      _sampleSize=sizeof(D_uint16);
      break;
    case DImage_u32:
      _sampleSize=sizeof(D_uint32);
      break;
    case DImage_RGB:
      _numChan = 3;
      _sampleSize=sizeof(D_uint8);
      break;
    case DImage_RGB_16:
      _numChan = 3;
      _sampleSize=sizeof(D_uint16);
      break;
    case DImage_flt_multi:
      _numChan = numChannels;
      _sampleSize=sizeof(float);
      _fInterleaved = false;
      break;
    case DImage_dbl_multi:
      _numChan = numChannels;
      _sampleSize=sizeof(double);
      _fInterleaved = false;
      break;
    case DImage_cmplx:
      _sampleSize=sizeof(std::complex<double>);
      break;
    default:
      fprintf(stderr, "DImage::create() unknown DImageType\n");
      abort();
  }
  bufSize = _numChan*w*h*_sampleSize;

  // allocate the data buffer.  Instead of normal malloc or new, ensure that
  // the buffer is aligned to fit the users needs (for example, if Altivec SIMD
  // instructions are used, data blocks must be aligned to 16-byte boundaries).
  deallocateBuffer();
  if(AllocationMethod_daligned == allocMeth)
    this->pData = (D_uint8*)daligned_malloc(bufSize, DImage::_data_alignment);
  else if(AllocationMethod_malloc == allocMeth)
    this->pData = (D_uint8*)malloc(bufSize);
  else if(AllocationMethod_new == allocMeth)
    this->pData = new D_uint8[bufSize];
  else{
    retVal = false;
    fprintf(stderr,"DImage::create() invalid allocation method (%d) used\n",
	    allocMeth);
  }
  _dataSize = bufSize;
  _allocMethod = allocMeth;

  // TODO: instead of just checking if it worked or not, maybe we should
  // (in the case of failure) also tell the memory manager to save some of
  // the managed images or buffers out to disk and then try allocating again.
  // if I change this, I should also change  _internal_DImage_assign() to match
  D_CHECKPTR(this->pData);
  
  _w = w;
  _h = h;
  _actualW = w;
  _actualH = h;

  return retVal;
}



/// load an image from a file into the DImage object
/** This function tries to automatically determine the file format
 *  from the file header.
 */
bool DImage::load(const char *stPath, DFileFormat fmt){
  FILE *fin;
  DFileFormat fmtHead;
  bool retVal = false;

  fin = fopen(stPath, "rb");
  if(!fin){
    fprintf(stderr, "DImage::load() can't open file '%s' for input!\n",stPath);
    return false;
  }
  fmtHead = getImageFileFormat(fin);
  if(DFileFormat_unknown == fmt){
    fmt = fmtHead;
  }
  if(DFileFormat_pnm == fmt){
    if((fmtHead < DFileFormat_pbm_plain) || (fmtHead > DFileFormat_ppm)){
      fprintf(stderr,
	      "DImage::load() file not of the type specified in parameter\n");
    }
    fmt = fmtHead;
  }
  if(DFileFormat_unknown == fmtHead){
    fmtHead = fmt;
  }
  if(fmt != fmtHead){
    fprintf(stderr,
	    "DImage::load() file is not of the type specified in parameter\n");
  }

  switch(fmt){
    case DFileFormat_pbm_plain:
      retVal=DImageIO::load_image_pbm_plain(this,fin);
      fclose(fin);
      break;
    case DFileFormat_pgm_plain:
      retVal=DImageIO::load_image_pgm_plain(this,fin);
      fclose(fin);
      break;
    case DFileFormat_ppm_plain:
      retVal=DImageIO::load_image_ppm_plain(this,fin);
      fclose(fin);
      break;
    case DFileFormat_pbm:
      retVal=DImageIO::load_image_pbm_raw(this,fin);
      fclose(fin);
      break;
    case DFileFormat_pgm:
      retVal=DImageIO::load_image_pgm_raw(this,fin);
      fclose(fin);
      break;
    case DFileFormat_ppm:
      retVal=DImageIO::load_image_ppm_raw(this,fin);
      fclose(fin);
      break;
    case DFileFormat_jpg:
      retVal=DImageIO::load_image_jpeg(this,fin);
      fclose(fin);
      break;
    case DFileFormat_png:
      retVal=DImageIO::load_image_png(this,fin);
      fclose(fin);
      break;
    case DFileFormat_tiff:
      fclose(fin);
      retVal=DImageIO::load_image_tiff(this,stPath);
      break;
    case DFileFormat_gif:
      retVal=DImageIO::load_image_gif(this,stPath);
      fclose(fin);
      break;
    case DFileFormat_bmp:
      retVal=DImageIO::load_image_bmp(this,stPath);
      fclose(fin);
      break;
    case DFileFormat_gnuplot:
      retVal=DImageIO::load_image_gnuplot(this,fin);
      fclose(fin);
      break;
    default:
      fprintf(stderr, "DImage::load() unknown image file type\n");
      fclose(fin);
  }
  return retVal;
}

///Saves the current image in the specified format (default=pnm).
/** If fStripComments is false, the image properties will be saved as
 *  comments (if supported by the file format) along with the image.
 *  Currently sipported formats are DFileFormat_pnm (saves as raw pgm or
 *  ppm of 8 or 16-bit depending on the image type), DFileFormat_gnuplot,
 *  DFileFormat_png, and DFileFormat_jpg (saves as grayscale or rgb jpeg
 *  using quality=75, optimized=false, progressive=false.  To use other
 *  jpeg parameter settings, call DImageIO::save_image_jpeg() directly).
 */
bool DImage::save(const char *stPath, DFileFormat fmt, bool fSaveProps,
		  bool fSaveComments){
  switch(fmt){
    case DFileFormat_pnm:
      return DImageIO::save_image_pnm(this, stPath, fSaveProps, fSaveComments);
      break;
    case DFileFormat_gnuplot:
      return DImageIO::save_image_gnuplot(this, stPath, fSaveProps,
					  fSaveComments);
      break;
    case DFileFormat_png:
      return DImageIO::save_image_png(this, stPath, fSaveProps,
				      fSaveComments);
      break;
    case DFileFormat_jpg:
      return DImageIO::save_image_jpeg(this, stPath, fSaveProps,
				       fSaveComments, 75, false, false);
      break;
    default:
      fprintf(stderr, "The only file formats currently supported are "
	      "pnm, png, jpg and gnuplot\n");
  }
  return false;
}


///Allows the DImage to take control (and responsibility for deleting) pBuf
/**The caller is responsible for ensuring that pBuf points to a valid
 * buffer of the correct size and type for the image data.  (The
 * buffer can actually be larger than the declared image size,
 * though).  Since DImage will now be responsible for managing the
 * memory, it is imperative that the correct D_AllocationMethod is
 * passed in and that the memory is not deallocated elsewhere.  If the
 * user wishes to deallocate the memory elsewhere, the data should be
 * copied instead of directly passing in a pointer to the memory.
 */
void DImage::setDataBuffer(void *pBuf, D_AllocationMethod allocMeth){
  deallocateBuffer();
  pData = (D_uint8*)pBuf;
  _allocMethod = allocMeth;
}

///releases the data buffer without deallocating it
/** A pointer to the data buffer is returned, and the image's internal
 *  pointer is set to NULL, as if no data buffer had been created.
 *  The user becomes responsible for making sure that the data buffer
 *  gets deallocated properly elsewhere when it is no longer needed.
 */
void* DImage::releaseDataBuffer(){
  void *retVal = (void*)pData;
  pData = NULL;
  return retVal;
}


///deallocates the buffer (if it has been allocated)
/** This function calls the appropriate deallocation function, depending on
 *  how the buffer was created (i.e. malloc/free, new/delete, etc.)
 */
void DImage::deallocateBuffer(){
  if(pData != NULL){
    if(AllocationMethod_malloc == _allocMethod)
      free(pData);
    else if(AllocationMethod_daligned == _allocMethod)
      daligned_free(pData);
    else if(AllocationMethod_new == _allocMethod)
      delete [] (D_uint8*)pData;
    else
      fprintf(stderr,"DImage::deallocateBuffer() unrecognized _allocMethod\n");
    pData = NULL;
  }
}

///Returns a DImage object that contains the sub-image specified by x,y,w,h
/** See copy_() for more details */
DImage DImage::copy(int x, int y, int w, int h) const{
  DImage dst(*this, false);
  copy_(dst,x,y,w,h);
  return dst;
}
///make a copy of a portion of this image and place it in imgDst
/**x,y is the upper-left corner of the region to be copied. w,h are the
 * width and height of the region to be copied.  imgDst will be set to
 * the same type of image and any previously held data buffer will be
 * deallocated.  Any image properties in imgDst will also be cleared, and
 * the image properties from this image will NOT be automatically copied over
 * to imgDst.
 */
void DImage::copy_(DImage &imgDst, int x, int y, int w, int h,
		   D_AllocationMethod allocMeth) const{
  size_t pxlLen;
  size_t rowLen;
  size_t srcRowLen;
  int srcChanOffs;
  int dstChanOffs;
  //TODO: compare how fast bcopy and memcpy are, and see if there is some
  // point (small row lengths) where we should just do a loop of assignments
  // instead of either.
#ifdef DEBUG
  if((x < 0)||(x >= this->_w)||(y < 0)||(y >= this->_h)){
    fprintf(stderr, "DImage::copy_() x,y==%d,%d out of range\n", x,y);
  }
  if((w < 0)||(x+w > this->_w)||(h < 0)||(y+h > this->_h)){
    fprintf(stderr, "DImage::copy_() w,h==%d,%d out of range\n", w,h);
  }
  if(this == &imgDst){
    fprintf(stderr, "DImage::copy_() source == dest!\n");
    return;
  }
  if(NULL == this->pData){
    fprintf(stderr, "DImage::copy_() this->pData == NULL\n");
    return;
  }
  // shouldn't do the next check! we're creating the buffer anyway
//   if(NULL == imgDst.pData){
//     fprintf(stderr, "DImage::copy_() imgDst.pData == NULL\n");
//     return;
//   }
#endif
  if(allocMeth == AllocationMethod_src)
    allocMeth = this->_allocMethod;
  imgDst.create(w,h,this->_imgType, this->_numChan, allocMeth);
  imgDst.mapProps.clear(); // clear any previous properties- they're not valid
  imgDst.vectComments.clear(); // clear any previous comments for same reason
  // copy the data buffer
  if(_fInterleaved){
    pxlLen = (this->_numChan) * (this->_sampleSize);
    rowLen = w * pxlLen;
    srcRowLen = this->_w * pxlLen;
    for(int yDst = 0; yDst < h; ++yDst){
      memcpy(&((imgDst.pData)[yDst*rowLen]),
	     &((this->pData)[(yDst+y)*srcRowLen + x * pxlLen]), rowLen);
    }
  }
  else{ /* multiple non-interleaved channels of data */
    pxlLen = (this->_sampleSize);
    rowLen = w * pxlLen;
    srcRowLen = this->_w * pxlLen;
    for(int chan = 0; chan < this->_numChan; ++chan){
      srcChanOffs = chan * srcRowLen * this->_h;
      dstChanOffs = chan * rowLen * h;
      for(int yDst = 0; yDst < h; ++yDst){
	memcpy(&((imgDst.pData)[srcChanOffs + yDst*rowLen]),
	       &((this->pData)[dstChanOffs + (yDst+y)*srcRowLen + x * pxlLen]),
	       rowLen);
      }/* end for yDst */
    }/* end for chan */
  } /* end else multiple non-interleaved channels of data */
}

/// paste a region of srcImg into this image
/** It is not an error to use this image as both source and dest.
 */
void DImage::pasteFromImage(int dstX, int dstY,
			    const DImage &imgSrc,int srcX,int srcY,int w,int h,
			    CompositeMode cmode){
  bool fOverlap = false; // if srcImg == this and regions overlap, use memmove
  size_t pxlLen; // number of bytes in a pixel (RGB=3, RGB16=6, etc)
  size_t rowLen; // number of bytes per row to be copied
  size_t srcRowLen; // # bytes in a row of the source image
  size_t dstRowLen; // # bytes in a row of the destination image
  size_t srcChanOffs;// #bytes to offset per channel in multi-channel src img
  size_t dstChanOffs;// #bytes to offset per channel in multi-channel dest img
  
#ifdef DEBUG
  if(imgSrc._imgType != this->_imgType){
    fprintf(stderr, "DImage::pasteFromImage() mismatched image types!\n");
    abort();
    return;
  }
  if((srcX < 0)||(srcX >= imgSrc._w)||(srcY<0)||(srcY >= imgSrc._h)){
    fprintf(stderr,"DImage::pasteFromImage() srcX,srcY (%d,%d) out of range w=%d h=%d\n", srcX, srcY, imgSrc._w, imgSrc._h);
    abort();
    return;
  }
  if(cmode != CompositeMode_Source)
    fprintf(stderr, "DImage::pasteFromImage() unsupported composite mode\n");
  if((w < 0)||(dstX+w > this->_w)||(h < 0)||(dstY+h > this->_h)){
    fprintf(stderr, "DImage::pasteFromImage() w,h==%d,%d out of range in dst\n", w,h);
    if(dstX+w > this->_w){
      w = this->_w - dstX;
      fprintf(stderr,"DImage::pasteFromImage() changing w to %d\n",w);
    }
    if(dstY+h > this->_h){
      h = this->_h - dstY;
      fprintf(stderr,"DImage::pasteFromImage() changing h to %d\n",h);
    }
  }
  if((w < 0)||(srcX+w > imgSrc._w)||(h < 0)||(srcY+h > imgSrc._h)){
    fprintf(stderr, "DImage::pasteFromImage() w,h==%d,%d out of range in src\n", w,h);
    if(srcX+w > imgSrc._w){
      w = imgSrc._w - srcX;
      fprintf(stderr,"DImage::pasteFromImage() changing w to %d\n",w);
    }
    if(srcY+h > imgSrc._h){
      h = imgSrc._h - srcY;
      fprintf(stderr,"DImage::pasteFromImage() changing h to %d\n",h);
    }
  }
  if((dstX < 0)||(dstX >= this->_w)||(dstY<0)||(dstY >= this->_h)){
    fprintf(stderr,"DImage::pasteFromImage() dstX,dstY (%d,%d) out of range\n",
	    dstX, dstY);
    abort();
    return;
  }
  if(NULL == this->pData){
    fprintf(stderr, "DImage::pasteFromImage() this->pData == NULL\n");
    abort();
    return;
  }
  if(NULL == imgSrc.pData){
    fprintf(stderr, "DImage::pasteFromImage() imgSrc.pData == NULL\n");
    abort();
    return;
  }
#endif
  
  if(this == &imgSrc){ // check if the region overlaps
    DRect *rsrc = new DRect(srcX, srcY, w, h);
    DRect *rdst = new DRect(dstX, dstY, w, h);
    fOverlap = rsrc->intersectsWith(*rdst);
    free(rsrc);
    free(rdst);
  }
  
  // copy the data buffer
  if(_fInterleaved){
    pxlLen = (this->_numChan) * (this->_sampleSize);
    rowLen = w * pxlLen;
    srcRowLen = imgSrc._w * pxlLen;
    dstRowLen = this->_w * pxlLen;
    if(!fOverlap){
      for(int yDst = 0; yDst < h; ++yDst){
	memcpy(&((this->pData)[(yDst+dstY)*dstRowLen + dstX * pxlLen]),
	       &((imgSrc.pData)[(yDst+srcY)*srcRowLen + srcX * pxlLen]),
	       rowLen);
      }
    } else{ // overlapped, use memmove instead of memcpy
      for(int yDst = 0; yDst < h; ++yDst){
	memmove(&((this->pData)[(yDst+dstY)*dstRowLen + dstX * pxlLen]),
		&((imgSrc.pData)[(yDst+srcY)*srcRowLen + srcX * pxlLen]),
		rowLen);
      }
    }
  }
  else{ /* multiple non-interleaved channels of data */
    pxlLen = (this->_sampleSize);
    rowLen = w * pxlLen;
    srcRowLen = imgSrc._w * pxlLen;
    dstRowLen = this->_w * pxlLen;
    for(int chan = 0; chan < this->_numChan; ++chan){
      srcChanOffs = chan * srcRowLen * imgSrc._h;
      dstChanOffs = chan * dstRowLen * this->_h;
      if(!fOverlap){
	for(int yDst = 0; yDst < h; ++yDst){
	  memcpy(&((this->pData)[dstChanOffs +
				 (yDst+dstY)*dstRowLen+dstX*pxlLen]),
		 &((imgSrc.pData)[srcChanOffs +
				  (yDst+srcY)*srcRowLen + srcX*pxlLen]),
		 rowLen);
	}/* end for yDst */
      } else{ // overlapped, use memmove instead of memcpy
	for(int yDst = 0; yDst < h; ++yDst){
	  memmove(&((this->pData)[dstChanOffs +
				  (yDst+dstY)*dstRowLen + dstX*pxlLen]),
		  &((imgSrc.pData)[srcChanOffs+
				   (yDst+srcY)*srcRowLen + srcX*pxlLen]),
		  rowLen);
	}/* end for yDgrepst */
      }
    }/* end for chan */
  } /* end else multiple non-interleaved channels of data */
}
  

///Linearly maps the image value range to [max..min]
/** perctMin and perctMax are the percent of pixels (as a fraction, so
 *  2% would be 0.02, for example) that map to max and min,
 *  respectively.  This is helpful in case there are extreme outliers.
 *  For example, the darkest 2% of an image can be mapped to 0(black)
 *  and the brightest 2% to 255(white).  if perctMin/perctMax are set
 *  to zero, then the actual extrema values will be used, even if they
 *  are outliers.  If there is multi-channel data (RGB, etc.) all
 *  channels will be used in calculating the old max/min, and the same
 *  linear scaling will be applied for all channels.  If you want them
 *  stretched independantly, you must do that yourself.  when
 *  calculating the source max/min values, 256 discrete bins are used,
 *  which limits how accurate perctMin/perctMax are with non-8-bit
 *  data.
 */
void DImage::setDataRange(double min, double max,
			  float perctMin, float perctMax){
  double oldMin;
  double oldMax;
  double oldRange;
  double newRange;
  double newVal;
  double scl;
  double offs;
  int len;
  int rgHist[256]; // histogram for figuring out min/max based on percentages
  int numPxls;
  double threshPxls;
  int idxVal;
  D_uint8 *pu8;
//   D_uint8 *pu16;
  D_uint32 *pu32;
  float *pflt;
  double *pdbl;

#ifdef DEBUG
  if((perctMin < 0.)||(perctMin > 1.0)||(perctMax<0.)||(perctMax>1.0)){
    fprintf(stderr,
	    "DImage::setDataRange() perctMin/perctMax must be [0..1]\n");
    return;
  }
#endif

  newRange = max - min;

  memset(rgHist, 0, 256*sizeof(int));
  switch(this->_imgType){
    case DImage_u8:
    case DImage_RGB:
      pu8 = this->pData; // temp pointer to iterate through the data
      len = _numChan * _w * _h;
      for(int idx = 0; idx < len; ++idx, ++pu8){
	++(rgHist[*pu8]); // histogram
      }
      threshPxls = perctMin * len;
      if(threshPxls == 0.)
	threshPxls = 1.;
      numPxls = 0;
      oldMin = 0;
      for(int i = 0; i < 256; ++i){
	numPxls += rgHist[i];
	if(numPxls >= threshPxls){
	  oldMin = i;
	  break;
	}
      }
      threshPxls = (perctMax) * len;
      if(threshPxls == 0.)
	threshPxls = 1.;
      numPxls = 0;
      oldMax = 255.;
      for(int i = 255; i >= 0; --i){
	numPxls += rgHist[i];
	if(numPxls >= threshPxls){
	  oldMax = i;
	  break;
	}
      }

      offs = 0 - oldMin;
      oldRange = oldMax - oldMin;
      if(oldRange < 1.)
	scl = 1.;
      else
	scl = newRange / oldRange;

#if 0
      printf("       DImage_u8 or DImage_RGB\n");
      printf("       oldMin=%f oldMax=%f oldRange=%f offs=%f\n",
	     oldMin, oldMax, oldRange, offs);
      printf("       min=%f max=%f newRange=%f scl=%f\n",
	     min, max, newRange, scl);
#endif      
      
      pu8 = this->pData;
      for(int idx = 0; idx < len; ++idx, ++pu8){
	newVal = ((*pu8) + offs) * scl + min;
	if(newVal < min)
	  newVal = min;
	if(newVal > max)
	  newVal = max;
	(*pu8) = (D_uint8)newVal;
      }
      break;
    case DImage_u32:
      pu32 = (D_uint32*)pData; //temp pointer to iterate through the data
      len = _numChan * _w * _h;
      oldMin = oldMax = (*pu32);
      for(int idx = 0; idx < len; ++idx, ++pu32){
	if((*pu32)<oldMin)
	  oldMin = (*pu32);
	else if((*pu32)>oldMax)
	  oldMax = (*pu32);
      }
      offs = 0 - oldMin;
      oldRange = oldMax - oldMin;
      if(oldRange == 0.)
	scl = 1.;
      else
	scl = newRange / oldRange;
//       printf("oldMin=%.1f oldMax=%.1f offs=%.1f scl=%.5f\n",
// 	     oldMin,oldMax,offs,scl);
      // if perctMin or perctMax not equal to 0., then bin the values
      // into a discrete histogram and do what we do for uint8 data
      if((perctMin != 0.) || (perctMax != 0.)){
	pu32 = (D_uint32*)pData; // reset the pointer to beginning
	for(int idx = 0; idx < len; ++idx, ++pu32){
	  idxVal = (int)( ((*pu32)+offs)*scl);
	  if(idxVal < 0)
	    idxVal = 0;
	  else if(idxVal > 255)
	    idxVal = 255;
	  ++(rgHist[idxVal]); // histogram
	}
	threshPxls = perctMin * len;
	if(threshPxls == 0.)
	  threshPxls = 1.;
	numPxls = 0;
	oldMin = 0;
	for(int i = 0; i < 256; ++i){
	  numPxls += rgHist[i];
	  if(numPxls >= threshPxls){
	    oldMin = i;
	    break;
	  }
	}
	threshPxls = (1.0-perctMax) * len;
	if(threshPxls == 0.)
	  threshPxls = 1.;
	numPxls = 0;
	oldMax = 255.;
	for(int i = 255; i >= 0; --i){
	  numPxls += rgHist[i];
	  if(numPxls >= threshPxls){
	    oldMax = i;
	    break;
	  }
	}
	
	offs = 0 - oldMin;
	oldRange = oldMax - oldMin;
	if(oldRange == 0.)
	  scl = 1.;
	else
	  scl = newRange / oldRange;
      }
      
      pu32 = (D_uint32*)pData; // reset the pointer to beginning
      for(int idx = 0; idx < len; ++idx, ++pu32){
	newVal = ((*pu32) + offs) * scl + min;
	if(newVal < min)
	  newVal = min;
	if(newVal > max)
	  newVal = max;
	(*pu32) = (D_uint32)newVal;
      }
      break;
    case DImage_flt_multi:
      pflt = (float*)pData; // temp pointer to iterate through the data
      len = _numChan * _w * _h;
      oldMin = oldMax = (*pflt);
      for(int idx = 0; idx < len; ++idx, ++pflt){
	if((*pflt)<oldMin)
	  oldMin = (*pflt);
	else if((*pflt)>oldMax)
	  oldMax = (*pflt);
      }
      offs = 0 - oldMin;
      oldRange = oldMax - oldMin;
      if(oldRange == 0.)
	scl = 1.;
      else
	scl = newRange / oldRange;

      // if perctMin or perctMax not equal to 0., then bin the values
      // into a discrete histogram and do what we do for uint8 data
      if((perctMin != 0.) || (perctMax != 0.)){
	pflt = (float*)pData; // reset the pointer to beginning
	for(int idx = 0; idx < len; ++idx, ++pflt){
	  idxVal = (int)( ((*pflt)+offs)*scl);
	  if(idxVal < 0)
	    idxVal = 0;
	  else if(idxVal > 255)
	    idxVal = 255;
	  ++(rgHist[idxVal]); // histogram
	}
	threshPxls = perctMin * len;
	if(threshPxls == 0.)
	  threshPxls = 1.;
	numPxls = 0;
	oldMin = 0;
	for(int i = 0; i < 256; ++i){
	  numPxls += rgHist[i];
	  if(numPxls >= threshPxls){
	    oldMin = i;
	    break;
	  }
	}
	threshPxls = (1.0-perctMax) * len;
	if(threshPxls == 0.)
	  threshPxls = 1.;
	numPxls = 0;
	oldMax = 255.;
	for(int i = 255; i >= 0; --i){
	  numPxls += rgHist[i];
	  if(numPxls >= threshPxls){
	    oldMax = i;
	    break;
	  }
	}
	
	offs = 0 - oldMin;
	oldRange = oldMax - oldMin;
	if(oldRange == 0.)
	  scl = 1.;
	else
	  scl = newRange / oldRange;
      }
#if 0
      printf("       DImage_flt_multi\n");
      printf("       oldMin=%f oldMax=%f oldRange=%f offs=%f\n",
	     oldMin, oldMax, oldRange, offs);
      printf("       min=%f max=%f newRange=%f scl=%f\n",
	     min, max, newRange, scl);
#endif      
      
      pflt = (float*)pData; // reset the pointer to beginning
      for(int idx = 0; idx < len; ++idx, ++pflt){
	newVal = ((*pflt) + offs) * scl + min;
	if(newVal < min)
	  newVal = min;
	if(newVal > max)
	  newVal = max;
	(*pflt) = (float)newVal;
      }
      break;
    case DImage_dbl_multi:
      pdbl = (double*)pData; // temp pointer to iterate through the data
      len = _numChan * _w * _h;
      oldMin = oldMax = (*pdbl);
      for(int idx = 0; idx < len; ++idx, ++pdbl){
	if((*pdbl)<oldMin)
	  oldMin = (*pdbl);
	else if((*pdbl)>oldMax)
	  oldMax = (*pdbl);
      }
      offs = 0. - oldMin;
      oldRange = oldMax - oldMin;
      if(oldRange == 0.)
	scl = 1.;
      else
	scl = newRange / oldRange;


      // if perctMin or perctMax not equal to 0., then bin the values
      // into a discrete histogram and do what we do for uint8 data
      if((perctMin != 0.) || (perctMax != 0.)){
	pdbl = (double*)pData; // reset the pointer to beginning
	for(int idx = 0; idx < len; ++idx, ++pdbl){
	  idxVal = (int)( ((*pdbl)+offs)*scl);
	  if(idxVal < 0)
	    idxVal = 0;
	  else if(idxVal > 255)
	    idxVal = 255;
	  ++(rgHist[idxVal]); // histogram
	}
	threshPxls = perctMin * len;
	if(threshPxls == 0.)
	  threshPxls = 1.;
	numPxls = 0;
	oldMin = 0;
	for(int i = 0; i < 256; ++i){
	  numPxls += rgHist[i];
	  if(numPxls >= threshPxls){
	    oldMin = i;
	    break;
	  }
	}
	threshPxls = (1.0-perctMax) * len;
	if(threshPxls == 0.)
	  threshPxls = 1.;
	numPxls = 0;
	oldMax = 255.;
	for(int i = 255; i >= 0; --i){
	  numPxls += rgHist[i];
	  if(numPxls >= threshPxls){
	    oldMax = i;
	    break;
	  }
	}
	
	offs = 0 - oldMin;
	oldRange = oldMax - oldMin;
	if(oldRange == 0.)
	  scl = 1.;
	else
	  scl = newRange / oldRange;
      }
#if 0
      printf("       DImage_dbl_multi\n");
      printf("       oldMin=%f oldMax=%f oldRange=%f offs=%f\n",
	     oldMin, oldMax, oldRange, offs);
      printf("       min=%f max=%f newRange=%f scl=%f\n",
	     min, max, newRange, scl);
#endif      
      
      pdbl = (double*)pData; // reset the pointer to beginning
      for(int idx = 0; idx < len; ++idx, ++pdbl){
	newVal = ((*pdbl) + offs) * scl + min;
// 	if(newVal < min)
// 	  newVal = min;
// 	if(newVal > max)
// 	  newVal = max;
	(*pdbl) = (double)newVal;
      }
      break;
    case DImage_cmplx:
    case DImage_u16:
    case DImage_RGB_16:
    default:
      fprintf(stderr,
	      "DImage::setDataRange() not yet implemented for specified "
	      "image type\n");
      abort();
    break;
  }
#if 0
  {
    double dmin, dmax;
    getDataRange(&dmin, &dmax);
    printf("DImage::setDataRange to min=%f max=%f   "
	   "resulted in min=%f max=%f\n", min, max, dmin, dmax);
  }
#endif
}

///Caps the image value range to [min..max]
/** Any value less than min is set to min and any value greater than max is
 *  set to max.
 */
void DImage::capDataRange(double min, double max){
  int len;

  if(max < min){
    fprintf(stderr, "DImage::capDataRange(min=%f, max=%f) min > max!\n",
	    min, max);
    abort();
  }

  switch(this->_imgType){
    case DImage_u8:
    case DImage_RGB:
      {
	D_uint8 min8, max8;
	D_uint8 *pu8;
	min8 = (D_uint8)min;
	max8 = (D_uint8)max;
	pu8 = this->pData; // temp pointer to iterate through the data
	len = _numChan * _w * _h;
	for(int idx = 0; idx < len; ++idx){
	  if(pu8[idx] < min)
	    pu8[idx] = min8;
	  else if(pu8[idx] > max)
	    pu8[idx] = max8;
	}
      }
      break;
    case DImage_u16:
    case DImage_RGB_16:
      {
	D_uint16 min16, max16;
	D_uint16 *pu16;
	min16 = (D_uint16)min;
	max16 = (D_uint16)max;
	pu16 = this->dataPointer_u16();
	len = _numChan * _w * _h;
	for(int idx = 0; idx < len; ++idx){
	  if(pu16[idx] < min)
	    pu16[idx] = min16;
	  else if(pu16[idx] > max)
	    pu16[idx] = max16;
	}
      }
      break;
    case DImage_u32:
      {
	D_uint32 min32, max32;
	D_uint32 *pu32;
	min32 = (D_uint32)min;
	max32 = (D_uint32)max;
	pu32 = this->dataPointer_u32();
	len = _numChan * _w * _h;
	for(int idx = 0; idx < len; ++idx){
	  if(pu32[idx] < min)
	    pu32[idx] = min32;
	  else if(pu32[idx] > max)
	    pu32[idx] = max32;
	}
      }
      break;
    case DImage_flt_multi:
      {
	float minflt, maxflt;
	float *pflt;
	minflt = (float)min;
	maxflt = (float)max;
	pflt = this->dataPointer_flt();
	len = _numChan * _w * _h;
	for(int idx = 0; idx < len; ++idx){
	  if(pflt[idx] < minflt)
	    pflt[idx] = minflt;
	  else if(pflt[idx] > maxflt)
	    pflt[idx] = maxflt;
	}
      }
      break;
    case DImage_dbl_multi:
      {
	double mindbl, maxdbl;
	double *pdbl;
	mindbl = (double)min;
	maxdbl = (double)max;
	pdbl = this->dataPointer_dbl();
	len = _numChan * _w * _h;
	for(int idx = 0; idx < len; ++idx){
	  if(pdbl[idx] < mindbl)
	    pdbl[idx] = mindbl;
	  else if(pdbl[idx] > maxdbl)
	    pdbl[idx] = maxdbl;
	}
      }
      break;
    case DImage_cmplx:
    default:
      fprintf(stderr,
	      "DImage::capDataRange() not yet implemented for specified "
	      "image type\n");
      abort();
    break;
  }
    
}




///Puts the data into a Log scale (only for float, double images)
/** This function makes each pixel value be the log of its previous
 *  value.  By default, it is the "log base 10", but you can change this
 *  with the base parameter.  If you are going to display the result,
 *  you probably want to perform a setDataRange() call after this
 *  method is used to scale to the appropriate range, such as [0..255].
 */
void DImage::setDataRangeLog(double base){
  double log_base;
  size_t len;

  log_base = log(base);
  if(_imgType == DImage_flt_multi){
    float *ptf =(float*)(this->pData);// temp pointer to iterate through data
    len = _numChan * _w * _h;
    for(unsigned int idx = 0; idx < len; ++idx, ++ptf){
      (*ptf) = (float)(log(*ptf)/log_base);
    }
  }
  else if(_imgType == DImage_dbl_multi){
    double *ptd = (double*)(this->pData); // temp pointer
    len = _numChan * _w * _h;
    for(unsigned int idx = 0; idx < len; ++idx, ++ptd){
      (*ptd) = (log(*ptd)/log_base);
    }
  }
  else{
    fprintf(stderr, "DImage::setDataRangeLog() only supports float and double"
	    " images.\n");
    abort();
  }
}

///returns the max and min values found in the image
void DImage::getDataRange(double *min, double *max){
  int len;
  switch(_imgType){
    case DImage_u8:
    case DImage_RGB:
      {
	D_uint8 *pu8;
	pu8 = dataPointer_u8();
	len = _numChan * _w * _h;
	(*min) = (*max) = (double)pu8[0];
	for(int idx = 0; idx < len; ++idx){
	  if(pu8[idx] < (*min))
	    (*min) = pu8[idx];
	  if(pu8[idx] > (*max))
	    (*max) = pu8[idx];
	}
      }
      break;
    case DImage_u16:
    case DImage_RGB_16:
      {
	D_uint16 *pu16;
	pu16 = dataPointer_u16();
	len = _numChan * _w * _h;
	(*min) = (*max) = (double)pu16[0];
	for(int idx = 0; idx < len; ++idx){
	  if(pu16[idx] < (*min))
	    (*min) = pu16[idx];
	  if(pu16[idx] > (*max))
	    (*max) = pu16[idx];
	}
      }
      break;
    case DImage_u32:
      {
	D_uint32 *pu32;
	pu32 = dataPointer_u32();
	len = _numChan * _w * _h;
	(*min) = (*max) = (double)pu32[0];
	for(int idx = 0; idx < len; ++idx){
	  if(pu32[idx] < (*min))
	    (*min) = pu32[idx];
	  if(pu32[idx] > (*max))
	    (*max) = pu32[idx];
	}
      }
      break;
    case DImage_flt_multi:
      {
	float *pflt;
	pflt = dataPointer_flt();
	len = _numChan * _w * _h;
	(*min) = (*max) = (double)pflt[0];
	for(int idx = 0; idx < len; ++idx){
	  if(pflt[idx] < (*min))
	    (*min) = pflt[idx];
	  if(pflt[idx] > (*max))
	    (*max) = pflt[idx];
	}
      }
      break;
    case DImage_dbl_multi:
      {
	double *pdbl;
	pdbl = dataPointer_dbl();
	len = _numChan * _w * _h;
	(*min) = (*max) = (double)pdbl[0];
	for(int idx = 0; idx < len; ++idx){
	  if(pdbl[idx] < (*min))
	    (*min) = pdbl[idx];
	  if(pdbl[idx] > (*max))
	    (*max) = pdbl[idx];
	}
      }
      break;
    default:
      fprintf(stderr, "DImage::getDataRange() not implemented for complex\n");
      abort();
  }
}

///Adds a scalar value to each pixel in the image (gives the image a bias)
/** out-of-range values are capped to the range of u8,16,32 data types
 *  so that they do not wrap.
 */
void DImage::addValueToPixels(int val){
  size_t len;
  double newVal;
  len = _numChan * _w * _h;
  switch(this->_imgType){
    case DImage_u8:
    case DImage_RGB:
      {
	D_uint8 *pu8 = pData;
	for(unsigned int idx = 0; idx < len; ++idx, ++pu8){
	  newVal = (*pu8)+val;
	  if(newVal < 0)
	    newVal = 0;
	  if(newVal > 255)
	    newVal = 255;
	  (*pu8) = (D_uint8)newVal;
	}
      }
      break;
    case DImage_u16:
    case DImage_RGB_16:
      {
	D_uint16 *pu16 = (D_uint16*)pData;
	for(unsigned int idx = 0; idx < len; ++idx, ++pu16){
	  newVal = (*pu16)+val;
	  if(newVal < 0)
	    newVal = 0;
	  if(newVal > 65535)
	    newVal = 65535;
	  (*pu16) = (D_uint16)newVal;
	}
      }
      break;
    case DImage_u32:
      {
	D_uint32 *pu32 = (D_uint32*)pData;
	for(unsigned int idx = 0; idx < len; ++idx, ++pu32){
	  newVal = (*pu32)+val;
	  if(newVal < 0)
	    newVal = 0;
	  if(newVal > (unsigned int)0xffffffff)
	    newVal = 0xffffffff;
	  (*pu32) = (D_uint32)newVal;
	}
      }
      break;
    case DImage_flt_multi:
    case DImage_dbl_multi:
    case DImage_cmplx:
      addValueToPixels((double)val);
      break;
    default:
      fprintf(stderr, "DImage::addValueToPixels() unknown image type\n");
      abort();
  }
}
void DImage::addValueToPixels(double val){
  size_t len;

  switch(this->_imgType){
    case DImage_u8:
    case DImage_u16:
    case DImage_u32:
    case DImage_RGB:
    case DImage_RGB_16:
#ifdef DEBUG
      if(val != (double)((int)val)){
	fprintf(stderr, "DImage::addValueToPixels() floating point value "
		"for integer-type image\n");
      }
#endif
      addValueToPixels((int)val);
      break;
    case DImage_flt_multi:
      {
	float *pflt = (float*)pData;
	len = _numChan * _w * _h;
	for(unsigned int idx = 0; idx < len; ++idx, ++pflt){
	  (*pflt) += (float)val;
	}
      }
      break;
    case DImage_dbl_multi:
      {
	double *pdbl = (double*)pData;
	len = _numChan * _w * _h;
	for(unsigned int idx = 0; idx < len; ++idx, ++pdbl){
	  (*pdbl) += val;
	}
      }
      break;
    case DImage_cmplx:
      fprintf(stderr, "DImage::addValueToPixels() doesn't support complex\n");
      abort();
      break;
    default:
      fprintf(stderr, "DImage::addValueToPixels() unsupported image type\n");
      abort();
  }
}
///Multiply a scalar value by each pixel in the image
/** out-of-range values are capped to the range of u8,16,32 data types
 *  so that they do not wrap.
 */
void DImage::multiplyPixelsByValue(double val){
  size_t len;
  double newVal;
  len = _numChan * _w * _h;
  switch(this->_imgType){
    case DImage_u8:
    case DImage_RGB:
      {
	D_uint8 *pu8 = pData;
#ifdef DEBUG
	if(val != (double)((int)val)){
	  fprintf(stderr, "DImage::multiplyPixelsByValue() floating point "
		  "value for integer-type image\n");
	}
#endif
	for(unsigned int idx = 0; idx < len; ++idx, ++pu8){
	  newVal = (*pu8)*val;
	  if(newVal < 0)
	    newVal = 0;
	  if(newVal > 255)
	    newVal = 255;
	  (*pu8) = (D_uint8)newVal;
	}
      }
      break;
    case DImage_u16:
    case DImage_RGB_16:
      {
	D_uint16 *pu16 = (D_uint16*)pData;
#ifdef DEBUG
	if(val != (double)((int)val)){
	  fprintf(stderr, "DImage::multiplyPixelsByValue() floating point "
		  "value for integer-type image\n");
	}
#endif
	for(unsigned int idx = 0; idx < len; ++idx, ++pu16){
	  newVal = (*pu16)*val;
	  if(newVal < 0)
	    newVal = 0;
	  if(newVal > 65535)
	    newVal = 65535;
	  (*pu16) = (D_uint16)newVal;
	}
      }
      break;
    case DImage_u32:
      {
	D_uint32 *pu32 = (D_uint32*)pData;
#ifdef DEBUG
	if(val != (double)((int)val)){
	  fprintf(stderr, "DImage::multiplyPixelsByValue() floating point "
		  "value for integer-type image\n");
	}
#endif
	for(unsigned int idx = 0; idx < len; ++idx, ++pu32){
	  newVal = (*pu32)*val;
	  if(newVal < 0)
	    newVal = 0;
	  if(newVal > (unsigned int)0xffffffff)
	    newVal = 0xffffffff;
	  (*pu32) = (D_uint32)newVal;
	}
      }
      break;
    case DImage_flt_multi:
      {
	float *pflt = (float*)pData;
	for(unsigned int idx = 0; idx < len; ++idx, ++pflt){
	  (*pflt) = (float)((*pflt)*val);
	}
      }
      break;
    case DImage_dbl_multi:
      {
	double *pdbl = (double*)pData;
	for(unsigned int idx = 0; idx < len; ++idx, ++pdbl){
	  (*pdbl) = ((*pdbl)*val);
	}
      }
      break;
    case DImage_cmplx:
      {
	std::complex<double> *pcmplx = (std::complex<double>*)pData;
	for(unsigned int idx = 0; idx < len; ++idx, ++pcmplx){
	  (*pcmplx) = ((*pcmplx)*val);
	}
      }
      break;
    default:
      fprintf(stderr, "DImage::multiplyPixelsByValue() unknown image type\n");
      abort();
  }
}

///Divide each pixel in the image by a scalar value
/** out-of-range values are capped to the range of u8,16,32 data types
 *  so that they do not wrap.
 */
void DImage::dividePixelsByValue(double val){
  size_t len;
  double newVal;
  len = _numChan * _w * _h;
  switch(this->_imgType){
    case DImage_u8:
    case DImage_RGB:
      {
	D_uint8 *pu8 = pData;
#ifdef DEBUG
	if(val != (double)((int)val)){
	  fprintf(stderr, "DImage::dividePixelsByValue() floating point value "
		  "for integer-type image\n");
	}
	if(0 == ((D_uint8)val)){
	  fprintf(stderr, "DImage::dividePixelsByValue() divide by zero!\n");
#ifdef DEBUG
	  abort();
#endif      
	}
#endif
	for(unsigned int idx = 0; idx < len; ++idx, ++pu8){
	  newVal = (*pu8)/val;
	  if(newVal < 0)
	    newVal = 0;
	  if(newVal > 255)
	    newVal = 255;
	  (*pu8) = (D_uint8)newVal;
	}
      }
      break;
    case DImage_u16:
    case DImage_RGB_16:
      {
	D_uint16 *pu16 = (D_uint16*)pData;
#ifdef DEBUG
	if(val != (double)((int)val)){
	  fprintf(stderr, "DImage::dividePixelsByValue() floating point value "
		  "for integer-type image\n");
	}
	if(0 == ((D_uint16)val)){
	  fprintf(stderr, "DImage::dividePixelsByValue() divide by zero!\n");
	  abort();
	}
#endif
	for(unsigned int idx = 0; idx < len; ++idx, ++pu16){
	  newVal = (*pu16)/val;
	  if(newVal < 0)
	    newVal = 0;
	  if(newVal > 65535)
	    newVal = 65535;
	  (*pu16) = (D_uint16)newVal;
	}
      }
      break;
    case DImage_u32:
      {
	D_uint32 *pu32 = (D_uint32*)pData;
#ifdef DEBUG
	if(val != (double)((int)val)){
	  fprintf(stderr, "DImage::dividePixelsByValue() floating point value "
		  "for integer-type image\n");
	}
	if(0 == ((D_uint32)val)){
	  fprintf(stderr, "DImage::dividePixelsByValue() divide by zero!\n");
	  abort();
	}
#endif
	for(unsigned int idx = 0; idx < len; ++idx, ++pu32){
	  newVal = (*pu32)/val;
	  if(newVal < 0)
	    newVal = 0;
	  if(newVal > (unsigned int)0xffffffff)
	    newVal = 0xffffffff;
	  (*pu32) = (D_uint32)newVal;
	}
      }
      break;
    case DImage_flt_multi:
      {
	float *pflt = (float*)pData;
#ifdef DEBUG
	if(0. == ((float)val)){
	  fprintf(stderr, "DImage::dividePixelsByValue() divide by zero!\n");
	  abort();
	}
#endif
	for(unsigned int idx = 0; idx < len; ++idx, ++pflt){
	  (*pflt) = (float)((*pflt)/val);
	}
      }
      break;
    case DImage_dbl_multi:
      {
	double *pdbl = (double*)pData;
#ifdef DEBUG
	if(0. == val){
	  fprintf(stderr, "DImage::dividePixelsByValue() divide by zero!\n");
	  abort();
	}
#endif
	for(unsigned int idx = 0; idx < len; ++idx, ++pdbl){
	  (*pdbl) = ((*pdbl)/val);
	}
      }
      break;
    case DImage_cmplx:
      {
	std::complex<double> *pcmplx = (std::complex<double>*)pData;
#ifdef DEBUG
	if(0. == val){
	  fprintf(stderr, "DImage::dividePixelsByValue() divide by zero!\n");
	  abort();
	}
#endif
	for(unsigned int idx = 0; idx < len; ++idx, ++pcmplx){
	  (*pcmplx) = ((*pcmplx)/val);
	}
      }
      break;
    default:
      fprintf(stderr, "DImage::dividePixelsByValue() unknown image type\n");
      abort();
  }
}

//draws a pixel into the grayscale image data buffer with transparency
/** If you don't need transparency, just use setPixel.  This
    function is like the drawLine function, but draws a pixel instead.*/
void DImage::drawPixel(int x, int y, int GS, float transparency){
  float opacity; // opacity fraction (1. - transparency)
  D_uint8 *pu8;
  D_uint16 *pu16;
  D_uint32 *pu32;
  int idx;
#ifdef DEBUG
  if((transparency < 0.) || (transparency > 1.0))
    fprintf(stderr, "DImage::drawPixel() transparency should be [0..1]\n");
#endif

  if((this->_imgType == DImage_RGB)||(this->_imgType == DImage_RGB_16)){
#ifdef DEBUG
    fprintf(stderr, "DImage::drawPixel() grayscale called with color image\n");
#endif
    drawPixel(x,y,GS,GS,GS,transparency);
    return;
  }
  if((_imgType == DImage_flt_multi)||(_imgType == DImage_dbl_multi)||
     (_imgType == DImage_cmplx)){
    fprintf(stderr, "DImage::drawPixel() unsupported image type\n");
    abort();
    return;
  }

  if( (x < 0) || (y < 0) || (x>=_w) || (y >=_h)){
    fprintf(stderr, "DImage::drawPixel(%d,%d) out of bounds\n",x,y);
    return;
#ifdef DEBUG
    //    abort();
#endif
    return;
  }
  opacity = 1. - transparency;

  switch(_imgType){
    case DImage_u8:
      pu8 = pData;
      idx = (y*_w+x);
      //I don't bother checking for value wraparound caused by rounding error
      pu8[idx] = (D_uint8)(GS * opacity + pu8[idx] * transparency);
      break;
    case DImage_u16:
      pu16 = (D_uint16*)pData;
      idx = (y*_w+x);
      pu16[idx] = (D_uint16)(GS * opacity + pu16[idx] * transparency);
      break;
    case DImage_u32:
      pu32 = (D_uint32*)pData;
      idx = (y*_w+x);
      pu32[idx] = (D_uint32)(GS * opacity + pu32[idx] * transparency);
      break;
    default:
      fprintf(stderr, "DImage::drawPixel() unsupported image format\n");
      abort();
  }

}
//draws a pixel into the RGB image data buffer with transparency
/** If you don't need transparency, just use setPixel.  This
    function is like the drawLine function, but draws a pixel instead.*/
void DImage::drawPixel(int x, int y, int R, int G, int B, float transparency){
  float opacity; // opacity fraction (1. - transparency)
  D_uint8 *pu8;
  D_uint16 *pu16;
  int idx = 0;

#ifdef DEBUG
  if((transparency < 0.) || (transparency > 1.0)){
    fprintf(stderr, "DImage::drawPixel() transparency should be [0..1]\n");
    abort();
  }
#endif

  if((this->_imgType != DImage_RGB)&&(this->_imgType != DImage_RGB_16)){
#ifdef DEBUG
    fprintf(stderr, "DImage::drawPixel() rgb called with non-RGB image\n");
    abort();
#endif
    return;
  }

  if( (x < 0) || (y < 0) || (x>=_w) || (y >=_h)){
    fprintf(stderr, "DImage::drawPixel(%d,%d) out of bounds\n",x,y);
    return;
    //    abort();
  }

  opacity = 1. - transparency;

  switch(_imgType){
    case DImage_RGB:
      pu8 = pData;
      idx = (y*_w+x)*3;
      pu8[idx] = (D_uint8)(R * opacity + pu8[idx] * transparency);
      pu8[idx+1] = (D_uint8)(G * opacity + pu8[idx+1] * transparency);
      pu8[idx+2] = (D_uint8)(B * opacity + pu8[idx+2] * transparency);
      break;
    case DImage_RGB_16:
      pu16 = (D_uint16*)pData;
      idx = (y*_w+x)*3;
      pu16[idx] = (D_uint16)(R * opacity + pu16[idx] * transparency);
      pu16[idx+1] = (D_uint16)(G * opacity + pu16[idx+1] * transparency);
      pu16[idx+2] = (D_uint16)(B * opacity + pu16[idx+2] * transparency);
      break;
    default:
      fprintf(stderr, "DImage::drawLine() unsupported image format\n");
      abort();
  }

}


// draws a line into the image data buffer.  (not anti-aliased).
/** Supports only grayscale or indexed (paletized) RGB images.  This
 *  function is meant as a convenience function for drawing rough
 *  lines into an image for debug or display purposes.  It is not
 *  implemented efficiently nor nicely.  Do not use this function for
 *  more advanced drawing needs. For completely opaque lines (the
 *  default), transparency=0. For perfectly transparent, transparency = 1.0.
 */
void DImage::drawLine(int x0, int y0, int x1, int y1, int GS,
		      float transparency, bool fIgnoreClipWarning){
  double dx, dy, incrX, incrY;
  double dist;
  double dblX, dblY;
  int x, y, i, idx = 0;
  int tmp;
  float opacity; // opacity fraction (1. - transparency)
  D_uint8 *pu8;
  D_uint16 *pu16;
  D_uint32 *pu32;

#ifdef DEBUG
  if((transparency < 0.) || (transparency > 1.0))
    fprintf(stderr, "DImage::drawLine() transparency should be [0..1]\n");
#endif

  if((this->_imgType == DImage_RGB)||(this->_imgType == DImage_RGB_16)){
#ifdef DEBUG
    fprintf(stderr, "DImage::drawLine() grayscale called with color image\n");
#endif
    drawLine(x0,y0,x1,y1,GS,GS,GS,transparency);
    return;
  }
  if((_imgType == DImage_flt_multi)||(_imgType == DImage_dbl_multi)||
     (_imgType == DImage_cmplx)){
    fprintf(stderr, "DImage::drawLine() unsupported image type\n");
    abort();
    return;
  }

  if( (x0 < 0) || (x1 < 0) ||  (y0 < 0) || (y1 < 0) ||
      (x0>=_w) || (x1>=_w) ||  (y0 >=_h) || (y1 >=_h)){
    if(!fIgnoreClipWarning)
      fprintf(stderr, "DImage::drawLine(%d,%d,%d,%d) line clipping enabled (but slow implementation)\n",x0,y0,x1,y1);
    // return;
#ifdef DEBUG
    //    abort();
#endif
  }

  opacity = 1. - transparency;

  if(x1 < x0){
    tmp = x1;
    x1 = x0;
    x0 = tmp;
    tmp = y1;
    y1 = y0;
    y0 = tmp;
  }
  
  dx = x1-x0;
  dy = y1-y0;
  dist = sqrt(dx*dx + dy*dy);
  incrX = dx / dist;
  incrY = dy / dist;
  dblX = x0;
  dblY = y0;

  switch(_imgType){
    case DImage_u8:
      pu8 = pData;
      for(i = 0; i <= dist; ++i){
	x = (int)(dblX+0.5); // round to nearest pixel
	y = (int)(dblY+0.5);
	if((x>=0)&&(x<_w)&&(y>=0)&&(y<_h)){
	  idx = (y*_w+x);
	  //I don't bother checking for value wraparound caused by rounding error
	  pu8[idx] = (D_uint8)(GS * opacity + pu8[idx] * transparency);
	}
	dblX += incrX;
	dblY += incrY;
      }
      if(idx != (y1*_w+x1)){//missed last point due to rounding
	if((x1>=0)&&(x1<_w)&&(y1>=0)&&(y1<_h)){
	  idx = (y1*_w+x1);
	  pu8[idx] = (D_uint8)(GS * opacity + pu8[idx] * transparency);
	}
      }
      break;
    case DImage_u16:
      pu16 = (D_uint16*)pData;
      for(i = 0; i <= dist; ++i){
	x = (int)(dblX+0.5);
	y = (int)(dblY+0.5);
	if((x>=0)&&(x<_w)&&(y>=0)&&(y<_h)){
	  idx = (y*_w+x);
	  pu16[idx] = (D_uint16)(GS * opacity + pu16[idx] * transparency);
	}
	dblX += incrX;
	dblY += incrY;
      }
      if(idx != (y1*_w+x1)){//missed last point due to rounding
	if((x1>=0)&&(x1<_w)&&(y1>=0)&&(y1<_h)){
	  idx = (y1*_w+x1);
	  pu16[idx] = (D_uint16)(GS * opacity + pu16[idx] * transparency);
	}
      }
      break;
    case DImage_u32:
      pu32 = (D_uint32*)pData;
      for(i = 0; i <= dist; ++i){
	x = (int)(dblX+0.5);
	y = (int)(dblY+0.5);
	if((x>=0)&&(x<_w)&&(y>=0)&&(y<_h)){
	  idx = (y*_w+x);
	  pu32[idx] = (D_uint32)(GS * opacity + pu32[idx] * transparency);
	}
	dblX += incrX;
	dblY += incrY;
      }
      if(idx != (y1*_w+x1)){//missed last point due to rounding
	if((x1>=0)&&(x1<_w)&&(y1>=0)&&(y1<_h)){
	  idx = (y1*_w+x1);
	  pu32[idx] = (D_uint32)(GS * opacity + pu32[idx] * transparency);
	}
      }
      break;
    default:
      fprintf(stderr, "DImage::drawLine() unsupported image format\n");
      abort();
  }

}
void DImage::drawLine(int x0, int y0, int x1, int y1, int R, int G, int B,
		      float transparency, bool fIgnoreClipWarning){
  double dx, dy, incrX, incrY;
  double dist;
  double dblX, dblY;
  int x, y, i, idx = 0;
  int tmp;
  float opacity; // opacity fraction (1. - transparency)
  D_uint8 *pu8;
  D_uint16 *pu16;

#ifdef DEBUG
  if((transparency < 0.) || (transparency > 1.0)){
    fprintf(stderr, "DImage::drawLine() transparency should be [0..1]\n");
    abort();
  }
#endif

  if((this->_imgType != DImage_RGB)&&(this->_imgType != DImage_RGB_16)){
#ifdef DEBUG
    fprintf(stderr, "DImage::drawLine() rgb called with non-RGB image\n");
    abort();
#endif
    return;
  }

  if( (x0 < 0) || (x1 < 0) ||  (y0 < 0) || (y1 < 0) ||
      (x0>=_w) || (x1>=_w) ||  (y0 >=_h) || (y1 >=_h)){
    if(!fIgnoreClipWarning)
      fprintf(stderr, "DImage::drawLine(%d,%d,%d,%d) line clipping enabled (but slow implementation)\n",x0,y0,x1,y1);
    // return;
    //    abort();
  }

  opacity = 1. - transparency;

  if(x1 < x0){
    tmp = x1;
    x1 = x0;
    x0 = tmp;
    tmp = y1;
    y1 = y0;
    y0 = tmp;
  }
  
  dx = x1-x0;
  dy = y1-y0;
  dist = sqrt(dx*dx + dy*dy);
  incrX = dx / dist;
  incrY = dy / dist;
  dblX = x0;
  dblY = y0;

  switch(_imgType){
    case DImage_RGB:
      pu8 = pData;
      for(i = 0; i <= dist; ++i){
	x = (int)(dblX+0.5);
	y = (int)(dblY+0.5);
	if((x>=0)&&(x<_w)&&(y>=0)&&(y<_h)){
	  idx = (y*_w+x)*3;
	  // I won't bother checking for wrap-around due to rounding error
	  pu8[idx] = (D_uint8)(R * opacity + pu8[idx] * transparency);
	  pu8[idx+1] = (D_uint8)(G * opacity + pu8[idx+1] * transparency);
	  pu8[idx+2] = (D_uint8)(B * opacity + pu8[idx+2] * transparency);
	}
	dblX += incrX;
	dblY += incrY;
      }
      if(idx != (y1*_w+x1)*3){//missed last point due to rounding
	if((x1>=0)&&(x1<_w)&&(y1>=0)&&(y1<_h)){
	  idx = (y1*_w+x1)*3;
	  pu8[idx] = (D_uint8)(R * opacity + pu8[idx] * transparency);
	  pu8[idx+1] = (D_uint8)(G * opacity + pu8[idx+1] * transparency);
	  pu8[idx+2] = (D_uint8)(B * opacity + pu8[idx+2] * transparency);
	}
      }
      break;
    case DImage_RGB_16:
      pu16 = (D_uint16*)pData;
      for(i = 0; i <= dist; ++i){
	x = (int)(dblX+0.5);
	y = (int)(dblY+0.5);
	if((x>=0)&&(x<_w)&&(y>=0)&&(y<_h)){
	  idx = (y*_w+x)*3;
	  pu16[idx] = (D_uint16)(R * opacity + pu16[idx] * transparency);
	  pu16[idx+1] = (D_uint16)(G * opacity + pu16[idx+1] * transparency);
	  pu16[idx+2] = (D_uint16)(B * opacity + pu16[idx+2] * transparency);
	}
	dblX += incrX;
	dblY += incrY;
      }
      if(idx != (y1*_w+x1)*3){//missed last point due to rounding
	if((x1>=0)&&(x1<_w)&&(y1>=0)&&(y1<_h)){
	  idx = (y1*_w+x1)*3;
	  pu16[idx] = (D_uint16)(R * opacity + pu16[idx] * transparency);
	  pu16[idx+1] = (D_uint16)(G * opacity + pu16[idx+1] * transparency);
	  pu16[idx+2] = (D_uint16)(B * opacity + pu16[idx+2] * transparency);
	}
      }
      break;
    default:
      fprintf(stderr, "DImage::drawLine() unsupported image format\n");
      abort();
  }

}

///Draws rectangle with opposite corners at x0,y0 and x1,y1 (not width,height)
void DImage::drawRect(int x0, int y0, int x1, int y1, int GS,
		      float transparency, bool fill){
  int tmp;
  // fix coordinates so x0,y0 is the upper-left and x1,y1 is bottom-right
  if(x1 < x0){
    tmp = x0;
    x0 = x1;
    x1 = tmp;
  }
  if(y1 < y0){
    tmp = y0;
    y0 = y1;
    y1 = tmp;
  }
  // now draw the rect (either filled or just the boundary pixels)
  if(fill){
    for(int y = y0; y <= y1; ++y){
      drawLine(x0, y, x1, y, GS, transparency);
    }
  }
  else{
    drawLine(x0, y0, x1, y0, GS, transparency);
    if(y1 != y0){
      drawLine(x1, y0+1, x1, y1-1, GS, transparency);
      drawLine(x1, y1, x0, y1, GS, transparency);
      drawLine(x0, y1-1, x0, y0+1, GS, transparency);
    }
  }
}
void DImage::drawRect(int x0, int y0, int x1, int y1, int R, int G, int B,
		      float transparency, bool fill){
  int tmp;
  // fix coordinates so x0,y0 is the upper-left and x1,y1 is bottom-right
  if(x1 < x0){
    tmp = x0;
    x0 = x1;
    x1 = tmp;
  }
  if(y1 < y0){
    tmp = y0;
    y0 = y1;
    y1 = tmp;
  }
  // now draw the rect (either filled or just the boundary pixels)
  if(fill){
    for(int y = y0; y <= y1; ++y){
      drawLine(x0, y, x1, y, R, G, B, transparency);
    }
  }
  else{
    drawLine(x0, y0, x1, y0, R, G, B, transparency);
    if(y1 != y0){
      drawLine(x1, y0+1, x1, y1-1, R, G, B, transparency);
      drawLine(x1, y1, x0, y1, R, G, B, transparency);
      drawLine(x0, y1-1, x0, y0+1, R, G, B, transparency);
    }
  }
}

///Changes w,h of image but only reallocates if the buffer needs to be bigger
/**The width and height of the image are changed to w and h, but the
 * data buffer, itself, is not reallocated unless w*h is greater than
 * it was when the buffer was originally allocated.  The contents of
 * the data buffer are left as is, so changing the size leaves the
 * image data in an undesirable state.  The image data type, number of
 * channels, etc. do not change.  If the image is copied after the
 * logical size has been changed, the copy will have the same actual
 * size as well as the same logical size as the original.
 */
void DImage::setLogicalSize(int w, int h){
  if( (w*h) > (_actualW*_actualH)){ // need to expand the buffer
    create(w,h,_imgType, _numChan,_allocMethod);
  }
  else{ // don't need to expand the buffer, just keep using the same one
    _w = w;
    _h = h;
  }
}


///Returns a resized version of this image, with any pad area set to pad_val
/**The width and height of the image are changed, but it is not scaled.  (The
 * image is clipped to the new width/height or padded to the new width/height
 * but the pixel data is otherwise unchanged.)
 */
DImage DImage::resized(int newW, int newH, double pad_val) const{
  DImage dst;
  resized_(dst, newW, newH, pad_val);
  return dst;
}

///Puts resized version of this image into imgDst, with pad area set to pad_val
/**The width and height of the image are changed, but it is not scaled.  (The
 * image is clipped to the new width/height or padded to the new width/height
 * but the pixel data is otherwise unchanged.)
 */
void DImage::resized_(DImage &imgDst, int newW, int newH, double pad_val)const{
  size_t dstRowLen;
  size_t cpyLen;
  int cpyChan;
  size_t padLen;
  int cpyH;
  size_t srcRowLen;
  D_uint8 *srcTmp;
  D_uint8 *dstTmp;
  

#ifdef DEBUG
  if((newW < 0)||(newH < 0)){
    fprintf(stderr, "DImage::resized_() w,h==%d,%d must be >= 0\n", newW,newH);
    abort();
  }
  if(this == &imgDst){
    fprintf(stderr, "DImage::resized_() source (this) == dest!\n");
    abort();
    return;
  }
  if(NULL == this->pData){
    fprintf(stderr, "DImage::resized_() this->pData == NULL\n");
    abort();
    return;
  }
#endif
  
  imgDst.create(newW,newH,this->_imgType, this->_numChan, this->_allocMethod);
  imgDst.mapProps.clear(); // clear any previous properties- they're not valid
  imgDst.vectComments.clear(); // clear any previous comments for same reason


  dstRowLen = newW * _sampleSize;
  srcRowLen = _w * _sampleSize;
  if(newW <= _w){
    cpyLen = newW * _sampleSize;
    padLen = 0;
  }
  else{
    cpyLen = _w * _sampleSize;
    padLen = (newW - _w); // number of samples, not number of bytes
  }
  cpyH = (newH <= _h) ? newH : _h;
  cpyChan = _numChan;


  if(_fInterleaved){
    cpyLen *= _numChan;
    padLen *= _numChan;
    dstRowLen *= _numChan;
    srcRowLen *= _numChan;
    cpyChan = 1;
  }
  

  for(int chan = 0; chan < cpyChan; ++chan){
    srcTmp = &(pData[chan * srcRowLen * _h]); // offset to channel being copied
    dstTmp = &(imgDst.pData[chan * dstRowLen * newH]); // offs to channel (dst)
    // copy main image data from this to destination image
    for(int y = 0; y < cpyH; ++y){
      memcpy(&(dstTmp[dstRowLen*y]), &(srcTmp[srcRowLen*y]), cpyLen);
    }
    // pad the side
    if(padLen > 0){
      if(pad_val == 0.){ // can use memset since it's zero
	for(int y = 0; y < cpyH; ++y){
	  memset(&(dstTmp[dstRowLen*y+cpyLen]), 0, padLen*_sampleSize);
	}
      }
      else{
	switch(_imgType){
          case DImage_u8:
          case DImage_RGB:
	    for(int y = 0; y < cpyH; ++y){
	      memset(&(dstTmp[dstRowLen*y+cpyLen]), (int)pad_val, padLen);
	    }
	    break;
          case DImage_u16:
          case DImage_RGB_16:
	    for(int y = 0; y < cpyH; ++y){
	      D_uint16 *pu16;
	      pu16 = (D_uint16*)&(dstTmp[dstRowLen*y+cpyLen]);
	      for(unsigned int x=0; x < padLen; ++x, ++pu16){
		(*pu16) =  (D_uint16)pad_val;
	      }
	    }
	    break;
          case DImage_u32:
	    for(int y = 0; y < cpyH; ++y){
	      D_uint32 *pu32;
	      pu32 = (D_uint32*)&(dstTmp[dstRowLen*y+cpyLen]);
	      for(unsigned int x=0; x < padLen; ++x, ++pu32){
		(*pu32) = (D_uint32)pad_val;
	      }
	    }
	    break;
          case DImage_cmplx:
	    for(int y = 0; y < cpyH; ++y){
	      std::complex<double> *pcmplx;
	      pcmplx = (std::complex<double>*)&(dstTmp[dstRowLen*y+cpyLen]);
	      for(unsigned int x=0; x < padLen; ++x, ++pcmplx){
		(*pcmplx) = (std::complex<double>)pad_val;
	      }
	    }
	    break;
          case DImage_flt_multi:
	    for(int y = 0; y < cpyH; ++y){
	      float *pflt;
	      pflt = (float*)&(dstTmp[dstRowLen*y+cpyLen]);
	      for(unsigned int x=0; x < padLen; ++x, ++pflt){
		(*pflt) = (float)pad_val;
	      }
	    }
	    break;
          case DImage_dbl_multi:
	    for(int y = 0; y < cpyH; ++y){
	      double *pdbl;
	      pdbl = (double*)&(dstTmp[dstRowLen*y+cpyLen]);
	      for(unsigned int x=0; x < padLen; ++x, ++pdbl){
		(*pdbl) = (double)pad_val;
	      }
	    }
	    break;
	  default:
	    fprintf(stderr, "DImage::resized() unsupported image type\n");
	    abort();
	}// end switch(_imgType)
      }// end else (pad_val != 0.)
    } // end if(padLen > 0)
    
    // pad the bottom
    if(newH > _h){
      padLen = (newH - _h) * newW;
      if(pad_val == 0.){ // can use memset since it's zero
	memset(&(dstTmp[_h*dstRowLen]), 0, padLen*_sampleSize);
      }
      else{ // pad_val != 0
	switch(_imgType){
          case DImage_u8:
          case DImage_RGB:
	    memset(&(dstTmp[_h*dstRowLen]), (int)pad_val, padLen*_sampleSize);
	    break;
          case DImage_u16:
          case DImage_RGB_16:
	    D_uint16 *pu16;
	    pu16 = (D_uint16*)&(dstTmp[_h*dstRowLen]);
	    for(unsigned int i = 0; i < padLen; ++i, ++pu16){
	      (*pu16) =  (D_uint16)pad_val;
	    }
	    break;
          case DImage_u32:
	    D_uint32 *pu32;
	    pu32 = (D_uint32*)&(dstTmp[_h*dstRowLen]);
	    for(unsigned int i = 0; i < padLen; ++i, ++pu32){
	      (*pu32) =  (D_uint32)pad_val;
	    }
	    break;
          case DImage_cmplx:
	    std::complex<double> *pcmplx;
	    pcmplx = (std::complex<double>*)&(dstTmp[_h*dstRowLen]);
	    for(unsigned int i = 0; i < padLen; ++i, ++pcmplx){
	      (*pcmplx) =  (std::complex<double>)pad_val;
	    }
	    break;
          case DImage_flt_multi:
	    float *pflt;
	    pflt = (float*)&(dstTmp[_h*dstRowLen]);
	    for(unsigned int i = 0; i < padLen; ++i, ++pflt){
	      (*pflt) =  (float)pad_val;
	    }
	    break;
          case DImage_dbl_multi:
	    double *pdbl;
	    pdbl = (double*)&(dstTmp[_h*dstRowLen]);
	    for(unsigned int i = 0; i < padLen; ++i, ++pdbl){
	      (*pdbl) =  (double)pad_val;
	    }
	    break;
	  default:
	    fprintf(stderr, "DImage::resized() unsupported image type\n");
	    abort();
	}// end switch(_imgType)
	
      }// end else pad_val != 0
    }// end if(newH > _h)
  }// end for(int chan=...)
}
///Returns a resized version of this image, with pad area set to given RGB vals
/**The width and height of the image are changed, but it is not scaled.  (The
 * image is clipped to the new width/height or padded to the new width/height
 * but the pixel data is otherwise unchanged.)
 */
DImage DImage::resized(int newW, int newH,
		       int pad_R, int pad_G, int pad_B) const{
  DImage dst;
  resized_(dst, newW, newH, pad_R, pad_G, pad_B);
  return dst;
}

///Puts resized version of this image into imgDst, with pad set as specified
/**The width and height of the image are changed, but it is not scaled.  (The
 * image is clipped to the new width/height or padded to the new width/height
 * but the pixel data is otherwise unchanged.)
 */
void DImage::resized_(DImage &imgDst, int newW, int newH,
		      int pad_R, int pad_G, int pad_B)const{
  size_t dstRowLen;
  size_t cpyLen;
  size_t padLen;
  int cpyH;
  size_t srcRowLen;
  D_uint8 *srcTmp;
  D_uint8 *dstTmp;
  

#ifdef DEBUG
  if((newW < 0)||(newH < 0)){
    fprintf(stderr, "DImage::resized_() w,h==%d,%d must be >= 0\n", newW,newH);
    abort();
  }
  if(this == &imgDst){
    fprintf(stderr, "DImage::resized_() source (this) == dest!\n");
    abort();
    return;
  }
  if(NULL == this->pData){
    fprintf(stderr, "DImage::resized_() this->pData == NULL\n");
    abort();
    return;
  }
  if((_imgType != DImage_RGB) && (_imgType != DImage_RGB_16)){
    fprintf(stderr, "DImage::resized_() RGB version used for non-RGB image\n");
    abort();
    return;
  }
#endif
  
  imgDst.create(newW,newH,this->_imgType, this->_numChan, this->_allocMethod);
  imgDst.mapProps.clear(); // clear any previous properties- they're not valid
  imgDst.vectComments.clear(); // clear any previous comments for same reason


  dstRowLen = newW * _sampleSize * _numChan;
  srcRowLen = _w * _sampleSize * _numChan;
  if(newW <= _w){
    cpyLen = newW * _sampleSize * _numChan;
    padLen = 0;
  }
  else{
    cpyLen = _w * _sampleSize * _numChan;
    padLen = (newW - _w); // number of pixels, not # of bytes nor samples
  }
  cpyH = (newH <= _h) ? newH : _h;

  srcTmp = pData; // offset to channel being copied
  dstTmp = imgDst.pData; // offset to channel (dst)
  // copy main image data from this to destination image
  for(int y = 0; y < cpyH; ++y){
    memcpy(&(dstTmp[dstRowLen*y]), &(srcTmp[srcRowLen*y]), cpyLen);
  }
  // pad the side
  if(padLen > 0){
    if((pad_R == 0)&&(pad_G==0)&&(pad_B==0)){// can use memset since all 0's
      for(int y = 0; y < cpyH; ++y){
	memset(&(dstTmp[dstRowLen*y+cpyLen]), 0, padLen*_sampleSize*_numChan);
      }
    }
    else{
      switch(_imgType){
        case DImage_RGB:
	  for(int y = 0; y < cpyH; ++y){
	    D_uint8 *pu8;
	    pu8 = &(dstTmp[dstRowLen*y+cpyLen]);
	    for(unsigned int x = 0; x < padLen; ++x, pu8 += 3){
	      (*pu8) = pad_R;
	      (*(pu8+1)) = pad_G;
	      (*(pu8+2)) = pad_B;
	    }
	    //TODO: instead of putting each value in for whole area, could do
	    //a single row, then use memcpy to replicate the row if it's faster
	  }
	  break;
        case DImage_RGB_16:
	  for(int y = 0; y < cpyH; ++y){
	    D_uint16 *pu16;
	    pu16 = (D_uint16*)&(dstTmp[dstRowLen*y+cpyLen]);
	    for(unsigned int x=0; x < padLen; ++x, pu16 += 3){
	      (*pu16) =  (D_uint16)pad_R;
	      (*(pu16+1)) =  (D_uint16)pad_G;
	      (*(pu16+2)) =  (D_uint16)pad_B;
	    }
	    //TODO: instead of putting each value in for whole area, could do
	    //a single row, then use memcpy to replicate the row if it's faster
	  }
	  break;
        default:
	  fprintf(stderr, "DImage::resized() unsupported image type\n");
	  abort();
      }// end switch(_imgType)
    }// end else (pad_val != 0.)
  } // end if(padLen > 0)
    
  // pad the bottom
  if(newH > _h){
    padLen = (newH - _h) * newW;
    if((pad_R == 0)&&(pad_G==0)&&(pad_B==0)){// can use memset since all 0's
      memset(&(dstTmp[_h*dstRowLen]), 0, padLen*_sampleSize);
    }
    else{ // pad_val != 0
      switch(_imgType){
        case DImage_RGB:
	  D_uint8 *pu8;
	  pu8 = (D_uint8*)&(dstTmp[_h*dstRowLen]);
	  for(unsigned int i = 0; i < padLen; ++i, pu8 += 3){
	    (*pu8) = pad_R;
	    (*(pu8+1)) = pad_G;
	    (*(pu8+2)) = pad_B;
	  }
	  //TODO: instead of putting each value in for whole area, could do
	  //a single row, then use memcpy to replicate the row if it's faster
	  break;
        case DImage_RGB_16:
	  D_uint16 *pu16;
	  pu16 = (D_uint16*)&(dstTmp[_h*dstRowLen]);
	  for(unsigned int i = 0; i < padLen; ++i, pu16 += 3){
	    (*pu16) =  (D_uint16)pad_R;
	    (*(pu16+1)) =  (D_uint16)pad_G;
	    (*(pu16+2)) =  (D_uint16)pad_B;
	  }
	  //TODO: instead of putting each value in for whole area, could do
	  //a single row, then use memcpy to replicate the row if it's faster
	  break;
        default:
	  fprintf(stderr, "DImage::resized() unsupported image type\n");
	  abort();
      }// end switch(_imgType)
    }// end else pad_val != 0
  }// end if(newH > _h)
}


///Returns a version of this image, with extra pad around the edges
DImage DImage::padEdges(int left, int right, int top, int bottom,
			DImagePadMethod padMethod, double val) const{
  DImage dst;
  padEdges_(dst, left, right, top, bottom, padMethod, val);
  return dst;
}

///Add extra pad around the edges of this image and store the result in imgDst
void DImage::padEdges_(DImage &imgDst,int left, int right, int top, int bottom,
		       DImagePadMethod padMethod, double val) const{
  unsigned int chanOffsDst;
  unsigned int chanOffsSrc;
  unsigned int lineLenSrc;
  unsigned int lineLenDst;
  unsigned int sampSz;
  unsigned int numCh;
  D_uint8 *pDst;
  D_uint8 *pSrc;
  imgDst.create(_w+left+right, _h+top+bottom,_imgType,_numChan,_allocMethod);
  imgDst.pasteFromImage(left, top, (*this), 0, 0, _w, _h);


  // now fill the pad according to the padMethod chosen
  switch(padMethod){
    case DImagePadValue:

      switch(_imgType){
        case DImage_u8:
	  {
	    D_uint8 u8Val;
	    u8Val = (D_uint8)val;
	    lineLenDst = imgDst._w;
	    if(top > 0)
	      memset(imgDst.pData, (int)val, lineLenDst * top);
	    pDst = &(imgDst.pData[lineLenDst * top]);
	    for(int y = 0; y < _h; ++y){
	      for(int x = 0; x < left; ++x){
		pDst[x] = u8Val;
	      }
	      pDst += _w+left;
	      for(int x = 0; x < right; ++x){
		pDst[x] = u8Val;
	      }
	      pDst += right;
	    }
	    if(bottom > 0)
	      memset(pDst, (int)val, lineLenDst * bottom);
	  }
	  break;
        case DImage_u16:
	  {
	    D_uint16 u16Val;
	    D_uint16 *pData16;
	    u16Val = (D_uint16)val;
	    pData16 = (D_uint16*)imgDst.pData;
	    lineLenDst = imgDst._w;
	    for(int idx = 0, len=lineLenDst*top; idx < len; ++idx){
	      pData16[idx] = u16Val;
	    }
	    pData16 =
	      (D_uint16*)&(imgDst.pData[lineLenDst * top * sizeof(D_uint16)]);
	    for(int y = 0; y < _h; ++y){
	      for(int x = 0; x < left; ++x){
		pData16[x] = u16Val;
	      }
	      pData16 += _w + left;
	      for(int x = 0; x < right; ++x){
		pData16[x] = u16Val;
	      }
	      pData16 += right;
	    }
	    for(int idx = 0, len=lineLenDst*bottom; idx < len; ++idx){
	      pData16[idx] = u16Val;
	    }
	  }
	  break;
        case DImage_u32:
	  {
	    D_uint32 u32Val;
	    D_uint32 *pData32;
	    u32Val = (D_uint32)val;
	    pData32 = (D_uint32*)imgDst.pData;
	    lineLenDst = imgDst._w;
	    for(int idx = 0, len=lineLenDst*top; idx < len; ++idx){
	      pData32[idx] = u32Val;
	    }
	    pData32 =
	      (D_uint32*)&(imgDst.pData[lineLenDst * top * sizeof(D_uint32)]);
	    for(int y = 0; y < _h; ++y){
	      for(int x = 0; x < left; ++x){
		pData32[x] = u32Val;
	      }
	      pData32 += _w + left;
	      for(int x = 0; x < right; ++x){
		pData32[x] = u32Val;
	      }
	      pData32 += right;
	    }
	    for(int idx = 0, len=lineLenDst*bottom; idx < len; ++idx){
	      pData32[idx] = u32Val;
	    }
	  }
	  break;
        case DImage_RGB:
	  {
	    D_uint8 u8Val;
	    u8Val = (D_uint8)val;
	    lineLenDst = imgDst._w*3;
	    if(top > 0)
	      memset(imgDst.pData, (int)val, lineLenDst * top);
	    pDst = &(imgDst.pData[lineLenDst * top]);
	    for(int y = 0; y < _h; ++y){
	      for(int x = 0; x < left*3; ++x){
		pDst[x] = u8Val;
	      }
	      pDst += (_w+left)*3;
	      for(int x = 0; x < right*3; ++x){
		pDst[x] = u8Val;
	      }
	      pDst += right*3;
	    }
	    if(bottom > 0)
	      memset(pDst, (int)val, lineLenDst * bottom);
	  }
	  break;
        case DImage_RGB_16:
	  {
	    D_uint16 u16Val;
	    D_uint16 *pData16;
	    u16Val = (D_uint16)val;
	    pData16 = (D_uint16*)imgDst.pData;
	    lineLenDst = imgDst._w * 3;
	    for(int idx = 0, len=lineLenDst*top; idx < len; ++idx){
	      pData16[idx] = u16Val;
	    }
	    pData16 =
	      (D_uint16*)&(imgDst.pData[lineLenDst * top * sizeof(D_uint16)]);
	    for(int y = 0; y < _h; ++y){
	      for(int x = 0; x < left*3; ++x){
		pData16[x] = u16Val;
	      }
	      pData16 += (_w+left)*3;
	      for(int x = 0; x < right*3; ++x){
		pData16[x] = u16Val;
	      }
	      pData16 += right;
	    }
	    for(int idx = 0, len=lineLenDst*bottom; idx < len; ++idx){
	      pData16[idx] = u16Val;
	    }
	  }
	  break;
        case DImage_flt_multi:
	  {
	    float fltVal;
	    float *pDataFlt;
	    fltVal = (float)val;
	    lineLenDst = imgDst._w;

	    for(int chan = 0; chan < _numChan; chan++){
	      pDataFlt = imgDst.dataPointer_flt(chan);
	      for(int idx = 0, len=lineLenDst*top; idx < len; ++idx){
		pDataFlt[idx] = fltVal;
	      }
	      pDataFlt =
		(float*)&(imgDst.pData[lineLenDst * top * sizeof(float)]);
	      for(int y = 0; y < _h; ++y){
		for(int x = 0; x < left; ++x){
		  pDataFlt[x] = fltVal;
		}
		pDataFlt += _w + left;
		for(int x = 0; x < right; ++x){
		  pDataFlt[x] = fltVal;
		}
		pDataFlt += right;
	      }
	      for(int idx = 0, len=lineLenDst*bottom; idx < len; ++idx){
		pDataFlt[idx] = fltVal;
	      }
	    }//end for chan
	  }
	  break;
        case DImage_dbl_multi:
	  {
	    double dblVal;
	    double *pDataDbl;
	    dblVal = (double)val;
	    lineLenDst = imgDst._w;

	    for(int chan = 0; chan < _numChan; chan++){
	      pDataDbl = imgDst.dataPointer_dbl(chan);
	      for(int idx = 0, len=lineLenDst*top; idx < len; ++idx){
		pDataDbl[idx] = dblVal;
	      }
	      pDataDbl =
		(double*)&(imgDst.pData[lineLenDst * top * sizeof(double)]);
	      for(int y = 0; y < _h; ++y){
		for(int x = 0; x < left; ++x){
		  pDataDbl[x] = dblVal;
		}
		pDataDbl += _w + left;
		for(int x = 0; x < right; ++x){
		  pDataDbl[x] = dblVal;
		}
		pDataDbl += right;
	      }
	      for(int idx = 0, len=lineLenDst*bottom; idx < len; ++idx){
		pDataDbl[idx] = dblVal;
	      }
	    }//end for chan
	  }
	  break;
        case DImage_cmplx:
	  {
	    std::complex<double> cmplxVal;
	    std::complex<double> *pDataCmplx;
	    cmplxVal = (std::complex<double>)val;
	    lineLenDst = imgDst._w;

	    pDataCmplx = imgDst.dataPointer_cmplx();
	    for(int idx = 0, len=lineLenDst*top; idx < len; ++idx){
	      pDataCmplx[idx] = cmplxVal;
	    }
	    pDataCmplx = (std::complex<double>*)
	      &(imgDst.pData[lineLenDst *top*sizeof(std::complex<double>)]);
	    for(int y = 0; y < _h; ++y){
	      for(int x = 0; x < left; ++x){
		pDataCmplx[x] = cmplxVal;
	      }
	      pDataCmplx += _w + left;
	      for(int x = 0; x < right; ++x){
		pDataCmplx[x] = cmplxVal;
	      }
	      pDataCmplx += right;
	    }
	    for(int idx = 0, len=lineLenDst*bottom; idx < len; ++idx){
	      pDataCmplx[idx] = cmplxVal;
	    }
	  }
	  break;
      
      }//end switch(_imgType)
      break;
    case DImagePadReplicate:
      numCh = _numChan;
      sampSz = _sampleSize;
      if((DImage_RGB==_imgType)||(DImage_RGB_16==_imgType)){
	numCh = 1;
	sampSz *= 3;
      }
      lineLenSrc = _w*sampSz;
      lineLenDst = imgDst._w*sampSz;
      for(unsigned int chan=0; chan < numCh; ++chan){
	chanOffsDst = chan*imgDst._w*imgDst._h*sampSz;
	//left and right
	pDst = &(imgDst.pData[chanOffsDst+lineLenDst*top]);
	for(int y = 0; y < _h; ++y){
	  pSrc = &(pDst[left*sampSz]);
	  for(int x = 0; x < left; ++x){
	    for(unsigned int j = 0; j < sampSz; ++j){
	      pDst[j] = pSrc[j];
	    }
	    pDst += sampSz;
	  }
	  pSrc = &(pDst[(_w-1)*sampSz]);
	  pDst += sampSz * _w;
	  for(int x = 0; x < right; ++x){
	    for(unsigned int j = 0; j < sampSz; ++j){
	      pDst[j] = pSrc[j];
	    }
	    pDst += sampSz;
	  }
	}//end for y
	//bottom
// 	pSrc = &(imgDst.pData[chanOffsDst+lineLenDst*(top+_w-1)]);
	pSrc = &(imgDst.pData[chanOffsDst+lineLenDst*(top+_h-1)]);
	for(int y = 0; y < bottom; ++y){
	  memcpy(pDst, pSrc, lineLenDst);
	  pDst += lineLenDst;
	}
	//top
	pSrc = &(imgDst.pData[chanOffsDst+lineLenDst*top]);
	pDst = &(imgDst.pData[chanOffsDst]);
	for(int y = 0; y < top; ++y){
	  memcpy(pDst, pSrc, lineLenDst);
	  pDst += lineLenDst;
	}
      }//end for chan
      break;
    case DImagePadMirror:
#ifdef DEBUG
      if((left > _w) || (right > _w) || (top > _h) || (bottom > _h)){
	fprintf(stderr,
		"DImage::padEdges_() requires that top,bottom < height and"
		" left,right < width when using DImagePadMirror or "
		"DImagePadWrap.  You need to pad multiple times or change "
		"padding method.\n");
	abort();
      }
#endif
      fprintf(stderr, "DImage::padEdges_() not yet implemented for "
	      "DImagePadMirror\n");
      abort();
      //TODO: need to code up DImagePadMirror in DImage::padEdges_()
      break;
    case DImagePadWrap:
#ifdef DEBUG
      if((left > _w) || (right > _w) || (top > _h) || (bottom > _h)){
	fprintf(stderr,
		"DImage::padEdges_() requires that top,bottom < height and"
		" left,right < width when using DImagePadMirror or "
		"DImagePadWrap.  You need to pad multiple times or change "
		"padding method.\n");
	abort();
      }
#endif
      numCh = _numChan;
      sampSz = _sampleSize;
      if((DImage_RGB==_imgType)||(DImage_RGB_16==_imgType)){
	numCh = 1;
	sampSz *= 3;
      }
      lineLenSrc = _w*sampSz;
      lineLenDst = imgDst._w*sampSz;
      for(unsigned int chan=0; chan < numCh; ++chan){
	chanOffsSrc = chan*_w*_h*sampSz;
	chanOffsDst = chan*imgDst._w*imgDst._h*sampSz;
	//top
// 	for(int y = 0; y < top; ++y){
// 	  // top-center
// 	  memcpy(&(imgDst.pData[chanOffsDst+y*lineLenDst + left*sampSz]),
// 		 &(pData[chanOffsSrc + (y+_h-top)*lineLenSrc]), lineLenSrc);
// 	  // top-left corner
// 	  if(left > 0){
// 	    memcpy(&(imgDst.pData[chanOffsDst+y*lineLenDst]),
// 		   &(pData[chanOffsSrc + (y+_h-top)*lineLenSrc +
// 			   (_w-left)*sampSz]), left*sampSz);
// 	  }
// 	  // top-right corner
// 	  if(right > 0){
// 	    memcpy(&(imgDst.pData[chanOffsDst+y*lineLenDst+
// 		     (left+_w)*sampSz]),
// 		   &(pData[chanOffsSrc + (y+_h-top)*lineLenSrc]),
// 		   right*sampSz);
// 	  }
// 	}// end for y

	//left and right
	pDst = &(imgDst.pData[chanOffsDst+top*lineLenDst]);
	for(int y = 0; y < _h; ++y){
	  if(left > 0){
	    memcpy(pDst, &(pDst[_w*sampSz]), left*sampSz);
	  }
	  if(right > 0){
	    memcpy(&(pDst[(left+_w)*sampSz]),
		   &(pDst[left*sampSz]), right*sampSz);
	  }
	  pDst += lineLenDst;
	}
	// bottom
	if(bottom > 0){
	  pSrc = &(imgDst.pData[chanOffsDst+top*lineLenDst]);
	  pDst = &(imgDst.pData[chanOffsDst+(top+_h)*lineLenDst]);
	  memcpy(pDst, pSrc, bottom*lineLenDst);
	}
	if(top > 0){
	  pSrc = &(imgDst.pData[chanOffsDst+_h*lineLenDst]);
	  pDst = &(imgDst.pData[chanOffsDst]);
	  memcpy(pDst, pSrc, top*lineLenDst);
	}
      }//end for chan
      break;
    case DImagePadNone:
      // do nothing.  Leave the memory uninitialized
      break;
    default:
      fprintf(stderr, "DImage::padEdges_() unknown padMethod\n");
      abort();
  }

}

//Invert the values of a grayscale image (white for black, black for white)
/**This function inverts the grayscale values of the image data.  For
 * an 8-bit image, any pixel with value 255 goes to 0, any pixel with
 * value 254 goes to 1, and so on. Similarly for 16-bit grayscale.
 * This function can only be used for DImage_u8 and DImage_u16
 * images. */
void DImage::invertGrayscale(){

  if(DImage_u8 == _imgType){
    D_uint8 *p8;
    p8 = pData;
    for(int y=0, idx=0; y < _h; ++y){
      for(int x=0; x < _w; ++x, ++idx){
	p8[idx] = (D_uint8)(255-p8[idx]);
      }
    }
  }
  else if(DImage_u16 == _imgType){
    D_uint16 *p16;
    p16 = (D_uint16*)pData;
    for(int y=0, idx=0; y < _h; ++y){
      for(int x=0; x < _w; ++x, ++idx){
	p16[idx] = (D_uint16)(0xffff-p16[idx]);
      }
    }
  }
  else{
    fprintf(stderr, "DImage::invertGrayscale() only supported for 8-bit and "
	    "16-bit grayscale images\n");
    abort();
  }
}

///Returns a scaled version of this image
DImage DImage::scaled(int dstW, int dstH,
		      DImageTransformMode mode) const{
  DImage dst;
  scaled_(dst, dstW, dstH, mode);
  return dst;
}
///Same as scaled() except it takes scale ratios instead of width, height
DImage DImage::scaled(double scaleX, double scaleY,
		      DImageTransformMode mode) const{
  DImage dst;
  scaled_(dst, scaleX, scaleY, mode);
  return dst;
}

///Special-case scaling when the image needs to be shrunk by a power of 2
/**Scales this image by a factor of 1/(2^numHalves) in each direction.
 * Puts the result in imgDst.  This is called automatically by the regular
 * scaled() and scaled_() functions if appropriate.
 */
void DImage::scaledDownPow2_(DImage &imgDst, int numHalves,
			     DImageTransformMode mode) const{
  int dstW, dstH;
  int dstWm1; // dstW minus 1
  int avgLen; // how many pixels to average in each direction
  int alm1; // avg len minus 1
  int numRows;
  float divisor, divisorRight, divisorBottom, divisorCorner;
  std::complex<double> zdivisor, zdivisorRight, zdivisorBottom, zdivisorCorner;
  int rightPixels;
  int bottomPixels;
  
  avgLen = 1 << numHalves;
  alm1 = avgLen - 1;
  dstW = (_w+alm1) >> numHalves;
  dstWm1 = dstW-1;
  dstH = (_h+alm1) >> numHalves;

  if(DImageTransSmooth == mode){
    divisor = avgLen * avgLen;
    rightPixels = _w % avgLen;
    if(0 == rightPixels)
      rightPixels = avgLen;
    bottomPixels = _h % avgLen;
    if(0 == bottomPixels)
      bottomPixels = avgLen;
    divisorRight = rightPixels * avgLen;
    divisorBottom = avgLen * bottomPixels;
    divisorCorner = rightPixels * bottomPixels;
    zdivisor = (std::complex<double>)divisor;
    zdivisorRight = (std::complex<double>)divisorRight;
    zdivisorBottom = (std::complex<double>)divisorBottom;
    zdivisorCorner = (std::complex<double>)divisorCorner;

    switch(this->_imgType){
      case DImage_u8:
	{
	  D_uint8 *pSrc;
	  float *pTmp;
	  float *rgTmp;
	  D_uint8 *pDst;
	  rgTmp = (float*)calloc(dstW*dstH, sizeof(float));
	  D_CHECKPTR(rgTmp);

	  imgDst.create(dstW, dstH, DImage_u8, _numChan, _allocMethod);
	  pSrc = pData;
	  pTmp = rgTmp;
	  pDst = imgDst.dataPointer_u8();
	  numRows = 0;
	  for(int y = 0; y < _h; ++y){
	    for(int x = 0; x < _w; ++x){
	      pTmp[x/avgLen] += (*pSrc);
	      ++pSrc;
	    }
	    ++numRows;
	    if(numRows >= avgLen){
	      for(int x = 0; x < dstWm1; ++x){
		pDst[x] = (D_uint8)(0.5 + pTmp[x] / divisor);
	      }
	      pDst[dstWm1] = (D_uint8)(0.5 +pTmp[dstWm1]/divisorRight);
	      pTmp += dstW;
	      pDst += dstW;
	      numRows = 0;
	    }
	  }
	  if(numRows > 0){ // last Row, so it didn't get divided yet
	    for(int x = 0; x < dstWm1; ++x){
	      pDst[x] = (D_uint8)(0.5 + pTmp[x] / divisorBottom);
	    }
	    pDst[dstWm1] = (D_uint8)(0.5 +pTmp[dstWm1]/divisorCorner);
	  }
	  free(rgTmp);
	}// end case
	break;
      case DImage_u16:
	{
	  D_uint16 *pSrc;
	  float *pTmp;
	  float *rgTmp;
	  D_uint16 *pDst;
	  rgTmp = (float*)calloc(dstW*dstH, sizeof(float));
	  D_CHECKPTR(rgTmp);

	  imgDst.create(dstW, dstH, DImage_u16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pTmp = rgTmp;
	  pDst = imgDst.dataPointer_u16();
	  numRows = 0;
	  for(int y = 0; y < _h; ++y){
	    for(int x = 0; x < _w; ++x){
	      pTmp[x/avgLen] += (*pSrc);
	      ++pSrc;
	    }
	    ++numRows;
	    if(numRows >= avgLen){
	      for(int x = 0; x < dstWm1; ++x){
		pDst[x] = (D_uint16)(0.5 + pTmp[x] / divisor);
	      }
	      pDst[dstWm1] = (D_uint16)(0.5 +pTmp[dstWm1]/divisorRight);
	      pTmp += dstW;
	      pDst += dstW;
	      numRows = 0;
	    }
	  }
	  if(numRows > 0){ // last Row, so it didn't get divided yet
	    for(int x = 0; x < dstWm1; ++x){
	      pDst[x] = (D_uint16)(0.5 + pTmp[x] / divisorBottom);
	    }
	    pDst[dstWm1] = (D_uint16)(0.5 +pTmp[dstWm1]/divisorCorner);
	  }
	  free(rgTmp);
	}// end case
	break;
      case DImage_u32:
	{
	  D_uint32 *pSrc;
	  double *rgTmp;
	  double *pTmp;
	  D_uint32 *pDst;
	  rgTmp = (double*)calloc(dstW*dstH, sizeof(double));
	  D_CHECKPTR(rgTmp);

	  imgDst.create(dstW, dstH, DImage_u32, _numChan, _allocMethod);
	  pSrc = (D_uint32*)pData;
	  pTmp = rgTmp;
	  pDst = imgDst.dataPointer_u32();
	  numRows = 0;
	  for(int y = 0; y < _h; ++y){
	    for(int x = 0; x < _w; ++x){
	      pTmp[x/avgLen] += (*pSrc);
	      ++pSrc;
	    }
	    ++numRows;
	    if(numRows >= avgLen){
	      for(int x = 0; x < dstWm1; ++x){
		pDst[x] = (D_uint32)(0.5 + pTmp[x] / divisor);
	      }
	      pDst[dstWm1] = (D_uint32)(0.5 +pTmp[dstWm1]/divisorRight);
	      pTmp += dstW;
	      pDst += dstW;
	      numRows = 0;
	    }
	  }
	  if(numRows > 0){ // last Row, so it didn't get divided yet
	    for(int x = 0; x < dstWm1; ++x){
	      pDst[x] = (D_uint32)(0.5 + pTmp[x] / divisorBottom);
	    }
	    pDst[dstWm1] = (D_uint32)(0.5 +pTmp[dstWm1]/divisorCorner);
	  }
	  free(rgTmp);
	}// end case
	break;
      case DImage_RGB:
	{
	  D_uint8 *pSrc;
	  float *rgTmp;
	  float *pTmp;
	  D_uint8 *pDst;
	  rgTmp = (float*)calloc(dstW*3*dstH, sizeof(float));
	  D_CHECKPTR(rgTmp);

	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrc = pData;
	  pTmp = rgTmp;
	  pDst = imgDst.dataPointer_u8();
	  numRows = 0;
	  for(int y = 0; y < _h; ++y){
	    for(int x = 0; x < _w; ++x){
	      pTmp[x/avgLen*3] += (*pSrc);
	      pTmp[x/avgLen*3+1] += pSrc[1];
	      pTmp[x/avgLen*3+2] += pSrc[2];
	      pSrc += 3;
	    }
	    ++numRows;
	    if(numRows >= avgLen){
	      for(int x = 0; x < dstWm1*3; ++x){
		pDst[x] = (D_uint8)(0.5 + pTmp[x] / divisor);
	      }
	      pDst[dstWm1*3] = (D_uint8)(0.5 +pTmp[dstWm1*3]/divisorRight);
	      pDst[dstWm1*3+1] = (D_uint8)(0.5 +pTmp[dstWm1*3+1]/divisorRight);
	      pDst[dstWm1*3+2] = (D_uint8)(0.5 +pTmp[dstWm1*3+2]/divisorRight);
	      pTmp += dstW*3;
	      pDst += dstW*3;
	      numRows = 0;
	    }
	  }
	  if(numRows > 0){ // last Row, so it didn't get divided yet
	    for(int x = 0; x < dstWm1*3; ++x){
	      pDst[x] = (D_uint8)(0.5 + pTmp[x] / divisorBottom);
	    }
	    pDst[dstWm1*3] = (D_uint8)(0.5 +pTmp[dstWm1*3]/divisorCorner);
	    pDst[dstWm1*3+1] = (D_uint8)(0.5 +pTmp[dstWm1*3+1]/divisorCorner);
	    pDst[dstWm1*3+2] = (D_uint8)(0.5 +pTmp[dstWm1*3+2]/divisorCorner);
	  }
	  free(rgTmp);
	}// end case
	break;
      case DImage_RGB_16:
	{
	  D_uint16 *pSrc;
	  float *rgTmp;
	  float *pTmp;
	  D_uint16 *pDst;
	  rgTmp = (float*)calloc(dstW*3*dstH, sizeof(float));
	  D_CHECKPTR(rgTmp);

	  imgDst.create(dstW, dstH, DImage_RGB_16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pTmp = rgTmp;
	  pDst = imgDst.dataPointer_u16();
	  numRows = 0;
	  for(int y = 0; y < _h; ++y){
	    for(int x = 0; x < _w; ++x){
	      pTmp[x/avgLen*3] += (*pSrc);
	      pTmp[x/avgLen*3+1] += pSrc[1];
	      pTmp[x/avgLen*3+2] += pSrc[2];
	      pSrc += 3;
	    }
	    ++numRows;
	    if(numRows >= avgLen){
	      for(int x = 0; x < dstWm1*3; ++x){
		pDst[x] = (D_uint16)(0.5 + pTmp[x] / divisor);
	      }
	      pDst[dstWm1*3] = (D_uint16)(0.5 +pTmp[dstWm1*3]/divisorRight);
	      pDst[dstWm1*3+1] = (D_uint16)(0.5+pTmp[dstWm1*3+1]/divisorRight);
	      pDst[dstWm1*3+2] = (D_uint16)(0.5+pTmp[dstWm1*3+2]/divisorRight);
	      pTmp += dstW*3;
	      pDst += dstW*3;
	      numRows = 0;
	    }
	  }
	  if(numRows > 0){ // last Row, so it didn't get divided yet
	    for(int x = 0; x < dstWm1*3; ++x){
	      pDst[x] = (D_uint16)(0.5 + pTmp[x] / divisorBottom);
	    }
	    pDst[dstWm1*3] = (D_uint16)(0.5 +pTmp[dstWm1*3]/divisorCorner);
	    pDst[dstWm1*3+1] = (D_uint16)(0.5 +pTmp[dstWm1*3+1]/divisorCorner);
	    pDst[dstWm1*3+2] = (D_uint16)(0.5 +pTmp[dstWm1*3+2]/divisorCorner);
	  }
	  free(rgTmp);
	}// end case
	break;
      case DImage_flt_multi:
	{
	  float *pSrc;
	  double *rgTmp;
	  double *pTmp;
	  float *pDst;
	  rgTmp = (double*)calloc(_numChan*dstW*dstH, sizeof(double));
	  D_CHECKPTR(rgTmp);

	  imgDst.create(dstW, dstH, DImage_flt_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_flt(chan);
	    pTmp = &rgTmp[chan*dstW*dstH];
	    pDst = imgDst.dataPointer_flt(chan);
	    numRows = 0;
	    for(int y = 0; y < _h; ++y){
	      for(int x = 0; x < _w; ++x){
		pTmp[x/avgLen] += (*pSrc);
		++pSrc;
	      }
	      ++numRows;
	      if(numRows >= avgLen){
		for(int x = 0; x < dstWm1; ++x){
		  pDst[x] = pTmp[x] / divisor;
		}
		pDst[dstWm1] = pTmp[dstWm1]/divisorRight;
		pTmp += dstW;
		pDst += dstW;
		numRows = 0;
	      }
	    }
	    if(numRows > 0){ // last Row, so it didn't get divided yet
	      for(int x = 0; x < dstWm1; ++x){
		pDst[x] = pTmp[x] / divisorBottom;
	      }
	      pDst[dstWm1] = pTmp[dstWm1]/divisorCorner;
	    }
	  }
	  free(rgTmp);
	}// end case
	break;
      case DImage_dbl_multi:
	{
	  double *pSrc;
	  double *rgTmp;
	  double *pTmp;
	  double *pDst;
	  rgTmp = (double*)calloc(_numChan*dstW*dstH, sizeof(double));
	  D_CHECKPTR(rgTmp);

	  imgDst.create(dstW, dstH, DImage_dbl_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_dbl(chan);
	    pTmp = &rgTmp[chan*dstW*dstH];
	    pDst = imgDst.dataPointer_dbl(chan);
	    numRows = 0;
	    for(int y = 0; y < _h; ++y){
	      for(int x = 0; x < _w; ++x){
		pTmp[x/avgLen] += (*pSrc);
		++pSrc;
	      }
	      ++numRows;
	      if(numRows >= avgLen){
		for(int x = 0; x < dstWm1; ++x){
		  pDst[x] = pTmp[x] / divisor;
		}
		pDst[dstWm1] = pTmp[dstWm1]/divisorRight;
		pTmp += dstW;
		pDst += dstW;
		numRows = 0;
	      }
	    }
	    if(numRows > 0){ // last Row, so it didn't get divided yet
	      for(int x = 0; x < dstWm1; ++x){
		pDst[x] = pTmp[x] / divisorBottom;
	      }
	      pDst[dstWm1] = pTmp[dstWm1]/divisorCorner;
	    }
	  }
	  free(rgTmp);
	}// end case
	break;
      case DImage_cmplx:
// 	fprintf(stderr, "scaledDownPow2_() doesn't support complex yet\n");

	{
	  std::complex<double> *pSrc;
	  std::complex<double> *rgTmp;
	  std::complex<double> *pTmp;
	  std::complex<double> *pDst;
	  rgTmp = (std::complex<double>*)
	    calloc(_numChan*dstW*dstH, sizeof(std::complex<double>));
	  D_CHECKPTR(rgTmp);

	  imgDst.create(dstW, dstH, DImage_dbl_multi, _numChan, _allocMethod);
	  
	  pSrc = this->dataPointer_cmplx();
	  pTmp = &rgTmp[dstW*dstH];
	  pDst = imgDst.dataPointer_cmplx();
	  numRows = 0;
	  for(int y = 0; y < _h; ++y){
	    for(int x = 0; x < _w; ++x){
	      pTmp[x/avgLen] += (*pSrc);
	      ++pSrc;
	    }
	    ++numRows;
	    if(numRows >= avgLen){
	      for(int x = 0; x < dstWm1; ++x){
		pDst[x] = pTmp[x] / zdivisor;
	      }
	      pDst[dstWm1] = pTmp[dstWm1]/zdivisorRight;
	      pTmp += dstW;
	      pDst += dstW;
	      numRows = 0;
	    }
	  }
	  if(numRows > 0){ // last Row, so it didn't get divided yet
	    for(int x = 0; x < dstWm1; ++x){
	      pDst[x] = pTmp[x] / zdivisorBottom;
	    }
	    pDst[dstWm1] = pTmp[dstWm1]/zdivisorCorner;
	  }
	  free(rgTmp);
	}// end case
	break;
	
    }// end switch
  }
  else if(DImageTransSample){

    switch(this->_imgType){
      case DImage_u8:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  D_uint8 *pDstRow;
	  D_uint8 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_u8, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y];
	    pSrcRow = &pSrc[_w*y*avgLen];
	    for(int x = 0; x < dstW; ++x){
	      pDstRow[x] = (*pSrcRow);
	      pSrcRow += avgLen;
	    }
	  }
	}// end case
	break;
      case DImage_u16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  D_uint16 *pDstRow;
	  D_uint16 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_u16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y];
	    pSrcRow = &pSrc[_w*y*avgLen];
	    for(int x = 0; x < dstW; ++x){
	      pDstRow[x] = (*pSrcRow);
	      pSrcRow += avgLen;
	    }
	  }
	}// end case
	break;
      case DImage_u32:
	{
	  D_uint32 *pSrc;
	  D_uint32 *pDst;
	  D_uint32 *pDstRow;
	  D_uint32 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_u32, _numChan, _allocMethod);
	  pSrc = (D_uint32*)pData;
	  pDst = imgDst.dataPointer_u32();
	  
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y];
	    pSrcRow = &pSrc[_w*y*avgLen];
	    for(int x = 0; x < dstW; ++x){
	      pDstRow[x] = (*pSrcRow);
	      pSrcRow += avgLen;
	    }
	  }
	}// end case
	break;
      case DImage_RGB:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  D_uint8 *pDstRow;
	  D_uint8 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y*3];
	    pSrcRow = &pSrc[_w*y*avgLen*3];
	    for(int x = 0; x < dstW*3; x+=3){
	      pDstRow[x] = (*pSrcRow);
	      pDstRow[x+1] = pSrcRow[1];
	      pDstRow[x+2] = pSrcRow[2];
	      pSrcRow += avgLen*3;
	    }
	  }
	}// end case
	break;
      case DImage_RGB_16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  D_uint16 *pDstRow;
	  D_uint16 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_RGB_16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y*3];
	    pSrcRow = &pSrc[_w*y*avgLen*3];
	    for(int x = 0; x < dstW*3; x+=3){
	      pDstRow[x] = (*pSrcRow);
	      pDstRow[x+1] = pSrcRow[1];
	      pDstRow[x+2] = pSrcRow[2];
	      pSrcRow += avgLen*3;
	    }
	  }
	}// end case
	break;
      case DImage_flt_multi:
	{
	  float *pSrc;
	  float *pDst;
	  float *pDstRow;
	  float *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_flt_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_flt(chan);
	    pDst = imgDst.dataPointer_flt(chan);
	    for(int y = 0; y < dstH; ++y){
	      pDstRow = &pDst[dstW*y];
	      pSrcRow = &pSrc[_w*y*avgLen];
	      for(int x = 0; x < dstW; ++x){
		pDstRow[x] = (*pSrcRow);
		pSrcRow += avgLen;
	      }
	    }
	  }
	}// end case
	break;
      case DImage_dbl_multi:
	{
	  double *pSrc;
	  double *pDst;
	  double *pDstRow;
	  double *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_dbl_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_dbl(chan);
	    pDst = imgDst.dataPointer_dbl(chan);
	    for(int y = 0; y < dstH; ++y){
	      pDstRow = &pDst[dstW*y];
	      pSrcRow = &pSrc[_w*y*avgLen];
	      for(int x = 0; x < dstW; ++x){
		pDstRow[x] = (*pSrcRow);
		pSrcRow += avgLen;
	      }
	    }
	  }
	}// end case
	break;
      case DImage_cmplx:
	{
	  std::complex<double> *pSrc;
	  std::complex<double> *pDst;
	  std::complex<double> *pDstRow;
	  std::complex<double> *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_cmplx, _numChan, _allocMethod);

	  pSrc = this->dataPointer_cmplx();
	  pDst = imgDst.dataPointer_cmplx();
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y];
	    pSrcRow = &pSrc[_w*y*avgLen];
	    for(int x = 0; x < dstW; ++x){
	      pDstRow[x] = (*pSrcRow);
	      pSrcRow += avgLen;
	    }
	  }
	}// end case
	break;
    }// end switch
  }
}


///Returns # of powers of 2 this image must be scaled down to fit constraints
/**Computes how many times this image must be scaled down by a power
 * of two (the numHalves parameter in scaledDownPow2_()) without
 * violating the constraints given.  The stopping criteria are that
 * the image should shrink enough times to be smaller or equal to maxW
 * by maxH, but only if that does not violate the minimum size (minW
 * by minH).  Once either dimension would go below the minimum, the
 * number of halves will not increase.  The return value may be zero
 * if the image is already smaller than maxW by maxH or if shrinking
 * the image by even a single power of 2 would break the minimum for
 * either dimension.
 */
int DImage::scaledDownNumHalves(int minW, int minH, int maxW, int maxH) const{
  int numHalves;
  int wLowRes;
  int hLowRes;

  numHalves = 0;
  wLowRes = _w;
  hLowRes = _h;
  while(((wLowRes > minW) && (hLowRes > minH)) &&
	((wLowRes > maxW) || (hLowRes > maxH))){
    wLowRes /= 2;
    hLowRes /= 2;
    ++numHalves;
  }
  return numHalves;
}

///Scale the image by scale factors scaleX and scaleY
/**Scales this image by a factor of scaleX and scaleY.  Puts the
 * result in imgDst.  When enlarging an image, you may want to use
 * DImageTransSample to get pixel replication.  When DImageTransSmooth
 * is used, each pixel in imgDst is backward-mapped to a location in
 * the source image (which may be between pixels), and bilinear
 * interpolation is used to calculate the value of the destination
 * pixel.  Since at most 4 source pixels contribute to any given
 * destination pixel's value, if you are scaling down by a factor of 2
 * or more withDImageTransSmooth , you probably would want to use
 * scaledDownPow2_() first before scaling to the final resolution so
 * that all data contributes to the final pixel values (otherwise you
 * are only doing a "smooth" scale for a subsampled version of the
 * image anyway).
 */
void DImage::scaled_(DImage &imgDst, double scaleX, double scaleY,
		     DImageTransformMode mode) const{
  int dstW, dstH;
  double invScaleX, invScaleY;
  double xTmp;

  double f00, f01, f10, f11, fx0, fx1, fxy; /* bilinear interpolation temps */
  int x0, y0;
  int x0p1, y0p1;
  double dx, dy; /* delta x and y for bilinear interpolation */
  double xp, yp; /* Xprime and Yprime for scaling */

  
  dstW = (int)(0.99999 + scaleX * _w);
  dstH = (int)(0.99999 + scaleY * _h);
  invScaleX = 1. / scaleX;
  invScaleY = 1. / scaleY;

  if(scaleX == scaleY){ // scale by same amount in both directions
    unsigned int uiInvScale;
    uiInvScale = (unsigned int)invScaleX;
    if( (invScaleX == (double)uiInvScale) &&
	(0 == (uiInvScale & (uiInvScale-1)))){
      // we are scaling down by a power of two, so special-case it for speed
      unsigned int numHalves;
      numHalves = 0;
      while(uiInvScale > 1){
	uiInvScale >>= 1;
	++numHalves;
      }
      scaledDownPow2_(imgDst, numHalves, mode);
      return;
    }
  }
//   if((scaleX >= 1.0) && (scaleY >= 1.0)){
//     // we always do pixel replication when we are scaling up
//     mode = DImageTransSample;
//   }

  if(DImageTransSmooth == mode){
    switch(this->_imgType){
      case DImage_u8:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  imgDst.create(dstW, dstH, DImage_u8, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = (x*invScaleX);
	      if(xp >= _w)
		xp = _w;
	      yp = (y*invScaleY);
	      if(yp >= _h)
		yp = _h;
	      x0 = (int)floor(xp);
	      y0 = (int)floor(yp);

	      dx = xp - x0;
	      dy = yp - y0;

	      x0p1 = x0+1;
	      if(x0p1 >= _w)
		x0p1 = x0;
	      y0p1 = y0+1;
	      if(y0p1 >= _h)
		y0p1 = y0;

	      f00 = pSrc[y0*_w+x0];
	      f10 = pSrc[y0*_w+(x0p1)];
	      f01 = pSrc[(y0p1)*_w+x0];
	      f11 = pSrc[(y0p1)*_w+(x0p1)];
	      
	      
	      fx0 = f00 + dx * (f10 - f00);
	      fx1 = f01 + dx * (f11 - f01);
	      fxy = fx0 + dy * (fx1 - fx0);
	      
	      if(fxy < 0.)
		fxy = 0;
	      if(fxy > 255.)
		fxy = 255.;
	      
	      pDst[y*dstW+x] = (D_uint8)(D_uint32)fxy;
	    } /* end for x */
	  } /* end for y */


	}// end case
	break;
      case DImage_RGB:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  imgDst.create(dstW, dstH, DImage_RGB, 3, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = (x*invScaleX);
	      if(xp >= _w)
		xp = _w;
	      yp = (y*invScaleY);
	      if(yp >= _h)
		yp = _h;
	      x0 = (int)floor(xp);
	      y0 = (int)floor(yp);

	      dx = xp - x0;
	      dy = yp - y0;

	      x0p1 = x0+1;
	      if(x0p1 >= _w)
		x0p1 = x0;
	      y0p1 = y0+1;
	      if(y0p1 >= _h)
		y0p1 = y0;

	      // Red
	      f00 = pSrc[3*(y0*_w+x0)];
	      f10 = pSrc[3*(y0*_w+(x0p1))];
	      f01 = pSrc[3*((y0p1)*_w+x0)];
	      f11 = pSrc[3*((y0p1)*_w+(x0p1))];
	      
	      fx0 = f00 + dx * (f10 - f00);
	      fx1 = f01 + dx * (f11 - f01);
	      fxy = fx0 + dy * (fx1 - fx0);
	      if(fxy < 0.)
		fxy = 0;
	      if(fxy > 255.)
		fxy = 255.;
	      pDst[3*(y*dstW+x)] = (D_uint8)(D_uint32)fxy;
	      // Green
	      f00 = pSrc[3*(y0*_w+x0)+1];
	      f10 = pSrc[3*(y0*_w+(x0p1))+1];
	      f01 = pSrc[3*((y0p1)*_w+x0)+1];
	      f11 = pSrc[3*((y0p1)*_w+(x0p1))+1];
	      
	      fx0 = f00 + dx * (f10 - f00);
	      fx1 = f01 + dx * (f11 - f01);
	      fxy = fx0 + dy * (fx1 - fx0);
	      if(fxy < 0.)
		fxy = 0;
	      if(fxy > 255.)
		fxy = 255.;
	      pDst[3*(y*dstW+x)+1] = (D_uint8)(D_uint32)fxy;
	      // Blue
	      f00 = pSrc[3*(y0*_w+x0)+2];
	      f10 = pSrc[3*(y0*_w+(x0p1))+2];
	      f01 = pSrc[3*((y0p1)*_w+x0)+2];
	      f11 = pSrc[3*((y0p1)*_w+(x0p1))+2];
	      
	      fx0 = f00 + dx * (f10 - f00);
	      fx1 = f01 + dx * (f11 - f01);
	      fxy = fx0 + dy * (fx1 - fx0);
	      if(fxy < 0.)
		fxy = 0;
	      if(fxy > 255.)
		fxy = 255.;
	      pDst[3*(y*dstW+x)+2] = (D_uint8)(D_uint32)fxy;
	    } /* end for x */
	  } /* end for y */
	}
	break;
      default:
	// TODO: implement smooth scaling for when width or height (or
	// both) are being shrunk.
	fprintf(stderr, "scaled_() Not Yet Implemented for DImageTransSmooth "
		"except when scaling down by power of 2 in both dimensions "
		"or for DImage_u8\n");
	abort();
    }
  }
  else if(DImageTransSample == mode){
    switch(this->_imgType){
      case DImage_u8:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  D_uint8 *pDstRow;
	  D_uint8 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_u8, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y];
	    pSrcRow = &pSrc[_w*(unsigned int)(y*invScaleY)];
	    xTmp = 0.;
	    for(int x = 0; x < dstW; ++x){
// 	      pDstRow[x] = pSrcRow[(unsigned int)(x*invScaleX)];
	      pDstRow[x] = pSrcRow[(unsigned int)xTmp];
	      xTmp += invScaleX;
	    }
	  }
	}// end case
	break;
      case DImage_u16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  D_uint16 *pDstRow;
	  D_uint16 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_u16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y];
	    pSrcRow = &pSrc[_w*(unsigned int)(y*invScaleY)];
	    xTmp = 0.;
	    for(int x = 0; x < dstW; ++x){
	      pDstRow[x] = pSrcRow[(unsigned int)xTmp];
	      xTmp += invScaleX;
	    }
	  }
	}// end case
	break;
      case DImage_u32:
	{
	  D_uint32 *pSrc;
	  D_uint32 *pDst;
	  D_uint32 *pDstRow;
	  D_uint32 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_u32, _numChan, _allocMethod);
	  pSrc = (D_uint32*)pData;
	  pDst = imgDst.dataPointer_u32();
	  
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y];
	    pSrcRow = &pSrc[_w*(unsigned int)(y*invScaleY)];
	    xTmp = 0.;
	    for(int x = 0; x < dstW; ++x){
	      pDstRow[x] = pSrcRow[(unsigned int)xTmp];
	      xTmp += invScaleX;
	    }
	  }
	}// end case
	break;
      case DImage_RGB:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  D_uint8 *pDstRow;
	  D_uint8 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y*3];
	    pSrcRow = &pSrc[3*_w*(unsigned int)(y*invScaleY)];
	    xTmp = 0.;
	    for(int x = 0; x < dstW*3; x+=3){
	      pDstRow[x] = pSrcRow[3*(unsigned int)xTmp];
	      pDstRow[x+1] = pSrcRow[1+3*(unsigned int)xTmp];
	      pDstRow[x+2] = pSrcRow[2+3*(unsigned int)xTmp];
	      xTmp += invScaleX;
	    }
	  }
	}// end case
	break;
      case DImage_RGB_16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  D_uint16 *pDstRow;
	  D_uint16 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_RGB_16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y*3];
	    pSrcRow = &pSrc[3*_w*(unsigned int)(y*invScaleY)];
	    xTmp = 0.;
	    for(int x = 0; x < dstW*3; x+=3){
	      pDstRow[x] = pSrcRow[3*(unsigned int)xTmp];
	      pDstRow[x+1] = pSrcRow[1+3*(unsigned int)xTmp];
	      pDstRow[x+2] = pSrcRow[2+3*(unsigned int)xTmp];
	      xTmp += invScaleX;
	    }
	  }
	}// end case
	break;
      case DImage_flt_multi:
	{
	  float *pSrc;
	  float *pDst;
	  float *pDstRow;
	  float *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_flt_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_flt(chan);
	    pDst = imgDst.dataPointer_flt(chan);
	    for(int y = 0; y < dstH; ++y){
	      pDstRow = &pDst[dstW*y];
	      pSrcRow = &pSrc[_w*(unsigned int)(y*invScaleY)];
	      xTmp = 0.;
	      for(int x = 0; x < dstW; ++x){
		pDstRow[x] = pSrcRow[(unsigned int)xTmp];
		xTmp += invScaleX;
	      }
	    }
	  }
	}// end case
	break;
      case DImage_dbl_multi:
	{
	  double *pSrc;
	  double *pDst;
	  double *pDstRow;
	  double *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_dbl_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_dbl(chan);
	    pDst = imgDst.dataPointer_dbl(chan);
	    for(int y = 0; y < dstH; ++y){
	      pDstRow = &pDst[dstW*y];
	      pSrcRow = &pSrc[_w*(unsigned int)(y*invScaleY)];
	      xTmp = 0.;
	      for(int x = 0; x < dstW; ++x){
		pDstRow[x] = pSrcRow[(unsigned int)xTmp];
		xTmp += invScaleX;
	      }
	    }
	  }
	}// end case
	break;
      case DImage_cmplx:
	{
	  std::complex<double> *pSrc;
	  std::complex<double> *pDst;
	  std::complex<double> *pDstRow;
	  std::complex<double> *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_cmplx, _numChan, _allocMethod);

	  pSrc = this->dataPointer_cmplx();
	  pDst = imgDst.dataPointer_cmplx();
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y];
	    pSrcRow = &pSrc[_w*(unsigned int)(y*invScaleY)];
	    xTmp = 0.;
	    for(int x = 0; x < dstW; ++x){
	      pDstRow[x] = pSrcRow[(unsigned int)xTmp];
	      xTmp += invScaleX;
	    }
	  }
	}// end case
	break;
    }// end switch
  } // end if sample transform mode
}
///Scales this image and stores the scaled image in destination image imgDst
/** mode specifies whether to use sampled values in the resulting image or
 * whether to use linear interpolation to choose the resulting values.
 */
void DImage::scaled_(DImage &imgDst, int dstW, int dstH,
		     DImageTransformMode mode) const{
  double scaleX, scaleY;
#ifdef DEBUG
  if((dstW < 1) || (dstH < 1)){
    fprintf(stderr, "DImage::scaled_() destination image "
	    "width, height must be positive\n");
    abort();
    return;
  }
#endif
  scaleX = dstW / (double)_w;
  scaleY = dstH / (double)_h;
  scaled_(imgDst, scaleX, scaleY, mode);
#ifdef DEBUG
  if((imgDst._w != dstW) || (imgDst._h != dstH)){
    fprintf(stderr,"DImage::scaled_() changed dstW,dstH from %d,%d to %d,%d\n",
	    dstW, dstH, imgDst._w, imgDst._h);
  }
#endif
}


///Returns a version of this image translated right and down by dx,dy pixels
/**The data will be clipped at the edge and the opposite side will be
 * filled with the pad_val grayscale value.
 */
DImage DImage::translated(double tx, double ty, double pad_val,
			  DImageTransformMode mode) const{
  DImage dst;
  translated_(dst, tx, ty, pad_val, mode);
  return dst;
}

///Returns a version of this image translated right and down by dx,dy pixels
/**The data will be clipped at the edge and the opposite side will be
 * filled with color pad_R,pad_G,pad_B.
 */
DImage DImage::translated(double tx, double ty, int pad_R,int pad_G,int pad_B,
			  DImageTransformMode mode) const{
  DImage dst;
  translated_(dst, tx, ty, pad_R, pad_G, pad_B, mode);
  return dst;
}


///Put a version of this image into imgDst, translated by dx,dy pixels
/**The data will be clipped at the edge and the opposite side will be
 * filled with grayscale value pad_val. Translation is positive right and down.
 */
void DImage::translated_(DImage &imgDst, double tx, double ty,
			 double pad_val, DImageTransformMode mode) const{
  double dx, dy; /* delta x and y for bilinear interpolation */
  double xp, yp; /* Xprime and Yprime for translation */
  int dstW, dstH; /* width and height of translated image */
  double f00, f01, f10, f11, fx0, fx1, fxy; /* bilinear interpolation temps */
  int x0, y0;

  dstW = _w;
  dstH = _h;

  if(DImageTransSmooth == mode){
    //TODO: make this faster
    switch(this->_imgType){
      case DImage_u8:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  imgDst.create(dstW, dstH, DImage_u8, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = (x-tx);
	      yp = (y-ty);
	      x0 = (int)floor(xp);
	      y0 = (int)floor(yp);
	      dx = xp - x0;
	      dy = yp - y0;
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h))
		f00 = pad_val;
	      else
		f00 = pSrc[y0*_w+x0];
	      if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h))
		f10 = pad_val;
	      else
		f10 = pSrc[y0*_w+(x0+1)];
	      if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h))
		f01 = pad_val;
	      else
		f01 = pSrc[(y0+1)*_w+x0];
	      if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1) >= _h))
		f11 = pad_val;
	      else
		f11 = pSrc[(y0+1)*_w+(x0+1)];
	      
	      fx0 = f00 + dx * (f10 - f00);
	      fx1 = f01 + dx * (f11 - f01);
	      fxy = fx0 + dy * (fx1 - fx0);
	      
	      if(fxy < 0.)
		fxy = 0;
	      if(fxy > 255.)
		fxy = 255.;
	      
	      pDst[y*dstW+x] = (D_uint8)(D_uint32)fxy;
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_u16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  imgDst.create(dstW, dstH, DImage_u16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = (x-tx);
	      yp = (y-ty);
	      x0 = (int)floor(xp);
	      y0 = (int)floor(yp);
	      dx = xp - x0;
	      dy = yp - y0;
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h))
		f00 = pad_val;
	      else
		f00 = pSrc[y0*_w+x0];
	      if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h))
		f10 = pad_val;
	      else
		f10 = pSrc[y0*_w+(x0+1)];
	      if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h))
		f01 = pad_val;
	      else
		f01 = pSrc[(y0+1)*_w+x0];
	      if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1) >= _h))
		f11 = pad_val;
	      else
		f11 = pSrc[(y0+1)*_w+(x0+1)];
	      
	      fx0 = f00 + dx * (f10 - f00);
	      fx1 = f01 + dx * (f11 - f01);
	      fxy = fx0 + dy * (fx1 - fx0);
	      
	      if(fxy < 0.)
		fxy = 0;
	      if(fxy > 255.)
		fxy = 255.;
	      
	      pDst[y*dstW+x] = (D_uint16)(D_uint32)fxy;
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_u32:
	{
	  D_uint32 *pSrc;
	  D_uint32 *pDst;
	  imgDst.create(dstW, dstH, DImage_u32, _numChan, _allocMethod);
	  pSrc = (D_uint32*)pData;
	  pDst = imgDst.dataPointer_u32();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = (x-tx);
	      yp = (y-ty);
	      x0 = (int)floor(xp);
	      y0 = (int)floor(yp);
	      dx = xp - x0;
	      dy = yp - y0;
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h))
		f00 = pad_val;
	      else
		f00 = pSrc[y0*_w+x0];
	      if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h))
		f10 = pad_val;
	      else
		f10 = pSrc[y0*_w+(x0+1)];
	      if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h))
		f01 = pad_val;
	      else
		f01 = pSrc[(y0+1)*_w+x0];
	      if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1) >= _h))
		f11 = pad_val;
	      else
		f11 = pSrc[(y0+1)*_w+(x0+1)];
	      
	      fx0 = f00 + dx * (f10 - f00);
	      fx1 = f01 + dx * (f11 - f01);
	      fxy = fx0 + dy * (fx1 - fx0);
	      
	      if(fxy < 0.)
		fxy = 0;
	      if(fxy > 255.)
		fxy = 255.;
	      
	      pDst[y*dstW+x] = (D_uint32)(D_uint32)fxy;
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_RGB:
      case DImage_RGB_16:
	fprintf(stderr,
		"translated_() called for RGB image with non-RGB pad value\n");
	abort();
	break;
      case DImage_flt_multi:
	{
	  float *pSrc;
	  float *pDst;
	  
	  imgDst.create(dstW, dstH, DImage_flt_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_flt(chan);
	    pDst = imgDst.dataPointer_flt(chan);
	    for(int y = 0; y < dstH; ++y){
	      for(int x = 0; x < dstW; ++x){
		xp = (x-tx);
		yp = (y-ty);
		x0 = (int)floor(xp);
		y0 = (int)floor(yp);
		dx = xp - x0;
		dy = yp - y0;
		if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h))
		  f00 = pad_val;
		else
		  f00 = pSrc[y0*_w+x0];
		if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h))
		  f10 = pad_val;
		else
		  f10 = pSrc[y0*_w+(x0+1)];
		if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h))
		  f01 = pad_val;
		else
		  f01 = pSrc[(y0+1)*_w+x0];
		if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1)>=_h))
		  f11 = pad_val;
		else
		  f11 = pSrc[(y0+1)*_w+(x0+1)];
		
		fx0 = f00 + dx * (f10 - f00);
		fx1 = f01 + dx * (f11 - f01);
		fxy = fx0 + dy * (fx1 - fx0);
		
		if(fxy < 0.)
		  fxy = 0;
		if(fxy > 255.)
		  fxy = 255.;
		
		pDst[y*dstW+x] = (float)fxy;
	      } /* end for x */
	    } /* end for y */
	  }
	}// end case
	break;
      case DImage_dbl_multi:
	{
	  double *pSrc;
	  double *pDst;
	  
	  imgDst.create(dstW, dstH, DImage_dbl_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_dbl(chan);
	    pDst = imgDst.dataPointer_dbl(chan);
	    for(int y = 0; y < dstH; ++y){
	      for(int x = 0; x < dstW; ++x){
		xp = (x-tx);
		yp = (y-ty);
		x0 = (int)floor(xp);
		y0 = (int)floor(yp);
		dx = xp - x0;
		dy = yp - y0;
		if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h))
		  f00 = pad_val;
		else
		  f00 = pSrc[y0*_w+x0];
		if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h))
		  f10 = pad_val;
		else
		  f10 = pSrc[y0*_w+(x0+1)];
		if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h))
		  f01 = pad_val;
		else
		  f01 = pSrc[(y0+1)*_w+x0];
		if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1)>=_h))
		  f11 = pad_val;
		else
		  f11 = pSrc[(y0+1)*_w+(x0+1)];
		
		fx0 = f00 + dx * (f10 - f00);
		fx1 = f01 + dx * (f11 - f01);
		fxy = fx0 + dy * (fx1 - fx0);
		
		if(fxy < 0.)
		  fxy = 0;
		if(fxy > 255.)
		  fxy = 255.;
		
		pDst[y*dstW+x] = fxy;
	      } /* end for x */
	    } /* end for y */
	  }
	}// end case
	break;
      case DImage_cmplx:
	fprintf(stderr,
		"translated_() not supported for "
		"smooth transform on complex images.\n");
	abort();
	break;
    }// end switch
  } // end if sample transform mode
  else if(DImageTransSample){
    switch(this->_imgType){
      case DImage_u8:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  imgDst.create(dstW, dstH, DImage_u8, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = (x-tx);
	      yp = (y-ty);
	      x0 = (int)rint(xp);
	      y0 = (int)rint(yp);
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		pDst[y*dstW+x] = (D_uint8)pad_val;
	      }
	      else{
		pDst[y*dstW+x] = pSrc[y0*_w+x0];
	      }
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_u16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  imgDst.create(dstW, dstH, DImage_u16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = (x-tx);
	      yp = (y-ty);
	      x0 = (int)rint(xp);
	      y0 = (int)rint(yp);
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		pDst[y*dstW+x] = (D_uint16)pad_val;
	      }
	      else{
		pDst[y*dstW+x] = pSrc[y0*_w+x0];
	      }
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_u32:
	{
	  D_uint32 *pSrc;
	  D_uint32 *pDst;
	  imgDst.create(dstW, dstH, DImage_u32, _numChan, _allocMethod);
	  pSrc = (D_uint32*)pData;
	  pDst = imgDst.dataPointer_u32();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = (x-tx);
	      yp = (y-ty);
	      x0 = (int)rint(xp);
	      y0 = (int)rint(yp);
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		pDst[y*dstW+x] = (D_uint32)pad_val;
	      }
	      else{
		pDst[y*dstW+x] = pSrc[y0*_w+x0];
	      }
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_RGB:
      case DImage_RGB_16:
	fprintf(stderr,
		"translated_() called for RGB image with non-RGB pad value\n");
	abort();
	break;
      case DImage_flt_multi:
	{
	  float *pSrc;
	  float *pDst;
	  
	  imgDst.create(dstW, dstH, DImage_flt_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_flt(chan);
	    pDst = imgDst.dataPointer_flt(chan);
	    for(int y = 0; y < dstH; ++y){
	      for(int x = 0; x < dstW; ++x){
		xp = (x-tx);
		yp = (y-ty);
		x0 = (int)rint(xp);
		y0 = (int)rint(yp);
		if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		  pDst[y*dstW+x] = (float)pad_val;
		}
		else{
		  pDst[y*dstW+x] = pSrc[y0*_w+x0];
		}
	      } /* end for x */
	    } /* end for y */
	  }
	}// end case
	break;
      case DImage_dbl_multi:
	{
	  double *pSrc;
	  double *pDst;
	  
	  imgDst.create(dstW, dstH, DImage_dbl_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_dbl(chan);
	    pDst = imgDst.dataPointer_dbl(chan);
	    for(int y = 0; y < dstH; ++y){
	      for(int x = 0; x < dstW; ++x){
		xp = (x-tx);
		yp = (y-ty);
		x0 = (int)rint(xp);
		y0 = (int)rint(yp);
		if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		  pDst[y*dstW+x] = (double)pad_val;
		}
		else{
		  pDst[y*dstW+x] = pSrc[y0*_w+x0];
		}
	      } /* end for x */
	    } /* end for y */
	  }
	}// end case
	break;
      case DImage_cmplx:
	{
	  std::complex<double> *pSrc;
	  std::complex<double> *pDst;
	  imgDst.create(dstW, dstH, DImage_cmplx, _numChan, _allocMethod);
	  pSrc = (std::complex<double>*)pData;
	  pDst = imgDst.dataPointer_cmplx();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = (x-tx);
	      yp = (y-ty);
	      x0 = (int)rint(xp);
	      y0 = (int)rint(yp);
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		pDst[y*dstW+x] = (std::complex<double>)pad_val;
	      }
	      else{
		pDst[y*dstW+x] = pSrc[y0*_w+x0];
	      }
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
    }// end switch
  }

}



///Put a version of this image into imgDst, translated by dx,dy pixels
/**The data will be clipped at the edge and the opposite side will be
 * filled with color pad_R,pad_G,pad_B. Translation is positive right and down.
 */
void DImage::translated_(DImage &imgDst, double tx, double ty,
			 int pad_R, int pad_G, int pad_B,
			 DImageTransformMode mode) const{
  double dx, dy; /* delta x and y for bilinear interpolation */
  double xp, yp; /* Xprime and Yprime for translation */
  int dstW, dstH; /* width and height of translated image */
  int x0, y0;
  /* temporary variables used for bilinear interpolation: */
  double f00_R, f01_R, f10_R, f11_R, fx0_R, fx1_R, fxy_R;
  double f00_G, f01_G, f10_G, f11_G, fx0_G, fx1_G, fxy_G;
  double f00_B, f01_B, f10_B, f11_B, fx0_B, fx1_B, fxy_B;

  dstW = _w;
  dstH = _h;


  if(DImageTransSmooth == mode){
    //TODO: make this faster
    switch(this->_imgType){
      case DImage_u8:
      case DImage_u16:
      case DImage_u32:
      case DImage_flt_multi:
      case DImage_dbl_multi:
      case DImage_cmplx:
	fprintf(stderr, "translated_() called for non-RGB image with RGB "
		"pad\n");
	abort();
	break;
      case DImage_RGB:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = (x-tx);
	      yp = (y-ty);
	      x0 = (int)floor(xp);
	      y0 = (int)floor(yp);
	      dx = xp - x0;
	      dy = yp - y0;
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		f00_R = pad_R;
		f00_G = pad_G;
		f00_B = pad_B;
	      }
	      else{
		f00_R = pSrc[(y0*_w+x0)*3];
		f00_G = pSrc[(y0*_w+x0)*3+1];
		f00_B = pSrc[(y0*_w+x0)*3+2];
	      }
	      if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h)){
		f10_R = pad_R;
		f10_G = pad_G;
		f10_B = pad_B;
	      }
	      else{
		f10_R = pSrc[(y0*_w+(x0+1))*3];
		f10_G = pSrc[(y0*_w+(x0+1))*3+1];
		f10_B = pSrc[(y0*_w+(x0+1))*3+2];
	      }
	      if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h)){
		f01_R = pad_R;
		f01_G = pad_G;
		f01_B = pad_B;
	      }
	      else{
		f01_R = pSrc[((y0+1)*_w+x0)*3];
		f01_G = pSrc[((y0+1)*_w+x0)*3+1];
		f01_B = pSrc[((y0+1)*_w+x0)*3+2];
	      }
	      if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1) >=_h)){
		f11_R = pad_R;
		f11_G = pad_G;
		f11_B = pad_B;
	      }
	      else{
		f11_R = pSrc[((y0+1)*_w+(x0+1))*3];
		f11_G = pSrc[((y0+1)*_w+(x0+1))*3+1];
		f11_B = pSrc[((y0+1)*_w+(x0+1))*3+2];
	      }
	      
	      fx0_R = f00_R + dx * (f10_R - f00_R);
	      fx0_G = f00_G + dx * (f10_G - f00_G);
	      fx0_B = f00_B + dx * (f10_B - f00_B);
	      fx1_R = f01_R + dx * (f11_R - f01_R);
	      fx1_G = f01_G + dx * (f11_G - f01_G);
	      fx1_B = f01_B + dx * (f11_B - f01_B);
	      fxy_R = fx0_R + dy * (fx1_R - fx0_R);
	      fxy_G = fx0_G + dy * (fx1_G - fx0_G);
	      fxy_B = fx0_B + dy * (fx1_B - fx0_B);
	      
	      if(fxy_R < 0.)
		fxy_R = 0;
	      if(fxy_R > 255.)
		fxy_R = 255.;
	      if(fxy_G < 0.)
		fxy_G = 0;
	      if(fxy_G > 255.)
		fxy_G = 255.;
	      if(fxy_B < 0.)
		fxy_B = 0;
	      if(fxy_B > 255.)
		fxy_B = 255.;
	      
	      pDst[(y*dstW+x)*3] = (D_uint8)(D_uint32)fxy_R;
	      pDst[(y*dstW+x)*3+1] = (D_uint8)(D_uint32)fxy_G;
	      pDst[(y*dstW+x)*3+2] = (D_uint8)(D_uint32)fxy_B;
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_RGB_16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = (x-tx);
	      yp = (y-ty);
	      x0 = (int)floor(xp);
	      y0 = (int)floor(yp);
	      dx = xp - x0;
	      dy = yp - y0;
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		f00_R = pad_R;
		f00_G = pad_G;
		f00_B = pad_B;
	      }
	      else{
		f00_R = pSrc[(y0*_w+x0)*3];
		f00_G = pSrc[(y0*_w+x0)*3+1];
		f00_B = pSrc[(y0*_w+x0)*3+2];
	      }
	      if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h)){
		f10_R = pad_R;
		f10_G = pad_G;
		f10_B = pad_B;
	      }
	      else{
		f10_R = pSrc[(y0*_w+(x0+1))*3];
		f10_G = pSrc[(y0*_w+(x0+1))*3+1];
		f10_B = pSrc[(y0*_w+(x0+1))*3+2];
	      }
	      if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h)){
		f01_R = pad_R;
		f01_G = pad_G;
		f01_B = pad_B;
	      }
	      else{
		f01_R = pSrc[((y0+1)*_w+x0)*3];
		f01_G = pSrc[((y0+1)*_w+x0)*3+1];
		f01_B = pSrc[((y0+1)*_w+x0)*3+2];
	      }
	      if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1) >=_h)){
		f11_R = pad_R;
		f11_G = pad_G;
		f11_B = pad_B;
	      }
	      else{
		f11_R = pSrc[((y0+1)*_w+(x0+1))*3];
		f11_G = pSrc[((y0+1)*_w+(x0+1))*3+1];
		f11_B = pSrc[((y0+1)*_w+(x0+1))*3+2];
	      }
	      
	      fx0_R = f00_R + dx * (f10_R - f00_R);
	      fx0_G = f00_G + dx * (f10_G - f00_G);
	      fx0_B = f00_B + dx * (f10_B - f00_B);
	      fx1_R = f01_R + dx * (f11_R - f01_R);
	      fx1_G = f01_G + dx * (f11_G - f01_G);
	      fx1_B = f01_B + dx * (f11_B - f01_B);
	      fxy_R = fx0_R + dy * (fx1_R - fx0_R);
	      fxy_G = fx0_G + dy * (fx1_G - fx0_G);
	      fxy_B = fx0_B + dy * (fx1_B - fx0_B);
	      
	      if(fxy_R < 0.)
		fxy_R = 0;
	      if(fxy_R > 255.)
		fxy_R = 255.;
	      if(fxy_G < 0.)
		fxy_G = 0;
	      if(fxy_G > 255.)
		fxy_G = 255.;
	      if(fxy_B < 0.)
		fxy_B = 0;
	      if(fxy_B > 255.)
		fxy_B = 255.;
	      
	      pDst[(y*dstW+x)*3] = (D_uint16)(D_uint32)fxy_R;
	      pDst[(y*dstW+x)*3+1] = (D_uint16)(D_uint32)fxy_G;
	      pDst[(y*dstW+x)*3+2] = (D_uint16)(D_uint32)fxy_B;
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
    }// end switch
  } // end if sample transform mode
  else if(DImageTransSample){
    switch(this->_imgType){
      case DImage_u8:
      case DImage_u16:
      case DImage_u32:
      case DImage_flt_multi:
      case DImage_dbl_multi:
      case DImage_cmplx:
	fprintf(stderr, "translated_() called for non-RGB image with RGB "
		"pad\n");
	abort();
	break;
      case DImage_RGB:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = (x-tx);
	      yp = (y-ty);
	      x0 = (int)rint(xp);
	      y0 = (int)rint(yp);
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		pDst[(y*dstW+x)*3] = (D_uint8)pad_R;
		pDst[(y*dstW+x)*3+1] = (D_uint8)pad_G;
		pDst[(y*dstW+x)*3+2] = (D_uint8)pad_B;
	      }
	      else{
		pDst[(y*dstW+x)*3] = pSrc[(y0*_w+x0)*3];
		pDst[(y*dstW+x)*3+1] = pSrc[(y0*_w+x0)*3+1];
		pDst[(y*dstW+x)*3+2] = pSrc[(y0*_w+x0)*3+2];
	      }
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_RGB_16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  imgDst.create(dstW, dstH, DImage_RGB_16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = (x-tx);
	      yp = (y-ty);
	      x0 = (int)rint(xp);
	      y0 = (int)rint(yp);
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		//TODO: move casts outside of loop by creating a temporary
		// variable of the correct type with the casted value stored.
		// Do the same for all other similar code
		pDst[(y*dstW+x)*3] = (D_uint16)pad_R;
		pDst[(y*dstW+x)*3+1] = (D_uint16)pad_G;
		pDst[(y*dstW+x)*3+2] = (D_uint16)pad_B;
	      }
	      else{
		pDst[(y*dstW+x)*3] = pSrc[(y0*_w+x0)*3];
		pDst[(y*dstW+x)*3+1] = pSrc[(y0*_w+x0)*3+1];
		pDst[(y*dstW+x)*3+2] = pSrc[(y0*_w+x0)*3+2];
	      }
	    } /* end for x */
	  } /* end for y */
	}
	break;
    }//end switch
  } // end if(mode==...
}




///Returns a version of this image rotated clockwise by ang degrees
/**Rotation is around the image center.  If fResize is true, then the
 * image will be resized so that the rotated image fits completely.
 * Otherwise the rotated image is cropped within the frame of the
 * original image.  Pixels that weren't part of the original image
 * will be set to pad_val.
 */
DImage DImage::rotated(double ang, double pad_val, bool fResize,
		       DImageTransformMode mode) const{
  DImage dst;
  rotated_(dst, ang, pad_val, fResize, mode);
  return dst;
}

///Returns a version of this image rotated clockwise by ang degrees
/**Rotation is around the image center.  If fResize is true, then the
 * image will be resized so that the rotated image fits completely.
 * Otherwise the rotated image is cropped within the frame of the
 * original image.  Pixels that weren't part of the original image
 * will be set to pad_R,pad_G,pad_B.
 */
DImage DImage::rotated(double ang, int pad_R,int pad_G,int pad_B,
		       bool fResize, DImageTransformMode mode) const{
  DImage dst;
  rotated_(dst, ang, pad_R, pad_G, pad_B, fResize, mode);
  return dst;
}

///Put a version of this image into imgDst, rotated clockwise ang degrees
/**The image is rotated around its center.  If fResize is true, then the image
 * is padded to make sure it all fits.  Pixels that aren't part of the
 * original image are filled with pad_val.
 */
void DImage::rotated_(DImage &imgDst, double ang,
		      double pad_val, bool fResize,
		      DImageTransformMode mode) const{
  double dx, dy; /* delta x and y for bilinear interpolation */
  double xp, yp; /* Xprime and Yprime for rotation */
  double xc, yc; /* position of center of image */
  double xOffs, yOffs;
  DSize rotSize; /* size (width, height) of rotated image */
  int dstW, dstH; /* width and height of rotated image */
  double cosTheta, sinTheta;
  double f00, f01, f10, f11, fx0, fx1, fxy; /* bilinear interpolation temps */
  int x0, y0;

  cosTheta = cos(DImage_radFromDeg(-1. * ang));
  sinTheta = sin(DImage_radFromDeg(-1. * ang));

  if(fResize){
    rotSize = getRotatedSize(ang);
    dstW = rotSize.w;
    dstH = rotSize.h;
  }
  else{
    dstW = _w;
    dstH = _h;
  }
  xOffs = (dstW - _w) / 2.;
  yOffs = (dstH - _h) / 2.;
  xc = _w / 2.;
  yc = _h / 2.;

  if(DImageTransSmooth == mode){
    //TODO: make this faster (check reference used by pnmrotate:
    // "A Fast Algorithm for General Raster Rotation" by Alan Paeth,
    // Graphics Interface '86, pp. 77-81.
    switch(this->_imgType){
      case DImage_u8:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  imgDst.create(dstW, dstH, DImage_u8, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) * sinTheta + xc;
	      yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) * cosTheta + yc;
	      x0 = (int)floor(xp);
	      y0 = (int)floor(yp);
	      dx = xp - x0;
	      dy = yp - y0;
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h))
		f00 = pad_val;
	      else
		f00 = pSrc[y0*_w+x0];
	      if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h))
		f10 = pad_val;
	      else
		f10 = pSrc[y0*_w+(x0+1)];
	      if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h))
		f01 = pad_val;
	      else
		f01 = pSrc[(y0+1)*_w+x0];
	      if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1) >= _h))
		f11 = pad_val;
	      else
		f11 = pSrc[(y0+1)*_w+(x0+1)];
	      
	      fx0 = f00 + dx * (f10 - f00);
	      fx1 = f01 + dx * (f11 - f01);
	      fxy = fx0 + dy * (fx1 - fx0);
	      
	      if(fxy < 0.)
		fxy = 0;
	      if(fxy > 255.)
		fxy = 255.;
	      
	      pDst[y*dstW+x] = (D_uint8)(D_uint32)fxy;
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_u16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  imgDst.create(dstW, dstH, DImage_u16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) * sinTheta + xc;
	      yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) * cosTheta + yc;
	      x0 = (int)floor(xp);
	      y0 = (int)floor(yp);
	      dx = xp - x0;
	      dy = yp - y0;
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h))
		f00 = pad_val;
	      else
		f00 = pSrc[y0*_w+x0];
	      if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h))
		f10 = pad_val;
	      else
		f10 = pSrc[y0*_w+(x0+1)];
	      if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h))
		f01 = pad_val;
	      else
		f01 = pSrc[(y0+1)*_w+x0];
	      if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1) >= _h))
		f11 = pad_val;
	      else
		f11 = pSrc[(y0+1)*_w+(x0+1)];
	      
	      fx0 = f00 + dx * (f10 - f00);
	      fx1 = f01 + dx * (f11 - f01);
	      fxy = fx0 + dy * (fx1 - fx0);
	      
	      if(fxy < 0.)
		fxy = 0;
	      if(fxy > 255.)
		fxy = 255.;
	      
	      pDst[y*dstW+x] = (D_uint16)(D_uint32)fxy;
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_u32:
	{
	  D_uint32 *pSrc;
	  D_uint32 *pDst;
	  imgDst.create(dstW, dstH, DImage_u32, _numChan, _allocMethod);
	  pSrc = (D_uint32*)pData;
	  pDst = imgDst.dataPointer_u32();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) * sinTheta + xc;
	      yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) * cosTheta + yc;
	      x0 = (int)floor(xp);
	      y0 = (int)floor(yp);
	      dx = xp - x0;
	      dy = yp - y0;
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h))
		f00 = pad_val;
	      else
		f00 = pSrc[y0*_w+x0];
	      if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h))
		f10 = pad_val;
	      else
		f10 = pSrc[y0*_w+(x0+1)];
	      if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h))
		f01 = pad_val;
	      else
		f01 = pSrc[(y0+1)*_w+x0];
	      if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1) >= _h))
		f11 = pad_val;
	      else
		f11 = pSrc[(y0+1)*_w+(x0+1)];
	      
	      fx0 = f00 + dx * (f10 - f00);
	      fx1 = f01 + dx * (f11 - f01);
	      fxy = fx0 + dy * (fx1 - fx0);
	      
	      if(fxy < 0.)
		fxy = 0;
	      if(fxy > 255.)
		fxy = 255.;
	      
	      pDst[y*dstW+x] = (D_uint32)(D_uint32)fxy;
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_RGB:
      case DImage_RGB_16:
	fprintf(stderr,
		"rotated_() called for RGB image with non-RGB pad value\n");
	abort();
	break;
      case DImage_flt_multi:
	{
	  float *pSrc;
	  float *pDst;
	  
	  imgDst.create(dstW, dstH, DImage_flt_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_flt(chan);
	    pDst = imgDst.dataPointer_flt(chan);
	    for(int y = 0; y < dstH; ++y){
	      for(int x = 0; x < dstW; ++x){
		xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) *sinTheta+ xc;
		yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) *cosTheta+ yc;
		x0 = (int)floor(xp);
		y0 = (int)floor(yp);
		dx = xp - x0;
		dy = yp - y0;
		if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h))
		  f00 = pad_val;
		else
		  f00 = pSrc[y0*_w+x0];
		if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h))
		  f10 = pad_val;
		else
		  f10 = pSrc[y0*_w+(x0+1)];
		if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h))
		  f01 = pad_val;
		else
		  f01 = pSrc[(y0+1)*_w+x0];
		if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1)>=_h))
		  f11 = pad_val;
		else
		  f11 = pSrc[(y0+1)*_w+(x0+1)];
		
		fx0 = f00 + dx * (f10 - f00);
		fx1 = f01 + dx * (f11 - f01);
		fxy = fx0 + dy * (fx1 - fx0);
		
		if(fxy < 0.)
		  fxy = 0;
		if(fxy > 255.)
		  fxy = 255.;
		
		pDst[y*dstW+x] = (float)fxy;
	      } /* end for x */
	    } /* end for y */
	  }
	}// end case
	break;
      case DImage_dbl_multi:
	{
	  double *pSrc;
	  double *pDst;
	  
	  imgDst.create(dstW, dstH, DImage_dbl_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_dbl(chan);
	    pDst = imgDst.dataPointer_dbl(chan);
	    for(int y = 0; y < dstH; ++y){
	      for(int x = 0; x < dstW; ++x){
		xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) *sinTheta+ xc;
		yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) *cosTheta+ yc;
		x0 = (int)floor(xp);
		y0 = (int)floor(yp);
		dx = xp - x0;
		dy = yp - y0;
		if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h))
		  f00 = pad_val;
		else
		  f00 = pSrc[y0*_w+x0];
		if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h))
		  f10 = pad_val;
		else
		  f10 = pSrc[y0*_w+(x0+1)];
		if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h))
		  f01 = pad_val;
		else
		  f01 = pSrc[(y0+1)*_w+x0];
		if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1)>=_h))
		  f11 = pad_val;
		else
		  f11 = pSrc[(y0+1)*_w+(x0+1)];
		
		fx0 = f00 + dx * (f10 - f00);
		fx1 = f01 + dx * (f11 - f01);
		fxy = fx0 + dy * (fx1 - fx0);
		
		if(fxy < 0.)
		  fxy = 0;
		if(fxy > 255.)
		  fxy = 255.;
		
		pDst[y*dstW+x] = fxy;
	      } /* end for x */
	    } /* end for y */
	  }
	}// end case
	break;
      case DImage_cmplx:
	fprintf(stderr,
		"rotated_() not supported for "
		"smooth transform on complex images.\n");
	abort();
	break;
    }// end switch
  } // end if sample transform mode
  else if(DImageTransSample){
    switch(this->_imgType){
      case DImage_u8:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  imgDst.create(dstW, dstH, DImage_u8, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) * sinTheta + xc;
	      yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) * cosTheta + yc;
	      x0 = (int)rint(xp);
	      y0 = (int)rint(yp);
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		pDst[y*dstW+x] = (D_uint8)pad_val;
	      }
	      else{
		pDst[y*dstW+x] = pSrc[y0*_w+x0];
	      }
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_u16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  imgDst.create(dstW, dstH, DImage_u16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) * sinTheta + xc;
	      yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) * cosTheta + yc;
	      x0 = (int)rint(xp);
	      y0 = (int)rint(yp);
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		pDst[y*dstW+x] = (D_uint16)pad_val;
	      }
	      else{
		pDst[y*dstW+x] = pSrc[y0*_w+x0];
	      }
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_u32:
	{
	  D_uint32 *pSrc;
	  D_uint32 *pDst;
	  imgDst.create(dstW, dstH, DImage_u32, _numChan, _allocMethod);
	  pSrc = (D_uint32*)pData;
	  pDst = imgDst.dataPointer_u32();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) * sinTheta + xc;
	      yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) * cosTheta + yc;
	      x0 = (int)rint(xp);
	      y0 = (int)rint(yp);
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		pDst[y*dstW+x] = (D_uint32)pad_val;
	      }
	      else{
		pDst[y*dstW+x] = pSrc[y0*_w+x0];
	      }
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_RGB:
      case DImage_RGB_16:
	fprintf(stderr,
		"rotated_() called for RGB image with non-RGB pad value\n");
	abort();
	break;
      case DImage_flt_multi:
	{
	  float *pSrc;
	  float *pDst;
	  
	  imgDst.create(dstW, dstH, DImage_flt_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_flt(chan);
	    pDst = imgDst.dataPointer_flt(chan);
	    for(int y = 0; y < dstH; ++y){
	      for(int x = 0; x < dstW; ++x){
		xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) *sinTheta+ xc;
		yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) *cosTheta+ yc;
		x0 = (int)rint(xp);
		y0 = (int)rint(yp);
		if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		  pDst[y*dstW+x] = (float)pad_val;
		}
		else{
		  pDst[y*dstW+x] = pSrc[y0*_w+x0];
		}
	      } /* end for x */
	    } /* end for y */
	  }
	}// end case
	break;
      case DImage_dbl_multi:
	{
	  double *pSrc;
	  double *pDst;
	  
	  imgDst.create(dstW, dstH, DImage_dbl_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_dbl(chan);
	    pDst = imgDst.dataPointer_dbl(chan);
	    for(int y = 0; y < dstH; ++y){
	      for(int x = 0; x < dstW; ++x){
		xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) *sinTheta+ xc;
		yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) *cosTheta+ yc;
		x0 = (int)rint(xp);
		y0 = (int)rint(yp);
		if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		  pDst[y*dstW+x] = (double)pad_val;
		}
		else{
		  pDst[y*dstW+x] = pSrc[y0*_w+x0];
		}
	      } /* end for x */
	    } /* end for y */
	  }
	}// end case
	break;
      case DImage_cmplx:
	{
	  std::complex<double> *pSrc;
	  std::complex<double> *pDst;
	  imgDst.create(dstW, dstH, DImage_cmplx, _numChan, _allocMethod);
	  pSrc = (std::complex<double>*)pData;
	  pDst = imgDst.dataPointer_cmplx();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) * sinTheta + xc;
	      yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) * cosTheta + yc;
	      x0 = (int)rint(xp);
	      y0 = (int)rint(yp);
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		pDst[y*dstW+x] = (std::complex<double>)pad_val;
	      }
	      else{
		pDst[y*dstW+x] = pSrc[y0*_w+x0];
	      }
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
    }// end switch
  }

}



///Put a version of this RGB image into imgDst, rotated clockwise ang degrees
/**The image is rotated around its center.  If fResize is true, then the image
 * is padded to make sure it all fits.  Pixels that aren't part of the
 * original image are filled with RGB value pad_R,pad_G,pad_B.
 */
void DImage::rotated_(DImage &imgDst, double ang,
		      int pad_R, int pad_G, int pad_B, bool fResize,
		      DImageTransformMode mode) const{
  double dx, dy; /* delta x and y for bilinear interpolation */
  double xp, yp; /* Xprime and Yprime for rotation */
  double xc, yc; /* position of center of image */
  double xOffs, yOffs;
  DSize rotSize; /* size (width, height) of rotated image */
  int dstW, dstH; /* width and height of rotated image */
  double cosTheta, sinTheta;
  int x0, y0;
  /* temporary variables used for bilinear interpolation: */
  double f00_R, f01_R, f10_R, f11_R, fx0_R, fx1_R, fxy_R;
  double f00_G, f01_G, f10_G, f11_G, fx0_G, fx1_G, fxy_G;
  double f00_B, f01_B, f10_B, f11_B, fx0_B, fx1_B, fxy_B;

  cosTheta = cos(DImage_radFromDeg(-1. * ang));
  sinTheta = sin(DImage_radFromDeg(-1. * ang));

  if(fResize){
    rotSize = getRotatedSize(ang);
    dstW = rotSize.w;
    dstH = rotSize.h;
  }
  else{
    dstW = _w;
    dstH = _h;
  }
  xOffs = (dstW - _w) / 2.;
  yOffs = (dstH - _h) / 2.;
  xc = _w / 2.;
  yc = _h / 2.;

  if(DImageTransSmooth == mode){
    //TODO: make this faster (check reference used by pnmrotate:
    // "A Fast Algorithm for General Raster Rotation" by Alan Paeth,
    // Graphics Interface '86, pp. 77-81.
    switch(this->_imgType){
      case DImage_u8:
      case DImage_u16:
      case DImage_u32:
      case DImage_flt_multi:
      case DImage_dbl_multi:
      case DImage_cmplx:
	fprintf(stderr, "rotated_() called for non-RGB image with RGB pad\n");
	abort();
	break;
      case DImage_RGB:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) * sinTheta + xc;
	      yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) * cosTheta + yc;
	      x0 = (int)floor(xp);
	      y0 = (int)floor(yp);
	      dx = xp - x0;
	      dy = yp - y0;
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		f00_R = pad_R;
		f00_G = pad_G;
		f00_B = pad_B;
	      }
	      else{
		f00_R = pSrc[(y0*_w+x0)*3];
		f00_G = pSrc[(y0*_w+x0)*3+1];
		f00_B = pSrc[(y0*_w+x0)*3+2];
	      }
	      if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h)){
		f10_R = pad_R;
		f10_G = pad_G;
		f10_B = pad_B;
	      }
	      else{
		f10_R = pSrc[(y0*_w+(x0+1))*3];
		f10_G = pSrc[(y0*_w+(x0+1))*3+1];
		f10_B = pSrc[(y0*_w+(x0+1))*3+2];
	      }
	      if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h)){
		f01_R = pad_R;
		f01_G = pad_G;
		f01_B = pad_B;
	      }
	      else{
		f01_R = pSrc[((y0+1)*_w+x0)*3];
		f01_G = pSrc[((y0+1)*_w+x0)*3+1];
		f01_B = pSrc[((y0+1)*_w+x0)*3+2];
	      }
	      if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1) >=_h)){
		f11_R = pad_R;
		f11_G = pad_G;
		f11_B = pad_B;
	      }
	      else{
		f11_R = pSrc[((y0+1)*_w+(x0+1))*3];
		f11_G = pSrc[((y0+1)*_w+(x0+1))*3+1];
		f11_B = pSrc[((y0+1)*_w+(x0+1))*3+2];
	      }
	      
	      fx0_R = f00_R + dx * (f10_R - f00_R);
	      fx0_G = f00_G + dx * (f10_G - f00_G);
	      fx0_B = f00_B + dx * (f10_B - f00_B);
	      fx1_R = f01_R + dx * (f11_R - f01_R);
	      fx1_G = f01_G + dx * (f11_G - f01_G);
	      fx1_B = f01_B + dx * (f11_B - f01_B);
	      fxy_R = fx0_R + dy * (fx1_R - fx0_R);
	      fxy_G = fx0_G + dy * (fx1_G - fx0_G);
	      fxy_B = fx0_B + dy * (fx1_B - fx0_B);
	      
	      if(fxy_R < 0.)
		fxy_R = 0;
	      if(fxy_R > 255.)
		fxy_R = 255.;
	      if(fxy_G < 0.)
		fxy_G = 0;
	      if(fxy_G > 255.)
		fxy_G = 255.;
	      if(fxy_B < 0.)
		fxy_B = 0;
	      if(fxy_B > 255.)
		fxy_B = 255.;
	      
	      pDst[(y*dstW+x)*3] = (D_uint8)(D_uint32)fxy_R;
	      pDst[(y*dstW+x)*3+1] = (D_uint8)(D_uint32)fxy_G;
	      pDst[(y*dstW+x)*3+2] = (D_uint8)(D_uint32)fxy_B;
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_RGB_16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) * sinTheta + xc;
	      yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) * cosTheta + yc;
	      x0 = (int)floor(xp);
	      y0 = (int)floor(yp);
	      dx = xp - x0;
	      dy = yp - y0;
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		f00_R = pad_R;
		f00_G = pad_G;
		f00_B = pad_B;
	      }
	      else{
		f00_R = pSrc[(y0*_w+x0)*3];
		f00_G = pSrc[(y0*_w+x0)*3+1];
		f00_B = pSrc[(y0*_w+x0)*3+2];
	      }
	      if(((x0+1) < 0) || (y0 < 0) || ((x0+1) >= _w) || (y0 >= _h)){
		f10_R = pad_R;
		f10_G = pad_G;
		f10_B = pad_B;
	      }
	      else{
		f10_R = pSrc[(y0*_w+(x0+1))*3];
		f10_G = pSrc[(y0*_w+(x0+1))*3+1];
		f10_B = pSrc[(y0*_w+(x0+1))*3+2];
	      }
	      if((x0 < 0) || ((y0+1) < 0) || (x0 >= _w) || ((y0+1) >= _h)){
		f01_R = pad_R;
		f01_G = pad_G;
		f01_B = pad_B;
	      }
	      else{
		f01_R = pSrc[((y0+1)*_w+x0)*3];
		f01_G = pSrc[((y0+1)*_w+x0)*3+1];
		f01_B = pSrc[((y0+1)*_w+x0)*3+2];
	      }
	      if(((x0+1) <0) || ((y0+1) <0) || ((x0+1) >= _w) ||((y0+1) >=_h)){
		f11_R = pad_R;
		f11_G = pad_G;
		f11_B = pad_B;
	      }
	      else{
		f11_R = pSrc[((y0+1)*_w+(x0+1))*3];
		f11_G = pSrc[((y0+1)*_w+(x0+1))*3+1];
		f11_B = pSrc[((y0+1)*_w+(x0+1))*3+2];
	      }
	      
	      fx0_R = f00_R + dx * (f10_R - f00_R);
	      fx0_G = f00_G + dx * (f10_G - f00_G);
	      fx0_B = f00_B + dx * (f10_B - f00_B);
	      fx1_R = f01_R + dx * (f11_R - f01_R);
	      fx1_G = f01_G + dx * (f11_G - f01_G);
	      fx1_B = f01_B + dx * (f11_B - f01_B);
	      fxy_R = fx0_R + dy * (fx1_R - fx0_R);
	      fxy_G = fx0_G + dy * (fx1_G - fx0_G);
	      fxy_B = fx0_B + dy * (fx1_B - fx0_B);
	      
	      if(fxy_R < 0.)
		fxy_R = 0;
	      if(fxy_R > 255.)
		fxy_R = 255.;
	      if(fxy_G < 0.)
		fxy_G = 0;
	      if(fxy_G > 255.)
		fxy_G = 255.;
	      if(fxy_B < 0.)
		fxy_B = 0;
	      if(fxy_B > 255.)
		fxy_B = 255.;
	      
	      pDst[(y*dstW+x)*3] = (D_uint16)(D_uint32)fxy_R;
	      pDst[(y*dstW+x)*3+1] = (D_uint16)(D_uint32)fxy_G;
	      pDst[(y*dstW+x)*3+2] = (D_uint16)(D_uint32)fxy_B;
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
    }// end switch
  } // end if sample transform mode
  else if(DImageTransSample){
    switch(this->_imgType){
      case DImage_u8:
      case DImage_u16:
      case DImage_u32:
      case DImage_flt_multi:
      case DImage_dbl_multi:
      case DImage_cmplx:
	fprintf(stderr, "rotated_() called for non-RGB image with RGB pad\n");
	abort();
	break;
      case DImage_RGB:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) * sinTheta + xc;
	      yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) * cosTheta + yc;
	      x0 = (int)rint(xp);
	      y0 = (int)rint(yp);
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		pDst[(y*dstW+x)*3] = (D_uint8)pad_R;
		pDst[(y*dstW+x)*3+1] = (D_uint8)pad_G;
		pDst[(y*dstW+x)*3+2] = (D_uint8)pad_B;
	      }
	      else{
		pDst[(y*dstW+x)*3] = pSrc[(y0*_w+x0)*3];
		pDst[(y*dstW+x)*3+1] = pSrc[(y0*_w+x0)*3+1];
		pDst[(y*dstW+x)*3+2] = pSrc[(y0*_w+x0)*3+2];
	      }
	    } /* end for x */
	  } /* end for y */
	}// end case
	break;
      case DImage_RGB_16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  imgDst.create(dstW, dstH, DImage_RGB_16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    for(int x = 0; x < dstW; ++x){
	      xp = ((x-xOffs)-xc) * cosTheta - ((y-yOffs)- yc) * sinTheta + xc;
	      yp = ((x-xOffs)-xc) * sinTheta + ((y-yOffs)- yc) * cosTheta + yc;
	      x0 = (int)rint(xp);
	      y0 = (int)rint(yp);
	      if((x0 < 0) || (y0 < 0) || (x0 >= _w) || (y0 >= _h)){
		//TODO: move casts outside of loop by creating a temporary
		// variable of the correct type with the casted value stored.
		// Do the same for all other similar code
		pDst[(y*dstW+x)*3] = (D_uint16)pad_R;
		pDst[(y*dstW+x)*3+1] = (D_uint16)pad_G;
		pDst[(y*dstW+x)*3+2] = (D_uint16)pad_B;
	      }
	      else{
		pDst[(y*dstW+x)*3] = pSrc[(y0*_w+x0)*3];
		pDst[(y*dstW+x)*3+1] = pSrc[(y0*_w+x0)*3+1];
		pDst[(y*dstW+x)*3+2] = pSrc[(y0*_w+x0)*3+2];
	      }
	    } /* end for x */
	  } /* end for y */
	}
	break;
    }//end switch
  } // end if(mode==...
}



///Returns the size this image would be if rotated clockwise by ang degrees
DSize DImage::getRotatedSize(double ang) const{
  //TODO: This is probably not the most efficient way to do this.
  //TODO: The size may be off by a pixel or so due to not being consistent
  // in how I round/truncate numbers here and in the actual rotation?
  double xc, yc; /* position of center of image */
  double xOffs, yOffs;
  double minx, miny, maxx, maxy;
  double tlxp, tlyp, trxp, tryp, blxp, blyp, brxp, bryp;
  double cosTheta, sinTheta;
  DSize rs;

  cosTheta= cos(DImage_radFromDeg(ang));
  sinTheta= sin(DImage_radFromDeg(ang));

  xc = _w / 2.;
  yc = _h / 2.;
  /* top left    x prime,yprime */
  tlxp = (0-xc) * cosTheta - (0 - yc) * sinTheta + xc;
  tlyp = (0-xc) * sinTheta + (0 - yc) * cosTheta + yc;
  /* top right    x prime,yprime */
  trxp = ((_w-1)-xc) * cosTheta - (0 - yc) * sinTheta + xc;
  tryp = ((_w-1)-xc) * sinTheta + (0 - yc) * cosTheta + yc;
  /* bottom left    x prime,yprime */
  blxp = (0-xc) * cosTheta - ((_h-1) - yc) * sinTheta + xc;
  blyp = (0-xc) * sinTheta + ((_h-1) - yc) * cosTheta + yc;
  /* bottom right    x prime,yprime */
  brxp = ((_w-1)-xc) * cosTheta - ((_h-1) - yc) * sinTheta + xc;
  bryp = ((_w-1)-xc) * sinTheta + ((_h-1) - yc) * cosTheta + yc;

  minx = DImage_minDbl(DImage_minDbl(tlxp,trxp),DImage_minDbl(blxp,brxp));
  miny = DImage_minDbl(DImage_minDbl(tlyp,tryp),DImage_minDbl(blyp,bryp));
  maxx = DImage_maxDbl(DImage_maxDbl(tlxp,trxp),DImage_maxDbl(blxp,brxp));
  maxy = DImage_maxDbl(DImage_maxDbl(tlyp,tryp),DImage_maxDbl(blyp,bryp));

  xOffs = -1*minx;
  yOffs = -1*miny;

  rs.w = _w + 2 * (int)(xOffs+0.999);
  rs.h = _h + 2 * (int)(yOffs+0.999);
  return rs;
}

///Put this image rotated clockwise by ang (multiple of 90) degrees into imgDst
/**Ang should be (either positive or negative) 90, 180, or 270  */
void DImage::rotate90_(DImage &imgDst, int ang) const{
  int w, h; 
  int sidx, didx; // src and dst idx
  //TODO: finish this function

  w = _w;
  h = _h;
  if(ang < 0)
    ang += 360;
  switch(ang){
    case 90:
      imgDst.create(_h, _w, _imgType, _numChan, _allocMethod);
      switch(_imgType){
        case DImage_u8:
	  {
	    D_uint8 *ps8;
	    D_uint8 *pd8;
	    ps8 = dataPointer_u8();
	    pd8 = imgDst.dataPointer_u8();
	    for(int sy = 0, sidx=0; sy < h; ++sy){
	      didx = h - 1 - sy;
	      for(int sx = 0; sx < w; ++sx, ++sidx, didx += h){
		pd8[didx] = ps8[sidx];
	      }
	    }
	  }
	  break;
        case DImage_u16:
	  {
	    D_uint16 *ps16;
	    D_uint16 *pd16;
	    ps16 = dataPointer_u16();
	    pd16 = imgDst.dataPointer_u16();
	    for(int sy = 0, sidx=0; sy < h; ++sy){
	      didx = h - 1 - sy;
	      for(int sx = 0; sx < w; ++sx, ++sidx, didx += h){
		pd16[didx] = ps16[sidx];
	      }
	    }
	  }
	  break;
        case DImage_u32:
	  {
	    D_uint32 *ps32;
	    D_uint32 *pd32;
	    ps32 = dataPointer_u32();
	    pd32 = imgDst.dataPointer_u32();
	    for(int sy = 0, sidx=0; sy < h; ++sy){
	      didx = h - 1 - sy;
	      for(int sx = 0; sx < w; ++sx, ++sidx, didx += h){
		pd32[didx] = ps32[sidx];
	      }
	    }
	  }
	  break;
        case DImage_RGB:
	  {
	    D_uint8 *ps8;
	    D_uint8 *pd8;
	    ps8 = dataPointer_u8();
	    pd8 = imgDst.dataPointer_u8();
	    for(int sy = 0, sidx=0; sy < h; ++sy){
	      didx = 3*(h - 1 - sy);
	      for(int sx = 0; sx < w; ++sx, sidx+=3, didx += h*3){
		pd8[didx] = ps8[sidx];
		pd8[didx+1] = ps8[sidx+1];
		pd8[didx+2] = ps8[sidx+2];
	      }
	    }
	  }
	  break;
        case DImage_RGB_16:
	  {
	    D_uint16 *ps16;
	    D_uint16 *pd16;
	    ps16 = dataPointer_u16();
	    pd16 = imgDst.dataPointer_u16();
	    for(int sy = 0, sidx=0; sy < h; ++sy){
	      didx = 3*(h - 1 - sy);
	      for(int sx = 0; sx < w; ++sx, sidx+=3, didx += h*3){
		pd16[didx] = ps16[sidx];
		pd16[didx+1] = ps16[sidx+1];
		pd16[didx+2] = ps16[sidx+2];
	      }
	    }
	  }
	  break;
        case DImage_flt_multi:
        case DImage_dbl_multi:
        case DImage_cmplx:
        default:
	  fprintf(stderr, "DImage::rotate90_() unsupported type\n");
	  abort();
      }
      break;
    case 270:
      imgDst.create(_h, _w, _imgType, _numChan, _allocMethod);
      switch(_imgType){
        case DImage_u8:
	  {
	    D_uint8 *ps8;
	    D_uint8 *pd8;
	    ps8 = dataPointer_u8();
	    pd8 = imgDst.dataPointer_u8();
	    for(int sy = 0, sidx=0; sy < h; ++sy){
	      didx = sy + h * (w-1);
	      for(int sx = 0; sx < w; ++sx, ++sidx, didx -= h){
		pd8[didx] = ps8[sidx];
	      }
	    }
	  }
	  break;
        case DImage_u16:
	  {
	    D_uint16 *ps16;
	    D_uint16 *pd16;
	    ps16 = dataPointer_u16();
	    pd16 = imgDst.dataPointer_u16();
	    for(int sy = 0, sidx=0; sy < h; ++sy){
	      didx = sy + h * (w-1);
	      for(int sx = 0; sx < w; ++sx, ++sidx, didx -= h){
		pd16[didx] = ps16[sidx];
	      }
	    }
	  }
	  break;
        case DImage_u32:
	  {
	    D_uint32 *ps32;
	    D_uint32 *pd32;
	    ps32 = dataPointer_u32();
	    pd32 = imgDst.dataPointer_u32();
	    for(int sy = 0, sidx=0; sy < h; ++sy){
	      didx = sy + h * (w-1);
	      for(int sx = 0; sx < w; ++sx, ++sidx, didx -= h){
		pd32[didx] = ps32[sidx];
	      }
	    }
	  }
	  break;
        case DImage_RGB:
	  {
	    D_uint8 *ps8;
	    D_uint8 *pd8;
	    ps8 = dataPointer_u8();
	    pd8 = imgDst.dataPointer_u8();
	    for(int sy = 0, sidx=0; sy < h; ++sy){
	      didx = 3*(sy + h * (w-1));
	      for(int sx = 0; sx < w; ++sx, sidx+=3, didx -= h*3){
		pd8[didx] = ps8[sidx];
		pd8[didx+1] = ps8[sidx+1];
		pd8[didx+2] = ps8[sidx+2];
	      }
	    }
	  }
	  break;
        case DImage_RGB_16:
	  {
	    D_uint16 *ps16;
	    D_uint16 *pd16;
	    ps16 = dataPointer_u16();
	    pd16 = imgDst.dataPointer_u16();
	    for(int sy = 0, sidx=0; sy < h; ++sy){
	      didx = 3*(sy + h * (w-1));
	      for(int sx = 0; sx < w; ++sx, sidx+=3, didx -= h*3){
		pd16[didx] = ps16[sidx];
		pd16[didx+1] = ps16[sidx+1];
		pd16[didx+2] = ps16[sidx+2];
	      }
	    }
	  }
        case DImage_flt_multi:
        case DImage_dbl_multi:
        case DImage_cmplx:
        default:
	  fprintf(stderr, "DImage::rotate90_() unsupported type\n");
	  abort();
      }
      break;
    case 180:
      imgDst.create(_w, _h, _imgType, _numChan, _allocMethod);
      switch(_imgType){
        case DImage_u8:
	  {
	    D_uint8 *ps8;
	    D_uint8 *pd8;
	    ps8 = dataPointer_u8();
	    pd8 = imgDst.dataPointer_u8();
	    for(int sy = 0, sidx=0, didx=w*h-1; sy < h; ++sy){
	      for(int sx = 0; sx < w; ++sx, ++sidx, --didx){
		pd8[didx] = ps8[sidx];
	      }
	    }
	  }
	  break;
        case DImage_u16:
	  {
	    D_uint16 *ps16;
	    D_uint16 *pd16;
	    ps16 = dataPointer_u16();
	    pd16 = imgDst.dataPointer_u16();
	    for(int sy = 0, sidx=0, didx=w*h-1; sy < h; ++sy){
	      for(int sx = 0; sx < w; ++sx, ++sidx, --didx){
		pd16[didx] = ps16[sidx];
	      }
	    }
	  }
	  break;
        case DImage_u32:
	  {
	    D_uint32 *ps32;
	    D_uint32 *pd32;
	    ps32 = dataPointer_u32();
	    pd32 = imgDst.dataPointer_u32();
	    for(int sy = 0, sidx=0, didx=w*h-1; sy < h; ++sy){
	      for(int sx = 0; sx < w; ++sx, ++sidx, --didx){
		pd32[didx] = ps32[sidx];
	      }
	    }
	  }
	  break;
        case DImage_RGB:
	  {
	    D_uint8 *ps8;
	    D_uint8 *pd8;
	    ps8 = dataPointer_u8();
	    pd8 = imgDst.dataPointer_u8();
	    for(int sy = 0, sidx=0, didx=3*(w*h-1); sy < h; ++sy){
	      for(int sx = 0; sx < w; ++sx, sidx+=3, didx -= 3){
		pd8[didx] = ps8[sidx];
		pd8[didx+1] = ps8[sidx+1];
		pd8[didx+2] = ps8[sidx+2];
	      }
	    }
	  }
	  break;
        case DImage_RGB_16:
	  {
	    D_uint16 *ps16;
	    D_uint16 *pd16;
	    ps16 = dataPointer_u16();
	    pd16 = imgDst.dataPointer_u16();
	    for(int sy = 0, sidx=0, didx=3*(w*h-1); sy < h; ++sy){
	      for(int sx = 0; sx < w; ++sx, sidx+=3, didx -= 3){
		pd16[didx] = ps16[sidx];
		pd16[didx+1] = ps16[sidx+1];
		pd16[didx+2] = ps16[sidx+2];
	      }
	    }
	  }
	  break;
        case DImage_flt_multi:
        case DImage_dbl_multi:
        case DImage_cmplx:
        default:
	  fprintf(stderr, "DImage::rotate90_() unsupported type\n");
	  abort();
      }
      break;
    default:
      fprintf(stderr, "DImage::rotate90_() only supports angles of + or - "
	      "90, 180, and 270 (%d was used)\n", ang);
      abort();
  }// end switch(ang)

}



///Returns a copy of this image that has been sheared ang degrees horizontally
DImage DImage::shearedH(double ang, double pad_val, bool fResize,
			DImageTransformMode mode) const{
  DImage dst;
  shearedH_(dst, ang, pad_val, fResize, mode);
  return dst;
}
///Returns a copy of this image that has been sheared ang degrees horizontally
DImage DImage::shearedH(double ang, int pad_R,int pad_G,int pad_B,
			bool fResize, DImageTransformMode mode) const{
  DImage dst;
  shearedH_(dst, ang, pad_R, pad_G, pad_B, fResize, mode);
  return dst;
}
///Returns the width that image would be if horizontally sheared ang degrees
int DImage::getShearedHWidth(double ang) const{
  if((ang >= 90.) || (ang <= -90.))
    return -1;
  return abs((int)(_h*tan(DImage_radFromDeg(ang))));
}


///Shear this image by ang degrees horizontally and put it in imgDst
void DImage::shearedH_(DImage &imgDst, double ang,
		       int pad_R, int pad_G, int pad_B,
		       bool fResize, DImageTransformMode mode) const{
  double xoffs; // x-offset within source image (relative to current x in dst)
  double dx;    // change in xoffs per scanline
  double srcX;
  int iSrcX, iSrcX_p1; // the left and right pixel index, respectively
  double hTanAng;
  double fracRight, fracLeft; // fraction belonging to right/left pixel
  int dstW, dstH;

#ifdef DEBUG
  if((ang >= 90.) || (ang <= -90.)){
    fprintf(stderr, "DImage::shearedH_() ang must be in range (-90 .. 90)\n");
    abort();
  }
#endif

  hTanAng = _h*tan(DImage_radFromDeg(ang));
  
  if(false == fResize){
    xoffs = hTanAng / 2.;
    dstW = _w;
  }
  else{
    xoffs = (ang >= 0.) ? 0. : hTanAng;
    dstW = _w + (int)(fabs(hTanAng) + 0.999);
  }
  dx = -1. * tan(DImage_radFromDeg(ang));
  dstH = _h;


  if(DImageTransSmooth == mode){
    switch(this->_imgType){
      case DImage_u8:
      case DImage_u16:
      case DImage_u32:
      case DImage_flt_multi:
      case DImage_dbl_multi:
      case DImage_cmplx:
	fprintf(stderr, "shearedH_() called for non-RGB image with RGB pad\n");
	abort();
	break;
      case DImage_RGB:
	{
	  D_uint8 pxlLeft_R, pxlRight_R;
	  D_uint8 pxlLeft_G, pxlRight_G;
	  D_uint8 pxlLeft_B, pxlRight_B;
	  D_uint8 *pSrcRow;
	  D_uint8 *pDst;
	  D_uint8 padR,padG,padB;

	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrcRow = pData;
	  pDst = imgDst.dataPointer_u8();
	  padR = (D_uint8)pad_R;
	  padG = (D_uint8)pad_G;
	  padB = (D_uint8)pad_B;

	  for(int y = 0; y < dstH; ++y){
	    srcX = xoffs + y*dx;
	    iSrcX = (int)floor(srcX); // the left pixel index
	    iSrcX_p1 = iSrcX+1; // right pixel index
	    fracRight = srcX-iSrcX;
	    fracLeft = 1.-fracRight;
	    for(int x = 0; x < dstW; ++x){
	      if( (iSrcX >= 0) && (iSrcX < _w)){
		pxlLeft_R = pSrcRow[iSrcX*3];
		pxlLeft_G = pSrcRow[iSrcX*3+1];
		pxlLeft_B = pSrcRow[iSrcX*3+2];
	      }
	      else{
		pxlLeft_R = padR;
		pxlLeft_G = padG;
		pxlLeft_B = padB;
	      }
	      if( (iSrcX_p1 >= 0) && (iSrcX_p1 < _w)){
		pxlRight_R = pSrcRow[iSrcX_p1*3];
		pxlRight_G = pSrcRow[iSrcX_p1*3+1];
		pxlRight_B = pSrcRow[iSrcX_p1*3+2];
	      }
	      else{
		pxlRight_R = padR;
		pxlRight_G = padG;
		pxlRight_B = padB;
	      }
	      (*pDst) = (D_uint8)(pxlRight_R*fracRight + pxlLeft_R*fracLeft);
	      pDst[1] = (D_uint8)(pxlRight_G*fracRight + pxlLeft_G*fracLeft);
	      pDst[2] = (D_uint8)(pxlRight_B*fracRight + pxlLeft_B*fracLeft);
	      pDst += 3;
	      ++iSrcX;
	      ++iSrcX_p1;
	    } // end for x
	    pSrcRow += _w*3;
	  }//end for y
	}//end case
	break;
      case DImage_RGB_16:
	{
	  D_uint16 pxlLeft_R, pxlRight_R;
	  D_uint16 pxlLeft_G, pxlRight_G;
	  D_uint16 pxlLeft_B, pxlRight_B;
	  D_uint16 *pSrcRow;
	  D_uint16 *pDst;
	  D_uint16 padR,padG,padB;

	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrcRow = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  padR = (D_uint16)pad_R;
	  padG = (D_uint16)pad_G;
	  padB = (D_uint16)pad_B;

	  for(int y = 0; y < dstH; ++y){
	    srcX = xoffs + y*dx;
	    iSrcX = (int)floor(srcX); // the left pixel index
	    iSrcX_p1 = iSrcX+1; // right pixel index
	    fracRight = srcX-iSrcX;
	    fracLeft = 1.-fracRight;
	    for(int x = 0; x < dstW; ++x){
	      if( (iSrcX >= 0) && (iSrcX < _w)){
		pxlLeft_R = pSrcRow[iSrcX*3];
		pxlLeft_G = pSrcRow[iSrcX*3+1];
		pxlLeft_B = pSrcRow[iSrcX*3+2];
	      }
	      else{
		pxlLeft_R = padR;
		pxlLeft_G = padG;
		pxlLeft_B = padB;
	      }
	      if( (iSrcX_p1 >= 0) && (iSrcX_p1 < _w)){
		pxlRight_R = pSrcRow[iSrcX_p1*3];
		pxlRight_G = pSrcRow[iSrcX_p1*3+1];
		pxlRight_B = pSrcRow[iSrcX_p1*3+2];
	      }
	      else{
		pxlRight_R = padR;
		pxlRight_G = padG;
		pxlRight_B = padB;
	      }
	      (*pDst) = (D_uint16)(pxlRight_R*fracRight + pxlLeft_R*fracLeft);
	      pDst[1] = (D_uint16)(pxlRight_G*fracRight + pxlLeft_G*fracLeft);
	      pDst[2] = (D_uint16)(pxlRight_B*fracRight + pxlLeft_B*fracLeft);
	      pDst += 3;
	      ++iSrcX;
	      ++iSrcX_p1;
	    } // end for x
	    pSrcRow += _w*3;
	  }//end for y
	}//end case
	break;
    }// end switch
  }//end if
  else if(DImageTransSample == mode){
    switch(this->_imgType){
      case DImage_u8:
      case DImage_u16:
      case DImage_u32:
      case DImage_flt_multi:
      case DImage_dbl_multi:
      case DImage_cmplx:
	fprintf(stderr, "shearedH_() called for non-RGB image with RGB pad\n");
	abort();
	break;
      case DImage_RGB:
	{
	  D_uint8 *pSrcRow;
	  D_uint8 *pDst;
	  D_uint8 padR,padG,padB;

	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrcRow = pData;
	  pDst = imgDst.dataPointer_u8();
	  padR = (D_uint8)pad_R;
	  padG = (D_uint8)pad_G;
	  padB = (D_uint8)pad_B;

	  for(int y = 0; y < dstH; ++y){
	    iSrcX = (int)round(xoffs + y*dx);
	    for(int x = 0; x < dstW; ++x){
	      if( (iSrcX >= 0) && (iSrcX < _w)){
		(*pDst) = pSrcRow[iSrcX*3];
		pDst[1] = pSrcRow[iSrcX*3+1];
		pDst[2] = pSrcRow[iSrcX*3+2];
	      }
	      else{
		(*pDst) = padR;
		pDst[1] = padG;
		pDst[2] = padB;
	      }
	      pDst += 3;
	      ++iSrcX;
	    } // end for x
	    pSrcRow += _w*3;
	  }//end for y
	}//end case
	break;
      case DImage_RGB_16:
	{
	  D_uint16 *pSrcRow;
	  D_uint16 *pDst;
	  D_uint16 padR,padG,padB;

	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrcRow = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  padR = (D_uint16)pad_R;
	  padG = (D_uint16)pad_G;
	  padB = (D_uint16)pad_B;

	  for(int y = 0; y < dstH; ++y){
	    iSrcX = (int)round(xoffs + y*dx);
	    for(int x = 0; x < dstW; ++x){
	      if( (iSrcX >= 0) && (iSrcX < _w)){
		(*pDst) = pSrcRow[iSrcX*3];
		pDst[1] = pSrcRow[iSrcX*3+1];
		pDst[2] = pSrcRow[iSrcX*3+2];
	      }
	      else{
		(*pDst) = padR;
		pDst[1] = padG;
		pDst[2] = padB;
	      }
	      pDst += 3;
	      ++iSrcX;
	    } // end for x
	    pSrcRow += _w*3;
	  }//end for y
	}//end case
	break;
    }// end switch
  }
}

///Shear this image by ang degrees horizontally and put it in imgDst
void DImage::shearedH_(DImage &imgDst, double ang, double pad_val,
		       bool fResize, DImageTransformMode mode) const{
  double xoffs; // x-offset within source image (relative to current x in dst)
  double dx;    // change in xoffs per scanline
  double srcX;
  int iSrcX, iSrcX_p1; // the left and right pixel index, respectively
  double hTanAng;
  double fracRight, fracLeft; // fraction belonging to right/left pixel
  int dstW, dstH;

#ifdef DEBUG
  if((ang >= 90.) || (ang <= -90.)){
    fprintf(stderr, "DImage::shearedH_() ang must be in range (-90 .. 90)\n");
    abort();
  }
#endif

  hTanAng = _h*tan(DImage_radFromDeg(ang));
  
  if(false == fResize){
    xoffs = hTanAng / 2.;
    dstW = _w;
  }
  else{
    xoffs = (ang >= 0.) ? 0. : hTanAng;
    dstW = _w + (int)(fabs(hTanAng) + 0.999);
  }
  dx = -1. * tan(DImage_radFromDeg(ang));
  dstH = _h;


  if(DImageTransSmooth == mode){
    switch(this->_imgType){
      case DImage_u8:
	{
	  D_uint8 pxlLeft, pxlRight;
	  D_uint8 *pSrcRow;
	  D_uint8 *pDst;
	  D_uint8 pad;
	  imgDst.create(dstW, dstH, DImage_u8, _numChan, _allocMethod);
	  pSrcRow = pData;
	  pDst = imgDst.dataPointer_u8();
	  pad = (D_uint8)pad_val;

	  for(int y = 0; y < dstH; ++y){
	    srcX = xoffs + y*dx;
	    iSrcX = (int)floor(srcX); // the left pixel index
	    iSrcX_p1 = iSrcX+1; // right pixel index
	    fracRight = srcX-iSrcX;
	    fracLeft = 1.-fracRight;
	    for(int x = 0; x < dstW; ++x){
	      if( (iSrcX >= 0) && (iSrcX < _w)){
		pxlLeft = pSrcRow[iSrcX];
	      }
	      else{
		pxlLeft = pad;
	      }
	      if( (iSrcX_p1 >= 0) && (iSrcX_p1 < _w)){
		pxlRight = pSrcRow[iSrcX_p1];
	      }
	      else{
		pxlRight = pad;
	      }
	      (*pDst) = (D_uint8)(pxlRight*fracRight + pxlLeft*fracLeft);
	      ++pDst;
	      ++iSrcX;
	      ++iSrcX_p1;
	    } // end for x
	    pSrcRow += _w;
	  }//end for y
	}//end case
	break;
      case DImage_u16:
	{
	  D_uint16 pxlLeft, pxlRight;
	  D_uint16 *pSrcRow;
	  D_uint16 *pDst;
	  D_uint16 pad;
	  imgDst.create(dstW, dstH, DImage_u16, _numChan, _allocMethod);
	  pSrcRow = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  pad = (D_uint16)pad_val;

	  for(int y = 0; y < dstH; ++y){
	    srcX = xoffs + y*dx;
	    iSrcX = (int)floor(srcX); // the left pixel index
	    iSrcX_p1 = iSrcX+1; // right pixel index
	    fracRight = srcX-iSrcX;
	    fracLeft = 1.-fracRight;
	    for(int x = 0; x < dstW; ++x){
	      if( (iSrcX >= 0) && (iSrcX < _w)){
		pxlLeft = pSrcRow[iSrcX];
	      }
	      else{
		pxlLeft = pad;
	      }
	      if( (iSrcX_p1 >= 0) && (iSrcX_p1 < _w)){
		pxlRight = pSrcRow[iSrcX_p1];
	      }
	      else{
		pxlRight = pad;
	      }
	      (*pDst) = (D_uint16)(pxlRight*fracRight + pxlLeft*fracLeft);
	      ++pDst;
	      ++iSrcX;
	      ++iSrcX_p1;
	    } // end for x
	    pSrcRow += _w;
	  }//end for y
	}//end case
	break;
      case DImage_u32:
	{
	  D_uint32 pxlLeft, pxlRight;
	  D_uint32 *pSrcRow;
	  D_uint32 *pDst;
	  D_uint32 pad;
	  imgDst.create(dstW, dstH, DImage_u32, _numChan, _allocMethod);
	  pSrcRow = (D_uint32*)pData;
	  pDst = imgDst.dataPointer_u32();
	  pad = (D_uint32)pad_val;

	  for(int y = 0; y < dstH; ++y){
	    srcX = xoffs + y*dx;
	    iSrcX = (int)floor(srcX); // the left pixel index
	    iSrcX_p1 = iSrcX+1; // right pixel index
	    fracRight = srcX-iSrcX;
	    fracLeft = 1.-fracRight;
	    for(int x = 0; x < dstW; ++x){
	      if( (iSrcX >= 0) && (iSrcX < _w)){
		pxlLeft = pSrcRow[iSrcX];
	      }
	      else{
		pxlLeft = pad;
	      }
	      if( (iSrcX_p1 >= 0) && (iSrcX_p1 < _w)){
		pxlRight = pSrcRow[iSrcX_p1];
	      }
	      else{
		pxlRight = pad;
	      }
	      (*pDst) = (D_uint32)(pxlRight*fracRight + pxlLeft*fracLeft);
	      ++pDst;
	      ++iSrcX;
	      ++iSrcX_p1;
	    } // end for x
	    pSrcRow += _w;
	  }//end for y
	}//end case
	break;
      case DImage_flt_multi:
	{
	  float pxlLeft, pxlRight;
	  float *pSrcRow;
	  float *pDst;
	  float pad;
	  imgDst.create(dstW, dstH, DImage_flt_multi, _numChan, _allocMethod);
	  pad = (float)pad_val;

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrcRow = this->dataPointer_flt(chan);
	    pDst = imgDst.dataPointer_flt(chan);


	    for(int y = 0; y < dstH; ++y){
	      srcX = xoffs + y*dx;
	      iSrcX = (int)floor(srcX); // the left pixel index
	      iSrcX_p1 = iSrcX+1; // right pixel index
	      fracRight = srcX-iSrcX;
	      fracLeft = 1.-fracRight;
	      for(int x = 0; x < dstW; ++x){
		if( (iSrcX >= 0) && (iSrcX < _w)){
		  pxlLeft = pSrcRow[iSrcX];
		}
		else{
		  pxlLeft = pad;
		}
		if( (iSrcX_p1 >= 0) && (iSrcX_p1 < _w)){
		  pxlRight = pSrcRow[iSrcX_p1];
		}
		else{
		  pxlRight = pad;
		}
		(*pDst) = (float)(pxlRight*fracRight + pxlLeft*fracLeft);
		++pDst;
		++iSrcX;
		++iSrcX_p1;
	      } // end for x
	      pSrcRow += _w;
	    }//end for y
	  } // end for chan
	}//end case
	break;
      case DImage_dbl_multi:
	{
	  double pxlLeft, pxlRight;
	  double *pSrcRow;
	  double *pDst;
	  double pad;
	  imgDst.create(dstW, dstH, DImage_dbl_multi, _numChan, _allocMethod);
	  pad = (double)pad_val;

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrcRow = this->dataPointer_dbl(chan);
	    pDst = imgDst.dataPointer_dbl(chan);


	    for(int y = 0; y < dstH; ++y){
	      srcX = xoffs + y*dx;
	      iSrcX = (int)floor(srcX); // the left pixel index
	      iSrcX_p1 = iSrcX+1; // right pixel index
	      fracRight = srcX-iSrcX;
	      fracLeft = 1.-fracRight;
	      for(int x = 0; x < dstW; ++x){
		if( (iSrcX >= 0) && (iSrcX < _w)){
		  pxlLeft = pSrcRow[iSrcX];
		}
		else{
		  pxlLeft = pad;
		}
		if( (iSrcX_p1 >= 0) && (iSrcX_p1 < _w)){
		  pxlRight = pSrcRow[iSrcX_p1];
		}
		else{
		  pxlRight = pad;
		}
		(*pDst) = (double)(pxlRight*fracRight + pxlLeft*fracLeft);
		++pDst;
		++iSrcX;
		++iSrcX_p1;
	      } // end for x
	      pSrcRow += _w;
	    }//end for y
	  } // end for chan
	}//end case
	break;
      case DImage_cmplx:
	{
	  std::complex<double> pxlLeft, pxlRight;
	  std::complex<double> *pSrcRow;
	  std::complex<double> *pDst;
	  std::complex<double> pad;
	  imgDst.create(dstW, dstH, DImage_cmplx, _numChan, _allocMethod);
	  pad = (std::complex<double>)pad_val;

	  pSrcRow = this->dataPointer_cmplx();
	  pDst = imgDst.dataPointer_cmplx();
	  
	  
	  for(int y = 0; y < dstH; ++y){
	    srcX = xoffs + y*dx;
	    iSrcX = (int)floor(srcX); // the left pixel index
	    iSrcX_p1 = iSrcX+1; // right pixel index
	    fracRight = srcX-iSrcX;
	    fracLeft = 1.-fracRight;
	    for(int x = 0; x < dstW; ++x){
	      if( (iSrcX >= 0) && (iSrcX < _w)){
		pxlLeft = pSrcRow[iSrcX];
	      }
	      else{
		pxlLeft = pad;
	      }
	      if( (iSrcX_p1 >= 0) && (iSrcX_p1 < _w)){
		pxlRight = pSrcRow[iSrcX_p1];
	      }
	      else{
		pxlRight = pad;
	      }
	      (*pDst) =
		(std::complex<double>)(pxlRight*fracRight +pxlLeft*fracLeft);
	      ++pDst;
	      ++iSrcX;
	      ++iSrcX_p1;
	    } // end for x
	    pSrcRow += _w;
	  }//end for y
	}//end case
	break;
      case DImage_RGB:
      case DImage_RGB_16:
	fprintf(stderr,
		"shearedH_() called for RGB image with non-RGB pad value\n");
	abort();
	break;
    }// end switch
  }//end if
  else if(DImageTransSample == mode){
    switch(this->_imgType){
      case DImage_u8:
	{
	  D_uint8 *pSrcRow;
	  D_uint8 *pDst;
	  D_uint8 pad;
	  imgDst.create(dstW, dstH, DImage_u8, _numChan, _allocMethod);
	  pDst = imgDst.dataPointer_u8();
	  pad = (D_uint8)pad_val;

	  pSrcRow = pData;
	  for(int y = 0; y < dstH; ++y){
	    iSrcX = (int)round(xoffs + y*dx); // the left pixel index
	    for(int x = 0; x < dstW; ++x){
	      if( (iSrcX >= 0) && (iSrcX < _w)){
		(*pDst) = pSrcRow[iSrcX];
	      }
	      else{
		(*pDst) = pad;
	      }
	      ++pDst;
	      ++iSrcX;
	    } // end for x
	    pSrcRow += _w;
	  }//end for y
	}//end case
	break;
      case DImage_u16:
	{
	  D_uint16 *pSrcRow;
	  D_uint16 *pDst;
	  D_uint16 pad;
	  imgDst.create(dstW, dstH, DImage_u16, _numChan, _allocMethod);
	  pDst = imgDst.dataPointer_u16();
	  pad = (D_uint16)pad_val;

	  pSrcRow = (D_uint16*)pData;
	  for(int y = 0; y < dstH; ++y){
	    iSrcX = (int)round(xoffs + y*dx); // the left pixel index
	    for(int x = 0; x < dstW; ++x){
	      if( (iSrcX >= 0) && (iSrcX < _w)){
		(*pDst) = pSrcRow[iSrcX];
	      }
	      else{
		(*pDst) = pad;
	      }
	      ++pDst;
	      ++iSrcX;
	    } // end for x
	    pSrcRow += _w;
	  }//end for y
	}//end case
	break;
      case DImage_u32:
	{
	  D_uint32 *pSrcRow;
	  D_uint32 *pDst;
	  D_uint32 pad;
	  imgDst.create(dstW, dstH, DImage_u32, _numChan, _allocMethod);
	  pDst = imgDst.dataPointer_u32();
	  pad = (D_uint32)pad_val;

	  pSrcRow = (D_uint32*)pData;
	  for(int y = 0; y < dstH; ++y){
	    iSrcX = (int)round(xoffs + y*dx); // the left pixel index
	    for(int x = 0; x < dstW; ++x){
	      if( (iSrcX >= 0) && (iSrcX < _w)){
		(*pDst) = pSrcRow[iSrcX];
	      }
	      else{
		(*pDst) = pad;
	      }
	      ++pDst;
	      ++iSrcX;
	    } // end for x
	    pSrcRow += _w;
	  }//end for y
	}//end case
	break;
      case DImage_flt_multi:
	{
	  float *pSrcRow;
	  float *pDst;
	  float pad;
	  imgDst.create(dstW, dstH, DImage_flt_multi, _numChan, _allocMethod);
	  pad = (float)pad_val;

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrcRow = this->dataPointer_flt(chan);
	    pDst = imgDst.dataPointer_flt(chan);
	    for(int y = 0; y < dstH; ++y){
	      iSrcX = (int)round(xoffs + y*dx); // the left pixel index
	      for(int x = 0; x < dstW; ++x){
		if( (iSrcX >= 0) && (iSrcX < _w)){
		  (*pDst) = pSrcRow[iSrcX];
		}
		else{
		  (*pDst) = pad;
		}
		++pDst;
		++iSrcX;
	      } // end for x
	      pSrcRow += _w;
	    }//end for y
	  }// end for chan
	}//end case
	break;
      case DImage_dbl_multi:
	{
	  double *pSrcRow;
	  double *pDst;
	  double pad;
	  imgDst.create(dstW, dstH, DImage_dbl_multi, _numChan, _allocMethod);
	  pad = (double)pad_val;

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrcRow = this->dataPointer_dbl(chan);
	    pDst = imgDst.dataPointer_dbl(chan);
	    for(int y = 0; y < dstH; ++y){
	      iSrcX = (int)round(xoffs + y*dx); // the left pixel index
	      for(int x = 0; x < dstW; ++x){
		if( (iSrcX >= 0) && (iSrcX < _w)){
		  (*pDst) = pSrcRow[iSrcX];
		}
		else{
		  (*pDst) = pad;
		}
		++pDst;
		++iSrcX;
	      } // end for x
	      pSrcRow += _w;
	    }//end for y
	  }// end for chan
	}//end case
	break;
      case DImage_cmplx:
	{
	  std::complex<double> *pSrcRow;
	  std::complex<double> *pDst;
	  std::complex<double> pad;
	  imgDst.create(dstW, dstH, DImage_cmplx, _numChan, _allocMethod);
	  pDst = imgDst.dataPointer_cmplx();
	  pad = (std::complex<double>)pad_val;

	  pSrcRow = (std::complex<double>*)pData;
	  for(int y = 0; y < dstH; ++y){
	    iSrcX = (int)round(xoffs + y*dx); // the left pixel index
	    for(int x = 0; x < dstW; ++x){
	      if( (iSrcX >= 0) && (iSrcX < _w)){
		(*pDst) = pSrcRow[iSrcX];
	      }
	      else{
		(*pDst) = pad;
	      }
	      ++pDst;
	      ++iSrcX;
	    } // end for x
	    pSrcRow += _w;
	  }//end for y
	}//end case
	break;
      case DImage_RGB:
      case DImage_RGB_16:
	fprintf(stderr,
		"shearedH_() called for RGB image with non-RGB pad value\n");
	abort();
	break;
    }// end switch
  }
 
}





///Returns a copy of this image that has been sheared ang degrees vertically
DImage DImage::shearedV(double ang, double pad_val, bool fResize,
			DImageTransformMode mode) const{
  DImage dst;
  shearedV_(dst, ang, pad_val, fResize, mode);
  return dst;
}
///Returns a copy of this image that has been sheared ang degrees vertically
DImage DImage::shearedV(double ang, int pad_R,int pad_G,int pad_B,
			bool fResize, DImageTransformMode mode) const{
  DImage dst;
  shearedV_(dst, ang, pad_R, pad_G, pad_B, fResize, mode);
  return dst;
}

void DImage::shearedV_(DImage &imgDst, double ang,
		       int pad_R, int pad_G, int pag_B,
		       bool fResize, DImageTransformMode mode) const{
}

///Shear this image by ang degrees horizontally and put it in imgDst
void DImage::shearedV_(DImage &imgDst, double ang, double pad_val,
		       bool fResize, DImageTransformMode mode) const{
}
///Returns the height that image would be if vertically sheared ang degrees
int DImage::getShearedVHeight(double ang) const{
  if((ang >= 90.) || (ang <= -90.))
    return -1;
  return abs((int)(_w*tan(DImage_radFromDeg(ang))));
}

DImage DImage::convertedImgType(DImageType imgType, int numChannels,
				D_uint32 srcChannelMask,
				D_AllocationMethod allocMeth) const{
  DImage dst;
  convertedImgType_(dst, imgType, numChannels, srcChannelMask, allocMeth);
  return dst;
}
///Convert this image's data into a different data type and put it into imgDst
/**If numChannels is -1, then the number of channels will be copied
 * from this image.  The range of the data will NOT be scaled to the
 * destination data range before conversion, so you may want to call
 * setDataRange first, for example, or scale the data values before
 * conversion.  For example, when converting from float to _u8, you
 * may want to first set the data range to [0..255].  When converting
 * from RGB_16 to RGB, you may want to dividePixelsByValue(256) first.
 * If you need to do that, and don't want your original data changed,
 * you will need to make a copy first.  If you are going from _u8 to
 * _u16, you may want to do the conversion, and then
 * multiplyPixelsByValue(256) after the conversion to use the full
 * range.
 *
 * \warning The destination image can't be the source image also.
 *
 * \warning Not all conversion combinations are supported, so you may
 * need to take an intermediary conversion step to get what you want.
 * For example, complex images must first be converted into double or
 * float images before they can be converted into anything
 * else. (Although removing only the real or imaginary part can be
 * accomplished in one step, as described below).
 *
 * srcChannelMask specifies which channels to copy.  The lowest order
 * bit in the mask corresponds to the first bitplane.  For example, to
 * extract the green channel from an rgb image called imgRGB, you
 * could do the following:
 *
 *\code imgRGB->convertedImgType_(imgDst, DImage::DImage_u8, 1, 0x00000002) \endcode
 *
 * Likewise, to extract the blue channel, you could do the following:
 *
 *\code imgRGB->convertedImgType_(imgDst, DImage::DImage_u8, 1, 0x00000004) \endcode
 *
 * To convert an RGB image to 8-bit grayscale, do:
 *
 *\code imgRGB->convertedImgType_(imgDst, DImage::DImage_u8, 1, 0x00000007) \endcode
 *
 * To remove only the imaginary part of a complex image imgCPLX, the
 * following would work:
 *
 *\code
 * imgCPLX->convertedImgType_(imgDst, DImage::DImage_dbl_multi, 1, 0x00000002)
 * \endcode
 *
 * To create a 2-channel double image with both real and imaginary parts:
 *\code
 * imgCPLX->convertedImgType_(imgDst, DImage::DImage_dbl_multi, 2, 0x00000003)\endcode
 * 
 */
void DImage::convertedImgType_(DImage &imgDst, DImageType imgType,
			       int numDstChannels, D_uint32 srcChannelMask,
			       D_AllocationMethod allocMeth) const{
  unsigned int len;

  if(this == &imgDst){
    fprintf(stderr,
	    "DImage::convertedImgType_() destination image is also source\n");
    abort();
  }

//   printf("convertedImgType_() imgType=%d numDstChannels=%d "
// 	 "srcChannelMask=0x%08x allocMeth=%d\n",
// 	 imgType, numDstChannels, srcChannelMask, allocMeth);
//   printf("this: _imgType=%d _numChan=%d _allocMethod=%d\n",
// 	 _imgType, _numChan, _allocMethod);

  if(AllocationMethod_src == allocMeth)
    allocMeth = _allocMethod;
  if(-1 == numDstChannels){
    switch(imgType){
      case DImage_u8:
      case DImage_u16:
      case DImage_u32:
	numDstChannels = 1;
	break;
      case DImage_RGB:
      case DImage_RGB_16:
	numDstChannels = 3;
	break;
      case DImage_flt_multi:
      case DImage_dbl_multi:
	numDstChannels = _numChan;
	break;
      case DImage_cmplx:
	numDstChannels = 1;
	break;
    }//end switch(imgType)
  }// end if -1==numDstChannels

  imgDst.create(_w, _h, imgType, numDstChannels, allocMeth);
  
  switch(imgType){
    case DImage_u8: // destination
      if(numDstChannels != 1){
	fprintf(stderr, "DImage::convertedImgType_() "
		"numDstChannels(%d)!=1\n", numDstChannels);
	abort();
      }
      {
	D_uint8 *pDst = imgDst.dataPointer_u8();
	switch(_imgType){ // src _imgType
          case DImage_u16: // src
	    {
	      D_uint16 *pSrc;
	      pSrc = dataPointer_u16();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDst) = (D_uint8)(*pSrc);
		++pDst;
		++pSrc;
	      }
	    }
	    break;
	  case DImage_u32: // src
	    {
	      D_uint32 *pSrc;
	      pSrc = dataPointer_u32();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDst) = (D_uint8)(*pSrc);
		++pDst;
		++pSrc;
	      }
	    }
	    break;
	  case DImage_RGB: // src
	    srcChannelMask &= 0x00000007;
	    {
	      D_uint8 *pSrc;
	      pSrc = dataPointer_u8();
	      len = _w * _h;
	      switch(srcChannelMask){
	        case 0x00000007: // convert from RGB to grayscale
		  for(unsigned int i = 0,j=0; i < len; ++i, j+=3){
		    pDst[i] = (D_uint8)(0.299*pSrc[j] + 0.587*pSrc[j+1] +
					0.114*pSrc[j+2]);
		  }
		  break;
	        case 0x00000001: // extract RED channel
		  for(unsigned int i = 0,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        case 0x00000002: // extract GREEN channel
		  for(unsigned int i = 1,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        case 0x00000004: // extract BLUE channel
		  for(unsigned int i = 2,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        default:
		  fprintf(stderr,
			  "DImage::convertedImgType_() bad srcChannelMask\n");
		  abort();
	      } // end switch(srcChannelMask)
	    }
	    break;
	  case DImage_RGB_16: // src
	    srcChannelMask &= 0x00000007;
	    {
	      D_uint16 *pSrc;
	      pSrc = dataPointer_u16();
	      len = _w * _h;
	      switch(srcChannelMask){
	        case 0x00000007: // convert from RGB to grayscale
		  for(unsigned int i = 0,j=0; i < len; ++i, j+=3){
		    pDst[i] = (D_uint8)(0.299*pSrc[j] + 0.587*pSrc[j+1] +
					0.114*pSrc[j+2]);
		  }
		  break;
	        case 0x00000001: // extract RED channel
		  for(unsigned int i = 0,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        case 0x00000002: // extract GREEN channel
		  for(unsigned int i = 1,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        case 0x00000004: // extract BLUE channel
		  for(unsigned int i = 2,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        default:
		  fprintf(stderr,
			  "DImage::convertedImgType_() bad srcChannelMask\n");
		  abort();
	      } // end switch(srcChannelMask)
	    }
	    break;
	  case DImage_flt_multi: // src
	    if(0 == srcChannelMask){
	      fprintf(stderr,
		      "DImage::convertedImgType_() srcChannelMask==0\n");
	      abort();
	    }
	    {
	      float *pSrc;
	      int whichChannel = 0;
	      while(0 == (srcChannelMask & 0x00000001)){
		srcChannelMask >>= 1;
		++whichChannel;
	      }
	      if(whichChannel >= _numChan){
		fprintf(stderr,
			"DImage::convertedImgType_() bad srcChannelMask\n");
		abort();
	      }
	      pSrc = dataPointer_flt(whichChannel);
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDst) = (D_uint8)(*pSrc);
		++pDst;
		++pSrc;
	      }
	    }
	    break;
	  case DImage_dbl_multi: // src
	    if(0 == srcChannelMask){
	      fprintf(stderr,
		      "DImage::convertedImgType_() srcChannelMask==0\n");
	      abort();
	    }
	    {
	      double *pSrc;
	      int whichChannel = 0;
	      while(0 == (srcChannelMask & 0x00000001)){
		srcChannelMask >>= 1;
		++whichChannel;
	      }
	      if(whichChannel >= _numChan){
		fprintf(stderr, "DImage::convertedImgType_() bad channel\n");
		abort();
	      }
	      pSrc = dataPointer_dbl(whichChannel);
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDst) = (D_uint8)(*pSrc);
		++pDst;
		++pSrc;
	      }
	    }
	    break;
	  default:
	    fprintf(stderr, "DImage::convertedImgType_() conversion "
		    "unsupported from type %d to type %d\n", _imgType,imgType);
	    abort();
	} // end switch(_imgType) (src image type)
      }//end case dst _u8
      break;
    case DImage_u16: // destination
      if(numDstChannels != 1){
	fprintf(stderr, "DImage::convertedImgType_() "
		"numDstChannels(%d)!=1\n", numDstChannels);
	abort();
      }
      {
	D_uint16 *pDst = imgDst.dataPointer_u16();
	switch(_imgType){ // src _imgType
          case DImage_u8: // src
	    {
	      D_uint8 *pSrc;
	      pSrc = dataPointer_u8();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDst) = (D_uint16)(*pSrc);
		++pDst;
		++pSrc;
	      }
	    }
	    break;
	  case DImage_u32: // src
	    {
	      D_uint32 *pSrc;
	      pSrc = dataPointer_u32();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDst) = (D_uint16)(*pSrc);
		++pDst;
		++pSrc;
	      }
	    }
	    break;
	  case DImage_RGB: // src
	    srcChannelMask &= 0x00000007;
	    {
	      D_uint8 *pSrc;
	      pSrc = dataPointer_u8();
	      len = _w * _h;
	      switch(srcChannelMask){
	        case 0x00000007: // convert from RGB to grayscale
		  for(unsigned int i = 0,j=0; i < len; ++i, j+=3){
		    pDst[i] = (D_uint16)(0.299*pSrc[j] + 0.587*pSrc[j+1] +
					 0.114*pSrc[j+2]);
		  }
		  break;
	        case 0x00000001: // extract RED channel
		  for(unsigned int i = 0,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        case 0x00000002: // extract GREEN channel
		  for(unsigned int i = 1,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        case 0x00000004: // extract BLUE channel
		  for(unsigned int i = 2,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        default:
		  fprintf(stderr,
			  "DImage::convertedImgType_() bad srcChannelMask\n");
		  abort();
	      } // end switch(srcChannelMask)
	    }
	    break;
	  case DImage_RGB_16: // src
	    srcChannelMask &= 0x00000007;
	    {
	      D_uint16 *pSrc;
	      pSrc = dataPointer_u16();
	      len = _w * _h;
	      switch(srcChannelMask){
	        case 0x00000007: // convert from RGB to grayscale
		  for(unsigned int i = 0,j=0; i < len; ++i, j+=3){
		    pDst[i] = (D_uint16)(0.299*pSrc[j] + 0.587*pSrc[j+1] +
					 0.114*pSrc[j+2]);
		  }
		  break;
	        case 0x00000001: // extract RED channel
		  for(unsigned int i = 0,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        case 0x00000002: // extract GREEN channel
		  for(unsigned int i = 1,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        case 0x00000004: // extract BLUE channel
		  for(unsigned int i = 2,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        default:
		  fprintf(stderr,
			  "DImage::convertedImgType_() bad srcChannelMask\n");
		  abort();
	      } // end switch(srcChannelMask)
	    }
	    break;
	  case DImage_flt_multi: // src
	    if(0 == srcChannelMask){
	      fprintf(stderr,
		      "DImage::convertedImgType_() srcChannelMask==0\n");
	      abort();
	    }
	    {
	      float *pSrc;
	      int whichChannel = 0;
	      while(0 == (srcChannelMask & 0x00000001)){
		srcChannelMask >>= 1;
		++whichChannel;
	      }
	      if(whichChannel >= _numChan){
		fprintf(stderr,
			"DImage::convertedImgType_() bad srcChannelMask\n");
		abort();
	      }
	      pSrc = dataPointer_flt(whichChannel);
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDst) = (D_uint16)(*pSrc);
		++pDst;
		++pSrc;
	      }
	    }
	    break;
	  case DImage_dbl_multi: // src
	    if(0 == srcChannelMask){
	      fprintf(stderr,
		      "DImage::convertedImgType_() srcChannelMask==0\n");
	      abort();
	    }
	    {
	      double *pSrc;
	      int whichChannel = 0;
	      while(0 == (srcChannelMask & 0x00000001)){
		srcChannelMask >>= 1;
		++whichChannel;
	      }
	      if(whichChannel >= _numChan){
		fprintf(stderr, "DImage::convertedImgType_() bad channel\n");
		abort();
	      }
	      pSrc = dataPointer_dbl(whichChannel);
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDst) = (D_uint16)(*pSrc);
		++pDst;
		++pSrc;
	      }
	    }
	    break;
	  default:
	    fprintf(stderr, "DImage::convertedImgType_() conversion "
		    "unsupported from type %d to type %d\n", _imgType,imgType);
	    abort();
	} // end switch(_imgType) (src image type)
      }//end case dst _u16
      break;
    case DImage_u32: // destination
      {
	if(numDstChannels != 1){
	  fprintf(stderr, "DImage::convertedImgType_() "
		  "numDstChannels(%d)!=1\n", numDstChannels);
	  abort();
	}
	D_uint32 *pDst = imgDst.dataPointer_u32();
	switch(_imgType){ // src _imgType
          case DImage_u8: // src
	    {
	      D_uint8 *pSrc;
	      pSrc = dataPointer_u8();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDst) = (D_uint32)(*pSrc);
		++pDst;
		++pSrc;
	      }
	    }
	    break;
	  case DImage_u16: // src
	    {
	      D_uint16 *pSrc;
	      pSrc = dataPointer_u16();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDst) = (D_uint32)(*pSrc);
		++pDst;
		++pSrc;
	      }
	    }
	    break;
	  case DImage_RGB: // src
	    srcChannelMask &= 0x00000007;
	    {
	      D_uint8 *pSrc;
	      pSrc = dataPointer_u8();
	      len = _w * _h;
	      switch(srcChannelMask){
	        case 0x00000007: // convert from RGB to grayscale
		  for(unsigned int i = 0,j=0; i < len; ++i, j+=3){
		    pDst[i] = (D_uint32)(0.299*pSrc[j] + 0.587*pSrc[j+1] +
					 0.114*pSrc[j+2]);
		  }
		  break;
	        case 0x00000001: // extract RED channel
		  for(unsigned int i = 0,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        case 0x00000002: // extract GREEN channel
		  for(unsigned int i = 1,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        case 0x00000004: // extract BLUE channel
		  for(unsigned int i = 2,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        default:
		  fprintf(stderr,
			  "DImage::convertedImgType_() bad srcChannelMask\n");
		  abort();
	      } // end switch(srcChannelMask)
	    }
	    break;
	  case DImage_RGB_16: // src
	    srcChannelMask &= 0x00000007;
	    {
	      D_uint16 *pSrc;
	      pSrc = dataPointer_u16();
	      len = _w * _h;
	      switch(srcChannelMask){
	        case 0x00000007: // convert from RGB to grayscale
		  for(unsigned int i = 0,j=0; i < len; ++i, j+=3){
		    pDst[i] = (D_uint32)(0.299*pSrc[j] + 0.587*pSrc[j+1] +
					 0.114*pSrc[j+2]);
		  }
		  break;
	        case 0x00000001: // extract RED channel
		  for(unsigned int i = 0,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        case 0x00000002: // extract GREEN channel
		  for(unsigned int i = 1,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        case 0x00000004: // extract BLUE channel
		  for(unsigned int i = 2,j=0; i < len; ++i, j+=3){
		    pDst[i] = pSrc[j];
		  }
		  break;
	        default:
		  fprintf(stderr,
			  "DImage::convertedImgType_() bad srcChannelMask\n");
		  abort();
	      } // end switch(srcChannelMask)
	    }
	    break;
	  case DImage_flt_multi: // src
	    if(0 == srcChannelMask){
	      fprintf(stderr,
		      "DImage::convertedImgType_() srcChannelMask==0\n");
	      abort();
	    }
	    {
	      float *pSrc;
	      int whichChannel = 0;
	      while(0 == (srcChannelMask & 0x00000001)){
		srcChannelMask >>= 1;
		++whichChannel;
	      }
	      if(whichChannel >= _numChan){
		fprintf(stderr,
			"DImage::convertedImgType_() bad srcChannelMask\n");
		abort();
	      }
	      pSrc = dataPointer_flt(whichChannel);
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDst) = (D_uint32)(*pSrc);
		++pDst;
		++pSrc;
	      }
	    }
	    break;
	  case DImage_dbl_multi: // src
	    if(0 == srcChannelMask){
	      fprintf(stderr,
		      "DImage::convertedImgType_() srcChannelMask==0\n");
	      abort();
	    }
	    {
	      double *pSrc;
	      int whichChannel = 0;
	      while(0 == (srcChannelMask & 0x00000001)){
		srcChannelMask >>= 1;
		++whichChannel;
	      }
	      if(whichChannel >= _numChan){
		fprintf(stderr, "DImage::convertedImgType_() bad channel\n");
		abort();
	      }
	      pSrc = dataPointer_dbl(whichChannel);
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDst) = (D_uint32)(*pSrc);
		++pDst;
		++pSrc;
	      }
	    }
	    break;
	  default:
	    fprintf(stderr, "DImage::convertedImgType_() conversion "
		    "unsupported from type %d to type %d\n", _imgType,imgType);
	    abort();
	} // end switch(_imgType) (src image type)
      }//end case dst _u32
      break;
    case DImage_RGB: // destination
      {
	if(numDstChannels != 3){
	  fprintf(stderr, "DImage::convertedImgType_() "
		  "numDstChannels(%d)!=3\n", numDstChannels);
	  abort();
	}
	D_uint8 *pDst = imgDst.dataPointer_u8();
	switch(_imgType){ // src _imgType
          case DImage_u8: // src
	    {
	      D_uint8 *pSrc;
	      pSrc = dataPointer_u8();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[2] = pDst[1] = (*pDst) = (D_uint8)(*pSrc);
		pDst += 3;
		++pSrc;
	      }
	    }
	    break;
          case DImage_u16: // src
	    {
	      D_uint16 *pSrc;
	      pSrc = dataPointer_u16();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[2] = pDst[1] = (*pDst) = (D_uint8)(*pSrc);
		pDst += 3;
		++pSrc;
	      }
	    }
	    break;
	  case DImage_u32: // src
	    {
	      D_uint32 *pSrc;
	      pSrc = dataPointer_u32();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[2] = pDst[1] = (*pDst) = (D_uint8)(*pSrc);
		pDst += 3;
		++pSrc;
	      }
	    }
	    break;
	  case DImage_RGB_16: // src
	    srcChannelMask &= 0x00000007;
	    if(0x00000007 != srcChannelMask){
	      fprintf(stderr,
		      "DImage::convertedImgType_() bad srcChannelMask\n");
	      abort();
	    }
	    {
	      D_uint16 *pSrc;
	      pSrc = dataPointer_u16();
	      len = _w * _h * 3;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[i] = (D_uint8)(pSrc[i]);
	      }
	    }
	    break;
	  case DImage_flt_multi: // src
	    if(3 != _numChan){
	      fprintf(stderr,
		      "DImage::convertedImgType_() "
		      "source image channels(%d)!=3\n", _numChan);
	      abort();
	    }
	    srcChannelMask &= 0x00000007;
	    if(0x00000007 != srcChannelMask){
	      fprintf(stderr,
		      "DImage::convertedImgType_() bad srcChannelMask\n");
	      abort();
	    }
	    {
	      float *pSrcR, *pSrcG, *pSrcB;
	      pSrcR = dataPointer_flt(0);
	      pSrcG = dataPointer_flt(1);
	      pSrcB = dataPointer_flt(2);
	      len = _w * _h;
	      for(unsigned int i = 0, j = 0; i < len; ++i, j+=3){
		pDst[j] = (D_uint8)(pSrcR[i]);
		pDst[j+1] = (D_uint8)(pSrcG[i]);
		pDst[j+2] = (D_uint8)(pSrcB[i]);
	      }
	    }
	    break;
	  case DImage_dbl_multi: // src
	    if(3 != _numChan){
	      fprintf(stderr,
		      "DImage::convertedImgType_() "
		      "source image channels(%d)!=3\n", _numChan);
	      abort();
	    }
	    srcChannelMask &= 0x00000007;
	    if(0x00000007 != srcChannelMask){
	      fprintf(stderr,
		      "DImage::convertedImgType_() bad srcChannelMask\n");
	      abort();
	    }
	    {
	      double *pSrcR, *pSrcG, *pSrcB;
	      pSrcR = dataPointer_dbl(0);
	      pSrcG = dataPointer_dbl(1);
	      pSrcB = dataPointer_dbl(2);
	      len = _w * _h;
	      for(unsigned int i = 0, j = 0; i < len; ++i, j+=3){
		pDst[j] = (D_uint8)(pSrcR[i]);
		pDst[j+1] = (D_uint8)(pSrcG[i]);
		pDst[j+2] = (D_uint8)(pSrcB[i]);
	      }
	    }
	    break;
	  default:
	    fprintf(stderr, "DImage::convertedImgType_() conversion "
		    "unsupported from type %d to type %d\n", _imgType,imgType);
	    abort();
	} // end switch(_imgType) (src image type)
      }//end case dst _RGB
      break;
    case DImage_RGB_16: // destination
      {
	if(numDstChannels != 3){
	  fprintf(stderr, "DImage::convertedImgType_() "
		  "numDstChannels(%d)!=3\n", numDstChannels);
	  abort();
	}
	D_uint16 *pDst = imgDst.dataPointer_u16();
	switch(_imgType){ // src _imgType
          case DImage_u8: // src
	    {
	      D_uint8 *pSrc;
	      pSrc = dataPointer_u8();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[2] = pDst[1] = (*pDst) = (D_uint16)(*pSrc);
		pDst += 3;
		++pSrc;
	      }
	    }
	    break;
          case DImage_u16: // src
	    {
	      D_uint16 *pSrc;
	      pSrc = dataPointer_u16();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[2] = pDst[1] = (*pDst) = (D_uint16)(*pSrc);
		pDst += 3;
		++pSrc;
	      }
	    }
	    break;
	  case DImage_u32: // src
	    {
	      D_uint32 *pSrc;
	      pSrc = dataPointer_u32();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[2] = pDst[1] = (*pDst) = (D_uint16)(*pSrc);
		pDst += 3;
		++pSrc;
	      }
	    }
	    break;
	  case DImage_RGB: // src
	    srcChannelMask &= 0x00000007;
	    if(0x00000007 != srcChannelMask){
	      fprintf(stderr,
		      "DImage::convertedImgType_() bad srcChannelMask\n");
	      abort();
	    }
	    {
	      D_uint8 *pSrc;
	      pSrc = dataPointer_u8();
	      len = _w * _h * 3;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[i] = (D_uint16)(pSrc[i]);
	      }
	    }
	    break;
	  case DImage_flt_multi: // src
	    if(3 != _numChan){
	      fprintf(stderr,
		      "DImage::convertedImgType_() "
		      "source image channels(%d)!=3\n", _numChan);
	      abort();
	    }
	    srcChannelMask &= 0x00000007;
	    if(0x00000007 != srcChannelMask){
	      fprintf(stderr,
		      "DImage::convertedImgType_() bad srcChannelMask\n");
	      abort();
	    }
	    {
	      float *pSrcR, *pSrcG, *pSrcB;
	      pSrcR = dataPointer_flt(0);
	      pSrcG = dataPointer_flt(1);
	      pSrcB = dataPointer_flt(2);
	      len = _w * _h;
	      for(unsigned int i = 0, j = 0; i < len; ++i, j+=3){
		pDst[j] = (D_uint16)(pSrcR[i]);
		pDst[j+1] = (D_uint16)(pSrcG[i]);
		pDst[j+2] = (D_uint16)(pSrcB[i]);
	      }
	    }
	    break;
	  case DImage_dbl_multi: // src
	    if(3 != _numChan){
	      fprintf(stderr,
		      "DImage::convertedImgType_() "
		      "source image channels(%d)!=3\n", _numChan);
	      abort();
	    }
	    srcChannelMask &= 0x00000007;
	    if(0x00000007 != srcChannelMask){
	      fprintf(stderr,
		      "DImage::convertedImgType_() bad srcChannelMask\n");
	      abort();
	    }
	    {
	      double *pSrcR, *pSrcG, *pSrcB;
	      pSrcR = dataPointer_dbl(0);
	      pSrcG = dataPointer_dbl(1);
	      pSrcB = dataPointer_dbl(2);
	      len = _w * _h;
	      for(unsigned int i = 0, j = 0; i < len; ++i, j+=3){
		pDst[j] = (D_uint16)(pSrcR[i]);
		pDst[j+1] = (D_uint16)(pSrcG[i]);
		pDst[j+2] = (D_uint16)(pSrcB[i]);
	      }
	    }
	    break;
	  default:
	    fprintf(stderr, "DImage::convertedImgType_() conversion "
		    "unsupported from type %d to type %d\n", _imgType,imgType);
	    abort();
	} // end switch(_imgType) (src image type)
      }//end case dst _RGB_16
      break;
    case DImage_flt_multi: // destination
      {
	float *pDst = imgDst.dataPointer_flt(0);
	switch(_imgType){ // src _imgType
          case DImage_u8: // src
	    {
	      D_uint8 *pSrc;
	      pSrc = dataPointer_u8();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[i] = (float)(pSrc[i]);
	      }
	    }
	    break;
          case DImage_u16: // src
	    {
	      D_uint16 *pSrc;
	      pSrc = dataPointer_u16();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[i] = (float)(pSrc[i]);
	      }
	    }
	    break;
	  case DImage_u32: // src
	    {
	      D_uint32 *pSrc;
	      pSrc = dataPointer_u32();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[i] = (float)(pSrc[i]);
	      }
	    }
	    break;
	  case DImage_RGB: // src
	    {
	      D_uint8 *pSrc;
	      float *pDstG;
	      float *pDstB;
	      pSrc = dataPointer_u8();
	      len = _w * _h;
	      pDstG = imgDst.dataPointer_flt(1);
	      pDstB = imgDst.dataPointer_flt(2);
	      for(unsigned int i = 0, j=0; i < len; ++i, j+=3){
		pDst[i] = (float)(pSrc[j]);
		pDstG[i] = (float)(pSrc[j+1]);
		pDstB[i] = (float)(pSrc[j+2]);
	      }
	    }
	    break;
	  case DImage_RGB_16: // src
	    {
	      D_uint16 *pSrc;
	      float *pDstG;
	      float *pDstB;
	      pSrc = dataPointer_u16();
	      len = _w * _h;
	      pDstG = imgDst.dataPointer_flt(1);
	      pDstB = imgDst.dataPointer_flt(2);
	      for(unsigned int i = 0, j=0; i < len; ++i, j+=3){
		pDst[i] = (float)(pSrc[j]);
		pDstG[i] = (float)(pSrc[j+1]);
		pDstB[i] = (float)(pSrc[j+2]);
	      }
	    }
	    break;
	  case DImage_flt_multi: // src
	    {
	      float *pSrc;
	      int chanNum = 0;

	      len = _w * _h;
	      for(int chn=0; chn < _numChan; ++chn){
		if(srcChannelMask & (0x00000001 << chn)){
		  pSrc = dataPointer_flt(chn);
		  pDst = imgDst.dataPointer_flt(chanNum);
		  for(unsigned int i = 0; i < len; ++i){
		    pDst[i] = (float)(pSrc[i]);
		  }
		  ++chanNum;
		  if(chanNum > _numChan){
		    fprintf(stderr, "DImage::convertedImgType_() bad channel "
			    "mask (too many channels)\n");
		    abort();
		  }
		}// end if
	      }//end for(chan...
	      
	      if(chanNum != numDstChannels){
		fprintf(stderr, "DImage::convertedImgType_() bad channel "
			"mask (doesn't match numDstChannels)\n");
		abort();
	      }
	    }
	    break;
	  case DImage_dbl_multi: // src
	    {
	      double *pSrc;
	      int chanNum = 0;

	      len = _w * _h;
	      for(int chn=0; chn < _numChan; ++chn){
		if(srcChannelMask & (0x00000001 << chn)){
		  pSrc = dataPointer_dbl(chn);
		  pDst = imgDst.dataPointer_flt(chanNum);
		  for(unsigned int i = 0; i < len; ++i){
		    pDst[i] = (float)(pSrc[i]);
		  }
		  ++chanNum;
		  if(chanNum > _numChan){
		    fprintf(stderr, "DImage::convertedImgType_() bad channel "
			    "mask (too many channels)\n");
		    abort();
		  }
		}// end if
	      }//end for(chan...

	      if(chanNum != numDstChannels){
		fprintf(stderr, "DImage::convertedImgType_() bad channel "
			"mask (doesn't match numDstChannels)\n");
		abort();
	      }
	    }
	    break;
	  case DImage_cmplx: // src
	    {
	      std::complex<double> *pSrc;
	      int getReal;// set to 1 iff extracting the real part (0 if not)
	      int getImag;// set to 1 iff extracting the imaginary part
	      float *pDstReal;
	      float *pDstImag;
	      float fltTmp; //dummy var to store value in when not extracting

	      srcChannelMask &= 0x00000003;

	      getReal = getImag = 0;
	      pDstReal = pDstImag = &fltTmp;
	      if(srcChannelMask & 0x00000001){
		getReal = 1;
		pDstReal = imgDst.dataPointer_flt(0);
	      }
	      if(srcChannelMask & 0x00000002){
		getImag = 1;
		pDstImag = imgDst.dataPointer_flt(getReal);
	      }
	      if(((int)(getReal + getImag)) != numDstChannels){
		fprintf(stderr, "DImage::convertedImgType_() bad channel "
			"mask (doesn't match numDstChannels)\n");
		abort();
	      }

	      pSrc = dataPointer_cmplx();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDstReal) =
		  (float)(getReal * pSrc[i].real());//times by 0 to ignore
		(*pDstImag) =
		  (float)(getImag * pSrc[i].imag());//times by 1 to use
		pDstReal += getReal; // dont advance ptr if ignore (use fltTmp)
		pDstImag += getImag; // dont advance ptr if ignore (use fltTmp)
	      }
	    }
	    break;
	  default:
	    fprintf(stderr, "DImage::convertedImgType_() conversion "
		    "unsupported from type %d to type %d\n", _imgType,imgType);
	    abort();
	} // end switch(_imgType) (src image type)
      }//end case _flt_multi
      break;
    case DImage_dbl_multi: // destination
      {
	double *pDst = imgDst.dataPointer_dbl(0);
	switch(_imgType){ // src _imgType
          case DImage_u8: // src
	    {
	      D_uint8 *pSrc;
	      pSrc = dataPointer_u8();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[i] = (double)(pSrc[i]);
	      }
	    }
	    break;
          case DImage_u16: // src
	    {
	      D_uint16 *pSrc;
	      pSrc = dataPointer_u16();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[i] = (double)(pSrc[i]);
	      }
	    }
	    break;
	  case DImage_u32: // src
	    {
	      D_uint32 *pSrc;
	      pSrc = dataPointer_u32();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[i] = (double)(pSrc[i]);
	      }
	    }
	    break;
	  case DImage_RGB: // src
	    {
	      D_uint8 *pSrc;
	      double *pDstG;
	      double *pDstB;
	      pSrc = dataPointer_u8();
	      len = _w * _h;
	      pDstG = imgDst.dataPointer_dbl(1);
	      pDstB = imgDst.dataPointer_dbl(2);
	      for(unsigned int i = 0, j=0; i < len; ++i, j+=3){
		pDst[i] = (double)(pSrc[j]);
		pDstG[i] = (double)(pSrc[j+1]);
		pDstB[i] = (double)(pSrc[j+2]);
	      }
	    }
	    break;
	  case DImage_RGB_16: // src
	    {
	      D_uint16 *pSrc;
	      double *pDstG;
	      double *pDstB;
	      pSrc = dataPointer_u16();
	      len = _w * _h;
	      pDstG = imgDst.dataPointer_dbl(1);
	      pDstB = imgDst.dataPointer_dbl(2);
	      for(unsigned int i = 0, j=0; i < len; ++i, j+=3){
		pDst[i] = (double)(pSrc[j]);
		pDstG[i] = (double)(pSrc[j+1]);
		pDstB[i] = (double)(pSrc[j+2]);
	      }
	    }
	    break;
	  case DImage_flt_multi: // src
	    {
	      float *pSrc;
	      int chanNum = 0;

	      len = _w * _h;
	      for(int chn=0; chn < _numChan; ++chn){
		if(srcChannelMask & (0x00000001 << chn)){
		  pSrc = dataPointer_flt(chn);
		  pDst = imgDst.dataPointer_dbl(chanNum);
		  for(unsigned int i = 0; i < len; ++i){
		    pDst[i] = (double)(pSrc[i]);
		  }
		  ++chanNum;
		  if(chanNum > _numChan){
		    fprintf(stderr, "DImage::convertedImgType_() bad channel "
			    "mask (too many channels)\n");
		    abort();
		  }
		}// end if
	      }//end for(chan...
	      
	      if(chanNum != numDstChannels){
		fprintf(stderr, "DImage::convertedImgType_() bad channel "
			"mask (doesn't match numDstChannels)\n");
		abort();
	      }
	    }
	    break;
	  case DImage_dbl_multi: // src
	    {
	      double *pSrc;
	      int chanNum = 0;

	      len = _w * _h;
	      for(int chn=0; chn < _numChan; ++chn){
		if(srcChannelMask & (0x00000001 << chn)){
		  pSrc = dataPointer_dbl(chn);
		  pDst = imgDst.dataPointer_dbl(chanNum);
		  for(unsigned int i = 0; i < len; ++i){
		    pDst[i] = (double)(pSrc[i]);
		  }
		  ++chanNum;
		  if(chanNum > _numChan){
		    fprintf(stderr, "DImage::convertedImgType_() bad channel "
			    "mask (too many channels)\n");
		    abort();
		  }
		}// end if
	      }//end for(chan...

	      if(chanNum != numDstChannels){
		fprintf(stderr, "DImage::convertedImgType_() bad channel "
			"mask (doesn't match numDstChannels)\n");
		abort();
	      }
	    }
	    break;
	  case DImage_cmplx: // src
	    {
	      std::complex<double> *pSrc;
	      int getReal;// set to 1 iff extracting the real part (0 if not)
	      int getImag;// set to 1 iff extracting the imaginary part
	      double *pDstReal;
	      double *pDstImag;
	      double dblTmp; //dummy var to store value in when not extracting

	      srcChannelMask &= 0x00000003;

	      getReal = getImag = 0;
	      pDstReal = pDstImag = &dblTmp;
	      if(srcChannelMask & 0x00000001){
		getReal = 1;
		pDstReal = imgDst.dataPointer_dbl(0);
	      }
	      if(srcChannelMask & 0x00000002){
		getImag = 1;
		// point to channel 0 if real was not also used, else channel 1
		pDstImag = imgDst.dataPointer_dbl(getReal);
	      }
	      if(((int)(getReal + getImag)) != numDstChannels){
		fprintf(stderr, "DImage::convertedImgType_() bad channel "
			"mask (doesn't match numDstChannels)\n");
		abort();
	      }

	      pSrc = dataPointer_cmplx();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDstReal) =
		  (double)(getReal * pSrc[i].real());//times by 0 to ignore
		(*pDstImag) =
		  (double)(getImag * pSrc[i].imag());//times by 1 to use
		pDstReal += getReal; // dont advance ptr if ignore (use dblTmp)
		pDstImag += getImag; // dont advance ptr if ignore (use dblTmp)
	      }
	    }
	    break;
	  default:
	    fprintf(stderr, "DImage::convertedImgType_() conversion "
		    "unsupported from type %d to type %d\n", _imgType,imgType);
	    abort();
	} // end switch(_imgType) (src image type)
      }//end case _dbl_multi
      break;
    case DImage_cmplx: // destination
      {
	std::complex<double> *pDst = imgDst.dataPointer_cmplx();
	switch(_imgType){ // src _imgType
          case DImage_u8: // src
	    {
	      D_uint8 *pSrc;
	      pSrc = dataPointer_u8();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[i] = std::complex<double>(pSrc[i],0.);
	      }
	    }
	    break;
          case DImage_u16: // src
	    {
	      D_uint16 *pSrc;
	      pSrc = dataPointer_u16();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[i] = std::complex<double>(pSrc[i],0.);
	      }
	    }
	    break;
	  case DImage_u32: // src
	    {
	      D_uint32 *pSrc;
	      pSrc = dataPointer_u32();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		pDst[i] = std::complex<double>(pSrc[i],0.);
	      }
	    }
	    break;
	  case DImage_RGB: // src
	  case DImage_RGB_16: // src
	    {
	      fprintf(stderr, "DImage::convertedImgType_() doesn't suport "
		      "conversion directly from RGB to complex\n");
	      abort();
	    }
	    break;
	  case DImage_flt_multi: // src
	    {
	      float *pSrcReal;
	      float *pSrcImag;
	      
	      if(numDstChannels != 1){
		fprintf(stderr, "DImage::convertedImgType_() "
			"numDstChannels(%d) should be 1 to convert to "
			"complex image\n", numDstChannels);
		abort();
	      }
	      if((_numChan < 1) || (_numChan > 2)){
		fprintf(stderr, "DImage::convertedImgType_() src image "
			"must have 1 or 2 channels to convert to complex\n");
		abort();
	      }
	      if( (srcChannelMask & 0x00000003) == 0x00000003){
		// copy two channels into real and imaginary parts
		pSrcReal = dataPointer_flt(0);
		pSrcImag = dataPointer_flt(1);
		len = _w * _h;
		for(unsigned int i = 0; i < len; ++i){
		  pDst[i] = std::complex<double>(pSrcReal[i],pSrcImag[i]);
		}
	      }
	      else if( (srcChannelMask & 0x00000003) == 0x00000001){
		// copy one channel into real parts
		pSrcReal = dataPointer_flt(0);
		len = _w * _h;
		for(unsigned int i = 0; i < len; ++i){
		  pDst[i] = std::complex<double>(pSrcReal[i],0.);
		}
	      }
	      else if( (srcChannelMask & 0x00000003) == 0x00000002){
		// copy two channels into real and imaginary parts
		pSrcImag = dataPointer_flt(0);
		len = _w * _h;
		for(unsigned int i = 0; i < len; ++i){
		  pDst[i] = std::complex<double>(0.,pSrcImag[i]);
		}
	      }
	      else{
		fprintf(stderr, "DImage::convertedImgType_() invalid "
			"srcChannelMask for complex conversion\n");
		abort();
	      }
	    }
	    break;
	  case DImage_dbl_multi: // src
	    {
	      double *pSrcReal;
	      double *pSrcImag;
	      
	      if(numDstChannels != 1){
		fprintf(stderr, "DImage::convertedImgType_() "
			"numDstChannels(%d) should be 1 to convert to "
			"complex image\n", numDstChannels);
		abort();
	      }
	      if((_numChan < 1) || (_numChan > 2)){
		fprintf(stderr, "DImage::convertedImgType_() src image "
			"must have 1 or 2 channels to convert to complex\n");
		abort();
	      }
	      if( (srcChannelMask & 0x00000003) == 0x00000003){
		// copy two channels into real and imaginary parts
		pSrcReal = dataPointer_dbl(0);
		pSrcImag = dataPointer_dbl(1);
		len = _w * _h;
		for(unsigned int i = 0; i < len; ++i){
		  pDst[i] = std::complex<double>(pSrcReal[i],pSrcImag[i]);
		}
	      }
	      else if( (srcChannelMask & 0x00000003) == 0x00000001){
		// copy one channel into real parts
		pSrcReal = dataPointer_dbl(0);
		len = _w * _h;
		for(unsigned int i = 0; i < len; ++i){
		  pDst[i] = std::complex<double>(pSrcReal[i],0.);
		}
	      }
	      else if( (srcChannelMask & 0x00000003) == 0x00000002){
		// copy two channels into real and imaginary parts
		pSrcImag = dataPointer_dbl(0);
		len = _w * _h;
		for(unsigned int i = 0; i < len; ++i){
		  pDst[i] = std::complex<double>(0.,pSrcImag[i]);
		}
	      }
	      else{
		fprintf(stderr, "DImage::convertedImgType_() invalid "
			"srcChannelMask for complex conversion\n");
		abort();
	      }
	    }
	    break;
	  case DImage_cmplx: // src
	    {
	      std::complex<double> *pSrc;
	      int getReal;// set to 1 iff extracting the real part (0 if not)
	      int getImag;// set to 1 iff extracting the imaginary part
	      double *pDstReal;
	      double *pDstImag;
	      double dblTmp; //dummy var to store value in when not extracting

	      srcChannelMask &= 0x00000003;

	      getReal = getImag = 0;
	      pDstReal = pDstImag = &dblTmp;
	      if(srcChannelMask & 0x00000001){
		getReal = 1;
		pDstReal = imgDst.dataPointer_dbl(0);
	      }
	      if(srcChannelMask & 0x00000002){
		getImag = 1;
		// point to channel 0 if real was not also used, else channel 1
		pDstImag = imgDst.dataPointer_dbl(getReal);
	      }
	      if(((int)(getReal + getImag)) != numDstChannels){
		fprintf(stderr, "DImage::convertedImgType_() bad channel "
			"mask (doesn't match numDstChannels)\n");
		abort();
	      }

	      pSrc = dataPointer_cmplx();
	      len = _w * _h;
	      for(unsigned int i = 0; i < len; ++i){
		(*pDstReal) =
		  (double)(getReal * pSrc[i].real());//times by 0 to ignore
		(*pDstImag) =
		  (double)(getImag * pSrc[i].imag());//times by 1 to use
		pDstReal += getReal; // dont advance ptr if ignore (use dblTmp)
		pDstImag += getImag; // dont advance ptr if ignore (use dblTmp)
	      }
	    }
	    break;
	  default:
	    fprintf(stderr, "DImage::convertedImgType_() conversion "
		    "unsupported from type %d to type %d\n", _imgType,imgType);
	    abort();
	} // end switch(_imgType) (src image type)
      }//end case _dbl_multi
      break;
    default:
      fprintf(stderr, "DImage::convertedImgType_() conversion "
	      "unsupported from type %d to type %d\n", _imgType,imgType);
      abort();
  }//end switch(dst imgType)
}

void DImage::setProperty(const std::string propName, std::string propVal){
  mapProps[propName] = propVal;
}
std::string DImage::getPropertyVal(std::string propName){
  std::string s;
  if(mapProps.end() != mapProps.find(propName)){
    s = mapProps[propName];
  }
  return s;
}
std::string DImage::getPropertyValByIndex(unsigned int propNum){
  std::string s;
  std::map<const std::string, std::string>::iterator iter = mapProps.begin();
#ifdef DEBUG
  if(propNum >= mapProps.size()){
    fprintf(stderr, "DImage::getPropertyValByIndex() propNum out of range\n");
    abort();
    return s;
  }
#endif
  for(unsigned int i = 0; i < propNum; ++i, ++iter);
  return (*iter).second;
}
int DImage::getNumProperties() const{
  return mapProps.size();
}
void DImage::clearProperties(){
  mapProps.clear();
}

void DImage::addComment(std::string stComment){
  vectComments.push_back(stComment);
}
std::string DImage::getCommentByIndex(unsigned int idx) const{
  std::string s;
#ifdef DEBUG
  if(idx >= vectComments.size()){
    fprintf(stderr, "DImage::getCommentByIndex() index out of range\n");
    abort();
    return s;
  }
#endif
  return vectComments[idx];
}
int DImage::getNumComments() const{
  return vectComments.size();
}
void DImage::clearComments(){
  vectComments.clear();
}

///Set the alignment boundary for newly allocated buffers
/**This alignment is only used when the allocation method is set to
 * AllocationMethod_daligned. */
void DImage::setAlignment(int alignment){
  int numBits = 0;
  int al = alignment;

  if(alignment < 4){
    fprintf(stderr, "DImage::setAlignment() alignment must be >=4 (was %d)\n",
	    alignment);
    abort();
    return;
  }
  while((al != 0)&&(numBits<=1)){
    numBits += (al & 0x00000001);
    al >>= 1;
  }
  if(numBits != 1){
    fprintf(stderr, "DImage::setAlignment() alignment must be a power of 2"
	    " (was %d)\n", alignment);
    abort();
    return;
  }
  if(0 != (alignment % (sizeof(void*)))){
    fprintf(stderr, "DImage::setAlignment() alignment must be a multple of "
	    " sizeof(void*), which is %d. (was %d)\n",
	    (int)sizeof(void*), alignment);
    abort();
    return;
  }
  DImage::_data_alignment = alignment;
}
int DImage::getAlignment(){
  return DImage::_data_alignment;
}


///Attempts to recognize (and return) the file format from the first few bytes
/**This function seeks to the beginning of the file upon entry, but
 * does not reset its position to the beginning before exiting.  It
 * also recognizes several types of images that may not actually be
 * supported yet
 */
DImage::DFileFormat DImage::getImageFileFormat(FILE *fin){
  DImage::DFileFormat retVal = DFileFormat_unknown;
  char rgHead[8];
  size_t bytesRead;

  rewind(fin); // rewind to the beginning of the file

  bytesRead = fread(rgHead, 1, 8,fin);
  if((bytesRead != 8) || (feof(fin)) || (ferror(fin))){
    fprintf(stderr,
	    "DImage::getImageFileFormat() error reading first 8 bytes\n");
    return DFileFormat_unknown;
  }

  switch((unsigned char)(rgHead[0])){
    case 'P': // may be a PNM
      if(isspace(rgHead[2])){// whitespace after magic number
	switch(rgHead[1]){
	  case '1':retVal = DFileFormat_pbm_plain; break;
	  case '2':retVal = DFileFormat_pgm_plain; break;
	  case '3':retVal = DFileFormat_ppm_plain; break;
	  case '4':retVal = DFileFormat_pbm; break;
	  case '5':retVal = DFileFormat_pgm; break;
	  case '6':retVal = DFileFormat_ppm; break;
	  default: break; // do nothing (leave retVal as DFileFormat_unknown)
	}
      }
      break;
    case 0xff: // may be a JPEG (ffd8ff hex)
      if(((char)0xd8==rgHead[1]) && ((char)0xff==rgHead[2]))
	retVal = DFileFormat_jpg;
      break;
    case 0x89: // may be a PNG (89504e470d0a1a0a hex)
      if((0x50==rgHead[1]) && (0x4e==rgHead[2]) && (0x47==rgHead[3]) &&
	 (0x0d==rgHead[4]) && (0x0a==rgHead[5]) && (0x1a==rgHead[6]) &&
	 (0x0a==rgHead[7]))
	retVal = DFileFormat_png;
      break;
    case 0x4d: // may be a big-endian TIFF 'M' ("MM" followed by (short)42)
      if((0x4d==rgHead[1]) && (0x00==rgHead[2]) && (0x2a==rgHead[3]))
	 retVal = DFileFormat_tiff;
      break;
    case 0x49: // may be a little-endian TIFF 'I' ("II" followed by (short)42)
      if((0x49==rgHead[1]) && (0x2a==rgHead[2]) && (0x00==rgHead[3]))
	 retVal = DFileFormat_tiff;
      break;
    case 'G': // may be a GIF ("GIF87a" or "GIF89a")
      if(('I'==rgHead[1]) && ('F'==rgHead[2]) && ('8'==rgHead[3]) &&
	 (('7'==rgHead[4])||('9'==rgHead[4]))
	 && ('a'==rgHead[5]))
	retVal = DFileFormat_gif;
      break;
    case 'B': // may be a BMP ("BM")
      if('M'==rgHead[1])
	retVal = DFileFormat_bmp;
      break;
    default: // can't tell type from file header - do nothing
      break;
  }
  return retVal;
}

///Attempts to recognize (and return) the file format from the first few bytes
/**This opens the file referenced in stPath.  It then calls the
   version of getImageFileFormat that takes a FILE* as a parameter and
   returns its result.
 */
DImage::DFileFormat DImage::getImageFileFormat(const char *stPath){
  DImage::DFileFormat retVal;
  FILE *fin;

  fin = fopen(stPath, "rb");
  if(!fin){
    fprintf(stderr, "DImage::getImageFileFormat() couldn't open '%s'\n",
	    stPath);
    return DFileFormat_unknown;
  }
  retVal = getImageFileFormat(fin);
  fclose(fin);

  return retVal;
}

///Split the RGB channels out into 3 distinct single-channel images
/**This only makes sense for images of type DImage_RGB or
   DImage_RGB_16.  For each channel image, if its width and height
   match this image and the number of channels and type are correct,
   then create() will not be called and the data will simply overwrite
   the buffer already allocated.*/
void DImage::splitRGB(DImage &imgR, DImage &imgG, DImage &imgB) const{
  if((_imgType != DImage_RGB) && (_imgType != DImage_RGB_16)){
    fprintf(stderr, "DImage::splitRGB() type must be RGB or RGB_16\n");
    abort();
    return;
  }

  switch(_imgType){
    case DImage_RGB:
      {
	D_uint8 *pR, *pG, *pB, *pSrc;
	unsigned int len;

	if((_w != imgR._w) || (_h != imgR._h) ||
	   (DImage_u8 != imgR._imgType) || (1 != imgR._numChan))
	  imgR.create(_w, _h, DImage_u8, 1, _allocMethod);
	if((_w != imgG._w) || (_h != imgG._h) ||
	   (DImage_u8 != imgG._imgType) || (1 != imgG._numChan))
	  imgG.create(_w, _h, DImage_u8, 1, _allocMethod);
	if((_w != imgB._w) || (_h != imgB._h) ||
	   (DImage_u8 != imgB._imgType) || (1 != imgB._numChan))
	  imgB.create(_w, _h, DImage_u8, 1, _allocMethod);

	pSrc = pData;
	pR = imgR.pData;
	pG = imgG.pData;
	pB = imgB.pData;

	len = _w * _h;
	for(unsigned int i = 0; i < len; ++i){
	  (*pR) = (*pSrc);
	  (*pG) = pSrc[1];
	  (*pB) = pSrc[2];
	  pSrc += 3;
	  ++pR;
	  ++pG;
	  ++pB;
	}
      }
      break;
    case DImage_RGB_16:
      {
	D_uint16 *pR, *pG, *pB, *pSrc;
	unsigned int len;
	if((_w != imgR._w) || (_h != imgR._h) ||
	   (DImage_u16 != imgR._imgType) || (1 != imgR._numChan))
	  imgR.create(_w, _h, DImage_u16, 1, _allocMethod);
	if((_w != imgG._w) || (_h != imgG._h) ||
	   (DImage_u16 != imgG._imgType) || (1 != imgG._numChan))
	  imgG.create(_w, _h, DImage_u16, 1, _allocMethod);
	if((_w != imgB._w) || (_h != imgB._h) ||
	   (DImage_u16 != imgB._imgType) || (1 != imgB._numChan))
	  imgB.create(_w, _h, DImage_u16, 1, _allocMethod);
	pSrc = (D_uint16*)pData;
	pR = (D_uint16*)imgR.pData;
	pG = (D_uint16*)imgG.pData;
	pB = (D_uint16*)imgB.pData;

	len = _w * _h;
	for(unsigned int i = 0; i < len; ++i){
	  (*pR) = (*pSrc);
	  (*pG) = pSrc[1];
	  (*pB) = pSrc[2];
	  pSrc += 3;
	  ++pR;
	  ++pG;
	  ++pB;
	}
      }
      break;
    default:
      return;
  }
}

///Combine three distinct single-channel images into one RGB image
/**The three images must be DImage_u8 or DImage_u16 and all of the
 *same image type, width, height, and one channel.  This image (the
 *destination image) buffer will not be reallocated if it is the same
 *width and height, and appropriate image type for imgR, imgG, and
 *imgB.  If any of these is different, create() will be called on this
 *to reallocate the buffer.
 */
void DImage::combineRGB(DImage &imgR, DImage &imgG, DImage &imgB){

  // check source images to make sure they all match
  if( (imgR._w != imgG._w) || (imgG._w != imgB._w) ||
      (imgR._h != imgG._h) || (imgG._h != imgB._h) ||
      (imgR._numChan != 1) || (imgG._numChan != 1) || (imgB._numChan != 1) ||
      (imgR._imgType != imgG._imgType) || (imgG._imgType != imgB._imgType)){
    fprintf(stderr, "DImage::combineRGB() imgR, imgG, and imgB must all "
	    "be the same size and type of RGB image\n");
    abort();
    return;
  }
  switch(imgR._imgType){
    case DImage_u8:
      {
	D_uint8 *pR, *pG, *pB, *pDst;
	unsigned int len;
	if((_w != imgR._w) || (_h != imgR._h) || (_imgType != DImage_RGB) ||
	   (_numChan != 3))
	  create(imgR._w, imgR._h, DImage_RGB, 3, _allocMethod);
	pDst = pData;
	pR = imgR.pData;
	pG = imgG.pData;
	pB = imgB.pData;
	len = _w * _h;
	for(unsigned int i = 0; i < len; ++i){
	  (*pDst) = (*pR);
	  pDst[1] = (*pG);
	  pDst[2] = (*pB);
	  pDst += 3;
	  ++pR;
	  ++pG;
	  ++pB;
	}
      }
      break;
    case DImage_u16:
      {
	D_uint16 *pR, *pG, *pB, *pDst;
	unsigned int len;
	if((_w != imgR._w) || (_h != imgR._h) || (_imgType != DImage_RGB_16) ||
	   (_numChan != 3))
	  create(imgR._w, imgR._h, DImage_RGB_16, 3, _allocMethod);
	pDst = (D_uint16*)pData;
	pR = (D_uint16*)imgR.pData;
	pG = (D_uint16*)imgG.pData;
	pB = (D_uint16*)imgB.pData;
	len = _w * _h;
	for(unsigned int i = 0; i < len; ++i){
	  (*pDst) = (*pR);
	  pDst[1] = (*pG);
	  pDst[2] = (*pB);
	  pDst += 3;
	  ++pR;
	  ++pG;
	  ++pB;
	}
      }
      break;
    default:
      fprintf(stderr, "DImage::combineRGB() only works for images of type "
	      "DImage_u8 or DImage_u16\n");
      abort();
  }
}


///Scale up the image by numPows powers of two
/**Scales this image by up by numPows powers of two, using pixel
 * replication.  Puts the result in imgDst.  dstW and dstH can be
 * specified to crop/pad the size of the destination image.  This
 * dstW/dstH feature is specifically included for convenience when you
 * need to first shrink the image by a certain number of powers of
 * two, then restore it to original size, without slightly modifying
 * the original size or pixel alignment
 */
void DImage::scaledUpPow2_(DImage &imgDst, int numPows, int dstW, int dstH,
			   DImageTransformMode mode) const{
  int numPowsDivisor;

  if(numPows < 0){
    fprintf(stderr, " DImage::scaledUpPow2_() numPows==%d (<0)\n", numPows);
    abort();
  }
  numPowsDivisor = 1 << numPows;

  if(-1 == dstW){
    dstW = _w << numPows;
  }
  if(-1 == dstH){
    dstH = _h << numPows;
  }
  if(DImageTransSample == mode){
    switch(this->_imgType){
      case DImage_u8:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  D_uint8 *pDstRow;
	  D_uint8 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_u8, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();

	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y];
	    pSrcRow = &pSrc[_w*(y/numPowsDivisor)];
	    for(int x = 0; x < dstW; ++x){
	      pDstRow[x] = pSrcRow[x/numPowsDivisor];
	    }
	  }
	}// end case
	break;
      case DImage_u16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  D_uint16 *pDstRow;
	  D_uint16 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_u16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y];
	    pSrcRow = &pSrc[_w*(y/numPowsDivisor)];
	    for(int x = 0; x < dstW; ++x){
	      pDstRow[x] = pSrcRow[x/numPowsDivisor];
	    }
	  }
	}// end case
	break;
      case DImage_u32:
	{
	  D_uint32 *pSrc;
	  D_uint32 *pDst;
	  D_uint32 *pDstRow;
	  D_uint32 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_u32, _numChan, _allocMethod);
	  pSrc = (D_uint32*)pData;
	  pDst = imgDst.dataPointer_u32();
	  
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y];
	    pSrcRow = &pSrc[_w*(y/numPowsDivisor)];
	    for(int x = 0; x < dstW; ++x){
	      pDstRow[x] = pSrcRow[x/numPowsDivisor];
	    }
	  }
	}// end case
	break;
      case DImage_RGB:
	{
	  D_uint8 *pSrc;
	  D_uint8 *pDst;
	  D_uint8 *pDstRow;
	  D_uint8 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_RGB, _numChan, _allocMethod);
	  pSrc = pData;
	  pDst = imgDst.dataPointer_u8();

	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y*3];
	    pSrcRow = &pSrc[3*_w*(y/numPowsDivisor)];
	    for(int x = 0; x < dstW; ++x){
	      int xTmp;
	      xTmp = 3*(x / numPowsDivisor);
	      pDstRow[3*x] = pSrcRow[xTmp];
	      pDstRow[3*x+1] = pSrcRow[1+xTmp];
	      pDstRow[3*x+2] = pSrcRow[2+xTmp];
	    }
	  }
	}// end case
	break;
      case DImage_RGB_16:
	{
	  D_uint16 *pSrc;
	  D_uint16 *pDst;
	  D_uint16 *pDstRow;
	  D_uint16 *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_RGB_16, _numChan, _allocMethod);
	  pSrc = (D_uint16*)pData;
	  pDst = imgDst.dataPointer_u16();
	  
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y*3];
	    pSrcRow = &pSrc[3*_w*(y/numPowsDivisor)];
	    for(int x = 0; x < dstW*3; x+=3){
	      int xTmp;
	      xTmp = x / numPowsDivisor;
	      pDstRow[x] = pSrcRow[xTmp];
	      pDstRow[x+1] = pSrcRow[1+xTmp];
	      pDstRow[x+2] = pSrcRow[2+xTmp];
	    }
	  }
	}// end case
	break;
      case DImage_flt_multi:
	{
	  float *pSrc;
	  float *pDst;
	  float *pDstRow;
	  float *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_flt_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_flt(chan);
	    pDst = imgDst.dataPointer_flt(chan);
	    for(int y = 0; y < dstH; ++y){
	      pDstRow = &pDst[dstW*y];
	      pSrcRow = &pSrc[_w*(y/numPowsDivisor)];
	      for(int x = 0; x < dstW; ++x){
		pDstRow[x] = pSrcRow[x/numPowsDivisor];
	      }
	    }
	  }
	}// end case
	break;
      case DImage_dbl_multi:
	{
	  double *pSrc;
	  double *pDst;
	  double *pDstRow;
	  double *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_dbl_multi, _numChan, _allocMethod);

	  for(int chan = 0; chan < _numChan; ++chan){
	    pSrc = this->dataPointer_dbl(chan);
	    pDst = imgDst.dataPointer_dbl(chan);
	    for(int y = 0; y < dstH; ++y){
	      pDstRow = &pDst[dstW*y];
	      pSrcRow = &pSrc[_w*(y/numPowsDivisor)];
	      for(int x = 0; x < dstW; ++x){
		pDstRow[x] = pSrcRow[x/numPowsDivisor];
	      }
	    }
	  }
	}// end case
	break;
      case DImage_cmplx:
	{
	  std::complex<double> *pSrc;
	  std::complex<double> *pDst;
	  std::complex<double> *pDstRow;
	  std::complex<double> *pSrcRow;
	  
	  imgDst.create(dstW, dstH, DImage_cmplx, _numChan, _allocMethod);

	  pSrc = this->dataPointer_cmplx();
	  pDst = imgDst.dataPointer_cmplx();
	  for(int y = 0; y < dstH; ++y){
	    pDstRow = &pDst[dstW*y];
	    pSrcRow = &pSrc[_w*(y/numPowsDivisor)];
	    for(int x = 0; x < dstW; ++x){
	      pDstRow[x] = pSrcRow[x/numPowsDivisor];
	    }
	  }
	}// end case
	break;
    }// end switch
  } // end if sample transform mode
  else if(DImageTransSmooth){
    // TODO: implement smooth scaling for when 
    fprintf(stderr, "scaledUpPow2_() Not Yet Implemented for "
	    "DImageTransSmooth\n");
    abort();
  }
}




///Fill the image with value val
/**Fills each pixel in the image with val (cast to the correct type
 * for image) The buffer should already be allocated.  If the size has
 * been changed without changing the actual allocated buffer size,
 * only width() by height() pixels will be filled.
 */
void DImage::fill(double val){
  int len;
  if((_w < 1) || (_h < 1)){
    fprintf(stderr, "DImage::fill() called for image with bad size (%dx%d)\n",
	    _w, _h);
    abort();
  }
  if(NULL == pData){
    fprintf(stderr, "DImage::fill() called for image with pData==NULL\n");
    abort();
  }
  switch(_imgType){
    case DImage_u8:
      {
	D_uint8 u8Val;
	u8Val = (D_uint8)val;
	memset(pData, (int)val, _w * _h * sizeof(D_uint8));
      }
      break;
    case DImage_u16:
      {
	D_uint16 u16Val;
	D_uint16 *pData16;
	u16Val = (D_uint16)val;
	pData16 = (D_uint16*)pData;
	len = _w * _h;
	for(int i = 0; i < len; ++i)
	  pData16[i] = u16Val;
      }
      break;
    case DImage_u32:
      {
	D_uint32 u32Val;
	D_uint32 *pData32;
	u32Val = (D_uint32)val;
	pData32 = (D_uint32*)pData;
	len = _w * _h;
	for(int i = 0; i < len; ++i)
	  pData32[i] = u32Val;
      }
      break;
    case DImage_RGB:
    case DImage_RGB_16:
      fprintf(stderr, "DImage::fill() grayscale called for RGB image\n");
      break;
    case DImage_flt_multi:
      {
	float fltVal;
	float *pDataFlt;
	fltVal = (float)val;
	len = _w * _h;
	for(int chan = 0; chan < _numChan; chan++){
	  pDataFlt = dataPointer_flt(chan);
	  for(int i = 0; i < len; ++i)
	    pDataFlt[i] = fltVal;
	}//end for chan
      }
      break;
    case DImage_dbl_multi:
      {
	double dblVal;
	double *pDataDbl;
	dblVal = (double)val;
	len = _w * _h;
	for(int chan = 0; chan < _numChan; chan++){
	  pDataDbl = dataPointer_dbl(chan);
	  for(int i = 0; i < len; ++i)
	    pDataDbl[i] = dblVal;
	}//end for chan
      }
      break;
    case DImage_cmplx:
      fprintf(stderr, "DImage::fill() grayscale called for complex image\n");
      break;
  } // end switch
}// end ::fill


///Fill the image with RGB value R,G,B
/**Fills each pixel in the image with R,G,B The buffer should already
 * be allocated.  If the size has been changed without changing the
 * actual allocated buffer size, only width() by height() pixels will
 * be filled.
 */
void DImage::fill(int R, int G, int B){
  int len;
  if((_w < 1) || (_h < 1)){
    fprintf(stderr, "DImage::fill() called for image with bad size (%dx%d)\n",
	    _w, _h);
    abort();
  }
  if(NULL == pData){
    fprintf(stderr, "DImage::fill() called for image with pData==NULL\n");
    abort();
  }

  switch(_imgType){
    case DImage_u8:
    case DImage_u16:
    case DImage_u32:
      fprintf(stderr, "DImage::fill() RGB called for grayscale image\n");
      break;
    case DImage_RGB:
      D_uint8 u8r, u8g, u8b;
      D_uint8 *pData8;
      u8r = (D_uint8)R;
      u8g = (D_uint8)G;
      u8b = (D_uint8)B;
      len = 3 * _w * _h;
      pData8 = (D_uint8*)pData;
      for(int i = 0; i < len; i += 3){
	pData8[i] = u8r;
	pData8[i+1] = u8g;
	pData8[i+2] = u8b;
      }
      break;
    case DImage_RGB_16:
      D_uint16 u16r, u16g, u16b;
      D_uint16 *pData16;
      u16r = (D_uint16)R;
      u16g = (D_uint16)G;
      u16b = (D_uint16)B;
      len = 3 * _w * _h;
      pData16 = (D_uint16*)pData;
      for(int i = 0; i < len; i += 3){
	pData16[i] = u16r;
	pData16[i+1] = u16g;
	pData16[i+2] = u16b;
      }
      break;
    case DImage_flt_multi:
      fprintf(stderr, "DImage::fill() RGB called for float image\n");
      break;
    case DImage_dbl_multi:
      fprintf(stderr, "DImage::fill() RGB called for double image\n");
      break;
    case DImage_cmplx:
      fprintf(stderr, "DImage::fill() RGB called for complex image\n");
      break;
  }//end switch(_imgType)

}


///Fill the image with complex value complexVal
/**Fills each pixel in the image with complexVal The buffer should
 * already be allocated.  If the size has been changed without
 * changing the actual allocated buffer size, only width() by height()
 * pixels will be filled.
 */
void DImage::fill(std::complex<double> complexVal){
  int len;
  if((_w < 1) || (_h < 1)){
    fprintf(stderr, "DImage::fill() called for image with bad size (%dx%d)\n",
	    _w, _h);
    abort();
  }
  if(NULL == pData){
    fprintf(stderr, "DImage::fill() called for image with pData==NULL\n");
    abort();
  }
  switch(_imgType){
    case DImage_u8:
    case DImage_u16:
    case DImage_u32:
      fprintf(stderr, "DImage::fill() complex called for grayscale image\n");
      break;
    case DImage_RGB:
      fprintf(stderr, "DImage::fill() complex called for RGB image\n");
      break;
    case DImage_RGB_16:
      fprintf(stderr, "DImage::fill() complex called for RGB16 image\n");
      break;
    case DImage_flt_multi:
      fprintf(stderr, "DImage::fill() complex called for float image\n");
      break;
    case DImage_dbl_multi:
      fprintf(stderr, "DImage::fill() complex called for double image\n");
      break;
    case DImage_cmplx:
      std::complex<double> *pDataCmplx;
      len = _w * _h;
      pDataCmplx = (std::complex<double> *)pData;
      len = _w * _h;
      for(int i = 0; i < len; ++i){
	pDataCmplx[i] = complexVal;
      }
      break;
  }//end switch(_imgType)
}
///sets all pixel values in data buffer to zero
void DImage::clear(){
  if(pData == NULL){
    fprintf(stderr, "DImage::clear() pData = NULL!\n");
    abort();
  }
  memset(pData, 0, _dataSize);
}


///Returns an image the (pixel-wise) product of this image and imgOther
/**Each pixel from this image is multiplied by each pixel of imgOther
 * and the product is stored in the result image, which is returned.
 * The two images being multiplied must be of the same image type and
 * have the same width and height.  The result image will be of the
 * same type, width, and height as this image and imgOther.  In the
 * case of 8-bit images, you may want to first convert to another type
 * to avoid overflow, and then scale appropriately and convert back
 * after the multiplication is done.  For RGB images, the R,G, and B
 * values in each pixel are multiplied independently with the
 * corresponding R, G, and B values of imgOther.
 */
DImage DImage::multiplyByImage(const DImage imgOther){
  DImage dst;
  multiplyByImage_(dst, imgOther);
  return dst;
}

///Multiply this image by imgOther and store the (pixel-wise) product in imgDst
/**Each pixel from this image is multiplied by each pixel of imgOther
 * and the product is stored in the result image, which is imgDst.
 * The two images being multiplied must be of the same image type and
 * have the same width and height.  The result image will be of the
 * same type, width, and height as this image and imgOther.  The
 * source images and destination image MAY be the same image.  In the
 * case of 8-bit images, you may want to first convert to another type
 * to avoid overflow, and then scale appropriately and convert back
 * after the multiplication is done.  For RGB images, the R,G, and B
 * values in each pixel are multiplied independently with the
 * corresponding R, G, and B values of imgOther.
 */
void DImage::multiplyByImage_(DImage &imgDst,const DImage imgOther){
  int len;
  if(this->getImageType() != imgOther.getImageType()){
    fprintf(stderr, "DImage::multiplyByImage_() requires both images to be "
	    "of the same type.\n");
    abort();
  }
  if((this->width() != imgOther.width()) ||
     (this->height() != imgOther.height()) ||
     (this->_numChan != imgOther._numChan)){
    fprintf(stderr, "DImage::multiplyByImage_() requires both images to be "
	    "of the same size and number of channels.\n");
    abort();
  }
  len = _numChan * _w * _h;
  if( ((&imgDst) != this) && ((&imgDst)!= (&imgOther))){
    imgDst.create(_w, _h, _imgType, _numChan, _allocMethod);
  }
  switch(_imgType){
    case DImage_u8:
    case DImage_RGB:
      {
	D_uint8 *s1u8;
	D_uint8 *s2u8;
	D_uint8 *dst8;
	s1u8 = dataPointer_u8();
	s2u8 = imgOther.dataPointer_u8();
	dst8 = imgDst.dataPointer_u8();
	for(int idx = 0; idx < len; ++idx){
	  dst8[idx] = s1u8[idx] * s2u8[idx];
	}
      }
      break;
    case DImage_u16:
    case DImage_RGB_16:
      {
	D_uint16 *s1u16;
	D_uint16 *s2u16;
	D_uint16 *dst16;
	s1u16 = dataPointer_u16();
	s2u16 = imgOther.dataPointer_u16();
	dst16 = imgDst.dataPointer_u16();
	for(int idx = 0; idx < len; ++idx){
	  dst16[idx] = s1u16[idx] * s2u16[idx];
	}
      }
      break;
    case DImage_u32:
      {
	D_uint32 *s1u32;
	D_uint32 *s2u32;
	D_uint32 *dst32;
	s1u32 = dataPointer_u32();
	s2u32 = imgOther.dataPointer_u32();
	dst32 = imgDst.dataPointer_u32();
	for(int idx = 0; idx < len; ++idx){
	  dst32[idx] = s1u32[idx] * s2u32[idx];
	}
      }
      break;
    case DImage_flt_multi:
      {
	float *s1flt;
	float *s2flt;
	float *dstflt;
	s1flt = dataPointer_flt();
	s2flt = imgOther.dataPointer_flt();
	dstflt = imgDst.dataPointer_flt();
	for(int idx = 0; idx < len; ++idx){
	  dstflt[idx] = s1flt[idx] * s2flt[idx];
	}
      }
      break;
    case DImage_dbl_multi:
      {
	double *s1dbl;
	double *s2dbl;
	double *dstdbl;
	s1dbl = dataPointer_dbl();
	s2dbl = imgOther.dataPointer_dbl();
	dstdbl = imgDst.dataPointer_dbl();
	for(int idx = 0; idx < len; ++idx){
	  dstdbl[idx] = s1dbl[idx] * s2dbl[idx];
	}
      }
      break;
    case DImage_cmplx:
      {
	std::complex<double> *s1cmplx;
	std::complex<double> *s2cmplx;
	std::complex<double> *dstcmplx;
	s1cmplx = dataPointer_cmplx();
	s2cmplx = imgOther.dataPointer_cmplx();
	dstcmplx = imgDst.dataPointer_cmplx();
	for(int idx = 0; idx < len; ++idx){
	  dstcmplx[idx] = s1cmplx[idx] * s2cmplx[idx];
	}
      }
      break;
  }
}



///Returns an image with the (pixel-wise) sum of this image and imgOther
/**Each pixel from this image is added to each pixel of imgOther and
 * the sum is stored in the result image, which is returned.  The two
 * images being added must be of the same image type and have the same
 * width and height.  The result image will be of the same type,
 * width, and height as this image and imgOther.  In the case of 8-bit
 * images, you may want to first convert to another type to avoid
 * overflow, and then scale appropriately and convert back after the
 * addition is done.  For RGB images, the R,G, and B values in each
 * pixel are added independently with the corresponding R, G, and B
 * values of imgOther.
 */
DImage DImage::addToImage(const DImage imgOther){
  DImage dst;
  addToImage_(dst, imgOther);
  return dst;
}

///Add this image to imgOther and store the (pixel-wise) sum in imgDst
/**Each pixel from this image is added to each pixel of imgOther and
 * the sum is stored in the result image, which is imgDst.  The two
 * images being added must be of the same image type and have the same
 * width and height.  The result image will be of the same type,
 * width, and height as this image and imgOther.  The source images
 * and destination image MAY be the same image.  In the case of 8-bit
 * images, you may want to first convert to another type to avoid
 * overflow, and then scale appropriately and convert back after the
 * addition is done.  For RGB images, the R,G, and B values in each
 * pixel are added independently with the corresponding R, G, and B
 * values of imgOther.
 */
void DImage::addToImage_(DImage &imgDst,const DImage imgOther){
  int len;
  if(this->getImageType() != imgOther.getImageType()){
    fprintf(stderr, "DImage::addToImage_() requires both images to be "
	    "of the same type.\n");
    abort();
  }
  if((this->width() != imgOther.width()) ||
     (this->height() != imgOther.height()) ||
     (this->_numChan != imgOther._numChan)){
    fprintf(stderr, "DImage::addToImage_() requires both images to be "
	    "of the same size and number of channels.\n");
    abort();
  }
  len = _numChan * _w * _h;
  if( ((&imgDst) != this) && ((&imgDst)!= (&imgOther))){
    imgDst.create(_w, _h, _imgType, _numChan, _allocMethod);
  }
  switch(_imgType){
    case DImage_u8:
    case DImage_RGB:
      {
	D_uint8 *s1u8;
	D_uint8 *s2u8;
	D_uint8 *dst8;
	s1u8 = dataPointer_u8();
	s2u8 = imgOther.dataPointer_u8();
	dst8 = imgDst.dataPointer_u8();
	for(int idx = 0; idx < len; ++idx){
	  dst8[idx] = s1u8[idx] + s2u8[idx];
	}
      }
      break;
    case DImage_u16:
    case DImage_RGB_16:
      {
	D_uint16 *s1u16;
	D_uint16 *s2u16;
	D_uint16 *dst16;
	s1u16 = dataPointer_u16();
	s2u16 = imgOther.dataPointer_u16();
	dst16 = imgDst.dataPointer_u16();
	for(int idx = 0; idx < len; ++idx){
	  dst16[idx] = s1u16[idx] + s2u16[idx];
	}
      }
      break;
    case DImage_u32:
      {
	D_uint32 *s1u32;
	D_uint32 *s2u32;
	D_uint32 *dst32;
	s1u32 = dataPointer_u32();
	s2u32 = imgOther.dataPointer_u32();
	dst32 = imgDst.dataPointer_u32();
	for(int idx = 0; idx < len; ++idx){
	  dst32[idx] = s1u32[idx] + s2u32[idx];
	}
      }
      break;
    case DImage_flt_multi:
      {
	float *s1flt;
	float *s2flt;
	float *dstflt;
	s1flt = dataPointer_flt();
	s2flt = imgOther.dataPointer_flt();
	dstflt = imgDst.dataPointer_flt();
	for(int idx = 0; idx < len; ++idx){
	  dstflt[idx] = s1flt[idx] + s2flt[idx];
	}
      }
      break;
    case DImage_dbl_multi:
      {
	double *s1dbl;
	double *s2dbl;
	double *dstdbl;
	s1dbl = dataPointer_dbl();
	s2dbl = imgOther.dataPointer_dbl();
	dstdbl = imgDst.dataPointer_dbl();
	for(int idx = 0; idx < len; ++idx){
	  dstdbl[idx] = s1dbl[idx] + s2dbl[idx];
	}
      }
      break;
    case DImage_cmplx:
      {
	std::complex<double> *s1cmplx;
	std::complex<double> *s2cmplx;
	std::complex<double> *dstcmplx;
	s1cmplx = dataPointer_cmplx();
	s2cmplx = imgOther.dataPointer_cmplx();
	dstcmplx = imgDst.dataPointer_cmplx();
	for(int idx = 0; idx < len; ++idx){
	  dstcmplx[idx] = s1cmplx[idx] + s2cmplx[idx];
	}
      }
      break;
  }
}

///blend two DImages and return the result
/**This is useful for generating crossfades between two images for
   example.  t is the time on the range 0.0 to 1.0.  The result is
   simply a weighted average.  So with t=0, the result equals img1,
   with t=1.0 the result equals img2. width and height will be the
   greater of the two images' widths and heights and coordinate 0,0 in
   each image will be aligned. */
DImage DImage::blendImages(DImage img1, DImage img2, double t){
  DImage imgResult;
  int w1,w2,w3, h1,h2,h3;

  if(img1._imgType != img2._imgType){
    fprintf(stderr,"DImage::blendImage() image types don't match\n");
    exit(1);
  }
  w1 = img1.width();
  w2 = img2.width();
  h1 = img1.height();
  h2 = img2.height();
  w3 = (w1 > w2) ? w1 : w2;
  h3 = (h1 > h2) ? h1 : h2;

  switch(img1._imgType){
    case DImage_u8:
      {
	D_uint8 *p1,*p2,*p3;
	p1 = img1.dataPointer_u8();
	p2 = img2.dataPointer_u8();
	imgResult.create(w3,h3,DImage_u8);
	p3 = imgResult.dataPointer_u8();
	for(int y=0; y < h3; ++y){
	  for(int x=0; x < w3; ++x){
	    double v1, v2;
	    if((x<w1)&&(y<h1))
	      v1 = (double)(unsigned int)(p1[y*w1+x]);
	    else
	      v1 = (double)(unsigned int)(p2[y*w2+x]);
	    if((x<w2)&&(y<h2))
	      v2 = (double)(unsigned int)(p2[y*w2+x]);
	    else
	      v2 = v1;
	    p3[y*w3+x] = (D_uint8)((1.-t)*v1 + t*v2);
	  }
	}
      }
    case DImage_RGB:
      {
	D_uint8 *p1,*p2,*p3;
	p1 = img1.dataPointer_u8();
	p2 = img2.dataPointer_u8();
	imgResult.create(w3,h3,DImage_RGB);
	p3 = imgResult.dataPointer_u8();
	for(int y=0; y < h3; ++y){
	  for(int x=0; x < w3; ++x){
	    double v1r, v2r;
	    double v1g, v2g;
	    double v1b, v2b;
	    if((x<w1)&&(y<h1)){
	      v1r = (double)(unsigned int)(p1[3*(y*w1+x)]);
	      v1g = (double)(unsigned int)(p1[3*(y*w1+x)+1]);
	      v1b = (double)(unsigned int)(p1[3*(y*w1+x)+2]);
	    }
	    else{
	      v1r = (double)(unsigned int)(p2[3*(y*w2+x)]);
	      v1g = (double)(unsigned int)(p2[3*(y*w2+x)+1]);
	      v1b = (double)(unsigned int)(p2[3*(y*w2+x)+2]);
	    }
	    if((x<w2)&&(y<h2)){
	      v2r = (double)(unsigned int)(p2[3*(y*w2+x)]);
	      v2g = (double)(unsigned int)(p2[3*(y*w2+x)+1]);
	      v2b = (double)(unsigned int)(p2[3*(y*w2+x)+2]);
	    }
	    else{
	      v2r = v1r;
	      v2g = v1g;
	      v2b = v1b;
	    }
	    p3[3*(y*w3+x)] = (D_uint8)((1.-t)*v1r + t*v2r);
	    p3[3*(y*w3+x)+1] = (D_uint8)((1.-t)*v1g + t*v2g);
	    p3[3*(y*w3+x)+2] = (D_uint8)((1.-t)*v1b + t*v2b);
	  }
	}
      }
      break;
    case DImage_u16:
    case DImage_RGB_16:
    case DImage_u32:
    case DImage_flt_multi:
    case DImage_dbl_multi:
    case DImage_cmplx:
      {
	fprintf(stderr,
		"DImage::blendImage() only implemented for u8 and RGB\n");
	exit(1);
      }
      break;
    default:
	fprintf(stderr,
		"DImage::blendImage() only implemented for u8 and RGB\n");
	exit(1);
  }
  return imgResult;
}
















#if 0



//draws a pixel into the grayscale image data buffer with transparency
/** If you don't need transparency, just use setPixel.  This
    function is like the drawLine function, but draws a pixel instead.*/
void DImage::drawPixelAntiAliased(double x, double y, int GS,
				  float transparency){
  float opacity; // opacity fraction (1. - transparency)
  D_uint8 *pu8;
  D_uint16 *pu16;
  D_uint32 *pu32;
  int idx;
#ifdef DEBUG
  if((transparency < 0.) || (transparency > 1.0))
    fprintf(stderr, "DImage::drawPixelAntiAliased() transparency should be [0..1]\n");
#endif

  if((this->_imgType == DImage_RGB)||(this->_imgType == DImage_RGB_16)){
#ifdef DEBUG
    fprintf(stderr, "DImage::drawPixelAntiAliased() grayscale called with color image\n");
#endif
    drawPixelAntiAliased(x,y,GS,GS,GS,transparency);
    return;
  }
  if((_imgType == DImage_flt_multi)||(_imgType == DImage_dbl_multi)||
     (_imgType == DImage_cmplx)){
    fprintf(stderr, "DImage::drawPixelAntiAliased() unsupported image type\n");
    abort();
    return;
  }

  if( (x < 0) || (y < 0) || (x>=_w) || (y >=_h)){
    fprintf(stderr, "DImage::drawPixelAntiAliased(%d,%d) out of bounds\n",x,y);
    return;
#ifdef DEBUG
    //    abort();
#endif
    return;
  }
  opacity = 1. - transparency;

  switch(_imgType){
    case DImage_u8:
      pu8 = pData;
      idx = (y*_w+x);
      //I don't bother checking for value wraparound caused by rounding error
      pu8[idx] = (D_uint8)(GS * opacity + pu8[idx] * transparency);
      break;
    case DImage_u16:
      pu16 = (D_uint16*)pData;
      idx = (y*_w+x);
      pu16[idx] = (D_uint16)(GS * opacity + pu16[idx] * transparency);
      break;
    case DImage_u32:
      pu32 = (D_uint32*)pData;
      idx = (y*_w+x);
      pu32[idx] = (D_uint32)(GS * opacity + pu32[idx] * transparency);
      break;
    default:
      fprintf(stderr, "DImage::drawPixelAntiAliased() unsupported image format\n");
      abort();
  }

}
//draws a pixel into the RGB image data buffer with transparency
/** If you don't need transparency, just use setPixel.  This
    function is like the drawLine function, but draws a pixel instead.*/
void DImage::drawPixelAntiAliased(double x, double y, int R, int G, int B, float transparency){
  float opacity; // opacity fraction (1. - transparency)
  D_uint8 *pu8;
  D_uint16 *pu16;
  int idx = 0;

#ifdef DEBUG
  if((transparency < 0.) || (transparency > 1.0)){
    fprintf(stderr, "DImage::drawPixelAntiAliased() transparency should be [0..1]\n");
    abort();
  }
#endif

  if((this->_imgType != DImage_RGB)&&(this->_imgType != DImage_RGB_16)){
#ifdef DEBUG
    fprintf(stderr, "DImage::drawPixelAntiAliased() rgb called with non-RGB image\n");
    abort();
#endif
    return;
  }

  if( (x < 0) || (y < 0) || (x>=_w) || (y >=_h)){
    fprintf(stderr, "DImage::drawPixelAntiAliased(%d,%d) out of bounds\n",x,y);
    return;
    //    abort();
  }

  opacity = 1. - transparency;

  switch(_imgType){
    case DImage_RGB:
      pu8 = pData;
      idx = (y*_w+x)*3;
      pu8[idx] = (D_uint8)(R * opacity + pu8[idx] * transparency);
      pu8[idx+1] = (D_uint8)(G * opacity + pu8[idx+1] * transparency);
      pu8[idx+2] = (D_uint8)(B * opacity + pu8[idx+2] * transparency);
      break;
    case DImage_RGB_16:
      pu16 = (D_uint16*)pData;
      idx = (y*_w+x)*3;
      pu16[idx] = (D_uint16)(R * opacity + pu16[idx] * transparency);
      pu16[idx+1] = (D_uint16)(G * opacity + pu16[idx+1] * transparency);
      pu16[idx+2] = (D_uint16)(B * opacity + pu16[idx+2] * transparency);
      break;
    default:
      fprintf(stderr, "DImage::drawLine() unsupported image format\n");
      abort();
  }

}


// draws a line into the image data buffer.  (not anti-aliased).
/** Supports only grayscale or indexed (paletized) RGB images.  This
 *  function is meant as a convenience function for drawing rough
 *  lines into an image for debug or display purposes.  It is not
 *  implemented efficiently nor nicely.  Do not use this function for
 *  more advanced drawing needs. For completely opaque lines (the
 *  default), transparency=0. For perfectly transparent, transparency = 1.0.
 */
void DImage::drawLineAntiAliased(int x0, int y0, int x1, int y1, int GS,
		      float transparency, bool fIgnoreClipWarning){
  double dx, dy, incrX, incrY;
  double dist;
  double dblX, dblY;
  int x, y, i, idx = 0;
  int tmp;
  float opacity; // opacity fraction (1. - transparency)
  D_uint8 *pu8;
  D_uint16 *pu16;
  D_uint32 *pu32;

#ifdef DEBUG
  if((transparency < 0.) || (transparency > 1.0))
    fprintf(stderr, "DImage::drawLineAntiAliased() transparency should be [0..1]\n");
#endif

  if((this->_imgType == DImage_RGB)||(this->_imgType == DImage_RGB_16)){
#ifdef DEBUG
    fprintf(stderr, "DImage::drawLineAntiAliased() grayscale called with color image\n");
#endif
    drawLineAntiAliased(x0,y0,x1,y1,GS,GS,GS,transparency);
    return;
  }
  if((_imgType == DImage_flt_multi)||(_imgType == DImage_dbl_multi)||
     (_imgType == DImage_cmplx)){
    fprintf(stderr, "DImage::drawLineAntiAliased() unsupported image type\n");
    abort();
    return;
  }

  if( (x0 < 0) || (x1 < 0) ||  (y0 < 0) || (y1 < 0) ||
      (x0>=_w) || (x1>=_w) ||  (y0 >=_h) || (y1 >=_h)){
    if(!fIgnoreClipWarning)
      fprintf(stderr, "DImage::drawLineAntiAliased(%d,%d,%d,%d) line clipping not yet implemented\n",x0,y0,x1,y1);
    return;
#ifdef DEBUG
    //    abort();
#endif
    return;
  }

  opacity = 1. - transparency;

  if(x1 < x0){
    tmp = x1;
    x1 = x0;
    x0 = tmp;
    tmp = y1;
    y1 = y0;
    y0 = tmp;
  }
  
  dx = x1-x0;
  dy = y1-y0;
  dist = sqrt(dx*dx + dy*dy);
  incrX = dx / dist;
  incrY = dy / dist;
  dblX = x0;
  dblY = y0;

  switch(_imgType){
    case DImage_u8:
      pu8 = pData;
      for(i = 0; i <= dist; ++i){
	x = (int)(dblX+0.5); // round to nearest pixel
	y = (int)(dblY+0.5);
	idx = (y*_w+x);
	//I don't bother checking for value wraparound caused by rounding error
	pu8[idx] = (D_uint8)(GS * opacity + pu8[idx] * transparency);
	dblX += incrX;
	dblY += incrY;
      }
      if(idx != (y1*_w+x1)){//missed last point due to rounding
	idx = (y1*_w+x1);
	pu8[idx] = (D_uint8)(GS * opacity + pu8[idx] * transparency);
      }
      break;
    case DImage_u16:
      pu16 = (D_uint16*)pData;
      for(i = 0; i <= dist; ++i){
	x = (int)(dblX+0.5);
	y = (int)(dblY+0.5);
	idx = (y*_w+x);
	pu16[idx] = (D_uint16)(GS * opacity + pu16[idx] * transparency);
	dblX += incrX;
	dblY += incrY;
      }
      if(idx != (y1*_w+x1)){//missed last point due to rounding
	idx = (y1*_w+x1);
	pu16[idx] = (D_uint16)(GS * opacity + pu16[idx] * transparency);
      }
      break;
    case DImage_u32:
      pu32 = (D_uint32*)pData;
      for(i = 0; i <= dist; ++i){
	x = (int)(dblX+0.5);
	y = (int)(dblY+0.5);
	idx = (y*_w+x);
	pu32[idx] = (D_uint32)(GS * opacity + pu32[idx] * transparency);
	dblX += incrX;
	dblY += incrY;
      }
      if(idx != (y1*_w+x1)){//missed last point due to rounding
	idx = (y1*_w+x1);
	pu32[idx] = (D_uint32)(GS * opacity + pu32[idx] * transparency);
      }
      break;
    default:
      fprintf(stderr, "DImage::drawLineAntiAliased() unsupported image format\n");
      abort();
  }

}
void DImage::drawLineAntiAliased(int x0, int y0, int x1, int y1, int R, int G, int B,
		      float transparency, bool fIgnoreClipWarning){
  double dx, dy, incrX, incrY;
  double dist;
  double dblX, dblY;
  int x, y, i, idx = 0;
  int tmp;
  float opacity; // opacity fraction (1. - transparency)
  D_uint8 *pu8;
  D_uint16 *pu16;

#ifdef DEBUG
  if((transparency < 0.) || (transparency > 1.0)){
    fprintf(stderr, "DImage::drawLineAntiAliased() transparency should be [0..1]\n");
    abort();
  }
#endif

  if((this->_imgType != DImage_RGB)&&(this->_imgType != DImage_RGB_16)){
#ifdef DEBUG
    fprintf(stderr, "DImage::drawLineAntiAliased() rgb called with non-RGB image\n");
    abort();
#endif
    return;
  }

  if( (x0 < 0) || (x1 < 0) ||  (y0 < 0) || (y1 < 0) ||
      (x0>=_w) || (x1>=_w) ||  (y0 >=_h) || (y1 >=_h)){
    if(!fIgnoreClipWarning)
      fprintf(stderr, "DImage::drawLineAntiAliased(%d,%d,%d,%d) line clipping not yet implemented\n",x0,y0,x1,y1);
    return;
    //    abort();
  }

  opacity = 1. - transparency;

  if(x1 < x0){
    tmp = x1;
    x1 = x0;
    x0 = tmp;
    tmp = y1;
    y1 = y0;
    y0 = tmp;
  }
  
  dx = x1-x0;
  dy = y1-y0;
  dist = sqrt(dx*dx + dy*dy);
  incrX = dx / dist;
  incrY = dy / dist;
  dblX = x0;
  dblY = y0;

  switch(_imgType){
    case DImage_RGB:
      pu8 = pData;
      for(i = 0; i <= dist; ++i){
	x = (int)(dblX+0.5);
	y = (int)(dblY+0.5);
	idx = (y*_w+x)*3;
	// I won't bother checking for wrap-around due to rounding error
	pu8[idx] = (D_uint8)(R * opacity + pu8[idx] * transparency);
	pu8[idx+1] = (D_uint8)(G * opacity + pu8[idx+1] * transparency);
	pu8[idx+2] = (D_uint8)(B * opacity + pu8[idx+2] * transparency);
	dblX += incrX;
	dblY += incrY;
      }
      if(idx != (y1*_w+x1)*3){//missed last point due to rounding
	idx = (y1*_w+x1)*3;
	pu8[idx] = (D_uint8)(R * opacity + pu8[idx] * transparency);
	pu8[idx+1] = (D_uint8)(G * opacity + pu8[idx+1] * transparency);
	pu8[idx+2] = (D_uint8)(B * opacity + pu8[idx+2] * transparency);
      }
      break;
    case DImage_RGB_16:
      pu16 = (D_uint16*)pData;
      for(i = 0; i <= dist; ++i){
	x = (int)(dblX+0.5);
	y = (int)(dblY+0.5);
	idx = (y*_w+x)*3;
	pu16[idx] = (D_uint16)(R * opacity + pu16[idx] * transparency);
	pu16[idx+1] = (D_uint16)(G * opacity + pu16[idx+1] * transparency);
	pu16[idx+2] = (D_uint16)(B * opacity + pu16[idx+2] * transparency);
	dblX += incrX;
	dblY += incrY;
      }
      if(idx != (y1*_w+x1)*3){//missed last point due to rounding
	idx = (y1*_w+x1)*3;
	pu16[idx] = (D_uint16)(R * opacity + pu16[idx] * transparency);
	pu16[idx+1] = (D_uint16)(G * opacity + pu16[idx+1] * transparency);
	pu16[idx+2] = (D_uint16)(B * opacity + pu16[idx+2] * transparency);
      }
      break;
    default:
      fprintf(stderr, "DImage::drawLineAntiAliased() unsupported image format\n");
      abort();
  }

}




#endif //antialiased draw functions

void DImage::getTextSize(const char *st, int pointSize, const char *stFont,
			 int *w, int *h){
  char stTmpFileName[2048];
  char stCmd[2048];
  // temp filename will be combination of the 'this' pointer address and 
  int fileNum = 0;

  while((fileNum+1) != DImage::_drawTextFileNum){
    fileNum = DImage::_drawTextFileNum;
    ++(DImage::_drawTextFileNum);
  }
  sprintf(stTmpFileName,"/tmp/dimage_drawtext_%p_%d.pgm",this,fileNum);
  
  sprintf(stCmd,"convert -background black -fill white -pointsize %d %s label:'%s' -depth 8 %s",pointSize, (NULL==stFont)?"":stFont, st, stTmpFileName);
  system(stCmd);

  DImage imgTxt;
  if(!imgTxt.load(stTmpFileName)){
    fprintf(stderr,"DImage::drawText() couldn't load temp image '%s'\n",
	    stTmpFileName);
    return;
  }
  sprintf(stCmd,"rm %s", stTmpFileName);
  system(stCmd);

  int tw, th;
  D_uint8 *p8txt;
  tw = imgTxt.width();
  th = imgTxt.height();
  if(NULL != w)
    (*w) = tw;
  if(NULL != h)
    (*h) = th;
}


///uses imagemagick 'convert' program to render text, then copies it into this
/**Warning: The temporary filename generated in this function is not
   thread-safe. st and stFont could also be exploited for shell access
   if they are passed in from user-supplied data instead of
   programmatically. Draws the text in st at position x,y using
   grayscale value gs_val and transparency 0..1 to blend with the
   image.  If the caller needs the width and height of the resulting
   text, w and h should be set to non-NULL pointing to integers.
   offsPctX and offsPctY can be used to automatically control the text
   offset around x,y (for example if offsPctX=-0.5 then the text will
   be centered horizontally around x). Eventually I will probably want
   to use freeType2 or some other library instead of slowly generating
   external images and copying them into this image. */
void DImage::drawText(const char *st, int x, int y, int pointSize,
		      int gs_val, double transparency,
		      const char *stFont, int *w, int *h,
		      double offsPctX, double offsPctY){
  char stTmpFileName[2048];
  char stCmd[2048];
  // temp filename will be combination of the 'this' pointer address and 
  int fileNum = 0;

  if(_imgType != DImage::DImage_u8){
    fprintf(stderr,"DImage::drawText() 8-bit Grayscale version called for other image type\n");
    if(_imgType == DImage::DImage_RGB){
      fprintf(stderr,"calling drawTextRGB() instead\n");
      drawTextRGB(st,x,y,pointSize,gs_val,gs_val,gs_val, transparency, stFont,
		  w, h, offsPctX, offsPctY);
      return;
    }
    exit(1);
  }

  while((fileNum+1) != DImage::_drawTextFileNum){
    fileNum = DImage::_drawTextFileNum;
    ++(DImage::_drawTextFileNum);
  }
  sprintf(stTmpFileName,"/tmp/dimage_drawtext_%p_%d.pgm",this,fileNum);
  
  sprintf(stCmd,"convert -background black -fill white -pointsize %d %s label:'%s' -depth 8 %s",pointSize, (NULL==stFont)?"":stFont, st, stTmpFileName);
  system(stCmd);

  DImage imgTxt;
  if(!imgTxt.load(stTmpFileName)){
    fprintf(stderr,"DImage::drawText() couldn't load temp image '%s'\n",
	    stTmpFileName);
    return;
  }
  sprintf(stCmd,"rm %s", stTmpFileName);
  system(stCmd);

  int tw, th;
  D_uint8 *p8txt;
  tw = imgTxt.width();
  th = imgTxt.height();
  if(NULL != w)
    (*w) = tw;
  if(NULL != h)
    (*h) = th;

  p8txt = imgTxt.dataPointer_u8();
  double opaqueNess;
  int offsX, offsY;
  offsX = tw * offsPctX + x;
  offsY = th * offsPctY + y;
  for(int ty=0, tidx=0; ty < th; ++ty){
    for(int tx=0; tx < tw; ++tx, ++tidx){
      if(((offsX+tx)>=0) && ((offsX+tx)<_w) && 
	 ((offsY+ty)>=0) && ((offsY+ty)<_h)){
	int dstIdx;
	int newVal;
	opaqueNess = (1.-transparency) * (p8txt[tidx]/255.);
	dstIdx = _w*(offsY+ty)+offsX+tx;
	newVal = (1.-opaqueNess)*this->pData[dstIdx*_numChan] +
	  opaqueNess*gs_val;
	if(newVal > 255)
	  newVal = 255;
	this->pData[dstIdx*_numChan] = (D_uint8)newVal;
      }
    }
  }
}

///uses imagemagick 'convert' program to render text, then copies it into this
/**Warning: The temporary filename generated in this function is not
   thread-safe. st and stFont could also be exploited for shell access
   if they are passed in from user-supplied data instead of
   programmatically. Draws the text in st at position x,y using
   grayscale value gs_val and transparency 0..1 to blend with the
   image.  If the caller needs the width and height of the resulting
   text, w and h should be set to non-NULL pointing to integers.
   offsPctX and offsPctY can be used to automatically control the text
   offset around x,y (for example if offsPctX=-0.5 then the text will
   be centered horizontally around x). Eventually I will probably want
   to use freeType2 or some other library instead of slowly generating
   external images and copying them into this image. */
void DImage::drawTextRGB(const char *st, int x, int y, int pointSize,
			 int R, int G, int B, double transparency,
			 const char *stFont, int *w, int *h,
			 double offsPctX, double offsPctY){
  char stTmpFileName[2048];
  char stCmd[2048];
  // temp filename will be combination of the 'this' pointer address and 
  int fileNum = 0;

  if(_imgType != DImage::DImage_RGB){
    fprintf(stderr,"DImage::drawText() RGB version called for other image type\n");
    exit(1);
  }

  while((fileNum+1) != DImage::_drawTextFileNum){
    fileNum = DImage::_drawTextFileNum;
    ++(DImage::_drawTextFileNum);
  }
  sprintf(stTmpFileName,"/tmp/dimage_drawtext_%p_%d.pgm",this,fileNum);
  
  sprintf(stCmd,"convert -background black -fill white -pointsize %d %s label:'%s' -depth 8 %s",pointSize, (NULL==stFont)?"":stFont, st, stTmpFileName);
  system(stCmd);

  DImage imgTxt;
  if(!imgTxt.load(stTmpFileName)){
    fprintf(stderr,"DImage::drawText() couldn't load temp image '%s'\n",
	    stTmpFileName);
    return;
  }
  sprintf(stCmd,"rm %s", stTmpFileName);
  system(stCmd);

  int tw, th;
  D_uint8 *p8txt;
  tw = imgTxt.width();
  th = imgTxt.height();
  if(NULL != w)
    (*w) = tw;
  if(NULL != h)
    (*h) = th;

  p8txt = imgTxt.dataPointer_u8();
  double opaqueNess;
  int offsX, offsY;
  offsX = tw * offsPctX + x;
  offsY = th * offsPctY + y;
  for(int ty=0, tidx=0; ty < th; ++ty){
    for(int tx=0; tx < tw; ++tx, ++tidx){
      if(((offsX+tx)>=0) && ((offsX+tx)<_w) && 
	 ((offsY+ty)>=0) && ((offsY+ty)<_h)){
	int dstIdx;
	int newR, newG, newB;
	opaqueNess = (1.-transparency) * (p8txt[tidx]/255.);
	dstIdx = _w*(offsY+ty)+offsX+tx;
	newR = (1.-opaqueNess)*this->pData[dstIdx*_numChan] +
	  opaqueNess*R;
	if(newR > 255)
	  newR = 255;
	this->pData[dstIdx*_numChan] = (D_uint8)newR;

	newG = (1.-opaqueNess)*this->pData[dstIdx*_numChan+1] +
	  opaqueNess*G;
	if(newG > 255)
	  newG = 255;
	this->pData[dstIdx*_numChan+1] = (D_uint8)newG;

	newB = (1.-opaqueNess)*this->pData[dstIdx*_numChan+2] +
	  opaqueNess*B;
	if(newB > 255)
	  newB = 255;
	this->pData[dstIdx*_numChan+2] = (D_uint8)newB;
      }
    }
  }
}


