#include "dbilateralfilter.h"
#include "dimage.h"
#include "dprogress.h"
#include "dinstancecounter.h"
#include <stdio.h>

#ifdef CHRONO
#undef CHRONO
#endif
#define NO_XML
#define USE_DPROGRESS
#include "paris_sylvain_bilateral_headers/fast_color_bf.h"
#include "paris_sylvain_bilateral_headers/fast_lbf.h"

#include "dtimer.h" // for debug timing

/// \cond
// doxygen will ignore anything between the cond and endcond commands

typedef Array_2D<double> gs_image_type;
typedef Array_2D<Geometry::Vec3<double> > rgb_image_type;



/// \endcond


DBilateralFilter::DBilateralFilter(){
  DInstanceCounter::addInstance("DBilateralFilter");
  _bilatFiltType = DBilatFilt_paris;
  _sigma_s = 16;
  _sigma_r = 0.1;
  _numThreads = 1;
}

DBilateralFilter::DBilateralFilter(double sigma_s, double sigma_r,
				   DBilatFiltType filtType){
  DInstanceCounter::addInstance("DBilateralFilter");
  _bilatFiltType = filtType;
  _sigma_s = sigma_s;
  _sigma_r = sigma_r;
  _numThreads = 1;
}

DBilateralFilter::~DBilateralFilter(){
  DInstanceCounter::removeInstance("DBilateralFilter");
}
void DBilateralFilter::setSigmas(double sigma_s, double sigma_r){
  _sigma_s = sigma_s;
  _sigma_r = sigma_r;
}
double DBilateralFilter::getSigmaS(){
  return _sigma_s;
}
double DBilateralFilter::getSigmaR(){
  return _sigma_r;
}
void DBilateralFilter::setType(DBilatFiltType filtType){
  _bilatFiltType = filtType;
}
DBilateralFilter::DBilatFiltType DBilateralFilter::getType(){
  return _bilatFiltType;
}
void DBilateralFilter::setNumThreads(int numThreads){
  if(numThreads > 1){
    fprintf(stderr, "DBilateralFilter::setNumThreads(%d) Warning! multiple threads not currently supported. Setting numThreads to 1.\n", numThreads);
  }
  _numThreads = 1;
}
int DBilateralFilter::getNumThreads(){
  return _numThreads;
}
///Bilateral filter imgSrc with current settings and return the result
DImage DBilateralFilter::filterImage(const DImage &imgSrc,
				     bool fAlreadyPadded, DProgress *pProg){
  DImage dst;
  filterImage_(dst, imgSrc, fAlreadyPadded, pProg);
  return dst;
}

///Static function provided for convenience
void DBilateralFilter::bilateralFilterImage(DImage &imgDst,
					    const DImage &imgSrc,
					    bool fAlreadyPadded,
					    double sigma_s, double sigma_r,
					    DBilatFiltType filtType,
					    DProgress *pProg,
					    int numThreads){
  DBilateralFilter bfilt(sigma_s, sigma_r, filtType);
  bfilt.filterImage_(imgDst, imgSrc, fAlreadyPadded, pProg);
}

///Bilateral filter imgSrc with current settings and store result in imgDst
void DBilateralFilter::filterImage_(DImage &imgDst, const DImage &imgSrc,
				    bool fAlreadyPadded, DProgress *pProg){
  DBilatFiltType filtType;
  filtType = _bilatFiltType;
  int isColor = 0;
  DImage::DImageType srcType;
  int width, height;
  D_uint8 *p8Tmp,*p8Dst;
  D_uint16 *p16Tmp,*p16Dst;
  double *pdblTmp,*pdblDst;
  float *pfltTmp,*pfltDst;
  double max, min, range;
  double dblT;
  unsigned int offs;
  DTimer t1;
  float split1,split2,split3;
  
  if(DBilatFilt_paris != filtType){
    fprintf(stderr, "DBilateralFilter::filterImage_() currently only supports "
	    "filter type DBilatFilt_paris.  Using that filter type.\n");
    filtType = DBilatFilt_paris;
  }
  
  width = imgSrc.width();
  height = imgSrc.height();
  srcType = imgSrc.getImageType();
  if((DImage::DImage_RGB == srcType) || (DImage::DImage_RGB_16 == srcType)){
    isColor = 1;
  }
  printf("DBilateralFilter::filterImage_() copying image to data structs\n");
  if(!isColor){ // grayscale
    if(imgSrc.numChannels() != 1){
      fprintf(stderr, "DBilateralFilter::filterImage_() only supports "
	      "grayscale/rgb images or single-channel float/double images\n");
      exit(1);
    }
    // copy the image data into Sylvain Paris's data structures
    gs_image_type gsimage(width, height);
    gs_image_type filtered_gsimage(width,height);
    imgDst.create(width,height,srcType,1,imgSrc.getAllocMethod());

    t1.start();
    switch(srcType){
      case DImage::DImage_u8:
	{
	  p8Tmp = imgSrc.dataPointer_u8();
	  max = 255.;
	  min = 0.;
	  range = 255.;
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, ++offs){
	      gsimage(x,y) = p8Tmp[offs] / 255.;
	    }
	  }
	}
	break;
      case DImage::DImage_u16:
	{
	  p16Tmp = imgSrc.dataPointer_u16();
	  max = 65535.;
	  min = 0.;
	  range = 65535.;
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, ++offs){
	      gsimage(x,y) = p16Tmp[offs] / 65535.;
	    }
	  }
	}
	break;
      case DImage::DImage_flt_multi:
	{
	  pfltTmp = imgSrc.dataPointer_flt();
	  max = min = pfltTmp[0];
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, ++offs){
	      if(pfltTmp[offs] > max)
		max = pfltTmp[offs];
	      if(pfltTmp[offs] < min)
		min = pfltTmp[offs];
	    }
	  }
	  range = max - min;
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, ++offs){
	      gsimage(x,y) = (pfltTmp[offs] - min) / range;
	    }
	  }
	}
	break;
      case DImage::DImage_dbl_multi:
	{
	  pdblTmp = imgSrc.dataPointer_dbl();
	  max = min = pdblTmp[0];
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, ++offs){
	      if(pdblTmp[offs] > max)
		max = pdblTmp[offs];
	      if(pdblTmp[offs] < min)
		min = pdblTmp[offs];
	    }
	  }
	  range = max - min;
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, ++offs){
	      gsimage(x,y) = (pdblTmp[offs] - min) / range;
	    }
	  }
	}
	break;
      default:
	fprintf(stderr, "DBilateralFilter::filterImage_() ERROR! NYI srcType=%d\n",srcType);
	exit(1);
    }
    t1.stop();
    split1 = t1.getSplit();
    // filter the image using Sylvain Paris's code
    printf("DBilateralFilter::filterImage_() filtering image\n");
    t1.resume();
    Image_filter::fast_LBF(gsimage,gsimage,
			   _sigma_s,_sigma_r,
			   false,
			   &filtered_gsimage,&filtered_gsimage,pProg);
    t1.stop();
    split2 = t1.getSplit();
    t1.resume();
    // copy the results to the destination image
    switch(srcType){
      case DImage::DImage_u8:
	{
	  p8Dst = imgDst.dataPointer_u8();
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, ++offs){
	      dblT = range * filtered_gsimage(x,y) + min;
	      if(dblT <= max){
		p8Dst[offs] = (D_uint8)dblT;
		if(dblT < 0.)
		  p8Dst[offs] = 0x00;
	      }
	      else{
		p8Dst[offs] = 0xff;
	      }
	    }
	  }
	}
	break;
      case DImage::DImage_u16:
	{
	  p16Dst = imgDst.dataPointer_u16();
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, ++offs){
	      dblT = range * filtered_gsimage(x,y) + min;
	      if(dblT <= max){
		p16Dst[offs] = (D_uint16)dblT;
		if(dblT < 0.)
		  p16Dst[offs] = 0x0000;
	      }
	      else{
		p16Dst[offs] = 0xffff;
	      }
	    }
	  }
	}
      case DImage::DImage_flt_multi:
	{
	  pfltDst = imgDst.dataPointer_flt();
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, ++offs){
	      dblT = range * filtered_gsimage(x,y) + min;
	      if(dblT <= max){
		pfltDst[offs] = (float)dblT;
		if(dblT < min)
		  pfltDst[offs] = (float)min;
	      }
	      else{
		pfltDst[offs] = (float)max;
	      }
	    }
	  }
	}
	break;
      case DImage::DImage_dbl_multi:
	{
	  pdblDst = imgDst.dataPointer_dbl();
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, ++offs){
	      dblT = range * filtered_gsimage(x,y) + min;
	      if(dblT <= max){
		pdblDst[offs] = (double)dblT;
		if(dblT < min)
		  pdblDst[offs] = min;
	      }
	      else{
		pdblDst[offs] = max;
	      }
	    }
	  }
	}
	break;
      default:
	fprintf(stderr, "DBilateralFilter::filterImage_() ERROR! Unsupported "
		"image type\n");
	exit(1);
    }
    t1.stop();
    split3 = t1.getSplit();
    
  }
  else{// color -------------------------------------------------------------
    // copy the image data into Sylvain Paris's data structures
    rgb_image_type rgbimage(width, height);
    rgb_image_type filtered_rgbimage(width,height);
    imgDst.create(width,height,srcType,3,imgSrc.getAllocMethod());

    switch(srcType){
      case DImage::DImage_RGB:
	{
	  p8Tmp = imgSrc.dataPointer_u8();
	  max = 255.;
	  min = 0.;
	  range = 255.;
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, offs+=3){
	      rgbimage(x,y)[0] = p8Tmp[offs] / 255.;
	      rgbimage(x,y)[1] = p8Tmp[offs+1] / 255.;
	      rgbimage(x,y)[2] = p8Tmp[offs+2] / 255.;
	    }
	  }
	}
	break;
      case DImage::DImage_RGB_16:
	{
	  p16Tmp = imgSrc.dataPointer_u16();
	  max = 65535.;
	  min = 0.;
	  range = 65535.;
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, offs+=3){
	      rgbimage(x,y)[0] = p16Tmp[offs] / 65535.;
	      rgbimage(x,y)[1] = p16Tmp[offs+1] / 65535.;
	      rgbimage(x,y)[2] = p16Tmp[offs+2] / 65535.;
	    }
	  }
	}
	break;
      default:
	fprintf(stderr, "DBilateralFilter::filterImage_() ERROR! (B) NYI srcType=%d\n",srcType);
	exit(1);
    }
    t1.stop();
    split1 = t1.getSplit();
    // filter the image using Sylvain Paris's code
    printf("DBilateralFilter::filterImage_() filtering image\n");
    t1.resume();
    Image_filter::fast_color_BF(rgbimage, _sigma_s, _sigma_r,
				&filtered_rgbimage,pProg);
    t1.stop();
    split2 = t1.getSplit();
    t1.resume();
    // copy the results to the destination image
    switch(srcType){
      case DImage::DImage_RGB:
	{
	  p8Dst = imgDst.dataPointer_u8();
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, offs+=3){
	      dblT = range * filtered_rgbimage(x,y)[0] + min;
	      if(dblT <= max){
		p8Dst[offs] = (D_uint8)dblT;
		if(dblT < 0.)
		  p8Dst[offs] = 0x00;
	      }
	      else{
		p8Dst[offs] = 0xff;
	      }
	      dblT = range * filtered_rgbimage(x,y)[1] + min;
	      if(dblT <= max){
		p8Dst[offs+1] = (D_uint8)dblT;
		if(dblT < 0.)
		  p8Dst[offs+1] = 0x00;
	      }
	      else{
		p8Dst[offs+1] = 0xff;
	      }
	      dblT = range * filtered_rgbimage(x,y)[2] + min;
	      if(dblT <= max){
		p8Dst[offs+2] = (D_uint8)dblT;
		if(dblT < 0.)
		  p8Dst[offs+2] = 0x00;
	      }
	      else{
		p8Dst[offs+2] = 0xff;
	      }
	    }
	  }
	}
	break;
      case DImage::DImage_RGB_16:
	{
	  p16Dst = imgDst.dataPointer_u16();
	  offs = 0;
	  for(int y = 0; y < height; ++y){
	    for(int x = 0; x < width; ++x, offs+=3){
	      dblT = range * filtered_rgbimage(x,y)[0] + min;
	      if(dblT <= max){
		p16Dst[offs] = (D_uint16)dblT;
		if(dblT < 0.)
		  p16Dst[offs] = 0x0000;
	      }
	      else{
		p16Dst[offs] = 0xffff;
	      }
	      dblT = range * filtered_rgbimage(x,y)[1] + min;
	      if(dblT <= max){
		p16Dst[offs+1] = (D_uint16)dblT;
		if(dblT < 0.)
		  p16Dst[offs+1] = 0x0000;
	      }
	      else{
		p16Dst[offs+1] = 0xffff;
	      }
	      dblT = range * filtered_rgbimage(x,y)[2] + min;
	      if(dblT <= max){
		p16Dst[offs+2] = (D_uint16)dblT;
		if(dblT < 0.)
		  p16Dst[offs+2] = 0x0000;
	      }
	      else{
		p16Dst[offs+2] = 0xffff;
	      }
	    }
	  }
	}
	break;
      default:
	fprintf(stderr, "DBilateralFilter::filterImage_() ERROR! (C) NYI srcType=%d\n",srcType);
	exit(1);
    }
    t1.stop();
    split3 = t1.getSplit();
    
  }
  if(NULL != pProg){
    pProg->reportStatus(width*12/10,0,width*12/10);
  }
  printf("DBilateralFilter::filterImage_() timings: copyIn:%.2fsec  filter:%.2fsec  copyOut:%.2fsec\n",split1,split2,split3);
}
