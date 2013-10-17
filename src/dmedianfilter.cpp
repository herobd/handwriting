#include "dmedianfilter.h"
#include "dimage.h"
#include "dprogress.h"
#include "dinstancecounter.h"
#include <string.h>

#ifndef D_NOTHREADS
#include "dthreads.h"

/*comparison function for qsort()ing doubles in nondecreasing order */
int compareDoubles(const void *p1, const void *p2){
  if ( (*(double*)(p1)) < (*(double*)(p2)) )
    return -1;
  else if ( (*(double*)(p1)) > (*(double*)(p2)) )
    return 1;
  return 0;
}

/*comparison function for qsort()ing floats in nondecreasing order */
int compareFloats(const void *p1, const void *p2){
  if ( (*(float*)(p1)) < (*(float*)(p2)) )
    return -1;
  else if ( (*(float*)(p1)) > (*(float*)(p2)) )
    return 1;
  return 0;
}


#if 1
int med_compareChars(const void *p1, const void *p2){
  if ( (*((unsigned char*)(p1))) < (*((unsigned char*)(p2))) )
    return -1;
  else if ( (*((unsigned char*)(p1))) > (*((unsigned char*)(p2))) )
    return 1;
  return 0;
}
#endif

//structure for passing parameters to the thread function
/// \cond
typedef struct{
  DImage *pImgDst;
  const DImage *pImgSrc;
  int radiusX;
  int radiusY;
  int wKern;
  int hKern;
  D_uint8 *rgKern;
  int numKernPxls;
  int *rgRightEdge;
  DProgress *pProg;
  int progStart;
  int progMax;
  int threadNumber;
  int numThreads;
} HUANG_8_THREAD_PARAMS_T;
/// \endcond

#define MAX_MEDFILT_THREADS 64

void* DMedianFilter::DMedianFilter_Huang8threadWrap(void* params){
  HUANG_8_THREAD_PARAMS_T *pParams;

  pParams = (HUANG_8_THREAD_PARAMS_T*)params;
  DMedianFilter::medFiltHuang_u8(*(pParams->pImgDst),
				 *(pParams->pImgSrc),
				 pParams->radiusX, pParams->radiusY,
				 pParams->wKern, pParams->hKern,
				 pParams->rgKern,
				 pParams->numKernPxls, pParams->rgRightEdge,
				 pParams->pProg, pParams->progStart,
				 pParams->progMax, pParams->threadNumber,
				 pParams->numThreads);
  return NULL;
}
#endif

///Default constructor
DMedianFilter::DMedianFilter(){
  DInstanceCounter::addInstance("DMedianFilter");
  _medFiltType = DMedFilt_default;
  rgKern = NULL;
  rgRightEdge = NULL;
  _radiusX = 1;
  _radiusY = 1;
  _numThreads = 1;
}

///Constructor with parameters
DMedianFilter::DMedianFilter(int radiusX, int radiusY, DMedFiltType filtType){
  DInstanceCounter::addInstance("DMedianFilter");
  _medFiltType = filtType;
  rgKern = NULL;
  rgRightEdge = NULL;
  _radiusX = radiusX;
  _radiusY = radiusY;
  _numThreads = 1;
}


///Destructor
DMedianFilter::~DMedianFilter(){
  DInstanceCounter::removeInstance("DMedianFilter");
  if(rgKern){
    free(rgKern);
    rgKern = NULL;
  }
  if(rgRightEdge){
    free(rgRightEdge);
    rgRightEdge = NULL;
  }
}

///Set the kernel (window) radii for the median filter
/**The actual width (and height) of the window will be 2*radius+1*/
void DMedianFilter::setRadius(int radiusX, int radiusY){
  if((radiusX < 0)||(radiusY < 0)){
    fprintf(stderr, "DMedianFilter::setRadius() negative value not allowed\n");
    return;
  }
  _radiusX = radiusX;
  _radiusY = radiusY;
  if(rgKern){
    free(rgKern);
    rgKern = NULL;
    free(rgRightEdge);
    rgRightEdge = NULL;
  }
}

///Get the horizontal radius of the kernel
int DMedianFilter::getRadiusX(){
  return _radiusX;
}
///Get the vertical radius of the kernel
int DMedianFilter::getRadiusY(){
  return _radiusY;
}

///Set the number of threads that will be used to do the median filter
/** The median filter may run faster if more than one thread is used,
 *  especially if the machine has more than one processor.  With
 *  numThreads threads specified, thread 0 will start with y=0, thread
 *  1 with y=1, and so on, and each will skip numThreads rows between
 *  each row that it operates on.  Only thread 0 will report status.
 *  Note: our experiments show that using multiple threads does not
 *  seem to speed up the median filter very much.  For example, on a
 *  very large image with radius=35, we saw a speedup from 33.6
 *  seconds to 29.1 seconds with 2 threads, and to 28.7 sec with 4
 *  threads.  When the code is not optimized, there is a more
 *  noticible speedup, but we assume you will optimize your code.
 */
void DMedianFilter::setNumThreads(int numThreads){
#ifndef D_NOTHREADS
  if((_numThreads > 0) && (_numThreads <= MAX_MEDFILT_THREADS)){
    _numThreads = numThreads;
  }
  else{
    fprintf(stderr, "DMedianFilter::setNumThreads() numThreads must be from "
	    "1 to %d\n", MAX_MEDFILT_THREADS);
  }
#else
  fprintf(stderr, "DMedianFilter::setNumThreads() thread support not "
	  "compiled");
#endif
}

///return the number of threads that the filter has been set to use
int DMedianFilter::getNumThreads(){
  return _numThreads;
}

///set the type of median filter algorithm
/**By default, RGB and _u8 images will use Huang_circle and all others
 * (including RGB_16 and _u16) will use the slow method */
void DMedianFilter::setType(DMedFiltType filtType){
  if(rgKern){
    free(rgKern);
    rgKern = NULL;
    free(rgRightEdge);
    rgRightEdge = NULL;
  }
  _medFiltType = filtType;
}

///return the type of the median filter
DMedianFilter::DMedFiltType DMedianFilter::getType(){
  return _medFiltType;
}

///private function that actually performs huang med filt on 8-bit images
/** imgDst will be 2*radius pixels less wide and high than imgSrc
 * because of the padding that is added before calling this function.
 * This function requires that imgDst.create() has already been called
 * with the proper w,h,imgType,etc.
 */
void DMedianFilter::medFiltHuang_u8(DImage &imgDst, const DImage &imgSrc,
				    int radiusX, int radiusY,
				    int wKern, int hKern, D_uint8 *rgKern,
				    int numKernPxls, int *rgRightEdge,
				    DProgress *pProg,
				    int progStart, int progMax,
				    int threadNumber, int numThreads){
  int th;
  int rgHist[256];
  int med;
  int lastMed;
  int lt_med;
  int numHistVals;
  unsigned char valTmp;
  int idxDst;
  int idx3;
  D_uint8 *pTmp; // pointer to padded image data
  int wTmp, hTmp; // width, height of imgSrc
  int w, h; // width, height of imgDst
  int histCount;
  D_uint8 *pDst;

  th = numKernPxls / 2;
  wTmp = imgSrc.width();
  hTmp = imgSrc.height();
  w = wTmp - radiusX*2;
  h = hTmp - radiusY*2;
  pDst = imgDst.dataPointer_u8();
  pTmp = imgSrc.dataPointer_u8();

  for(int y = threadNumber; y < h; y += numThreads){
    // update progress report and check if user cancelled the operation
    if((NULL != pProg) && (0 == (y & 0x000000ff))){
      if(0 != pProg->reportStatus(progStart + y, 0, progMax)){
	// the operation has been cancelled
	pProg->reportStatus(-1, 0, progMax); // report cancel acknowledged
	return;
      }
    }

    // position window at the beginning of a new row and fill the kernel, hist
    memset(rgHist, 0, sizeof(int)*256);
    for(int kr = 0, kidx =0; kr < hKern; ++kr){
      for(int kc = 0; kc < wKern; ++kc, ++kidx){
	if(rgKern[kidx]){ // pixel is part of the kernel mask
	  ++(rgHist[pTmp[(y+kr)*wTmp+kc]]);//add pixel val to histogram
	}
      }
    }
    // calculate median for first spot
    med = 0;
    numHistVals = 0;
    lt_med = 0;
    while(lt_med <= th){
      lt_med += rgHist[med];
      ++med;
    }
    --med;
    lastMed = med;

    // put the median in the spot we're at
    idxDst = y*w;
    pDst[idxDst] = (unsigned char)med;
    
    // remove pixels from leftmost column
    idx3 = y*wTmp+radiusX;
    for(int ky = 0; ky < hKern; ++ky){
      valTmp = pTmp[idx3 - rgRightEdge[ky]];
      --(rgHist[valTmp]);
      idx3 += wTmp;
    }

    for(int x=1;  x < w;  ++x){
      ++idxDst;
      // add pixels from the right-hand side of kernel (after moving over one)
      idx3 = y*wTmp+x+radiusX;
      for(int ky = 0; ky < hKern; ++ky){
	valTmp = pTmp[idx3 + rgRightEdge[ky]];
	++(rgHist[valTmp]);
	idx3 += wTmp;
      }

      // find median from the histogram.
      // count high-to-low if lastMed was high, otherwise low-to-high
      histCount = 0;
      if(lastMed >= 128){
	med = 255;
	while(histCount <= th){
	  histCount += rgHist[med];
	  --med;
	}
	++med;
      }
      else{
	med = 0;
	while(histCount <= th){
	  histCount += rgHist[med];
	  ++med;
	}
	--med;
      }
      
      // put the median value in the destination pixel
      pDst[idxDst] = (unsigned char)med;

      // remove pixels from leftmost column for next time through loop
      if(x < (w-1)){//don't need to remove left edge if going to a new row
	idx3 = y*wTmp+x+radiusX;
	for(int ky = 0; ky < hKern; ++ky){
	  valTmp = pTmp[idx3 - rgRightEdge[ky]];
	  --(rgHist[valTmp]);
	  idx3 += wTmp;
	} // end for(ky...
      } // end if

    } // end for (x=1; ...

  }// end for(y=0...
  // report progress
  if(NULL != pProg){
    pProg->reportStatus(progStart + h, 0, progMax);
  }
}


///private function that actually performs huang med filt on 8-bit images
/** imgDst will be 2*radius pixels less wide and high than imgSrc
 * because of the padding that is added before calling this function.
 * This function requires that imgDst.create() has already been called
 * with the proper w,h,imgType,etc.
 */
void DMedianFilter::medFiltHuang_u8_square(DImage &imgDst,
					   const DImage &imgSrc,
					   int radiusX, int radiusY,
					   int wKern, int hKern,
					   D_uint8 *rgKern,
					   int numKernPxls,
					   DProgress *pProg,
					   int progStart, int progMax,
					   int threadNumber, int numThreads){
  int th;
  int rgHist[256];
  int med;
  int lastMed;
  int lt_med;
  int numHistVals;
  unsigned char valTmp;
  int idxDst;
  int idx3;
  D_uint8 *pTmp; // pointer to padded image data
  int wTmp, hTmp; // width, height of imgSrc
  int w, h; // width, height of imgDst
  int histCount;
  D_uint8 *pDst;

  th = numKernPxls / 2;
  wTmp = imgSrc.width();
  hTmp = imgSrc.height();
  w = wTmp - radiusX*2;
  h = hTmp - radiusY*2;
  pDst = imgDst.dataPointer_u8();
  pTmp = imgSrc.dataPointer_u8();

  for(int y = threadNumber; y < h; y += numThreads){
    // update progress report and check if user cancelled the operation
    if((NULL != pProg) && (0 == (y & 0x000000ff))){
      if(0 != pProg->reportStatus(progStart + y, 0, progMax)){
	// the operation has been cancelled
	pProg->reportStatus(-1, 0, progMax); // report cancel acknowledged
	return;
      }
    }

    // position window at the beginning of a new row and fill the kernel, hist
    memset(rgHist, 0, sizeof(int)*256);
    for(int kr = 0, kidx =0; kr < hKern; ++kr){
      for(int kc = 0; kc < wKern; ++kc, ++kidx){
	if(rgKern[kidx]){ // pixel is part of the kernel mask
	  ++(rgHist[pTmp[(y+kr)*wTmp+kc]]);//add pixel val to histogram
	}
      }
    }
    // calculate median for first spot
    med = 0;
    numHistVals = 0;
    lt_med = 0;
    while(lt_med <= th){
      lt_med += rgHist[med];
      ++med;
    }
    --med;
    lastMed = med;

    // put the median in the spot we're at
    idxDst = y*w;
    pDst[idxDst] = (unsigned char)med;
    
    // remove pixels from leftmost column
    idx3 = y*wTmp+radiusX;
    for(int ky = 0; ky < hKern; ++ky){
      valTmp = pTmp[idx3 - wKern];
      --(rgHist[valTmp]);
      idx3 += wTmp;
    }

    for(int x=1;  x < w;  ++x){
      ++idxDst;
      // add pixels from the right-hand side of kernel (after moving over one)
      idx3 = y*wTmp+x+radiusX;
      for(int ky = 0; ky < hKern; ++ky){
	valTmp = pTmp[idx3 + wKern];
	++(rgHist[valTmp]);
	idx3 += wTmp;
      }

      // find median from the histogram.
      // count high-to-low if lastMed was high, otherwise low-to-high
      histCount = 0;
      if(lastMed >= 128){
	med = 255;
	while(histCount <= th){
	  histCount += rgHist[med];
	  --med;
	}
	++med;
      }
      else{
	med = 0;
	while(histCount <= th){
	  histCount += rgHist[med];
	  ++med;
	}
	--med;
      }
      
      // put the median value in the destination pixel
      pDst[idxDst] = (unsigned char)med;

      // remove pixels from leftmost column for next time through loop
      if(x < (w-1)){//don't need to remove left edge if going to a new row
	idx3 = y*wTmp+x+radiusX;
	for(int ky = 0; ky < hKern; ++ky){
	  valTmp = pTmp[idx3 - wKern];
	  --(rgHist[valTmp]);
	  idx3 += wTmp;
	} // end for(ky...
      } // end if

    } // end for (x=1; ...

  }// end for(y=0...
  // report progress
  if(NULL != pProg){
    pProg->reportStatus(progStart + h, 0, progMax);
  }
}



///Median filter imgSrc with current settings and return the result
DImage DMedianFilter::filterImage(const DImage &imgSrc,
				  bool fAlreadyPadded, DProgress *pProg){
  DImage dst;
  filterImage_(dst, imgSrc, fAlreadyPadded, pProg);
  return dst;
}

///Median filter imgSrc with current settings and store result in imgDst
void DMedianFilter::filterImage_(DImage &imgDst, const DImage &imgSrc,
				 bool fAlreadyPadded, DProgress *pProg){
  DMedFiltType filtType;
  DImage *pImgPad;
  int wKern, hKern;
  int numKernPxls;
  int wUnpad, hUnpad;
#ifndef D_NOTHREADS
  HUANG_8_THREAD_PARAMS_T rgParms[MAX_MEDFILT_THREADS];
  pthread_t rgThreadID[MAX_MEDFILT_THREADS];
#endif


  filtType = _medFiltType;
  if(DMedFilt_default == filtType){
    switch(imgSrc.getImageType()){
      case DImage::DImage_u8:
      case DImage::DImage_RGB:
	filtType = DMedFilt_Huang_circle;
	break;
      default:
	filtType = DMedFilt_slow;
	break;
    }
  }

  pImgPad = (DImage*)&imgSrc;
  if(!fAlreadyPadded){
    pImgPad = new DImage();
    imgSrc.padEdges_(*pImgPad, _radiusX, _radiusX, _radiusY, _radiusY,
		     DImage::DImagePadReplicate);
  }

  //  printf("_radiusX=%d _radiusY=%d\n",_radiusX,_radiusY);

  wUnpad = pImgPad->width()-(_radiusX*2);
  hUnpad = pImgPad->height()-(_radiusY*2);

  wKern = _radiusX * 2 + 1;
  hKern = _radiusY * 2 + 1;
  if(NULL == rgKern){
    rgKern = (unsigned char*)malloc(sizeof(unsigned char) * wKern * hKern);
    if(!rgKern){
      fprintf(stderr, "DMedianFilter::filterImage_() out of memory\n");
      exit(1);
    }
    rgRightEdge = (int*)malloc(sizeof(int)*hKern);
    if(!rgRightEdge){
      fprintf(stderr, "DMedianFilter::filterImage_() out of memory\n");
      exit(1);
    }
    if(DMedFilt_Huang_circle == filtType){
      fill_circle_kern_offsets(_radiusX, _radiusY, rgKern,
			       rgRightEdge, &numKernPxls);
    }
    else{
      fill_square_kern_offsets(_radiusX, _radiusY, rgKern,
			       rgRightEdge, &numKernPxls);
    }
  }
  
  switch(imgSrc.getImageType()){
    case DImage::DImage_u8:
      {
	imgDst.create(wUnpad, hUnpad, DImage::DImage_u8, 1,
		      imgSrc.getAllocMethod());
#ifndef D_NOTHREADS
	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  rgParms[tnum].pImgDst = &imgDst;
	  rgParms[tnum].pImgSrc = pImgPad;
	  rgParms[tnum].radiusX = _radiusX;
	  rgParms[tnum].radiusY = _radiusY;
	  rgParms[tnum].wKern = wKern;
	  rgParms[tnum].hKern = hKern;
	  rgParms[tnum].rgKern = rgKern;
	  rgParms[tnum].numKernPxls = numKernPxls;
	  rgParms[tnum].rgRightEdge = rgRightEdge;
	  rgParms[tnum].pProg = NULL;
	  rgParms[tnum].progStart = 0;
	  rgParms[tnum].progMax = 1;
	  rgParms[tnum].threadNumber = tnum;
	  rgParms[tnum].numThreads = _numThreads;

	  if(0 != pthread_create(&rgThreadID[tnum], NULL,
				 DMedianFilter::DMedianFilter_Huang8threadWrap,
				 &rgParms[tnum])){
	    fprintf(stderr, "DMedianFilter::filterImage_() failed to spawn "
		    "thread #%d. Exiting.\n", tnum);
	    exit(1);
	  }
	}
#endif
	medFiltHuang_u8(imgDst, *pImgPad, _radiusX, _radiusY,
			wKern, hKern, rgKern, numKernPxls,
			rgRightEdge, pProg, 0, hUnpad+1, 0, _numThreads);
#ifndef D_NOTHREADS
	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  if(pthread_join(rgThreadID[tnum],NULL))
	    fprintf(stderr, "DMedianFilter::filterImage_() failed to join "
		    "thread %d\n", tnum);
	}
#endif
	if(NULL != pProg){
	  pProg->reportStatus(hUnpad+1, 0, hUnpad+1);//report progress complete
	}	  
      }
      break;
    case DImage::DImage_RGB:
      {
	DImage imgR, imgG, imgB;
	DImage imgRDst, imgGDst, imgBDst;

	imgRDst.create(wUnpad, hUnpad, DImage::DImage_u8, 1,
		       imgSrc.getAllocMethod());
	imgGDst.create(wUnpad, hUnpad, DImage::DImage_u8, 1,
		       imgSrc.getAllocMethod());
	imgBDst.create(wUnpad, hUnpad, DImage::DImage_u8, 1,
		       imgSrc.getAllocMethod());

	pImgPad->splitRGB(imgR, imgG, imgB);

#ifndef D_NOTHREADS
	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  rgParms[tnum].pImgDst = &imgRDst;
	  rgParms[tnum].pImgSrc = &imgR;
	  rgParms[tnum].radiusX = _radiusX;
	  rgParms[tnum].radiusY = _radiusY;
	  rgParms[tnum].wKern = wKern;
	  rgParms[tnum].hKern = hKern;
	  rgParms[tnum].rgKern = rgKern;
	  rgParms[tnum].numKernPxls = numKernPxls;
	  rgParms[tnum].rgRightEdge = rgRightEdge;
	  rgParms[tnum].pProg = NULL;
	  rgParms[tnum].progStart = 0;
	  rgParms[tnum].progMax = 1;
	  rgParms[tnum].threadNumber = tnum;
	  rgParms[tnum].numThreads = _numThreads;

	  if(0 != pthread_create(&rgThreadID[tnum], NULL,
				 DMedianFilter::DMedianFilter_Huang8threadWrap,
				 &rgParms[tnum])){
	    fprintf(stderr, "DMedianFilter::filterImage_() failed to spawn "
		    "thread #%d. Exiting.\n",tnum);
	    exit(1);
	  }
	}
#endif
	medFiltHuang_u8(imgRDst, imgR, _radiusX, _radiusY,
			wKern, hKern, rgKern, numKernPxls,
			rgRightEdge, pProg, 0, 3 * hUnpad);
#ifndef D_NOTHREADS
	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  if(pthread_join(rgThreadID[tnum],NULL))
	    fprintf(stderr, "DMedianFilter::filterImage_() failed to join "
		    "thread %d\n", tnum);
	}
	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  rgParms[tnum].pImgDst = &imgGDst;
	  rgParms[tnum].pImgSrc = &imgG;
	  if(0 != pthread_create(&rgThreadID[tnum], NULL,
				 DMedianFilter::DMedianFilter_Huang8threadWrap,
				 &rgParms[tnum])){
	    fprintf(stderr, "DMedianFilter::filterImage_() failed to spawn "
		    "thread #%d. Exiting.\n",tnum);
	    exit(1);
	  }
	}
#endif
	medFiltHuang_u8(imgGDst, imgG, _radiusX, _radiusY,
			wKern, hKern, rgKern, numKernPxls,
			rgRightEdge, pProg, hUnpad, 3 * hUnpad);
#ifndef D_NOTHREADS
	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  if(pthread_join(rgThreadID[tnum],NULL))
	    fprintf(stderr, "DMedianFilter::filterImage_() failed to join "
		    "thread %d\n", tnum);
	}

	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  rgParms[tnum].pImgDst = &imgBDst;
	  rgParms[tnum].pImgSrc = &imgB;
	  if(0 != pthread_create(&rgThreadID[tnum], NULL,
				 DMedianFilter::DMedianFilter_Huang8threadWrap,
				 &rgParms[tnum])){
	    fprintf(stderr, "DMedianFilter::filterImage_() failed to spawn "
		    "thread #%d. Exiting.\n",tnum);
	    exit(1);
	  }
	}
#endif
	medFiltHuang_u8(imgBDst, imgB, _radiusX, _radiusY,
			wKern, hKern, rgKern, numKernPxls,
			rgRightEdge, pProg, 2 * hUnpad, 3 * hUnpad+1);
#ifndef D_NOTHREADS
	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  if(pthread_join(rgThreadID[tnum],NULL))
	    fprintf(stderr, "DMedianFilter::filterImage_() failed to join "
		    "thread %d\n", tnum);
	}
#endif
	if(NULL != pProg){
	  pProg->reportStatus(3*hUnpad+1,0,3*hUnpad+1);//report complete
	}	  
	imgDst.combineRGB(imgRDst, imgGDst, imgBDst);
      }
      break;
    case DImage::DImage_dbl_multi:
      {
	int w, h, wm1, hm1; // width, height of padded image
	double *pDst;
	double *pPad;
	double *rgWindowBuff;
	
	fprintf(stderr, "DMedianFilter::filterImage_() performing brute-force "
		"(slow) median filter on double image...\n");
	w = pImgPad->width();
	h = pImgPad->height();
	wm1=w-1;
	hm1 = h-1;
	imgDst.create(wUnpad, hUnpad, DImage::DImage_dbl_multi,
		      imgSrc.numChannels(), imgSrc.getAllocMethod());
	rgWindowBuff = (double*)malloc(sizeof(double)*wKern*hKern);
	D_CHECKPTR(rgWindowBuff);
	for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
	  pDst = imgDst.dataPointer_dbl(chan);
	  pPad = pImgPad->dataPointer_dbl(chan);
	  for(int y = 0, idxDst = 0; y < hUnpad; ++y){
	    int idxPad;
	    idxPad = (y+_radiusY)*w+_radiusX;
	    for(int x=0; x < wUnpad; ++x, ++idxDst, ++idxPad){
	      int count;
	      count = 0;
	      for(int dy = -_radiusY; dy <= _radiusY; ++dy){
		for(int dx = -_radiusX; dx <= _radiusX; ++dx){
		  rgWindowBuff[count] = pPad[idxPad+(dy*w)+dx];
		  ++count;
		}
	      }
	      // find median
	      qsort((void*)rgWindowBuff, count, sizeof(double),compareDoubles);
	      pDst[idxDst] = rgWindowBuff[count / 2];
	    }//end for(x...
	  }//end for(y...
	}//end for(chan...
	free(rgWindowBuff);
      }//end block in case DImage_dbl_multi
      break;
    default:
      //TODO: finish implementing this
      fprintf(stderr, "DMedianFilter::filterImage_() not yet implemented for "
	      "some image types\n");
      exit(1);
      break;
  }
  if(!fAlreadyPadded){
    delete pImgPad;
    pImgPad = NULL;
  }

}




///Fills rgKernFull with a circle and sets rgRightEdge and numKernPxls
/**Pixels that are part of the circle will be 0xff and the rest will
 * be 0x00.  rgRightEdge is an array that says how many pixels to the
 * right the leading edge of the kernel is (the right-most pixel), and
 * it is assumed that the trailing edge is the same number of pixels
 * to the left. numKernPxls is a count of how many pixels in the
 * kernel are set to 0xff.
 */
void DMedianFilter::fill_circle_kern_offsets(int radiusX, int radiusY,
					     unsigned char *rgKernFull,
					     int *rgRightEdge,
					     int *numKernPxls){
  double rx2,ry2;

  (*numKernPxls)= 0;

  rx2=(0.5+radiusX)*(0.5+radiusX);
  ry2=(0.5+radiusY)*(0.5+radiusY);
  for(int ky = -radiusY, idx = 0; ky <= radiusY; ++ky){
    for(int kx = -radiusX; kx <= radiusX; ++kx, ++idx){
      if( ((kx*kx/rx2)+(ky*ky/ry2)) <= 1.){
	rgKernFull[idx] = 0xff;
	rgRightEdge[ky+radiusY] = kx;
	++(*numKernPxls);
      }
      else
	rgKernFull[idx] = 0;
    }
  }

#if 0
  for(int ky = -radiusY, idx = 0; ky <= radiusY; ++ky){
    printf("%02d:  ", ky+radiusY);
    for(int kx = -radiusX; kx <= radiusX; ++kx, ++idx){
      printf("%c", rgKernFull[idx] ? 'X' : '.');
    }
    printf("  %d\n", rgRightEdge[ky+radiusY]);
  }
#endif
}

///Fills rgKernFull with a square and sets rgRightEdge and numKernPxls
/**Similar to fill_circle_kern_offsets() but all of the pixels are set
 * to 0xff, in effect, forming a square kernel that is (radius*2+1)
 * wide and high.
 */
void DMedianFilter::fill_square_kern_offsets(int radiusX, int radiusY,
					     unsigned char *rgKernFull,
					     int *rgRightEdge,
					     int *numKernPxls){
  int kernWidth, kernHeight;
  kernWidth = radiusX*2+1;
  kernHeight = radiusY*2+1;
  
  memset(rgKernFull, 1, kernWidth*kernHeight);
  for(int ky = 0; ky < kernHeight; ++ky){
    rgRightEdge[ky] = radiusX;
  }
  (*numKernPxls) = kernWidth*kernHeight;
}

///Static function provided for convenience
void DMedianFilter::medianFilterImage(DImage &imgDst, const DImage &imgSrc,
				      bool fAlreadyPadded,
				      int radiusX, int radiusY,
				      DMedFiltType filtType, DProgress *pProg,
				      int numThreads){
  DMedianFilter mfilt(radiusX, radiusY, filtType);
  mfilt.setNumThreads(numThreads);
  mfilt.filterImage_(imgDst, imgSrc, fAlreadyPadded, pProg);
}

///this is part-implemented only for some debug/comparison work.  do not use it
// void DMedianFilter::medianFilter_separable_u8(DImage &imgDst,
// 					      const DImage &imgSrc,
// 					      bool fAlreadyPadded,
// 					      int radius){
//   DImage imgPad;
//   DImage *pImgPad;
//   DImage imgMed; // temporary image to store horizontal medians
//   unsigned int rgHist[256];
//   D_uint8 *p8pad, *p8med;
//   bool fSearchFromZero = true;
//   int w, h;

//   //pad the image first
//   if(fAlreadyPadded)
//     pImgPad = (DImage*)&imgSrc;
//   else{
//     pImgPad = &imgPad;
//     imgSrc.padEdges_(*pImgPad, radius, radius, radius, radius,
// 		     DImage::DImagePadReplicate);
//   }
  
//   printf("doing horizontal medians...\n");
//   // now take one-dimensional medians in horizontal direction
//   w = pImgPad->width();
//   h = pImgPad->height();
//   imgMed.create(w,h,DImage::DImage_u8);
//   p8med = imgMed.dataPointer_u8();
//   p8pad = pImgPad->dataPointer_u8();
//   for(int y=0; y < h; ++y){
//     D_uint8 *p8;
//     D_uint8 *p8_dst;
//     int midx;//median idx in rgHist;
//     //initialize hist for first real pixel
//     p8 = &(p8pad[y*w]);
//     p8_dst = &(p8med[y*w]);
//     memset(rgHist,0,256*sizeof(unsigned int));
//     for(int x=0; x < 2*radius+1; ++x)
//       ++(rgHist[p8[x]]);
//     for(int x=radius; x < (w-radius); ++x){
//       //find median in histogram
//       int hcount;
//       hcount = 0;
//       if(fSearchFromZero){//start searching from 0
// 	for(midx=0; midx < 256; ++midx){
// 	  hcount += rgHist[midx];
// 	  if(hcount>=radius)
// 	    break;
// 	}
//       }
//       else{
// 	for(midx=255; midx >=0; --midx){
// 	  hcount += rgHist[midx];
// 	  if(hcount>=radius)
// 	    break;
// 	}
//       }
//       if(hcount < radius){//sanity check
// 	fprintf(stderr,"error!(%s:%d)\n",__FILE__,__LINE__);
// 	exit(1);
//       }
//       //assign median value
//       p8_dst[x] = (D_uint8)midx;
//       //update the histogram for next pixel (drop left edge, add new right edge)
//       if((x+1)<(w-radius)){
// 	--(rgHist[p8[x-radius]]);
// 	++(rgHist[p8[x+1+radius]]);
//       }
//       fSearchFromZero = (midx < 128) ? true : false;
//     }
//   }
//   printf("done with horizontal medians\n");
//   printf("the rest of this function isn't written yet\n");
//   exit(1);
// }






///this is part-implemented only for some debug/comparison work.  do not use it
/** This is actually a brute-force approach.  The commented code is
    the partially-implemented version of the separable approximation
    (median of 1-D medians).**/
void DMedianFilter::medianFilter_separable_u8(DImage &imgDst,
					      const DImage &imgSrc,
					      bool fAlreadyPadded,
					      int radius){
  DImage imgPad;
  DImage *pImgPad;
  DImage imgMed; // temporary image to store horizontal medians
  unsigned int rgHist[256];
  D_uint8 *p8pad, *p8med;
  bool fSearchFromZero = true;
  int w, h;
  //pad the image first
  if(fAlreadyPadded)
    pImgPad = (DImage*)&imgSrc;
  else{
    pImgPad = &imgPad;
    imgSrc.padEdges_(*pImgPad, radius, radius, radius, radius,
		     DImage::DImagePadReplicate);
  }
  
  printf("doing brute-force median...\n");
  // now take one-dimensional medians in horizontal direction
  w = pImgPad->width();
  h = pImgPad->height();
  imgMed.create(w,h,DImage::DImage_u8);
  p8med = imgMed.dataPointer_u8();
  p8pad = pImgPad->dataPointer_u8();

#if 1
  D_uint8 *rgSortArray;
  rgSortArray=new D_uint8[radius];
#endif

  D_uint8 *p8;
  D_uint8 *p8_dst;

  p8 = p8pad;
  p8_dst = p8med;

  for(int y=radius; y < (h-radius-1); ++y){
    int midx;//median idx in rgHist;
    //initialize hist for first real pixel
    // p8 = &(p8pad[y*w]);
    // p8_dst = &(p8med[y*w]);
    memset(rgHist,0,256*sizeof(unsigned int));
    // for(int x=0; x < 2*radius+1; ++x)
    //   ++(rgHist[p8[x]]);
    for(int x=radius; x < (w-radius-1); ++x){

#if 0
      int countp;
      countp = 0;
      for(int yp=y-radius; yp <=(y+radius); ++yp){
	for(int xp=x-radius; xp <=(x+radius); ++xp){
	  // if((yp<0)||(xp<0)||(yp>=hm1)||(xp>=wm1))
	  //   continue;
	  rgSortArray[countp] = p8[yp*w+xp];
	  ++countp;
	}
      }
      
      qsort((void*)rgSortArray, countp,
	    sizeof(D_uint8),med_compareChars);
      p8_dst[y*w+x] = (int)(unsigned int)rgSortArray[countp/2];
#else
      //do historgram median for each pixel
      int countp;
      countp = 0;
      memset(rgHist,0,sizeof(unsigned int)*256);
      for(int yp=y-radius; yp <=(y+radius); ++yp){
	for(int xp=x-radius; xp <=(x+radius); ++xp){
	  if((yp<0)||(xp<0)||(yp>=h)||(xp>=w))
	  	      continue;
	  ++(rgHist[p8[yp*w+xp]]);
	  ++countp;
	}
      }
      //find median
      int hcount;
      int hidx = 0;
      hcount = 0;
      if(fSearchFromZero){//start searching from 0
	for(hidx=0; hidx<256; ++hidx){
	  hcount += rgHist[hidx];
	  if(hcount >= countp/2)
	    break;
	}
      }
      else{
	for(hidx=255; hidx>=0; --hidx){
	  hcount += rgHist[hidx];
	  if(hcount >= countp/2)
	    break;
	}
      }
      if(hcount == 0)
	fprintf(stderr,"error!\n");
      p8_dst[y*w+x] = (D_uint8)hidx;

      fSearchFromZero = (hidx < 128) ? true : false;

#endif

    }
  }
  printf("done with brute-force median\n");
  printf("the rest of this function isn't written yet\n");
  imgDst = imgMed;
#if 1
  delete [] rgSortArray;
#endif
  //  exit(1);
}



