#include "dminfilter.h"
#include "dimage.h"
#include "dprogress.h"
#include "dinstancecounter.h"
#include <string.h>

#ifndef D_NOTHREADS
#include "dthreads.h"

/*comparison function for qsort()ing doubles in nondecreasing order */
// int compareDoubles(const void *p1, const void *p2){
//   if ( (*(double*)(p1)) < (*(double*)(p2)) )
//     return -1;
//   else if ( (*(double*)(p1)) > (*(double*)(p2)) )
//     return 1;
//   return 0;
// }

// /*comparison function for qsort()ing floats in nondecreasing order */
// int compareFloats(const void *p1, const void *p2){
//   if ( (*(float*)(p1)) < (*(float*)(p2)) )
//     return -1;
//   else if ( (*(float*)(p1)) > (*(float*)(p2)) )
//     return 1;
//   return 0;
// }


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
} MIN_HUANG_8_THREAD_PARAMS_T;
/// \endcond

#define MAX_MINFILT_THREADS 64

void* DMinFilter::DMinFilter_Huang8threadWrap(void* params){
  MIN_HUANG_8_THREAD_PARAMS_T *pParams;

  pParams = (MIN_HUANG_8_THREAD_PARAMS_T*)params;
  DMinFilter::minFiltHuang_u8(*(pParams->pImgDst),
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
DMinFilter::DMinFilter(){
  DInstanceCounter::addInstance("DMinFilter");
  _minFiltType = DMinFilt_circle;
  rgKern = NULL;
  rgRightEdge = NULL;
  _radiusX = 1;
  _radiusY = 1;
  _numThreads = 1;
}

///Constructor with parameters
DMinFilter::DMinFilter(int radiusX, int radiusY, DMinFiltType filtType){
  DInstanceCounter::addInstance("DMinFilter");
  _minFiltType = filtType;
  rgKern = NULL;
  rgRightEdge = NULL;
  _radiusX = radiusX;
  _radiusY = radiusY;
  _numThreads = 1;
}


///Destructor
DMinFilter::~DMinFilter(){
  DInstanceCounter::removeInstance("DMinFilter");
  if(rgKern){
    free(rgKern);
    rgKern = NULL;
  }
  if(rgRightEdge){
    free(rgRightEdge);
    rgRightEdge = NULL;
  }
}

///Set the kernel (window) radii for the min filter
/**The actual width (and height) of the window will be 2*radius+1*/
void DMinFilter::setRadius(int radiusX, int radiusY){
  if((radiusX < 0)||(radiusY < 0)){
    fprintf(stderr, "DMinFilter::setRadius() negative value not allowed\n");
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
int DMinFilter::getRadiusX(){
  return _radiusX;
}
///Get the vertical radius of the kernel
int DMinFilter::getRadiusY(){
  return _radiusY;
}

///Set the number of threads that will be used to do the min filter
/** The min filter may run faster if more than one thread is used,
 *  especially if the machine has more than one processor.  With
 *  numThreads threads specified, thread 0 will start with y=0, thread
 *  1 with y=1, and so on, and each will skip numThreads rows between
 *  each row that it operates on.  Only thread 0 will report status.
 *  Note: our experiments show that using multiple threads does not
 *  seem to speed up the min filter very much.  For example, on a
 *  very large image with radius=35, we saw a speedup from 33.6
 *  seconds to 29.1 seconds with 2 threads, and to 28.7 sec with 4
 *  threads.  When the code is not optimized, there is a more
 *  noticible speedup, but we assume you will optimize your code.
 */
void DMinFilter::setNumThreads(int numThreads){
#ifndef D_NOTHREADS
  if((_numThreads > 0) && (_numThreads <= MAX_MINFILT_THREADS)){
    _numThreads = numThreads;
  }
  else{
    fprintf(stderr, "DMinFilter::setNumThreads() numThreads must be from "
	    "1 to %d\n", MAX_MINFILT_THREADS);
  }
#else
  fprintf(stderr, "DMinFilter::setNumThreads() thread support not "
	  "compiled");
#endif
}

///return the number of threads that the filter has been set to use
int DMinFilter::getNumThreads(){
  return _numThreads;
}

///set the type of min filter algorithm
/**By default, RGB and _u8 images will use Huang_circle and all others
 * (including RGB_16 and _u16) will use the slow method */
void DMinFilter::setType(DMinFiltType filtType){
  if(rgKern){
    free(rgKern);
    rgKern = NULL;
    free(rgRightEdge);
    rgRightEdge = NULL;
  }
  _minFiltType = filtType;
}

///return the type of the min filter
DMinFilter::DMinFiltType DMinFilter::getType(){
  return _minFiltType;
}

///private function that actually performs huang min filt on 8-bit images
/** imgDst will be 2*radius pixels less wide and high than imgSrc
 * because of the padding that is added before calling this function.
 * This function requires that imgDst.create() has already been called
 * with the proper w,h,imgType,etc.
 */
void DMinFilter::minFiltHuang_u8(DImage &imgDst, const DImage &imgSrc,
				    int radiusX, int radiusY,
				    int wKern, int hKern, D_uint8 *rgKern,
				    int numKernPxls, int *rgRightEdge,
				    DProgress *pProg,
				    int progStart, int progMax,
				    int threadNumber, int numThreads){
  int rgHist[256];
  int min;
  int lastMin;
  unsigned char valTmp;
  int idxDst;
  int idx3;
  D_uint8 *pTmp; // pointer to padded image data
  int wTmp, hTmp; // width, height of imgSrc
  int w, h; // width, height of imgDst
  D_uint8 *pDst;

  wTmp = imgSrc.width();
  hTmp = imgSrc.height();
  w = wTmp - radiusX*2;
  h = hTmp - radiusY*2;
  pDst = imgDst.dataPointer_u8();
  pTmp = imgSrc.dataPointer_u8();

  for(int y = threadNumber; y < h; y += numThreads){
    // update progress report and check if user cancelled the operation
    if((NULL != pProg) && (0 == (y & 0x0000003f))){
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
    // calculate min for first spot
    // for(min = 255; (min > 0) && (0==rgHist[min]); --min){
    for(min = 0; (min < 255) && (0==rgHist[min]); ++min){
      // do nothing
    }
    lastMin = min;

    // put the min in the spot we're at
    idxDst = y*w;
    pDst[idxDst] = (unsigned char)min;
    
    // remove pixels from leftmost column
    idx3 = y*wTmp+radiusX;
    for(int ky = 0; ky < hKern; ++ky){
      valTmp = pTmp[idx3 - rgRightEdge[ky]];
      --(rgHist[valTmp]);
      if((valTmp==min)&&(0 == rgHist[valTmp])){//update the min
	for(;(min<255)&&(0==rgHist[min]); ++min){
	  //do nothing
	}
      }
      idx3 += wTmp;
    }

    for(int x=1;  x < w;  ++x){
      ++idxDst;
      // add pixels from the right-hand side of kernel (after moving over one)
      idx3 = y*wTmp+x+radiusX;
      for(int ky = 0; ky < hKern; ++ky){
	valTmp = pTmp[idx3 + rgRightEdge[ky]];
	if(valTmp < min)//update the min
	  min = valTmp;
	++(rgHist[valTmp]);
	idx3 += wTmp;
      }

      // find min from the histogram for DEBUG purposes (remove this later)
      for(lastMin = 0; (lastMin < 255) && (0==rgHist[lastMin]); ++lastMin){
	// do nothing
      }
      if(lastMin != min){
	fprintf(stderr, "in dminfilter.cpp, lastMin!=min!\n");
	exit(1);
      }
      
      // put the min value in the destination pixel
      pDst[idxDst] = (unsigned char)min;

      // remove pixels from leftmost column for next time through loop
      if(x < (w-1)){//don't need to remove left edge if going to a new row
	idx3 = y*wTmp+x+radiusX;
	for(int ky = 0; ky < hKern; ++ky){
	  valTmp = pTmp[idx3 - rgRightEdge[ky]];
	  --(rgHist[valTmp]);
	  if((valTmp==min)&&(0 == rgHist[valTmp])){//update the min
	    for(;(min<255)&&(0==rgHist[min]); ++min){
	      //do nothing
	    }
	  }
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


///private function that actually performs huang min filt on 8-bit images
/** imgDst will be 2*radius pixels less wide and high than imgSrc
 * because of the padding that is added before calling this function.
 * This function requires that imgDst.create() has already been called
 * with the proper w,h,imgType,etc.
 */
void DMinFilter::minFiltHuang_u8_square(DImage &imgDst,
					   const DImage &imgSrc,
					   int radiusX, int radiusY,
					   int wKern, int hKern,
					   D_uint8 *rgKern,
					   int numKernPxls,
					   DProgress *pProg,
					   int progStart, int progMax,
					   int threadNumber, int numThreads){
  int rgHist[256];
  int min;
  unsigned char valTmp;
  int idxDst;
  int idx3;
  D_uint8 *pTmp; // pointer to padded image data
  int wTmp, hTmp; // width, height of imgSrc
  int w, h; // width, height of imgDst
  D_uint8 *pDst;

  wTmp = imgSrc.width();
  hTmp = imgSrc.height();
  w = wTmp - radiusX*2;
  h = hTmp - radiusY*2;
  pDst = imgDst.dataPointer_u8();
  pTmp = imgSrc.dataPointer_u8();

  for(int y = threadNumber; y < h; y += numThreads){
    // update progress report and check if user cancelled the operation
    if((NULL != pProg) && (0 == (y & 0x0000003f))){
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
    // calculate min for first spot
    // for(min = 255; (min > 0) && (0==rgHist[min]); --min){
    for(min = 0; (min < 255) && (0==rgHist[min]); ++min){
      // do nothing
    }

    // put the min in the spot we're at
    idxDst = y*w;
    pDst[idxDst] = (unsigned char)min;
    
    // remove pixels from leftmost column
    idx3 = y*wTmp+radiusX;
    for(int ky = 0; ky < hKern; ++ky){
      valTmp = pTmp[idx3 - wKern];
      --(rgHist[valTmp]);
      if((valTmp==min)&&(0 == rgHist[valTmp])){//update the min
	for(;(min<255)&&(0==rgHist[min]); ++min){
	  //do nothing
	}
      }
      idx3 += wTmp;
    }

    for(int x=1;  x < w;  ++x){
      ++idxDst;
      // add pixels from the right-hand side of kernel (after moving over one)
      idx3 = y*wTmp+x+radiusX;
      for(int ky = 0; ky < hKern; ++ky){
	valTmp = pTmp[idx3 + wKern];
	if(valTmp < min)//update the min
	  min = valTmp;
	++(rgHist[valTmp]);
	idx3 += wTmp;
      }
      
      // put the min value in the destination pixel
      pDst[idxDst] = (unsigned char)min;

      // remove pixels from leftmost column for next time through loop
      if(x < (w-1)){//don't need to remove left edge if going to a new row
	idx3 = y*wTmp+x+radiusX;
	for(int ky = 0; ky < hKern; ++ky){
	  valTmp = pTmp[idx3 - wKern];
	  --(rgHist[valTmp]);
	  if((valTmp==min)&&(0 == rgHist[valTmp])){//update the min
	    for(;(min<255)&&(0==rgHist[min]); ++min){
	      //do nothing
	    }
	  }
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



///Min filter imgSrc with current settings and return the result
DImage DMinFilter::filterImage(const DImage &imgSrc,
				  bool fAlreadyPadded, DProgress *pProg){
  DImage dst;
  filterImage_(dst, imgSrc, fAlreadyPadded, pProg);
  return dst;
}

///Min filter imgSrc with current settings and store result in imgDst
void DMinFilter::filterImage_(DImage &imgDst, const DImage &imgSrc,
				 bool fAlreadyPadded, DProgress *pProg){
  DMinFiltType filtType;
  DImage *pImgPad;
  int wKern, hKern;
  int numKernPxls;
  int wUnpad, hUnpad;
#ifndef D_NOTHREADS
  MIN_HUANG_8_THREAD_PARAMS_T rgParms[MAX_MINFILT_THREADS];
  pthread_t rgThreadID[MAX_MINFILT_THREADS];
#endif


  filtType = _minFiltType;

  pImgPad = (DImage*)&imgSrc;
  if(!fAlreadyPadded){
    pImgPad = new DImage();
    imgSrc.padEdges_(*pImgPad, _radiusX, _radiusX, _radiusY, _radiusY,
		     DImage::DImagePadReplicate);
  }

  wUnpad = pImgPad->width()-(_radiusX*2);
  hUnpad = pImgPad->height()-(_radiusY*2);

  wKern = _radiusX * 2 + 1;
  hKern = _radiusY * 2 + 1;
  if(NULL == rgKern){
    rgKern = (unsigned char*)malloc(sizeof(unsigned char) * wKern * hKern);
    if(!rgKern){
      fprintf(stderr, "DMinFilter::filterImage_() out of memory\n");
      exit(1);
    }
    rgRightEdge = (int*)malloc(sizeof(int)*hKern);
    if(!rgRightEdge){
      fprintf(stderr, "DMinFilter::filterImage_() out of memory\n");
      exit(1);
    }
    if(DMinFilt_circle == filtType){
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
				 DMinFilter::DMinFilter_Huang8threadWrap,
				 &rgParms[tnum])){
	    fprintf(stderr, "DMinFilter::filterImage_() failed to spawn "
		    "thread #%d. Exiting.\n", tnum);
	    exit(1);
	  }
	}
#endif
	minFiltHuang_u8(imgDst, *pImgPad, _radiusX, _radiusY,
			wKern, hKern, rgKern, numKernPxls,
			rgRightEdge, pProg, 0, hUnpad+1, 0, _numThreads);
#ifndef D_NOTHREADS
	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  if(pthread_join(rgThreadID[tnum],NULL))
	    fprintf(stderr, "DMinFilter::filterImage_() failed to join "
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
				 DMinFilter::DMinFilter_Huang8threadWrap,
				 &rgParms[tnum])){
	    fprintf(stderr, "DMinFilter::filterImage_() failed to spawn "
		    "thread #%d. Exiting.\n",tnum);
	    exit(1);
	  }
	}
#endif
	minFiltHuang_u8(imgRDst, imgR, _radiusX, _radiusY,
			wKern, hKern, rgKern, numKernPxls,
			rgRightEdge, pProg, 0, 3 * hUnpad);
#ifndef D_NOTHREADS
	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  if(pthread_join(rgThreadID[tnum],NULL))
	    fprintf(stderr, "DMinFilter::filterImage_() failed to join "
		    "thread %d\n", tnum);
	}
	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  rgParms[tnum].pImgDst = &imgGDst;
	  rgParms[tnum].pImgSrc = &imgG;
	  if(0 != pthread_create(&rgThreadID[tnum], NULL,
				 DMinFilter::DMinFilter_Huang8threadWrap,
				 &rgParms[tnum])){
	    fprintf(stderr, "DMinFilter::filterImage_() failed to spawn "
		    "thread #%d. Exiting.\n",tnum);
	    exit(1);
	  }
	}
#endif
	minFiltHuang_u8(imgGDst, imgG, _radiusX, _radiusY,
			wKern, hKern, rgKern, numKernPxls,
			rgRightEdge, pProg, hUnpad, 3 * hUnpad);
#ifndef D_NOTHREADS
	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  if(pthread_join(rgThreadID[tnum],NULL))
	    fprintf(stderr, "DMinFilter::filterImage_() failed to join "
		    "thread %d\n", tnum);
	}

	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  rgParms[tnum].pImgDst = &imgBDst;
	  rgParms[tnum].pImgSrc = &imgB;
	  if(0 != pthread_create(&rgThreadID[tnum], NULL,
				 DMinFilter::DMinFilter_Huang8threadWrap,
				 &rgParms[tnum])){
	    fprintf(stderr, "DMinFilter::filterImage_() failed to spawn "
		    "thread #%d. Exiting.\n",tnum);
	    exit(1);
	  }
	}
#endif
	minFiltHuang_u8(imgBDst, imgB, _radiusX, _radiusY,
			wKern, hKern, rgKern, numKernPxls,
			rgRightEdge, pProg, 2 * hUnpad, 3 * hUnpad+1);
#ifndef D_NOTHREADS
	for(int tnum = 1; tnum < _numThreads; ++tnum){
	  if(pthread_join(rgThreadID[tnum],NULL))
	    fprintf(stderr, "DMinFilter::filterImage_() failed to join "
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
	
	fprintf(stderr, "DMinFilter::filterImage_() performing brute-force "
		"(slow) min filter on double image... NOT YET IMPLEMENTED!\n");
	exit(1);
	// w = pImgPad->width();
	// h = pImgPad->height();
	// wm1=w-1;
	// hm1 = h-1;
	// imgDst.create(wUnpad, hUnpad, DImage::DImage_dbl_multi,
	// 	      imgSrc.numChannels(), imgSrc.getAllocMethod());
	// rgWindowBuff = (double*)malloc(sizeof(double)*wKern*hKern);
	// D_CHECKPTR(rgWindowBuff);
	// for(int chan = 0; chan < imgSrc.numChannels(); ++chan){
	//   pDst = imgDst.dataPointer_dbl(chan);
	//   pPad = pImgPad->dataPointer_dbl(chan);
	//   for(int y = 0, idxDst = 0; y < hUnpad; ++y){
	//     int idxPad;
	//     idxPad = (y+_radiusY)*w+_radiusX;
	//     for(int x=0; x < wUnpad; ++x, ++idxDst, ++idxPad){
	//       int count;
	//       count = 0;
	//       for(int dy = -_radiusY; dy <= _radiusY; ++dy){
	// 	for(int dx = -_radiusX; dx <= _radiusX; ++dx){
	// 	  rgWindowBuff[count] = pPad[idxPad+(dy*w)+dx];
	// 	  ++count;
	// 	}
	//       }
	//       // find min
	//       qsort((void*)rgWindowBuff, count, sizeof(double),compareDoubles);
	//       pDst[idxDst] = rgWindowBuff[count / 2];
	//     }//end for(x...
	//   }//end for(y...
	// }//end for(chan...
	// free(rgWindowBuff);
      }//end block in case DImage_dbl_multi
      break;
    default:
      //TODO: finish implementing this
      fprintf(stderr, "DMinFilter::filterImage_() not yet implemented for "
	      "some image types\n");
      exit(1);
      break;
  }
  if(!fAlreadyPadded){
    delete pImgPad;
    pImgPad = NULL;
  }

}



#if 0
// ///Fills rgKernFull with a circle and sets rgRightEdge and numKernPxls
// /**Pixels that are part of the circle will be 0xff and the rest will
//  * be 0x00.  rgRightEdge is an array that says how many pixels to the
//  * right the leading edge of the kernel is (the right-most pixel), and
//  * it is assumed that the trailing edge is the same number of pixels
//  * to the left. numKernPxls is a count of how many pixels in the
//  * kernel are set to 0xff.
//  */
// void DMinFilter::fill_circle_kern_offsets(int radiusX, int radiusY,
// 					     unsigned char *rgKernFull,
// 					     int *rgRightEdge,
// 					     int *numKernPxls){
//   double rx2,ry2;

//   (*numKernPxls)= 0;

//   rx2=(0.5+radiusX)*(0.5+radiusX);
//   ry2=(0.5+radiusY)*(0.5+radiusY);
//   for(int ky = -radiusY, idx = 0; ky <= radiusY; ++ky){
//     for(int kx = -radiusX; kx <= radiusX; ++kx, ++idx){
//       if( ((kx*kx/rx2)+(ky*ky/ry2)) <= 1.){
// 	rgKernFull[idx] = 0xff;
// 	rgRightEdge[ky+radiusY] = kx;
// 	++(*numKernPxls);
//       }
//       else
// 	rgKernFull[idx] = 0;
//     }
//   }

// #if 0
//   for(int ky = -radiusY, idx = 0; ky <= radiusY; ++ky){
//     printf("%02d:  ", ky+radiusY);
//     for(int kx = -radiusX; kx <= radiusX; ++kx, ++idx){
//       printf("%c", rgKernFull[idx] ? 'X' : '.');
//     }
//     printf("  %d\n", rgRightEdge[ky+radiusY]);
//   }
// #endif
// }
#endif

///Fills rgKernFull with a circle and sets rgRightEdge and numKernPxls
/**Pixels that are part of the circle will be 0xff and the rest will
 * be 0x00.  rgRightEdge is an array that says how many pixels to the
 * right the leading edge of the kernel is (the right-most pixel), and
 * it is assumed that the trailing edge is the same number of pixels
 * to the left. numKernPxls is a count of how many pixels in the
 * kernel are set to 0xff.
 */
void DMinFilter::fill_circle_kern_offsets(int radiusX, int radiusY,
					     unsigned char *rgKernFull,
					     int *rgRightEdge,
					     int *numKernPxls){
  double rx2,ry2;

  int r2;


  if(radiusX != radiusY){
    fprintf(stderr, "DMinFilter::fill_circle_kern_offsets() radiusX should equal radiusY for now.\n");
    exit(1);
  }

  r2 = radiusX*radiusX;

  (*numKernPxls)= 0;

  for(int ky = -radiusY, idx = 0; ky <= radiusY; ++ky){
    for(int kx = -radiusX; kx <= radiusX; ++kx, ++idx){
      if( ((kx*kx)+(ky*ky)) <= r2){
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
void DMinFilter::fill_square_kern_offsets(int radiusX, int radiusY,
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
void DMinFilter::minFilterImage(DImage &imgDst, const DImage &imgSrc,
				      bool fAlreadyPadded,
				      int radiusX, int radiusY,
				      DMinFiltType filtType, DProgress *pProg,
				      int numThreads){
  DMinFilter mfilt(radiusX, radiusY, filtType);
  mfilt.setNumThreads(numThreads);
  mfilt.filterImage_(imgDst, imgSrc, fAlreadyPadded, pProg);
}
