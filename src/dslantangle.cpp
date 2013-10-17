#include "dslantangle.h"
#include "dmath.h"
#include "dprofile.h"
#include "dthreads.h"
#include <string.h>

///Calculate the slant angle for this black and white textline image
/**The slant angle is assumed to be between 60 and -45 degrees (0 deg=vertical,
 * negative values are left-slanted, positive values right-slanted).
 * To determine slant: at each x-position, the longest runlength at each angle
 * is found and its squared value is added into the accumulator for that angle.
 * The histogram is smoothed, and the angle corresponding to the highest value 
 * in the histogram is the returned angle (in degrees).
 *
 * Runlengths of less than rlThresh pixels are ignored.
 *
 * The image should be black(0) and white(255).  The portion of the image
 * specified by x0,y0 - x1,y1 is considered to be the textline of interest.
 * If no coordinates are specified, then the entire image is used as the 
 * textline.
 *
 * If weight is not NULL, it will be the sum of max runlengths (not squared) at
 * all 120 angles.  Weights are used in determination of weighted average angle
 * for all textlines in getAllTextlinesSlantAngleDeg() before adjusting angles.
 *
 * If rgSlantHist is not NULL, the squared max RL values in the angle histogram
 * will be copied into the rgSlantHist array. It must already be allocated to
 * 120*sizeof(unsigned int).
 * 
 * if imgAngleHist is not NULL, then the image is set to w=120 and h=y1-y0+1.
 * It is a (grayscale) graphical representation of what is in rgSlantHist.
 */
double DSlantAngle::getTextlineSlantAngleDeg(DImage &imgBW,
					     int rlThresh,
					     int x0,int y0,int x1,int y1,
					     double *weight,
					     unsigned int *rgSlantHist,
					     DImage *imgAngleHist){
  int *rgLineSlantAngles;
  int lineH;
  int slantOffset, slantAngle, angle;
  unsigned int rgSlantSums[120];
  unsigned int rgSlantSumsTmp[120];
  int runlen, maxrl; /* maximum slant runlen */
  double slantDx;
  int w, h;
  D_uint8 *p8;
  double dblWeight = 0;
  
  w = imgBW.width();
  h = imgBW.height();
  p8 = imgBW.dataPointer_u8();
  if(-1 == x1)
    x1 = w-1;
  if(-1 == y1)
    y1 = h-1;
  
  lineH = y1-y0+1;

  /* estimate the predominant slant angle (0=vertical, +right, -left) */
  slantOffset = (int)(0.5+ (lineH / 2.0) / tan(DMath::degreesToRadians(30.)));
  for(int j = 0; j < 120; ++j){
    rgSlantSums[j] = 0;
    rgSlantSumsTmp[j] = 0;
  }
  for(angle = -45; angle <= 60; angle += 1){
    /* at each x-position, sum the maximum run length at that angle into the
       accumulator */
    if(0 == angle) /* vertical, so tangent is infinity */
      slantDx = 0.;
    else
      slantDx = -1.0 / tan(DMath::degreesToRadians(90-angle));
    //       for(j = slantOffset; j < (hdr.w-slantOffset); ++j){
    for(int j = x0; j <= x1; ++j){
      maxrl = 0;
      runlen = 0;
      for(int y = 0; y < lineH; ++y){
	int x;
	x = (int)(0.5+ j + y * slantDx);
	if( (x>=x0) && (x <= x1)){ /* make sure we are within bounds */
	  int idxtmp;
	  idxtmp = (y+y0)*w+x;
// 	    imgCoded[idxtmp*3] = 0;
	  if(0 == p8[idxtmp]){
	    ++runlen;
	    if(runlen > maxrl){
	      maxrl = runlen;
	    }
	  }
	  else
	    runlen = 0;
	} /* end if in bounds */
	else{
	  runlen = 0; /* ignore runs that go off edge of image */
	}
      }
      if(maxrl > rlThresh){
	rgSlantSums[angle+45] += maxrl*maxrl;
	dblWeight += maxrl;
      }
    } /* end for j */
  } /* end for angle */

  //smooth the histogram
  rgSlantSumsTmp[0] = (rgSlantSums[0] + rgSlantSums[1]) / 2;
  for(int aa = 1; aa < 119; ++aa){
    rgSlantSumsTmp[aa]=(rgSlantSums[aa-1]+rgSlantSums[aa]+rgSlantSums[aa+1])/3;
  }
  for(int aa = 0; aa < 120; ++aa){
    rgSlantSums[aa] = rgSlantSumsTmp[aa];
  }

  //use the histogram peak as the slant angle
  slantAngle = 0;
  for(angle = -45; angle <= 60; angle += 1){
    if(rgSlantSums[angle+45] > rgSlantSums[slantAngle+45]){
      slantAngle = angle;
    }
  } /* end for angle */

  if(NULL != weight)
    (*weight) = dblWeight;

  if(NULL != rgSlantHist){
    for(int aa = 0; aa < 120; ++aa){
      rgSlantHist[aa] = rgSlantSums[aa];
    }
  }

  if(NULL != imgAngleHist){//debug tool- return an image of the angle histogram
    DProfile prof;
    int max = 0;
    int htmp;
    D_uint8 *p8ang;
    
    htmp = y1-y0+1;
    imgAngleHist->create(120,htmp,DImage::DImage_u8);
    imgAngleHist->clear();
    p8ang = imgAngleHist->dataPointer_u8();
    for(int i=0; i < 120; ++i){
      if(rgSlantSums[i] > max)
	max = rgSlantSums[i];
    }
    if(0==max)
      max = 1;
    // for(int y=0, idx=0; y < htmp; ++y){
    //   for(int x=0; x < 120; ++x, ++idx){
    // 	if((rgSlantSums[x]/(double)max) >= ((htmp-y)/(double)htmp))
    // 	  p8ang[idx] = 0xee;
    // 	else
    // 	  p8ang[idx] = 0x88;
    //   }
    // }
    // printf("htmp=%d\n", htmp);
    for(int x=0; x < 120; ++x){
      double pct;
      pct = 1.-rgSlantSums[x] / (double)max;
      imgAngleHist->drawLine(x,htmp-1,x,(int)(pct*(htmp-1)), 128);
    }


  }

  return (double)slantAngle;
}

#ifndef D_NOTHREADS
#include "dthreads.h"
/// \cond
typedef struct{
  DImage *pImgSrc;//pointer to the image (shared by all threads)
  int numThreads;//how many threads are processing
  int threadNum;//which thread number this is (0..numThreads-1)
  unsigned int *rgSlantSums;//each thread points to different array.
  double *rgAngles;//all threads point to same array of double[numTextlines]
  double *rgWeights;//all threads point to same array of double[numTextlines]
  DRect *rgRects;
  int numTextlines;
  int rlThresh; // runlength thresh -- how long an rl needs to be to count as rl
} SLANT_THREAD_PARMS;
/// \endcond
#endif /* D_NOTHREADS not defined */


void* DSlantAngle::getSlant_thread_func(void *params){
  SLANT_THREAD_PARMS *pparms;
  int numThreads;
  int w, h;
  D_uint8 *p8;

  int runlen, maxrl; /* maximum slant runlen */
  double slantDx;
  int lineH;
  int slantOffset, slantAngle, angle;
  double dblWeight;
  DImage *pimg;
  int rlThresh;

  pparms = (SLANT_THREAD_PARMS*)params;
  
  numThreads = pparms->numThreads;
  pimg = pparms->pImgSrc;
  rlThresh = pparms->rlThresh;
  w = pimg->width();
  h = pimg->height();
  p8 = pimg->dataPointer_u8();
  for(int i=0; i < 120; ++i)
    pparms->rgSlantSums[i] = 0;
  for(int tl=pparms->threadNum; tl < (pparms->numTextlines); tl+=numThreads){
    int x0, y0, x1, y1;
    unsigned int rgSlantSums[120];
    x0 = pparms->rgRects[tl].x;
    y0 = pparms->rgRects[tl].y;
    x1 = pparms->rgRects[tl].x + pparms->rgRects[tl].w - 1;
    y1 = pparms->rgRects[tl].y + pparms->rgRects[tl].h - 1;
    lineH = y1-y0+1;
    memset(rgSlantSums, 0, sizeof(int)*120);
    dblWeight = 0.;
    
    for(angle = -45; angle <= 60; angle += 1){
      /* at each x-position, sum the maximum run length at that angle into the
	 accumulator */
      if(0 == angle) /* vertical, so tangent is infinity */
	slantDx = 0.;
      else
	slantDx = -1.0 / tan(DMath::degreesToRadians(90-angle));
      //       for(j = slantOffset; j < (hdr.w-slantOffset); ++j){
      for(int j = x0; j <= x1; ++j){
	maxrl = 0;
	runlen = 0;
	for(int y = 0; y < lineH; ++y){
	  int x;
	  x = (int)(0.5+ j + y * slantDx);
	  if( (x>=x0) && (x <= x1)){ /* make sure we are within bounds */
	    int idxtmp;
	    idxtmp = (y+y0)*w+x;
	    // 	    imgCoded[idxtmp*3] = 0;
	    if(0 == p8[idxtmp]){
	      ++runlen;
	      if(runlen > maxrl){
		maxrl = runlen;
	      }
	    }
	    else
	      runlen = 0;
	  } /* end if in bounds */
	  else{
	    runlen = 0; /* ignore runs that go off edge of image */
	  }
	}
	if(maxrl > rlThresh){
	  rgSlantSums[angle+45] += maxrl*maxrl;
	  dblWeight += maxrl;
	}
      } /* end for j */
    } /* end for angle */
    for(int i=0; i < 120; ++i)
      pparms->rgSlantSums[i] += rgSlantSums[i];
    if(NULL != (pparms->rgWeights)){
      pparms->rgWeights[tl] = dblWeight;
    }
    if(NULL != (pparms->rgAngles)){
      // need to independently figure out the angle for this particular textline
      unsigned int rgSlantSumsTmp[120];
      
      //smooth the histogram
      rgSlantSumsTmp[0] = (rgSlantSums[0] + rgSlantSums[1]) / 2;
      for(int aa = 1; aa < 119; ++aa){
	rgSlantSumsTmp[aa]=(rgSlantSums[aa-1]+rgSlantSums[aa]+rgSlantSums[aa+1])/3;
      }
      // for(int aa = 0; aa < 120; ++aa){
      // 	rgSlantSums[aa] = rgSlantSumsTmp[aa];
      // }
      
      //use the smoothed histogram peak as the slant angle
      slantAngle = 0;
      for(angle = -45; angle <= 60; angle += 1){
	if(rgSlantSumsTmp[angle+45] > rgSlantSumsTmp[slantAngle+45]){
	  slantAngle = angle;
	}
      } /* end for angle */
      pparms->rgAngles[tl] = slantAngle;
    }
  }

}

/** If rgAngles or rgWeights are non-NULL, they must be already allocated
    with enough space for numTextlines doubles.
 */
double DSlantAngle::getAllTextlinesSlantAngleDeg(DImage &imgBW, int rlThresh,
						 DRect *rgTextlineRects,
						 int numRects,
						 int numThreads,
						 double *rgAngles,
						 double *rgWeights){
					     // int rlThresh,
					     // int x0,int y0,int x1,int y1,
					     // double *weight,
					     // unsigned int *rgSlantHist,
					     // DImage *imgAngleHist){
  int *rgLineSlantAngles;
  int lineH;
  int slantOffset, slantAngle, angle;
  unsigned int rgSlantSums[120];
  unsigned int rgSlantSumsTmp[120];
  int runlen, maxrl; /* maximum slant runlen */
  double slantDx;
  int w, h;
  D_uint8 *p8;
  double dblWeight = 0;
  SLANT_THREAD_PARMS *rgThreadParms;
#ifndef D_NOTHREADS
  pthread_t rgThreadID[256];
#else
  numThreads = 1;
#endif

  rgThreadParms = (SLANT_THREAD_PARMS*)malloc(sizeof(SLANT_THREAD_PARMS) *
					      numThreads);
  D_CHECKPTR(rgThreadParms);
  
  

  for(int tnum = numThreads-1; tnum >= 0; --tnum){//so other threads launch
    //before thread 0 calls the function so they don't have to wait for it
    rgThreadParms[tnum].pImgSrc = &imgBW;
    rgThreadParms[tnum].rlThresh = rlThresh;
    rgThreadParms[tnum].numThreads = numThreads;
    rgThreadParms[tnum].threadNum = tnum;
    rgThreadParms[tnum].rgSlantSums = new unsigned int[120];
    rgThreadParms[tnum].rgAngles = rgAngles;
    rgThreadParms[tnum].rgWeights = rgWeights;
    rgThreadParms[tnum].numTextlines = numRects;
    rgThreadParms[tnum].rgRects = rgTextlineRects;

    D_CHECKPTR(rgThreadParms[tnum].rgSlantSums);
#ifdef D_NOTHREADS
    getSlant_thread_func(rgThreadParms);
#else
    if(0 == tnum){//don't spawn thread zero. Use the current thread
      DSlantAngle::getSlant_thread_func(rgThreadParms);
    }
    else{//spawn all other threads besides thread zero
      if(0 != pthread_create(&rgThreadID[tnum], NULL,
			     DSlantAngle::getSlant_thread_func,
			     &rgThreadParms[tnum])){
	fprintf(stderr,"DSlantAngle::getAllTextlinesSlantAngleDeg() "
		"failed to spawn thread "
		"#%d. Exiting.\n",tnum);
	exit(1);
      }
    }
#endif

  }
#ifndef D_NOTHREADS
  // wait for all threads to finish
  for(int tnum = 1; tnum < numThreads; ++tnum){
    if(pthread_join(rgThreadID[tnum],NULL)){
      fprintf(stderr, "DSlantAngle::getAllTextlinesSlantAngleDeg() "
	      "failed to join thread %d. Exiting.\n", tnum);
      exit(1);
    }
  }
#endif
  for(int i=0; i < 120; ++i){
    rgSlantSums[i] = 0;
    rgSlantSumsTmp[i] = 0;
  }
  for(int tnum=0; tnum < numThreads; ++tnum){
    for(int i=0; i < 120; ++i){
      rgSlantSums[i] += rgThreadParms[tnum].rgSlantSums[i];
    }
  }
  for(int tnum=0; tnum < numThreads; ++tnum){
    delete [] (rgThreadParms[tnum].rgSlantSums);
  }
  



  //smooth the histogram
  rgSlantSumsTmp[0] = (rgSlantSums[0] + rgSlantSums[1]) / 2;
  for(int aa = 1; aa < 119; ++aa){
    rgSlantSumsTmp[aa]=(rgSlantSums[aa-1]+rgSlantSums[aa]+rgSlantSums[aa+1])/3;
  }
  for(int aa = 0; aa < 120; ++aa){
    rgSlantSums[aa] = rgSlantSumsTmp[aa];
  }

  //use the histogram peak as the slant angle
  slantAngle = 0;
  for(angle = -45; angle <= 60; angle += 1){
    if(rgSlantSums[angle+45] > rgSlantSums[slantAngle+45]){
      slantAngle = angle;
    }
  } /* end for angle */


#if 0
  if(NULL != imgAngleHist){//debug tool- return an image of the angle histogram
    DProfile prof;
    int max = 0;
    int htmp;
    D_uint8 *p8ang;
    
    htmp = y1-y0+1;
    imgAngleHist->create(120,htmp,DImage::DImage_u8);
    imgAngleHist->clear();
    p8ang = imgAngleHist->dataPointer_u8();
    for(int i=0; i < 120; ++i){
      if(rgSlantSums[i] > max)
	max = rgSlantSums[i];
    }
    if(0==max)
      max = 1;
   for(int x=0; x < 120; ++x){
      double pct;
      pct = 1.-rgSlantSums[x] / (double)max;
      imgAngleHist->drawLine(x,htmp-1,x,(int)(pct*(htmp-1)), 128);
    }
  }
#endif

  free(rgThreadParms);

  return (double)slantAngle;
}


