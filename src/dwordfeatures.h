#ifndef DWORDFEATURES_H
#define DWORDFEATURES_H

#include "dimage.h"
#include "dfeaturevector.h"

class DWordFeatures{
public:
  static DFeatureVector extractWordFeatures(DImage &img,
					    bool fUseGrayscaleProf=true,
					    bool fUseMyTransitions=false,
					    bool fInkIsBlack=true,
					    bool fRangeIs255=true,
					    D_uint8 tval=127,
					    double weight_prof=1.,
					    double weight_upper=1.,
					    double weight_lower=1.,
					    double weight_trans=1.);
  static DFeatureVector extractWordFeatures_flt(DImage &img,
						bool fUseGrayscaleProf=false,
						bool fUseMyTransitions=true,
						bool fInkIsBlack=true,
						bool fRangeIs255=true,
						float weight_prof=1.,
						float weight_upper=1.,
						float weight_lower=1.,
						float weight_trans=1.);
  static DFeatureVector extractWordFourierFeatures_flt();


};



#endif
