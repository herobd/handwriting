#ifndef DCONNECTEDCOMPONENTLABELER_H
#define DCONNECTEDCOMPONENTLABELER_H

#include "dimage.h"

class DConnectedComponentInfo; // forward declaration

class DConnectedComponentLabeler{
public:
  static DImage getCCimage(const DImage &imgSrc, int *numCCs,
			   D_uint32 BGval = (D_uint32)-1,
			   bool f8connected=true,
			   bool fSequentialIDs = true);
  static void getCCimage_(DImage &imgDst, const DImage &imgSrc, int *numCCs,
			  D_uint32 BGval = (D_uint32)-1,
			  bool f8connected=true,
			  bool fSequentialIDs = true);
  static void getCCimageForVal_(DImage &imgDst, const DImage &imgSrc,
				int *numCCs,
				D_uint32 val,
				bool f8connected=true,
				bool fSequentialIDs = true);
  static void getCCInfoFromCCimage(DImage &imgCCMap,
				   DConnectedComponentInfo *rgCCInfo,
				   int numComponentsPlusBG);

  static DImage getRGBImageFromCCImage(DImage &imgSrc, bool fKeepNumbers=false);
};


#endif
