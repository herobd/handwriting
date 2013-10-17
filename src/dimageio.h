#ifndef DIMAGEIO
#define DIMAGEIO

#include <stdio.h>
#include "ddefs.h"

class DImage;

#define MAX_PNM_COMMENTS 100 /* this can be increased if necessary */

class DImageIO{
public:
  static bool load_image_pnm(DImage *pImg, const char *stPath);
  static bool load_image_pbm_plain(DImage *pImg, FILE *fin);
  static bool load_image_pgm_plain(DImage *pImg, FILE *fin);
  static bool load_image_ppm_plain(DImage *pImg, FILE *fin);
  static bool load_image_pbm_raw(DImage *pImg, FILE *fin);
  static bool load_image_pgm_raw(DImage *pImg, FILE *fin);
  static bool load_image_ppm_raw(DImage *pImg, FILE *fin);
  static bool load_image_png(DImage *pImg, FILE *fin);
  static bool load_image_jpeg(DImage *pImg, FILE *fin);
  static bool load_image_tiff(DImage *pImg, const char *stPath);
  static bool load_image_gif(DImage *pImg, const char *stPath);
  static bool load_image_bmp(DImage *pImg, const char *stPath);
  static bool load_image_gnuplot(DImage *pImg, FILE *fin);

  static bool save_image_pnm(DImage *pImg, const char *stPath,
			     bool fSaveProps=true, bool fSaveComments=true);
  static bool save_image_ppm_raw(DImage *pImg, const char *stPath,
				 bool fSaveProps=true,bool fSaveComments=true);
  static bool save_image_pgm_raw(DImage *pImg, const char *stPath,
				 bool fSaveProps=true,bool fSaveComments=true);
  static bool save_image_pbm_raw(DImage *pImg, const char *stPath,
				 bool fSaveProps=true,bool fSaveComments=true);
  static bool save_image_gnuplot(DImage *pImg, const char *stPath,
				 bool fSaveProps=true,bool fSaveComments=true,
				 int everyNthPixel=1);
  static bool save_image_png(DImage *pImg, const char *stPath,
			     bool fSaveProps=true,bool fSaveComments=true,
			     int compression=6);
  static bool save_image_jpeg(DImage *pImg, const char *stPath,
			      bool fSaveProps=true,bool fSaveComments=true,
			      int quality=75, bool fProgressive=false,
			      bool fOptimize=false);

  static void set_alloc_method(D_AllocationMethod allocMeth);
  static bool get_image_width_height_chan(const char *stPath,
					  int *width, int *height, int *chan);
  
private:
  typedef struct {
    int type; // 1-6 (P1...P6)
    int w;
    int h;
    int max; // will usually be 255
    long dataOffs; // offset at which first data byte can be read
    int numComments; // number of comments in the file
    char *rgComments[MAX_PNM_COMMENTS];
  } D_PNM_HEADER_S;

  static bool read_pnm_header(D_PNM_HEADER_S *hdr, FILE *fin);
  static void clear_pnm_header_comments(D_PNM_HEADER_S *hdr);
  static void extractCommentProps(DImage *pImg, D_PNM_HEADER_S *hdr);
  static bool readDataBlock(FILE *fin, D_uint8 *rgBuff, size_t len);
  static bool writeDataBlock(FILE *fout, D_uint8 *rgBuff, size_t len);
  static int readComments(D_PNM_HEADER_S *hdr, FILE *fin);
  static int write_img_props_pnm(DImage *pImg, FILE *fout);
  static int write_img_comments_pnm(DImage *pImg, FILE *fout);
  static D_AllocationMethod _allocMethod;
};


#endif
