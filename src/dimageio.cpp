#include "dimageio.h"
#include "dimage.h"
#include <string>
#include <cstring>

#ifndef D_NOJPEG

#ifdef WIN32
#ifndef XMD_H
#define ADDED_XMD_H
#define XMD_H
#endif /*ndef XMD_H*/
#endif /* WIN32 */
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
/* stdio.h needs to be included before jpeglib.h.  It is included indirectly
   because dimageio.h includes it (above) */
#include <jpeglib.h>
#include <jerror.h>
#ifdef __cplusplus
}
#endif /* __cplusplus */
#ifdef ADDED_XMD_H
#undef XMD_H
#undef ADDED_XMD_H
#endif /* ADDED_XMD_H */
#endif /*D_NOJPEG*/

#ifndef D_NOTIFF
#include "tiffio.h"
#endif /*D_NOTIFF*/

#ifndef D_NOPNG
#include "png.h"
#endif /*D_NOPNG*/

#ifndef WIN32
#include <unistd.h>
#include <time.h>
#endif

D_AllocationMethod DImageIO::_allocMethod = AllocationMethod_malloc;


int DImageIO::readComments(D_PNM_HEADER_S *hdr, FILE *fin){
  char stTmp[1025];// not static for purpose of being thread-safe
  int peeker;
  int fCommentSkipped;
  int fSpaceSkipped;
  do{
    fCommentSkipped = 0;
    fSpaceSkipped = 0;
    peeker = fgetc(fin);
    if('#' == peeker){
      fgets(stTmp, 1024, fin); /* ignore comment */
      if(hdr->numComments < MAX_PNM_COMMENTS){
	hdr->rgComments[hdr->numComments]=(char*)malloc(1+strnlen(stTmp,1024));
	D_CHECKPTR(hdr->rgComments[hdr->numComments]);
	strncpy(hdr->rgComments[hdr->numComments], stTmp, strnlen(stTmp,1024));
	hdr->rgComments[hdr->numComments][strnlen(stTmp,1024)] = '\0';
	// get rid of newline if it is in the string
	if(hdr->rgComments[hdr->numComments][strnlen(stTmp,1024)-1]=='\n')
	  hdr->rgComments[hdr->numComments][strnlen(stTmp,1024)-1]='\0';
	++(hdr->numComments);
      }
      fCommentSkipped = 1;
    }
    else if(isspace(peeker)){
      fSpaceSkipped = 1;
    }
    else if(EOF == peeker){
      return -1;
    }
    else{
      ungetc(peeker, fin);
    }
  } while(fCommentSkipped || fSpaceSkipped);
  return 0;
}

bool DImageIO::read_pnm_header(D_PNM_HEADER_S *hdr, FILE *fin){
  char cTmp, cTmp2;
  if(0 != fseek(fin, 0L, SEEK_SET)){
    fprintf(stderr, "DImageIO:read_pnm_header() fseek failed\n");
    return false;
  }
  (hdr->type) = 0;
  (hdr->w) = 0;
  (hdr->h) = 0;
  (hdr->max) = 0;
  (hdr->dataOffs) = 0;
  (hdr->numComments) = 0;
  if(0 != readComments(hdr, fin))
    return false;
  fscanf(fin, "%c%c", &cTmp, &cTmp2);
  if(('P' != cTmp) ||
      (cTmp2 < '1') ||
      (cTmp2 > '6')){
    return false;
  }
  (hdr->type) = ((int)cTmp2) - ((int)'0');
  if(0 != readComments(hdr, fin))
    return false;
  fscanf(fin, "%d", &(hdr->w));
  if(0 != readComments(hdr, fin))
    return false;
  fscanf(fin, "%d", &(hdr->h));
  if(((hdr->type) != 1) &&((hdr->type)!=4)){
    if(0 != readComments(hdr, fin))
      return false;
    fscanf(fin, "%d", &(hdr->max));
  }
  else
    (hdr->max) = 255;
  fscanf(fin, "%*c"); /* last whitespace before data */
  (hdr->dataOffs) = ftell(fin);
  if(-1 == hdr->dataOffs)
    return false;
  /* sanity checks */
  if((hdr->max < 1) || (hdr->max > 65535)){
    fprintf(stderr, "DImageIO::read_pnm_header() max out of range\n");
    DImageIO::clear_pnm_header_comments(hdr);
    return false;
  }


  return true;
}

void DImageIO::clear_pnm_header_comments(D_PNM_HEADER_S *hdr){
  for(int i = 0; i < hdr->numComments; ++i){
    free(hdr->rgComments[i]);
  }
  hdr->numComments = 0;
}

void DImageIO::extractCommentProps(DImage *pImg, D_PNM_HEADER_S *hdr){
  pImg->clearProperties();
  pImg->clearComments();
  for(int i = 0; i < hdr->numComments; ++i){
    if(0 == strncmp("DProp:", hdr->rgComments[i],6)){
      for(unsigned int j = 6; j < strlen(hdr->rgComments[i]); ++j){
	if('=' == hdr->rgComments[i][j]){
	  pImg->mapProps[std::string(&(hdr->rgComments[i][6]), j-6)] =
	    std::string(&(hdr->rgComments[i][j+1]));
	  break;
	}
      }
    }
    else{ // not a property, a regular comment or blank line
      pImg->vectComments.push_back(std::string(hdr->rgComments[i]));
    }
  }
}

bool DImageIO::readDataBlock(FILE *fin, D_uint8 *rgBuff, size_t len){
  size_t bytesRead = 0;
  size_t readRet;

  while((!feof(fin)) && (!ferror(fin)) && (bytesRead < len)){
    readRet = fread(&(rgBuff[bytesRead]), sizeof(D_uint8), len-bytesRead, fin);
    if((int)readRet >= 0)
      bytesRead += readRet;
  }
  if(bytesRead < len){ // error or EOF
#ifdef DEBUG
    fprintf(stderr,"only read %d bytes from file (tried to read %d bytes)\n",
	    bytesRead, len);
    //    memset(&rgBuff[bytesRead],0,len-bytesRead);
#endif
    return false;
  }
  return true;
}
bool DImageIO::writeDataBlock(FILE *fout, D_uint8 *rgBuff, size_t len){
  size_t bytesWrote;
  size_t bytesWroteTotal = 0;

  while((bytesWroteTotal < len) && (!ferror(fout))){
    bytesWrote = fwrite(&(rgBuff[bytesWroteTotal]), sizeof(unsigned char),
			len-bytesWroteTotal, fout);
    if(bytesWrote < 1){
      return false;
    }
    bytesWroteTotal += bytesWrote;
  }
  if(bytesWroteTotal < len) // error or EOF
    return false;
  return true;
}

bool DImageIO::load_image_pbm_plain(DImage *pImg, FILE *fin){
  fprintf(stderr,"NYI! (%s:%d)\n", __FILE__, __LINE__);
  pImg->clearProperties();
  pImg->clearComments();
  return false;
}
bool DImageIO::load_image_pgm_plain(DImage *pImg, FILE *fin){
  D_PNM_HEADER_S hdr;
  bool fHeaderOK;
  size_t numPxls = 0;
  int iTmp;
  D_uint8 *pData8;
  D_uint16 *pData16;

  fHeaderOK = DImageIO::read_pnm_header(&hdr, fin);
  if(fHeaderOK){
    pImg->clearProperties();
    pImg->clearComments();
    DImageIO::extractCommentProps(pImg, &hdr);
  }
  DImageIO::clear_pnm_header_comments(&hdr);
  
  /* sanity checks */
  if(!fHeaderOK){
    fprintf(stderr, "DImageIO::load_image_pgm_plain() failed to read "
	    "header\n");
    return false;
  }
  if(2 != hdr.type){ // not a raw ppm
    fprintf(stderr, "DImageIO::load_image_pgm_plain() image format not P2!\n");
    return false;
  }
  
  numPxls = hdr.w * hdr.h;
  if(hdr.max < 256){
    pImg->create(hdr.w, hdr.h, DImage::DImage_u8,
		 1, DImageIO::_allocMethod);
    pData8 = pImg->dataPointer_u8();
    for(unsigned int i = 0; i < numPxls; ++i){
      if(1 != fscanf(fin, "%d", &iTmp)){
	fprintf(stderr, "DImageIO::load_image_pgm_plain() premature EOF!\n");
	return false;
      }
      pData8[i] = (D_uint8)iTmp;
    }
  }
  else{
    pImg->create(hdr.w, hdr.h, DImage::DImage_u16, 1, DImageIO::_allocMethod);
    pData16 = pImg->dataPointer_u16();
    for(unsigned int i = 0; i < numPxls; ++i){
      if(1 != fscanf(fin, "%d", &iTmp)){
	fprintf(stderr, "DImageIO::load_image_pgm_plain() premature EOF!\n");
	return false;
      }
      pData16[i] = (D_uint16)iTmp;
    }
  }
  
  return true;
}
bool DImageIO::load_image_ppm_plain(DImage *pImg, FILE *fin){
  D_PNM_HEADER_S hdr;
  bool fHeaderOK;
  size_t numPxls = 0;
  int iTmp;
  D_uint8 *pData8;
  D_uint16 *pData16;

  fHeaderOK = DImageIO::read_pnm_header(&hdr, fin);
  if(fHeaderOK){
    pImg->clearProperties();
    pImg->clearComments();
    DImageIO::extractCommentProps(pImg, &hdr);
  }
  DImageIO::clear_pnm_header_comments(&hdr);
  
  /* sanity checks */
  if(!fHeaderOK){
    fprintf(stderr, "DImageIO::load_image_pgm_plain() failed to read "
	    "header\n");
    return false;
  }
  if(3 != hdr.type){ // not a plain ppm
    fprintf(stderr, "DImageIO::load_image_pgm_plain() image format not P3!\n");
    return false;
  }
  
  numPxls = 3 * hdr.w * hdr.h;
  if(hdr.max < 256){
    pImg->create(hdr.w, hdr.h, DImage::DImage_RGB,
		 3, DImageIO::_allocMethod);
    pData8 = pImg->dataPointer_u8();
    for(unsigned int i = 0; i < numPxls; ++i){
      if(1 != fscanf(fin, "%d", &iTmp)){
	fprintf(stderr, "DImageIO::load_image_pgm_plain() premature EOF!\n");
	return false;
      }
      pData8[i] = (D_uint8)iTmp;
    }
  }
  else{
    pImg->create(hdr.w, hdr.h, DImage::DImage_RGB_16,
		 3, DImageIO::_allocMethod);
    pData16 = pImg->dataPointer_u16();
    for(unsigned int i = 0; i < numPxls; ++i){
      if(1 != fscanf(fin, "%d", &iTmp)){
	fprintf(stderr, "DImageIO::load_image_pgm_plain() premature EOF!\n");
	return false;
      }
      pData16[i] = (D_uint16)iTmp;
    }
  }
  return true;
}
bool DImageIO::load_image_pbm_raw(DImage *pImg, FILE *fin){
  fprintf(stderr,"NYI! (%s:%d)\n", __FILE__, __LINE__);
  pImg->clearProperties();
  pImg->clearComments();
  return false;
}
bool DImageIO::load_image_pgm_raw(DImage *pImg, FILE *fin){
  D_PNM_HEADER_S hdr;
  bool fHeaderOK;
  size_t numBytes = 0;

  fHeaderOK = DImageIO::read_pnm_header(&hdr, fin);
  if(fHeaderOK){
    pImg->clearProperties();
    pImg->clearComments();
    DImageIO::extractCommentProps(pImg, &hdr);
  }
  DImageIO::clear_pnm_header_comments(&hdr);
  
  /* sanity checks */
  if(!fHeaderOK){
    fprintf(stderr, "DImageIO::load_image_pgm_raw() failed to read header\n");
    return false;
  }
  if(5 != hdr.type){ // not a raw ppm
    fprintf(stderr, "DImageIO::load_image_pgm_raw() image format not P5!\n");
    return false;
  }
  
  if(hdr.max < 256){
    pImg->create(hdr.w, hdr.h, DImage::DImage_u8,
		 1, DImageIO::_allocMethod);
    numBytes = hdr.w * hdr.h;
  }
  else{
    pImg->create(hdr.w, hdr.h, DImage::DImage_u16, 1, DImageIO::_allocMethod);
    numBytes = hdr.w * hdr.h * 2; // "*2" is because 16-bit per pixel
  }
//   if(0 != fseek(fin, hdr.dataOffs, SEEK_SET)){
//     fprintf(stderr, "DImageIO:load_image_pgm_raw() fseek failed\n");
//     return false;
//   }
  if(!DImageIO::readDataBlock(fin, pImg->pData, numBytes)){
    fprintf(stderr, "DImageIO:load_image_pgm_raw() failed reading data\n");
    return false;
  }
  // TODO: if two-bytes per sample, swap bytes if necessary
  
  return true;
}
bool DImageIO::load_image_ppm_raw(DImage *pImg, FILE *fin){
  D_PNM_HEADER_S hdr;
  bool fHeaderOK;
  size_t numBytes = 0;

  fHeaderOK = DImageIO::read_pnm_header(&hdr, fin);
  if(fHeaderOK){
    pImg->clearProperties();
    pImg->clearComments();
    DImageIO::extractCommentProps(pImg, &hdr);
  }
  DImageIO::clear_pnm_header_comments(&hdr);
    
  /* sanity checks */
  if(!fHeaderOK){
    fprintf(stderr, "DImageIO::load_image_ppm_raw() failed to read header\n");
    return false;
  }
  if(6 != hdr.type){ // not a raw ppm
    fprintf(stderr, "DImageIO::load_image_ppm_raw() image format not P6!\n");
    return false;
  }

  if(hdr.max < 256){
    pImg->create(hdr.w, hdr.h, DImage::DImage_RGB, 3, DImageIO::_allocMethod);
    numBytes = hdr.w * hdr.h * 3; // "*3" is because RGB
  }
  else{
    pImg->create(hdr.w, hdr.h, DImage::DImage_RGB_16,
		 3, DImageIO::_allocMethod);
    numBytes = hdr.w * hdr.h * 6; // "*6" is because RRGGBB
  }
//   if(0 != fseek(fin, hdr.dataOffs, SEEK_SET)){
//     fprintf(stderr, "DImageIO:load_image_ppm_raw() fseek failed\n");
//     return false;
//   }
  if(!DImageIO::readDataBlock(fin, pImg->pData, numBytes)){
    fprintf(stderr, "DImageIO:load_image_ppm_raw() failed reading data\n");
    return false;
  }
  // TODO: if two-bytes per sample, swap bytes if necessary
  
  return true;
}

// ///This function is a callback used by DImageIO::load_image_png()
// void png_unknown_chunk_callback(png_ptr ptr, png_unknown_chunkp chunk){
//   // this function is called when an unrecognized PNG chuck is found.
//   //TODO: implement functionality of this function so that my own 
//   //image comments and properties can be read back out of PNG files that
//   //I create.
//   printf("png_unknown_chunk_callback()\n");
//   return 0;
// }


///Load a PNG image into pImg.
/**This function uses the PNG Reference Library (libpng), and also
 *follows the code examples provided with that library for its
 *implementation
 */
bool DImageIO::load_image_png(DImage *pImg, FILE *fin){
#ifdef D_NOPNG
  fprintf(stderr, "PNG image support is not compiled.  Need to make sure that"
	  "libpng is available on your system, and recompile without the "
	  "D_NOPNG definition.\n");
  return false;
#else
  png_uint_32 width, height;
  D_uint8 *pDst;
  int bit_depth;
  int color_type;
  png_bytep *row_pointers;
  png_bytep rCur;

  png_structp png_ptr;
  png_infop info_ptr;
  png_infop end_info;
//   png_unknown_chunkp user_chunk_ptr;

  pImg->clearProperties();
  pImg->clearComments();


//   png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
// 				   (png_voidp)user_error_ptr,
// 				   user_error_fn, user_warning_fn);
  png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL,NULL);
// 				   (png_voidp)user_error_ptr,
// 				   user_error_fn, user_warning_fn);
  if(!png_ptr)
    return false;

  info_ptr = png_create_info_struct(png_ptr);
  if(!info_ptr){
    png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
    return false;
  }

  end_info = png_create_info_struct(png_ptr);
  if(!end_info){
    png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);
    return false;
  }

  rewind(fin);
  png_init_io(png_ptr, fin);


//   user_chunk_ptr = png_get_ser_chunk_ptr(png_ptr);
//   if(!user_chunk_ptr){
//     png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
//     return false;
//   }

//   png_set_read_user_chunk_fn(png_ptr, user_chunk_ptr,
// 			     png_unknown_chunk_callback);
  
  //   png_set_keep_unknown_chunks(png_ptr, info_ptr, keep,
  //         chunk_list, num_chunks);
  
  png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_STRIP_ALPHA |
	       PNG_TRANSFORM_PACKING, NULL);

  png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
	       NULL, NULL, NULL);

  
//   printf("PNG info: w=%d h=%d bit_depth=%d color_type=%d\n",
// 	 (int)width, (int)height, bit_depth, color_type);

  row_pointers = png_get_rows(png_ptr, info_ptr);

  switch(color_type){
    case PNG_COLOR_TYPE_GRAY:
    case PNG_COLOR_TYPE_GRAY_ALPHA:
      if(16 == bit_depth){
	pImg->create(width, height, DImage::DImage_u16, 1,
		     DImageIO::_allocMethod);
	pDst = pImg->pData;
	for(unsigned int y = 0; y < height; ++y){
	  rCur = row_pointers[y];
	  for(unsigned int x = 0; x < width*2; ++x){
	    (*pDst) = rCur[x];
	    ++pDst;
	  }
	}
      }
      else{
	pImg->create(width, height, DImage::DImage_u8, 1,
		     DImageIO::_allocMethod);
	pDst = pImg->pData;
	for(unsigned int y = 0; y < height; ++y){
	  rCur = row_pointers[y];
	  for(unsigned int x = 0; x < width; ++x){
	    (*pDst) = rCur[x];
	    ++pDst;
	  }
	}
      }
      break;
    case PNG_COLOR_TYPE_RGB:
    case PNG_COLOR_TYPE_RGB_ALPHA:
      if(16 == bit_depth){
	pImg->create(width, height, DImage::DImage_RGB_16, 3,
		     DImageIO::_allocMethod);
	pDst = pImg->pData;
	for(unsigned int y = 0; y < height; ++y){
	  rCur = row_pointers[y];
	  for(unsigned int x = 0; x < width*6; ++x){
	    (*pDst) = rCur[x];
	    ++pDst;
	  }
	}
      }
      else{
	pImg->create(width, height, DImage::DImage_RGB, 3,
		     DImageIO::_allocMethod);
	pDst = pImg->pData;
	for(unsigned int y = 0; y < height; ++y){
	  rCur = row_pointers[y];
	  for(unsigned int x = 0; x < width*3; ++x){
	    (*pDst) = rCur[x];
	    ++pDst;
	  }
	}
      }
      break;
    case PNG_COLOR_TYPE_PALETTE:
      pImg->create(width, height, DImage::DImage_RGB, 3,
		   DImageIO::_allocMethod);
      fprintf(stderr, "NYI!(%s:%d)\n", __FILE__, __LINE__);
      break;
    default:
      fprintf(stderr, "DImageIO::load_image_png() unknown PNG color type\n");
      break;
  }


  // grab the comments out here before destroying the structures
  // grab the comments out here before destroying the structures
  // grab the comments out here before destroying the structures


  png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);

  return true;
#endif
}

///load a JPEG-encoded image from fin into pImg
/**"This software is based in part on the work of the Independent JPEG
 * Group".  Code in this function was created by following
 * instructions in libjpeg.doc, as well as copying examples found in
 * djpeg.c, example.c, and jdatasrc.c in the IJG's distribution of
 * libjpeg.
*/
bool DImageIO::load_image_jpeg(DImage *pImg, FILE *fin){
#ifdef D_NOJPEG
  fprintf(stderr, "JPEG image support is not compiled.  Need to make sure that"
	  "libjpeg is available on your system, and recompile without the "
	  "D_NOJPEG definition.\n");
  return false;
#else
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPARRAY output_pointers;
  int linesRead;

  if(1 != sizeof(JSAMPLE)){
    fprintf(stderr,
	    "DImageIO::load_image_jpeg() expected 8-bit JSAMPLE size."
	    " Other sizes not supported\n");
    return false;
  }

  if(0 != fseek(fin, 0L, SEEK_SET)){
    fprintf(stderr, "DImageIO:load_image_jpeg() fseek failed\n");
    return false;
  }
  pImg->clearProperties();
  pImg->clearComments();
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fin);
  jpeg_save_markers(&cinfo, JPEG_COM, 2048);//keep the comments in memory
  if(JPEG_HEADER_OK != jpeg_read_header(&cinfo, boolean(true))){
    fprintf(stderr, "DImageIO:load_image_jpeg() error in jpeg header\n");
    jpeg_abort_decompress(&cinfo);
    return false;
  }
  //follow the comment marker chain and save them into the dimage comments
  jpeg_saved_marker_ptr nextMarker;
  char markerData[2049];
  nextMarker = cinfo.marker_list;
  while(nextMarker != NULL){
    if(nextMarker->marker == JPEG_COM){//comment
      int len = nextMarker->data_length;
      if(len > 2048)
	len = 2048;
      strncpy(markerData, (const char*)nextMarker->data,
	      nextMarker->data_length);
      markerData[len]='\0';
      if('\n'==markerData[len-1])
	markerData[len-1]='\0';
      std::string strTmp;
      strTmp = std::string(markerData);
      //check if it is a property or just a regular comment
      if(0 == strncasecmp(markerData,"DProp:",6)){
	int eqidx=-1;
	for(int i=6; i<len; ++i){
	  if('=' == markerData[i]){
	    eqidx = i;
	    break;
	  }
	}
	if(eqidx > -1){//found equals sign
	  pImg->setProperty(strTmp.substr(6,eqidx-6),
			    strTmp.substr(eqidx+1));
	}	
      }
      else{
	pImg->addComment(strTmp);
      }
    }
    nextMarker = nextMarker->next;
  }



  // set the parameter for floating point calculation (default is
  // integer, which is slightly less accurate)
  cinfo.dct_method = JDCT_FLOAT;
  jpeg_start_decompress(&cinfo); // reads the header, etc.
  // now we have width and height, so create the image buffer, and pointers
  // to each row in the buffer
  output_pointers = (JSAMPARRAY)malloc(sizeof(JSAMPROW) * cinfo.output_height);
  D_CHECKPTR(output_pointers);
  if(3 == cinfo.out_color_components){ // RGB
    pImg->create(cinfo.output_width, cinfo.output_height,
		 DImage::DImage_RGB, 3, DImageIO::_allocMethod);
    for(int li = 0; li < pImg->_h; ++li){
      output_pointers[li] = (JSAMPROW)&(pImg->pData[li * pImg->_w * 3]);
    }
  }
  else if(1==cinfo.out_color_components){ // GrayScale
    pImg->create(cinfo.output_width, cinfo.output_height,
		 DImage::DImage_u8, 1, DImageIO::_allocMethod);
    for(int li = 0; li < pImg->_h; ++li){
      output_pointers[li] = (JSAMPROW)&(pImg->pData[li * pImg->_w]);
    }
  }
  else{
    fprintf(stderr,
	    "DImageIO:load_image_jpeg() unexpected number of channels\n");
    jpeg_destroy_decompress(&cinfo);
    return false;
  }
  linesRead = 0;
  while((int)cinfo.output_scanline < (pImg->_h)){
    linesRead += jpeg_read_scanlines(&cinfo,
				     &output_pointers[cinfo.output_scanline],
				     pImg->_h);
  }
  free(output_pointers);
  output_pointers = NULL;
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  return true;
#endif
}
///Loads a TIFF image into pImg (given the file path)
/**Although this function can be used to load 16-bit RGB or Grayscale TIFFs,
 *they will actually be converted into 8-bit images when loaded.
 *Images are always converted into 8-bit RGBA images, and then converted into
 *either 8-bit grayscale or 8-bit RGB images, so this function could use some
 *work to be made more efficient.
 *TIFF images are loaded using libtiff, which has the following notice:
 *
 * Copyright (c) 1988-1997 Sam Leffler
 * Copyright (c) 1991-1997 Silicon Graphics, Inc.
 *
 * Permission to use, copy, modify, distribute, and sell this software and 
 * its documentation for any purpose is hereby granted without fee, provided
 * that (i) the above copyright notices and this permission notice appear in
 * all copies of the software and related documentation, and (ii) the names of
 * Sam Leffler and Silicon Graphics may not be used in any advertising or
 * publicity relating to the software without the specific, prior written
 * permission of Sam Leffler and Silicon Graphics.
 * 
 * THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY 
 * WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 * IN NO EVENT SHALL SAM LEFFLER OR SILICON GRAPHICS BE LIABLE FOR
 * ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND,
 * OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
 * WHETHER OR NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF 
 * LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE 
 * OF THIS SOFTWARE.
 */
bool DImageIO::load_image_tiff(DImage *pImg, const char *stPath){
#ifdef D_NOTIFF
  fprintf(stderr, "TIFF image support is not compiled.  Need to make sure that"
	  "libtiff is available on your system, and recompile without the "
	  "D_NOTIFF definition.\n");
  return false;
#else
  //TODO: I'm currently just loading as if it were an RGBA image, and then
  // repacking the buffer.  Would be better to actually use a put method to
  // pack the data in the proper manner to begin with.  Also, I haven't
  // tried to do anything with 16-bit images at all.  I don't know if this
  // would automatically reduce the bits per pixel and lose data, or if it
  // would fail on data that is more than 8-bits per pixel.
  TIFF *tif;
//   char stErr[1025];
  uint32 w, h;
  uint16 sampsPerPxl, bitsPerSamp; // # samples per pixel, #bits per sample
  uint32 *pRGBA;
  bool fErr = false;
  D_uint8 *pDst;
  uint32 *pSrc;
  unsigned int len;

  pImg->clearProperties();
  pImg->clearComments();
  tif = TIFFOpen(stPath, "r");
  if(tif == NULL){
    fprintf(stderr,"DImageIO::load_image_tiff() failed to open '%s'\n",stPath);
    return false;
  }
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
  TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &sampsPerPxl);
  TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitsPerSamp);
//   printf("w=%d h=%d sampsPerPxl=%d(0x%02x) bitsPerSamp=%d(0x%02x)\n",
// 	 (int)w,(int)h,sampsPerPxl,sampsPerPxl,bitsPerSamp,bitsPerSamp);
  pRGBA = (uint32*)_TIFFmalloc(w*h*sizeof(uint32));
  if(!pRGBA){
    fprintf(stderr, "DImageIO::load_image_tiff() no memory!\n");
    TIFFClose(tif);
  }
  if(TIFFReadRGBAImageOriented(tif, w, h, pRGBA, ORIENTATION_TOPLEFT, 0)){
    // image was read ok, so pack it correctly
    pSrc = pRGBA;
    len = w*h;
    if(sampsPerPxl == 1){ // grayscale or black and white
      pImg->create(w,h,DImage::DImage_u8,1, DImageIO::_allocMethod);
      pDst = pImg->pData;
      for(unsigned int i = 0; i < len; ++i){
	(*pDst) = TIFFGetR(pRGBA[i]);
	++pDst;
	++pSrc;
      }
    }
    else if((sampsPerPxl == 3) || (sampsPerPxl == 4)){
      // RGB image
      pImg->create(w,h,DImage::DImage_RGB, 3, DImageIO::_allocMethod);
      pDst = pImg->pData;
      for(unsigned int i = 0; i < len; ++i){
	(*pDst) = TIFFGetR(*pSrc);
	pDst[1] = TIFFGetG(*pSrc);
	pDst[2] = TIFFGetB(*pSrc);
	pDst += 3;
	++pSrc;
      }
    }
    else{
      fprintf(stderr, "DImageIO::load_image_tiff() sampsPerPxl=%d\n",
	      sampsPerPxl);
      fErr = true;
    }
  }
  else{
    fprintf(stderr, "DImageIO::load_image_tiff() error reading TIFF image\n");
    fErr = true;
  }
  _TIFFfree(pRGBA);
  TIFFClose(tif);
  return !fErr;
#endif
}

bool DImageIO::load_image_pnm(DImage *pImg, const char *stPath){
  DImage::DFileFormat fmt;
  bool retVal = false;
  FILE *fin;

  fin = fopen(stPath, "rb");
  if(!fin){
    fprintf(stderr, "DImageIO::load_image_pnm() could't open '%s'\n",stPath);
    return false;
  }
  fmt = DImage::getImageFileFormat(stPath);
  switch(fmt){
    case DImage::DFileFormat_pbm_plain:
      retVal = load_image_pbm_plain(pImg, fin);
      break;
    case DImage::DFileFormat_pgm_plain:
      retVal = load_image_pgm_plain(pImg, fin);
      break;
    case DImage::DFileFormat_ppm_plain:
      retVal = load_image_ppm_plain(pImg, fin);
      break;
    case DImage::DFileFormat_pbm:
      retVal = load_image_pbm_raw(pImg, fin);
      break;
    case DImage::DFileFormat_pgm:
      retVal = load_image_pgm_raw(pImg, fin);
      break;
    case DImage::DFileFormat_ppm:
      retVal = load_image_ppm_raw(pImg, fin);
      break;
    default:
      fprintf(stderr, "DImageIO::load_image_pnm() '%s' isn't a PNM file\n",
	      stPath);
      retVal = false;
      break;
  }
  fclose(fin);
  return retVal;
}

bool DImageIO::load_image_gif(DImage *pImg,  const char *stPath){
#ifdef WIN32
  fprintf(stderr,"NYI! (%s:%d)\n", __FILE__, __LINE__);
  pImg->clearProperties();
  pImg->clearComments();
  return false;
#else
  char stTempPnmFile[1025];
  char stBase[1025];
  char stCmd[1025];
  pid_t procid;
  time_t curtime;
  int jj;
  bool retVal;

  // just make a system call to giftopnm and then load it as a pnm file
  //TODO: use a more robust temp file name
  jj = strlen(stPath)-1;
  while((jj > 0) && (stPath[jj]!='/'))
    --jj;
  strcpy(stBase, &stPath[jj+1]);
  procid = getpid();
  curtime = time(NULL);
  sprintf(stTempPnmFile, "/tmp/tempgif_%s_%d_%d.pnm",
	  stBase, procid,(int)curtime);
  sprintf(stCmd, "giftopnm %s > %s", stPath, stTempPnmFile);
  system(stCmd);
  retVal = load_image_pnm(pImg, stTempPnmFile);
  unlink(stTempPnmFile);
  return retVal;
#endif
}
bool DImageIO::load_image_bmp(DImage *pImg, const char *stPath){
#ifdef WIN32
  fprintf(stderr,"NYI! (%s:%d)\n", __FILE__, __LINE__);
  pImg->clearProperties();
  pImg->clearComments();
  return false;
#else
  char stTempPnmFile[1025];
  char stBase[1025];
  char stCmd[1025];
  pid_t procid;
  time_t curtime;
  int jj;
  bool retVal;
  // just make a system call to giftopnm and then load it as a pnm file
  //TODO: use a more robust temp file name
  jj = strlen(stPath)-1;
  while((jj > 0) && (stPath[jj]!='/'))
    --jj;
  strcpy(stBase, &stPath[jj+1]);
  procid = getpid();
  curtime = time(NULL);
  sprintf(stTempPnmFile, "/tmp/tempgif_%s_%d_%d.pnm",
	  stBase, procid, (int)curtime);
  sprintf(stCmd, "bmptopnm %s > %s", stPath, stTempPnmFile);
  system(stCmd);
  retVal = load_image_pnm(pImg, stTempPnmFile);
  unlink(stTempPnmFile);
  return retVal;
#endif
}
bool DImageIO::load_image_gnuplot(DImage *pImg, FILE *fin){
  fprintf(stderr,"NYI! (%s:%d)\n", __FILE__, __LINE__);
  pImg->clearProperties();
  pImg->clearComments();
  return false;
}

bool DImageIO::save_image_pnm(DImage *pImg, const char *stPath,
			      bool fSaveProps, bool fSaveComments){
  switch(pImg->_imgType){
    case DImage::DImage_u8:
    case DImage::DImage_u16:
    case DImage::DImage_u32:
      return DImageIO::save_image_pgm_raw(pImg, stPath,
					  fSaveProps, fSaveComments);
    case DImage::DImage_RGB:
    case DImage::DImage_RGB_16:
      return DImageIO::save_image_ppm_raw(pImg, stPath, fSaveProps);
    default:
      fprintf(stderr,
	      "DImageIO::save_image_pnm() unsupported image type for PNM\n");
#ifdef DEBUG
      exit(1);
#endif
      return false;
  }
  return false;
}

bool DImageIO::save_image_ppm_raw(DImage *pImg, const char *stPath,
				  bool fSaveProps, bool fSaveComments){
  FILE *fout;
  DImage imgTmp;
  size_t bufLen = 0;

  switch(pImg->_imgType){
    case DImage::DImage_u8:
      pImg->convertedImgType_(imgTmp, DImage::DImage_RGB, 3);
      pImg = &imgTmp; // pointer is only changed inside this function
      break;
    case DImage::DImage_u16:
      pImg->convertedImgType_(imgTmp, DImage::DImage_RGB_16, 3);
      pImg = &imgTmp; // pointer is only changed inside this function
      break;
    case DImage::DImage_RGB:
    case DImage::DImage_RGB_16:
      // do nothing because already in the correct format
      break;
    default:
      fprintf(stderr,"DImageIO::save_image_ppm_raw() wrong image type. use "
	      "convertedImgType_() to convert it to RGB or RGB_16\n");
#ifdef DEBUG
      exit(1);
#endif
      return false;
  } // end switch


  // open the output file
  fout = fopen(stPath, "wb");
  if(!fout){
    fprintf(stderr,
	    "DImageIO::save_image_ppm_raw() couldn't open '%s' for writing\n",
	    stPath);
#ifdef DEBUG
      exit(1);
#endif
    return false;
  }
  // write the header
  fprintf(fout,"P6\n"); // magic number
  if(fSaveProps){
    if(0 != write_img_props_pnm(pImg, fout)){
      fprintf(stderr,
	      "DImageIO::save_image_ppm_raw() err writing properties\n");
      fclose(fout);
#ifdef DEBUG
      exit(1);
#endif
      return false;
    }
  }
  if(fSaveComments){
    if(0 != write_img_comments_pnm(pImg, fout)){
      fprintf(stderr,
	      "DImageIO::save_image_ppm_raw() err writing comments\n");
#ifdef DEBUG
      exit(1);
#endif
      fclose(fout);
      return false;
    }
  }
  fprintf(fout,"%d %d\n", pImg->_w, pImg->_h);
  if(DImage::DImage_RGB == pImg->_imgType){
    fprintf(fout,"255\n");
    bufLen = pImg->_w * pImg->_h * 3;
  }
  else if(DImage::DImage_RGB_16 == pImg->_imgType){
    fprintf(fout,"65535\n");
    bufLen = pImg->_w * pImg->_h * 6;
    fprintf(stderr,
	    "WARNING: DImageIO::save_image_ppm_raw(%s) not swapping bytes for "
	    "16-bit image.  Image may not be read right on other machines.\n",
	    stPath);
    //TODO: if two-byte data and needs swapped bytes, then do the swap and
    // swap back afterwards (unless temp image, then don't bother to swap back)
  }
  // write the block of data
  if(!writeDataBlock(fout, pImg->pData, bufLen)){
    fprintf(stderr, "DImageIO::save_image_ppm_raw() couldn't write data "
	    "for '%s'\n", stPath);
    fclose(fout);
#ifdef DEBUG
      exit(1);
#endif
    return false;
  }
  fclose(fout);
  return true;
}
bool DImageIO::save_image_pgm_raw(DImage *pImg, const char *stPath,
				  bool fSaveProps, bool fSaveComments){
  FILE *fout;
  DImage imgTmp;
  size_t bufLen = 0;

  switch(pImg->_imgType){
    case DImage::DImage_u8:
    case DImage::DImage_u16:
      // do nothing because already in the correct format
      break;
    case DImage::DImage_u32:
#ifdef DEBUG
      fprintf(stderr, "WARNINIG: DImageIO::save_image_pgm_raw(%s) converting "
	      "32-bit image to 16-bit image for output.  Image may not be "
	      "accurate\n",stPath);
#endif
      pImg->convertedImgType_(imgTmp, DImage::DImage_u16, 1);
      pImg = &imgTmp; // pointer is only changed inside this function
      break;
    case DImage::DImage_RGB:
      pImg->convertedImgType_(imgTmp, DImage::DImage_u8, 1);
      pImg = &imgTmp; // pointer is only changed inside this function
      break;
    case DImage::DImage_RGB_16:
      pImg->convertedImgType_(imgTmp, DImage::DImage_u16, 1);
      pImg = &imgTmp; // pointer is only changed inside this function
      break;
    default:
      fprintf(stderr,"DImageIO::save_image_pgm_raw() wrong image type. use "
	      "convertedImgType_() to convert it to DImage_u8 or _u16\n");
#ifdef DEBUG
      exit(1);
#endif
      return false;
  } // end switch

  // open the output file
  fout = fopen(stPath, "wb");
  if(!fout){
    fprintf(stderr,
	    "DImageIO::save_image_pgm_raw() couldn't open '%s' for writing\n",
	    stPath);
#ifdef DEBUG
      exit(1);
#endif
    return false;
  }
  // write the header
  fprintf(fout,"P5\n"); // magic number
  if(fSaveProps){
    if(0 != write_img_props_pnm(pImg, fout)){
      fprintf(stderr,
	      "DImageIO::save_image_pgm_raw() err writing properties\n");
      fclose(fout);
#ifdef DEBUG
      exit(1);
#endif
      return false;
    }
  }
  if(fSaveComments){
    if(0 != write_img_comments_pnm(pImg, fout)){
      fprintf(stderr,
	      "DImageIO::save_image_pgm_raw() err writing comments\n");
#ifdef DEBUG
      exit(1);
#endif
      fclose(fout);
      return false;
    }
  }
  fprintf(fout,"%d %d\n", pImg->_w, pImg->_h);
  if(DImage::DImage_u8 == pImg->_imgType){
    fprintf(fout,"255\n");
    bufLen = pImg->_w * pImg->_h;
  }
  else if(DImage::DImage_u16 == pImg->_imgType){
    fprintf(fout,"65535\n");
    bufLen = pImg->_w * pImg->_h *sizeof(D_uint16);
    fprintf(stderr,
	    "WARNING: DImageIO::save_image_ppm_raw() not swapping bytes for "
	    "16-bit image.  Image may not be read right on other machines.\n");
    //TODO: if two-byte data and needs swapped bytes, then do the swap and
    // swap back afterwards (unless temp image, then don't bother to swap back)
  }
  else{
    fprintf(stderr, "DImageIO::save_image_pgm_raw() err: pImg->_imgType=%d\n",
	    pImg->_imgType);
  }
  // write the block of data
  if(!writeDataBlock(fout, pImg->pData, bufLen)){
    fprintf(stderr, "DImageIO::save_image_pgm_raw() couldn't write data "
	    "for '%s'\n", stPath);
    fclose(fout);
#ifdef DEBUG
      exit(1);
#endif
    return false;
  }
  fclose(fout);
  return true;
}
bool DImageIO::save_image_pbm_raw(DImage *pImg, const char *stPath,
				  bool fSaveProps, bool fSaveComments){
  fprintf(stderr,"NYI! (%s:%d)\n", __FILE__, __LINE__);
#ifdef DEBUG
      exit(1);
#endif
  return false;
}

int DImageIO::write_img_props_pnm(DImage *pImg, FILE *fout){
  std::map<const std::string, std::string>::iterator iter;

  iter = pImg->mapProps.begin();
  while(iter != pImg->mapProps.end()){
    fprintf(fout, "#DProp:%s=%s\n",
	    (*iter).first.c_str(),(*iter).second.c_str());
    ++iter;
  }
  return 0;
}
int DImageIO::write_img_comments_pnm(DImage *pImg, FILE *fout){
  std::vector<std::string>::iterator iter;

  iter = pImg->vectComments.begin();
  while(iter != pImg->vectComments.end()){
    fprintf(fout, "#%s\n", (*iter).c_str());
    ++iter;
  }
  return 0;
}

///Set the allocation method to be used to create the data buffers when loading
/** By default, AllocationMethod_malloc is used.  If you need your
 *  images to be aligned to boundaries larger than those guaranteed by
 *  malloc(), allocMeth should be called with
 *  AllocationMethod_daligned, and the DImage::setAlignment() method
 *  should be used to set the number of bytes to align to.  This may
 *  be the case when using code that manipulates your images by using
 *  SIMD instructions (AltiVec, SSE, SSE2, etc.), for example.
 */
void DImageIO::set_alloc_method(D_AllocationMethod allocMeth){
  if((allocMeth >= AllocationMethod_daligned) &&
     (allocMeth <= AllocationMethod_new)){
    DImageIO::_allocMethod = allocMeth;
  }
  else
    fprintf(stderr, "DImageIO::set_alloc_method() invalid allocMethod\n");
}

///save the image in a format that can be viewed with gnuplot's splot command
bool DImageIO::save_image_gnuplot(DImage *pImg, const char *stPath,
				  bool fSaveProps, bool fSaveComments,
				  int everyNthPixel){
  FILE *fout;
  DImage imgTmp;
  double *pDbl;
  int w, h;

  if(DImage::DImage_dbl_multi != pImg->_imgType){
    switch(pImg->_imgType){
      case DImage::DImage_u8:
      case DImage::DImage_u16:
      case DImage::DImage_u32:
      case DImage::DImage_flt_multi:
	pImg->convertedImgType_(imgTmp, DImage::DImage_dbl_multi, 1);
	pImg = &imgTmp; // pointer is only changed inside this function
	break;
      default:
	fprintf(stderr,"DImageIO::save_image_gnuplot() wrong image type. use "
		"convertedImgType_() to convert it to single-channel image\n");
#ifdef DEBUG
      exit(1);
#endif
	return false;
    } // end switch
  }

  // open the output file
  fout = fopen(stPath, "wb");
  if(!fout){
    fprintf(stderr,
	    "DImageIO::save_image_gnuplot() couldn't open '%s' for writing\n",
	    stPath);
#ifdef DEBUG
      exit(1);
#endif
    return false;
  }
  // write the header
  fprintf(fout,"#Document Project gnuplot output\n");
  fprintf(fout,"#to view, use the gnuplot splot command.  for example:\n");
  fprintf(fout,"#splot \"<filename>\"\n");
  fprintf(fout,"#or: splot \"<filename>\" with image\n");
  fprintf(fout,"#for large files, use only every 10th point (in the file, not both dimensions): "
	  "splot \"<file>\" every 10 with image\n");
  fprintf(fout,"#or, you could save the file using DImageIO::save_image_gnuplot() with everyNthPixel set to something other than 1 to do both dimensions\n");
  fprintf(fout,"#\n");
  if(fSaveProps){
    if(0 != write_img_props_pnm(pImg, fout)){
      fprintf(stderr,
	      "DImageIO::save_image_gnuplot() err writing properties\n");
      fclose(fout);
#ifdef DEBUG
      exit(1);
#endif
      return false;
    }
  }
  if(fSaveComments){
    if(0 != write_img_comments_pnm(pImg, fout)){
      fprintf(stderr,
	      "DImageIO::save_image_gnuplot() err writing comments\n");
      fclose(fout);
#ifdef DEBUG
      exit(1);
#endif
      return false;
    }
  }
  w = pImg->_w;
  h = pImg->_h;
  fprintf(fout,"#width=%d height=%d\n", w, h);
  pDbl = pImg->dataPointer_dbl(0);
  for(int y = 0; y < h; y+=everyNthPixel){
    for(int x = 0; x < w; x+=everyNthPixel){
      fprintf(fout, "%d %d %.12f\n", x, y, pDbl[y*w+x]);
      //      ++pDbl;
    }
    fprintf(fout, "\n");
  }

  fclose(fout);
  return true;
}



///save an image in PNG format
/**This function uses the PNG Reference Library (libpng), and also
 *follows the code examples provided with that library for its
 *implementation
 */
bool DImageIO::save_image_png(DImage *pImg, const char *stPath,
			      bool fSaveProps, bool fSaveComments,
			      int compression){
#ifdef D_NOPNG
  fprintf(stderr, "PNG image support is not compiled.  Need to make sure that"
	  "libpng is available on your system, and recompile without the "
	  "D_NOPNG definition.\n");
  return false;
#else
  FILE *fout;
  png_uint_32 width, height;
  int bit_depth;
  int color_type;
  png_bytep *row_pointers;
  png_structp png_ptr;
  png_infop info_ptr;
  int numComments;
  png_text *rgText=NULL;
  char stKey[80];
  int w, h; // image width, height
  D_uint8 *p8;
  int bytesPerSample;

  sprintf(stKey, "Comment");
//   png_unknown_chunkp user_chunk_ptr;

  switch(pImg->getImageType()){
    case DImage::DImage_u8:
      bit_depth = 8;
      color_type = PNG_COLOR_TYPE_GRAY;
      bytesPerSample = 1;
      break;
    case DImage::DImage_u16:
      bit_depth = 16;
      color_type = PNG_COLOR_TYPE_GRAY;
      bytesPerSample = 2;
      break;
    case DImage::DImage_RGB:
      bit_depth = 8;
      color_type = PNG_COLOR_TYPE_RGB;
      bytesPerSample = 3;
      break;
    case DImage::DImage_RGB_16:
      bit_depth = 16;
      color_type = PNG_COLOR_TYPE_RGB;
      bytesPerSample = 6;
      break;
    default:
      fprintf(stderr, "DImageIO::save_image_png() does not support all image "
	      "types\n");
      return false;
  }

  fout = fopen(stPath, "wb");
  if(!fout){
    fprintf(stderr,
	    "DImageIO::save_image_png() couldn't open '%s' for writing\n",
	    stPath);
#ifdef DEBUG
    //      exit(1);
#endif
    return false;
  }

  //set up row pointers
  w = pImg->width();
  h = pImg->height();
  if((h < 1) || (w < 1)){
    fprintf(stderr, "DImageIO::save_image_png() width or height is zero\n");
    fclose(fout);
#ifdef DEBUG
    //    exit(1);
#endif
    return false;
  }
  
  row_pointers = new png_bytep[h];
  if(!row_pointers){
    fprintf(stderr,
	    "DImageIO::save_image_png() could not allocate row pointers\n");
    fclose(fout);
    return false;
  }
  p8 = pImg->dataPointer_u8();
  for(int r = 0; r < h; ++r){
    row_pointers[r] = (png_bytep)(&(p8[bytesPerSample*w*r]));
  }

  // set up PNG structures and error handler
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if(!png_ptr){
    delete [] row_pointers;
    fclose(fout);
    return false;
  }

  info_ptr = png_create_info_struct(png_ptr);
  if(!info_ptr){
    png_destroy_write_struct(&png_ptr, &info_ptr);
    delete [] row_pointers;
    fclose(fout);
    return false;
  }

//attempt to fix undef long_jmp problem
#if (PNG_LIBPNG_VER < 10400 || PNG_LIBPNG_VER >= 10500)
  if(setjmp(png_jmpbuf(png_ptr))){
#else
  if (setjmp(png_ptr->jmpbuf)){
#endif
    png_destroy_write_struct(&png_ptr, &info_ptr);
    delete [] row_pointers;
    fclose(fout);
    return false;
  }

  //tell the png library what file to write to
  png_init_io(png_ptr, fout);

  //set options and header information
  png_set_compression_level(png_ptr, compression);
  
  width = w;
  height = h;

  png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth, color_type,
	       PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
	       PNG_FILTER_TYPE_DEFAULT);

  // set the comment text
  numComments = pImg->getNumComments();
  if(fSaveComments){
    if(numComments > 0){
      rgText = new png_text[numComments];
      if(rgText){
	for(int cnum = 0; cnum < numComments; ++cnum){
	  std::string strTmp;
	  rgText[cnum].compression = PNG_TEXT_COMPRESSION_NONE;
	  rgText[cnum].key = stKey;
	  strTmp = pImg->getCommentByIndex(cnum);
	  rgText[cnum].text_length = strTmp.length();
	  rgText[cnum].text = new char[1+rgText[cnum].text_length];
	  if(NULL != rgText[cnum].text){
	    strncpy(rgText[cnum].text, strTmp.c_str(),
		    rgText[cnum].text_length);
	    rgText[cnum].text[rgText[cnum].text_length] = '\0';
	  }
#ifdef PNG_iTXt_SUPPORTED
	  rgText[cnum].itxt_length = 0;
	  rgText[cnum].lang = NULL;
	  rgText[cnum].lang_key = NULL;
#endif
	}
	png_set_text(png_ptr, info_ptr, rgText, numComments);
      }
      else{
	numComments=0;
      }
    }
  }

  //write the header information
  png_write_info(png_ptr, info_ptr);

  // for little-endian machines, make sure data bytes are swapped correctly:
  if(D_LITTLE_ENDIAN && (bit_depth > 8))
    png_set_swap(png_ptr);
  

  // now write the data
  png_write_image(png_ptr, row_pointers);

  // finish the sequential write and close the output file
  png_write_end(png_ptr, info_ptr);
  fclose(fout);

  //clean up the memory that has been allocated
  if(fSaveComments){
    if(numComments > 0){
      if(rgText){
	for(int cnum = 0; cnum < numComments; ++cnum){
	  rgText[cnum].key= NULL; //it's pointing to stKey on the stack
	  if(NULL != rgText[cnum].text)
	    delete [] rgText[cnum].text;
	}
	delete [] rgText;
      }
    }
  }
  delete [] row_pointers;
  png_destroy_write_struct(&png_ptr, &info_ptr);

  return true;
#endif
}













///save an image in JPEG format
/**"This software is based in part on the work of the Independent JPEG
 * Group".  Code in this function was created by following
 * instructions in libjpeg.doc, as well as copying examples found in
 * cjpeg.c, example.c, and jdatasrc.c in the IJG's distribution of
 * libjpeg.
*/
bool DImageIO::save_image_jpeg(DImage *pImg, const char *stPath,
			       bool fSaveProps, bool fSaveComments,
			       int quality, bool fProgressive, bool fOptimize){
#ifdef D_NOJPEG
  fprintf(stderr, "JPEG image support is not compiled.  Need to make sure that"
	  "libjpeg is available on your system, and recompile without the "
	  "D_NOJPEG definition.\n");
  return false;
#else
  FILE *fout;
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPROW row_pointer[1];	/* pointer to JSAMPLE row[s] */
  int row_stride;		/* physical row width in image buffer */
  int idx;
  D_uint8 *p8;
  J_COLOR_SPACE in_color_space = JCS_GRAYSCALE;
  int input_components = 0;
  int numComments;
  int numProps;

  if(DImage::DImage_u8 == pImg->getImageType()){
    input_components = 1;
    in_color_space = JCS_GRAYSCALE;
  }
  else if(DImage::DImage_RGB == pImg->getImageType()){
    input_components = 3;
    in_color_space = JCS_RGB;
  }
  else{
    fprintf(stderr, "DImageIO::save_image_jpeg() only handles "
	    "8-bit RGB or 8-bit grayscale images\n");
    return false;
  }
  if((quality < 0) || (quality > 100)){
    fprintf(stderr, "DImageIO::save_image_jpeg() quality must be 0..100 (was "
	    "%d)\n", quality);
    return false;
  }
  // use default jpeg error handler (if error, print to stderr then exit)
  cinfo.err = jpeg_std_error(&jerr);
  // create compressor
  jpeg_create_compress(&cinfo);
  //open output file
  fout = fopen(stPath, "wb");
  if(!fout){
    fprintf(stderr, "DImageIO::save_image_jpeg() can't open '%s'\n",stPath);
    return false;
  }
  //set output file for compressor
  jpeg_stdio_dest(&cinfo, fout);
  //set compression parameters
  cinfo.image_width = pImg->width(); //required
  cinfo.image_height = pImg->height(); //required
  cinfo.input_components = input_components; //required
  cinfo.in_color_space = in_color_space; //required
  //set default parms (cinfo.in_color_space must already be set)
  jpeg_set_defaults(&cinfo);
  //now set non-default parameters
  jpeg_set_quality(&cinfo, quality, TRUE);
  if(fProgressive){//progressive image
#ifdef C_PROGRESSIVE_SUPPORTED
    jpeg_simple_progression(&cinfo);
#else
    fprintf(stderr, "DImageIO::save_image_jpeg() C_PROGRESSIVE_SUPPORTED is "
	    "not enabled in the JPEG library\n");
#endif
  }
  if(fOptimize){
    cinfo.optimize_coding = TRUE;
  }
  // start compressor
  jpeg_start_compress(&cinfo, TRUE);

  // add comments, if desired
  // (do it after jpeg_start_compress(), but before jpeg_write_scanlines() )
  numComments = pImg->getNumComments();
  if(fSaveComments){
    if(numComments > 0){
      for(int cnum = 0; cnum < numComments; ++cnum){
	  std::string strTmp;
	  strTmp = pImg->getCommentByIndex(cnum);
	  strTmp += "\n";
	  jpeg_write_marker(&cinfo, JPEG_COM,
			    (const unsigned char*)strTmp.c_str(),
			    strlen(strTmp.c_str()));
      }
    }
  }
  //add props
  numProps = pImg->getNumProperties();
  if(fSaveProps){
    std::map<const std::string, std::string>::iterator iter;
    if(numProps > 0){
      iter = pImg->mapProps.begin();
      while(iter != pImg->mapProps.end()){
	std::string strTmp;
	strTmp = std::string("DProp:");
	strTmp += (*iter).first + std::string("=") + (*iter).second +
	  std::string("\n");
	jpeg_write_marker(&cinfo, JPEG_COM,
			  (const unsigned char*)strTmp.c_str(),
			  strlen(strTmp.c_str()));
	++iter;
      }
    }
  }

  //compress the image
  row_stride=pImg->width() * input_components;
  idx = 0;
  p8 = pImg->dataPointer_u8();
  while(cinfo.next_scanline < cinfo.image_height){
    row_pointer[0] = &p8[row_stride*cinfo.next_scanline];
    (void)jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
  // finish the compression
  jpeg_finish_compress(&cinfo);
  //close the output file
  fclose(fout);
  //release the compressor
  jpeg_destroy_compress(&cinfo);

  return true;
#endif
}









//check beginning of file for header info about w,h, num channels in image
/**If you do not want any particular value, set the pointer for that parameter
 * to NULL.  Otherwise, the variable pointed to will be filled with the 
 * corresponding information.  The same copyright/license notices apply to
 * the code used in the implementation of this function as apply for the 
 * various image loading functions for PNG, TIFF, and JPEG images.
 */
bool DImageIO::get_image_width_height_chan(const char *stPath,
					   int *width, int *height, int *chan){
  DImage::DFileFormat fmt;
  D_PNM_HEADER_S pnmhdr;
  bool fHeaderOK;
  FILE *fin;

  fmt = DImage::getImageFileFormat(stPath);
  switch(fmt){
    case DImage::DFileFormat_pbm_plain:
    case DImage::DFileFormat_pgm_plain:
    case DImage::DFileFormat_ppm_plain:
    case DImage::DFileFormat_pbm:
    case DImage::DFileFormat_pgm:
    case DImage::DFileFormat_ppm:
    case DImage::DFileFormat_pnm:
      fin = fopen(stPath, "rb");
      if(!fin){
	fprintf(stderr, "DImageIO::get_image_width_height_chan() couldn't open"
		" '%s'\n", stPath);
	return false;
      }
      fHeaderOK = DImageIO::read_pnm_header(&pnmhdr, fin);
      DImageIO::clear_pnm_header_comments(&pnmhdr);
      if(NULL != width)
	(*width) = pnmhdr.w;
      if(NULL != height)
	(*height) = pnmhdr.h;
      if(NULL != chan){
	if((fmt == DImage::DFileFormat_ppm_plain) ||
	   (fmt == DImage::DFileFormat_ppm))
	  (*chan) = 3;
	else
	  (*chan) = 1;
      }
      fclose(fin);
      return true;
      break;
    case DImage::DFileFormat_jpg:
#ifdef D_NOJPEG
      fprintf(stderr,
	      "JPEG image support is not compiled.  Need to make sure that"
	      "libjpeg is available on your system, and recompile without the "
	      "D_NOJPEG definition.\n");
      return false;
#else
      struct jpeg_decompress_struct cinfo;
      struct jpeg_error_mgr jerr;

      if(1 != sizeof(JSAMPLE)){
	fprintf(stderr,
		"DImageIO::get_image_width_height_chan() expected 8-bit "
		"JSAMPLE size. Other sizes not supported\n");
	return false;
      }
      fin = fopen(stPath, "rb");
      if(!fin){
	fprintf(stderr, "DImageIO::get_image_width_height_chan() couldn't open"
		" '%s'\n", stPath);
	return false;
      }
      cinfo.err = jpeg_std_error(&jerr);
      jpeg_create_decompress(&cinfo);
      jpeg_stdio_src(&cinfo, fin);
      if(JPEG_HEADER_OK != jpeg_read_header(&cinfo, boolean(true))){
	fprintf(stderr, "DImageIO::get_image_width_height_chan() error in jpeg"
		" header\n");
	jpeg_abort_decompress(&cinfo);
	return false;
      }
      jpeg_start_decompress(&cinfo); // reads the header, etc.
      // now we have width and height, so create the image buffer, and pointers
      // to each row in the buffer
      if(NULL != width)
	(*width) = cinfo.output_width;
      if(NULL != height)
	(*height) = cinfo.output_height;
      if(NULL != chan){
	if(3 == cinfo.out_color_components)
	  (*chan) = 3;
	else if(1 == cinfo.out_color_components)
	  (*chan) = 1;
	else{
	  fprintf(stderr, "DImageIO::get_image_width_height_chan() "
		  "unexpected number of channels\n");
	  (*chan) = 0;
	  jpeg_abort_decompress(&cinfo);
	  jpeg_destroy_decompress(&cinfo);
	  fclose(fin);
	  return false;
	}
      }
      jpeg_abort_decompress(&cinfo);
      jpeg_destroy_decompress(&cinfo);
      fclose(fin);
      return true;
#endif
      break;
    case DImage::DFileFormat_png:
#ifdef D_NOPNG
      fprintf(stderr,
	      "PNG image support is not compiled.  Need to make sure that"
	      "libpng is available on your system, and recompile without the "
	      "D_NOPNG definition.\n");
      return false;
#else
      png_structp png_ptr;
      png_infop info_ptr;

      png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL,NULL);
      if(!png_ptr)
	return false;

      info_ptr = png_create_info_struct(png_ptr);
      if(!info_ptr){
	png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
	return false;
      }


      fin = fopen(stPath, "rb");
      if(!fin){
	fprintf(stderr, "DImageIO::get_image_width_height_chan() couldn't open"
		" '%s'\n", stPath);
	return false;
      }
    
      png_init_io(png_ptr, fin);
      png_read_info(png_ptr, info_ptr);
      
      if(NULL != width)
	(*width) = png_get_image_width(png_ptr, info_ptr);
      if(NULL != height)
	(*height) = png_get_image_height(png_ptr, info_ptr);
      if(NULL != chan){
	(*chan) = png_get_channels(png_ptr, info_ptr);
	if(2==(*chan)){ // ignore alpha channel
	  (*chan) = 1;
	}
	if(4==(*chan)){ // ignore alpha channel
	  (*chan) = 3;
	}
      }
      fclose(fin);
      png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);

      return true;
#endif
      break;
    case DImage::DFileFormat_tiff:
#ifdef D_NOTIFF
      fprintf(stderr,
	      "TIFF image support is not compiled.  Need to make sure that"
	      "libtiff is available on your system, and recompile without the "
	      "D_NOTIFF definition.\n");
      return false;
#else
      TIFF *tif;
      uint32 w, h;
      uint16 sampsPerPxl, bitsPerSamp; // # samples per pixel, #bits per sample

      tif = TIFFOpen(stPath, "r");
      if(tif == NULL){
	fprintf(stderr, "DImageIO::get_image_width_height_chan() couldn't open"
		" '%s'\n", stPath);
	return false;
      }
      TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
      TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
      TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &sampsPerPxl);
      TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitsPerSamp);
      TIFFClose(tif);
      if(NULL != width)
	(*width) = w;
      if(NULL != height)
	(*height) = h;
      if(NULL != chan){
	(*chan) = sampsPerPxl;
      }
      return true;
#endif
      break;
    default:
      fprintf(stderr, "DImageIO::get_image_width_height_chan() currently only "
	      "supports pnm, jpg, png, and tiff images.\n");
      return false;
  }
  fprintf(stderr, "DImageIO::get_image_width_height_chan() unexpected end of "
	  "function.\n");
  exit(1);
  return false;
}
