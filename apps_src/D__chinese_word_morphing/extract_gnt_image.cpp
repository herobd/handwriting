//g++ -Wall -O3 -g -fPIC extract_gnt_image.cpp -o extract_gnt_image -I../../../familysearch_documentproject_2013.06.21/src -L../../../familysearch_documentproject_2013.06.21/lib/ -ldocumentproj_2013.06.21 -ljpeg -ltiff -lpng -lpthread -lm

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <stdint.h>
#include "dimage.h"
#include "dthresholder.h"

void loadImageFromGnt(const char *stGntName, DImage &imgGS, DImage &imgThresh,
		     int imgNumToRead, std::string &strLabel){
  FILE *fin;
  fin = fopen(stGntName, "rb");
  if(!fin){
    fprintf(stderr,"file '%s' could not be opened for reading\n",stGntName);
    exit(1);
  }
  size_t imgNum;
  imgNum = 0;

  while(!feof(fin)){
    uint32_t blockSize32;
    uint16_t tagCode16;
    uint16_t width16;
    uint16_t height16;
    unsigned char *p8;

    if(1!=fread(&blockSize32, sizeof(uint32_t), 1, fin)){
      //we must have reached end of file
      break;
    }
    if(1!=fread(&tagCode16, sizeof(uint16_t), 1, fin)){
      fprintf(stderr,"an error occurred reading the tagcode for image #%ld\n",
	      imgNum);
      exit(1);
    }
    if(1!=fread(&width16, sizeof(uint16_t), 1, fin)){
      fprintf(stderr,"an error occurred reading the width for image #%ld\n",
	      imgNum);
      exit(1);
    }
    if(1!=fread(&height16, sizeof(uint16_t), 1, fin)){
      fprintf(stderr,"an error occurred reading the height for image #%ld\n",
	      imgNum);
      exit(1);
    }
    //printf("image %ld w=%d h=%d\n",imgNum, (int)width16, (int)height16);
    size_t datalen;
    datalen = (unsigned int)width16 * (unsigned int)height16;

    DImage img;
    img.create((int)width16,(int)height16, DImage::DImage_u8, 1);
    p8 = img.dataPointer_u8();
    if(1!=fread(p8, datalen, 1, fin)){
      fprintf(stderr,"error reading image data\n");
      exit(1);
    }
    if((int)imgNum==imgNumToRead){
      imgGS = img;
      //find threshold using otsu
      DThresholder::otsuThreshImage_(imgThresh, img);
      char stLabelTmp[20];
      sprintf(stLabelTmp,"%03u%03u",
	      (unsigned char)tagCode16, (unsigned char)(tagCode16>>8));
      strLabel = std::string(stLabelTmp);
      imgGS.setProperty(std::string("tagcode"),strLabel);
      imgThresh.setProperty(std::string("tagcode"),strLabel);
    }
    ++imgNum;
    if((int)imgNum>imgNumToRead)
      break;
  }
  fclose(fin);
}



int main(int argc, char **argv){
  DImage img, imgT;
  int imgNum;
  std::string strLabel;

  if(argc!=5){
    fprintf(stderr,
	    "usage: %s <input.gnt> <imgNum> <output.pgm> <outputThresh.pgm>\n",
	    argv[0]);
    fprintf(stderr,"  imgNum is 0-based and is the index of the image within "
	    "the gnt file\n  output.pgm is the extracted grayscale image\n"
	    "  outputThresh.pgm is the non-255 pixels are black\n"
	    " *To not output one of the files, use /dev/null as its name\n");
    exit(1);
  }
  imgNum = atoi(argv[2]);
  if((imgNum<0)||(imgNum>5000)){
    fprintf(stderr,"expected imgNum to be from 0 to 5000 (was %d)\n",imgNum);
    exit(1);
  }
  loadImageFromGnt(argv[1], img, imgT, imgNum, strLabel);
  img.save(argv[3]);
  imgT.save(argv[4]);
  return 0;
}
