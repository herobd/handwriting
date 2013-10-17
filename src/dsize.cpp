#include "dsize.h"
#include "dinstancecounter.h"

DSize::DSize(){
  DInstanceCounter::addInstance("DSize");
  w = h = 0;
}
DSize::DSize(int w, int h){
  DInstanceCounter::addInstance("DSize");
  this->w = w;
  this->h = h;
}
DSize::DSize(const DSize &src){
  DInstanceCounter::addInstance("DSize");
  w = src.w;
  h = src.h;
}
DSize::~DSize(){
  DInstanceCounter::removeInstance("DSize");
}



