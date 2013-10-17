#include "drect.h"
#include "dinstancecounter.h"
#include "dpoint.h"

DRect::DRect(){
  DInstanceCounter::addInstance("DRect");
  x = y = w = h = 0;
}
DRect::DRect(int x, int y, int w, int h){
  DInstanceCounter::addInstance("DRect");
  this->x = x;
  this->y = y;
  this->w = w;
  this->h = h;
}
DRect::DRect(const DPoint &p0, DPoint &p1){
  if(p0.x <= p1.x){
    x = p0.x;
    w = p1.x-p0.x+1;
  }
  else{
    x = p1.x;
    w = p0.x-p1.x+1;
  }
  if(p0.y <= p1.y){
    y = p0.y;
    h = p1.y-p0.y+1;
  }
  else{
    y = p1.y;
    h = p0.y-p1.y+1;
  }
}
DRect::DRect(const DRect &src){
  DInstanceCounter::addInstance("DRect");
  x = src.x;
  y = src.y;
  w = src.w;
  h = src.h;
}
DRect::~DRect(){
  DInstanceCounter::removeInstance("DRect");
}

