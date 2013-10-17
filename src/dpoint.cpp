#include "dpoint.h"
#include "dinstancecounter.h"

DPoint::DPoint(){
  DInstanceCounter::addInstance("DPoint");
  x = y = 0;
}
DPoint::DPoint(int x, int y){
  DInstanceCounter::addInstance("DPoint");
  this->x = x;
  this->y = y;
}
DPoint::DPoint(const DPoint &src){
  DInstanceCounter::addInstance("DPoint");
  x = src.x;
  y = src.y;
}
DPoint::~DPoint(){
  DInstanceCounter::removeInstance("DPoint");
}



