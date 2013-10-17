#include "dconnectedcompinfo.h"
#include "dinstancecounter.h"

DConnectedComponentInfo::DConnectedComponentInfo(){
  DInstanceCounter::addInstance("DConnectedComponentInfo");
  label=0xffffffff;
  bbLeft = -1;
  bbTop = -1;
  bbRight = -1;
  bbBottom = -1;
  startX = -1;
  startY = -1;
  pixels = 0;
  centroidX = -1;
  centroidY = -1;
}

DConnectedComponentInfo::~DConnectedComponentInfo(){
  DInstanceCounter::removeInstance("DConnectedComponentInfo");
}
