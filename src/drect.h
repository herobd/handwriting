#ifndef DRECT_H
#define DRECT_H

#include "dpoint.h"

class DRect{
public:

  DRect();
  DRect(int x, int y, int w, int h);
  DRect(const DPoint &p0, DPoint &p1);
  DRect(const DRect &src);
  ~DRect();
  bool containsPoint(int x, int y) const;
  bool containsPoint(const DPoint &p) const;
  bool intersectsWith(const DRect &r2) const;
  int x;
  int y;
  int w;
  int h;
};


inline bool DRect::containsPoint(const int x, const int y) const{
  return ( (x>=this->x) && (x<(this->x + w)) &&
	   (y>=this->y) && (y<(this->y + h)) );
}

inline bool DRect::containsPoint(const DPoint &p) const{
  return containsPoint(p.x, p.y);
}

inline bool DRect::intersectsWith(const DRect &r2) const{
  return ( (r2.containsPoint(this->x,this->y)) ||
	   (this->containsPoint(r2.x,r2.y)) ||
	   (r2.containsPoint(this->x+this->w-1,this->y)) ||
	   (r2.containsPoint(this->x,this->y+this->h-1)) ||
	   (r2.containsPoint(this->x+this->w-1,this->y+this->h-1)) ||
	   (this->containsPoint(r2.x+r2.w-1, r2.y)) ||
	   (this->containsPoint(r2.x, r2.y+r2.h-1)) ||
	   (this->containsPoint(r2.x+r2.w-1, r2.y+r2.h-1)));
}

#endif
