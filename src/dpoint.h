#ifndef DPOINT_H
#define DPOINT_H

class DPoint{
public:

  DPoint();
  DPoint(int x, int y);
  DPoint(const DPoint &src);
  ~DPoint();
  int x;
  int y;
};

#endif
