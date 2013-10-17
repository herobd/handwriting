#ifndef DSIZE_H
#define DSIZE_H

class DSize{
public:

  DSize();
  DSize(int w, int h);
  DSize(const DSize &src);
  ~DSize();
  int w;
  int h;
};

#endif
