#ifndef DPROGRESS_H
#define DPROGRESS_H

///Provides a mechanism for reporting progress of long operations
/**Some long operations permit a pointer to a DProgress object to be
 * passed into them.  If the pointer is NULL, then no status will be
 * reported.  If the pointer is not NULL, the long operation will
 * periodically call the reportStatus() function with the current
 * status.  The long operations should also check the return value of
 * reoprtStatus() to see if they should continue processing or not (a
 * non-zero return value from reportStatus() tells the long operation
 * that it should cancel itself).  The reportStatus() function in this
 * base class just prints the current progress status in numbers and
 * percentage to stdout, but the DProgress class can be sub-classed
 * and the reportStatus() method can be changed to do whatever you
 * need it to do when status is reported.
 */
class DProgress {
public:
  DProgress();
  virtual ~DProgress();
  virtual int reportStatus(int cur, int min, int max);
  int getMin() const;
  int getMax() const;
  int getCur() const;
  void cancelOperation();
  void reset();
protected:
  int min, max, cur;
  bool fCancelRequested;
};

#endif
