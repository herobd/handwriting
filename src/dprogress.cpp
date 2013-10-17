#include <stdio.h>
#include "dprogress.h"
#include "dinstancecounter.h"

DProgress::DProgress(){
  DInstanceCounter::addInstance("DProgress");
  reset();
}

DProgress::~DProgress(){
  DInstanceCounter::removeInstance("DProgress");
  // do nothing
}


///The callback that long operations call to update the current progress status
/**This function by default prints the status to stdout.  You can
 * create a child class of DProgress that overrides this callback to
 * do whatever else you want.  However, it should be quick so the
 * processing can continue (whatever you do will happen in the same
 * thread as the processing). Functions that call reportStatus()
 * should check the return value, to see if they should continue to
 * run or not.  If the return value is 0, the long operation can
 * continue to run, but if the return value is non-zero, the long
 * operation should stop processing, as this means the operation has
 * been cancelled.  It should call reportStatus one last time with
 * cur=min-1 to notify the monitoring function that the cancel has
 * occured.  When reportStatus() is called with cur>=max, the
 * operation is complete.  When reportStatus() is called with cur<min,
 * the operation was cancelled.
 */
int DProgress::reportStatus(int cur, int min, int max){
  printf("%d/%d (%.1f%%)\n", cur, max-min, 100.*cur/(double)(max-min));

  if(cur == min-1){
    fprintf(stderr, "--- operation cancelled ---\n");
  }
  else if(cur >= max){
    fprintf(stderr, "--- operation complete ---\n");
  }
  return fCancelRequested ? 1 : 0;
}

///Accessor function for applications to check the min (as last reported)
inline int DProgress::getMin() const{
  return min;
}
///Accessor function for applications to check the max (as last reported)
inline int DProgress::getMax() const{
  return max;
}
///Accessor function for applications to check the last reported progress value
inline int DProgress::getCur() const{
  return cur;
}

///Requests that a long operation be cancelled.
/**The next time the long operation calls reportStatus(), reportStatus
 * will return a non-zero number to the long operation to signal it
 * that it should stop processing.
 */
void DProgress::cancelOperation(){
  fCancelRequested = true;
}

///Clear any cancel request and reset the min, max, and cur values to defaults.
void DProgress::reset(){
  min = 0;
  max = 100;
  cur = 0;
  fCancelRequested = false;
}
