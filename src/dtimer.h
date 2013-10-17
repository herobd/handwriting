#ifndef DTIMER_H
#define DTIMER_H

/// This class provides simple functionality for timing code, etc.
/** Timing code can be as simple as creating a new timer (or calling
    start() on an existing timer) just before the section of code to
    be timed, and then calling getAccumulated() just after the section
    of code te return the number of elapsed seconds.  For eample:

    \code
    DTimer dt;

    doSomeLengthyCalculations(...);
    cout << "elapsed time = " << dt.getAccumulated() << " seconds\n";
    \endcode

    \warning Resolution of timers is not perfect.
    Also, the timers, themselves can throw off timing if they are being
    used frequently, such as inside of a tight loop.  Timing should be done
    AROUND (not within) the code that executes frequently.
*/

//#ifdef DEBUG
#include <stdio.h>
//#endif

#ifndef _WIN32
#include <sys/time.h>
#endif

class DTimer{
public:
  DTimer();
  ///<Default Constructor: Start the new timer now with accumulated time = 0.0.
  DTimer(const DTimer &src);///<Copy Constructor: The state is copied from src.
  ~DTimer();///< Destructor

  void start();///<Start timing beginning now with accumulated time = 0.0.
  void start(double stimeSecs);
  ///<Same as start(), but use stimeSecs instead of now
  void stop();///<Stop accumulating now, and don't clear the accumulated time.
  void stop(double stimeSecs);
  ///<Same as stop(), but use stimeSecs instead of now
  void resume();///<Continue accumulating time starting now (without resetting)
  void resume(double stimeSecs);
  ///<Same as resume(), but use stimeSecs instead of now

  double getAccumulated();
  ///<Total secs accumulated, including current split if timer is running.
  double getSplit();///<Secs since most recent start(), resume(), or creation
  static double getNow();///<Current system time in seconds
  static void sleep(double secs);///< sleep for at least secs seconds.
  
private:
  double startTime; // time of most recent of: creation, start, or resume
  double accumulatedTime; // total time accumulated (not including the current
                          // amount of time since most recent start/resume if
                          // not currently stopped)
  bool fStopped;
};


#ifdef _WIN32
#include <windows.h> // for timeGetTime() function and Sleep() function
inline double DTimer::getNow(){
  return ((double)timeGetTime()) / 1000.;
}
#else
inline double DTimer::getNow(){
  struct timeval tv;

  if(0 != gettimeofday(&tv, NULL)){
#ifdef DEBUG
    fprintf(stderr, "DTimer::getNow() gettimeofday failed\n");
#endif
  }
  return ( (double)(tv.tv_sec) + tv.tv_usec * 1.0E-6);
}
#endif

inline void DTimer::start(){
  this->start(DTimer::getNow());
}

inline void DTimer::start(double stimeSecs){
  fStopped = false;
  accumulatedTime = 0.;
  startTime = stimeSecs;
}

inline void DTimer::stop(){
  this->stop(DTimer::getNow());
}

inline void DTimer::stop(double stimeSecs){
#ifdef DEBUG
  if(fStopped)
    fprintf(stderr,"DTimer::stop() called while timer not running\n");
#endif
  if(!fStopped){
    accumulatedTime += stimeSecs - startTime;
    fStopped = true;
  }
}

inline void DTimer::resume(){
  this->resume(DTimer::getNow());
}

inline void DTimer::resume(double stimeSecs){
#ifdef DEBUG
  if(!fStopped)
    fprintf(stderr,"DTimer::resume() called while timer running\n");
#endif
  fStopped = false;
  startTime = stimeSecs;
}
inline double DTimer::getAccumulated(){
  if(!fStopped)
    return (accumulatedTime + getSplit());
  return accumulatedTime;
}
inline double DTimer::getSplit(){
  return (getNow() - startTime);
}

#endif
