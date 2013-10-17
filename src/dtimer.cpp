#include "dtimer.h"
#include "dinstancecounter.h"
#include <time.h>
#include <errno.h>
#include <stdio.h>

DTimer::DTimer(){
  DInstanceCounter::addInstance("DTimer");
  this->start();
}

DTimer::DTimer(const DTimer &src){
  DInstanceCounter::addInstance("DTimer");
  this->startTime = src.startTime;
  this->accumulatedTime = src.accumulatedTime;
  this->fStopped = src.fStopped;
}

DTimer::~DTimer(){
  DInstanceCounter::removeInstance("DTimer");
}

#ifdef _WIN32
void DTimer::sleep(double secs){
  Sleep((DWORD)(secs*1000));
}
#else
void DTimer::sleep(double secs){
  struct timespec tsSleep;
  struct timespec tsRem;
  int retVal = 0;
  int tries = 0;

  tsSleep.tv_sec = (int)secs;
  tsSleep.tv_nsec = (long)(( secs - (int)secs) * 1.0E9);

  retVal = nanosleep(&tsSleep, &tsRem);
  while((-1==retVal) && (EINTR==errno) && (tries < 10)){
    ++tries;
#ifdef DEBUGLIBRARY
    fprintf(stderr, "DTimer::sleep() "
	    "try #%d failed due to interrupt. Retry.\n", tries);
#endif
    tsSleep = tsRem;
    retVal = nanosleep(&tsSleep, &tsRem);
  }
#ifdef DEBUG
  if(tries == 10){
    fprintf(stderr, "DTimer::sleep() "
	    "failed due to 10 or more interrupts while trying\n");
  }
#endif
}
#endif
