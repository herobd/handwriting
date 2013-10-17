#ifndef DINSTANCECOUNTER_H
#define DINSTANCECOUNTER_H

//***THIS CLASS IS NOT THREAD-SAFE, SO I HAVE DISABLED IT BELOW!!!****

#include <map>
#include <stdio.h>

class DInstanceCounter{
public:
  static void addInstance(const char *stClassName);///< Increment count of constructed instances of class stClassName (does nothing if DEBUG not defined)
  static void removeInstance(const char *stClassName);///< Increment count of destroyed instances of class stClassName (does nothing if DEBUG not defined)
  static void report();///< Output a list of all tracked-classes' instantiation counts to stderr (DEBUG must be defined) 
  static void cleanup();
  static int numUnfreedInstances();///< current number of non-destroyed instances of all tracked objects (DEBUG must be defined)
private:
//   /* the lstr struct is copied from an example on SGI's stl documentation page
//      URL: http://www.sgi.com/tech/stl/Map.html.  It provides the comparison.*/
//   struct lstr{
//     bool operator()(const char* s1, const char* s2) const{
//       return strcmp(s1,s2) < 0;
//     }
//   };

  static std::map<const char*, int/*, lstr*/> mapNumConstructed;
  static std::map<const char*, int/*, lstr*/> mapNumDestructed;
};


// #ifdef DEBUG
#if 0 //I disabled this class because it is not thread-safe
inline void DInstanceCounter::addInstance(const char *stClassName){
  std::map<const char*, int/*, lstr*/>::iterator iter;
  int numInstances = 0;

  iter = mapNumConstructed.find(stClassName);
  if(iter != mapNumConstructed.end()){
    numInstances = (*iter).second;
  }
  else{
    mapNumDestructed[stClassName] = 0;
  }
  
  mapNumConstructed[stClassName] = numInstances+1;
}

inline void DInstanceCounter::removeInstance(const char *stClassName){
  std::map<const char*, int/*, DInstanceCounter::lstr*/>::iterator iter;
  int numInstances = 0;
  iter = mapNumDestructed.find(stClassName);
  if(iter == mapNumDestructed.end()){
    fprintf(stderr, "DInstanceCounter::removeInstance() couldn't find '%s'\n",
	    stClassName);
    return;
  }
  numInstances = (*iter).second;
  mapNumDestructed[stClassName] = numInstances+1;
}

inline void DInstanceCounter::report(){
  std::map<const char*, int/*, lstr*/>::iterator iter1;
  std::map<const char*, int/*, lstr*/>::iterator iter2;

  iter1 = mapNumConstructed.begin();
  iter2 = mapNumDestructed.begin();

  fprintf(stderr, "------DInstanceCounter::report():\n");
  while(iter1 != mapNumConstructed.end()){
    fprintf(stderr,"  %s: constructed=%d destructed=%d current=%d\n",
	    (*iter1).first, (*iter1).second, (*iter2).second,
	    ((*iter1).second-(*iter2).second));
    ++iter1;
    ++iter2;
  }
}

inline int DInstanceCounter::numUnfreedInstances(){
  int retVal = 0;
  std::map<const char*, int/*, DInstanceCounter::lstr*/>::iterator iter1;
  std::map<const char*, int/*, DInstanceCounter::lstr*/>::iterator iter2;

  iter1 = mapNumConstructed.begin();
  iter2 = mapNumDestructed.begin();

  while(iter1 != mapNumConstructed.end()){
    retVal += ((*iter1).second-(*iter2).second);
    ++iter1;
    ++iter2;
  }
  return retVal;
}
#else /* make the funtions NO-OPS when not compiling debug so it's not slow*/
inline void DInstanceCounter::addInstance(const char *stClassName){
  (void)stClassName; // NOOP - avoid compiler warning for unused parameter
}
inline void DInstanceCounter::removeInstance(const char *stClassName){
  (void)stClassName; // NOOP - avoid compiler warning for unused parameter
}
inline void DInstanceCounter::report(){
  fprintf(stderr, "DInstanceCounter::report() called from non-debug code!\n");
}
inline int DInstanceCounter::numUnfreedInstances(){
  fprintf(stderr, "DInstanceCounter::numUnfreedInstances() called from "
	  "non-debug code!\n");
  return 0;
}
#endif /* DEBUG */

#endif
