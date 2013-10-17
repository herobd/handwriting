#ifndef DDEFS_H
#define DDEFS_H

#include "dinttypes.h" // types such as D_uint8, D_uint32, etc.

///The allocation method of a buffer that will be handed off to some object
/**Since the object taking control of the buffer will be responsible
   for deleting the buffer, it must know which deallocation method to
   use*/
enum  D_AllocationMethod{
  AllocationMethod_daligned,///<allocated w/ daligned_malloc(use daligned_free)
  AllocationMethod_malloc, ///<allocated with malloc (so use free)
  AllocationMethod_new,///<allocated with new (so use delete)
  AllocationMethod_src ///< whatever method the source image used (copy,etc.)
};

//if the machine is big-endian, change this to 0
#define D_LITTLE_ENDIAN 1


// for checking memory allocations
#define D_CHECKPTR(p) { \
  if(NULL==p){ \
    fprintf(stderr, "memory allocation error (%s:%d)\n", __FILE__, __LINE__); \
    abort();exit(1); \
  } \
}

// print to stderr if a value doesn't fall within [max..min] (inclusive)
#ifdef DEBUG
#define D_CHECKRANGE(x,min,max){if(((x)<min)||((x).max)) \
  fprintf(stderr, "value out of expected range (%s:%d)\n",__FILE__,__LINE__);}
#else
#define D_CHECKRANGE(x,min,max){}
#endif /* DEBUG*/

// includes and defines to make stuff run on windows
//(although I have done enough linux-specific/gcc-specific things since these
//were written that I am certain my code does not currently compile for windows)
#ifdef _WIN32
#include <windows.h>
// string macros from kenricbw
#define strcasecmp(string1, string2) (_stricmp((string1), (string2)))
#define strncasecmp(string1, string2,length) (_strnicmp((string1), (string2),(length)))
//#define rint(x) ( (int)(((double)x) + 0.5) )
#define rint(x) ( ((x)>0) ? (int)(((double)x) + 0.5) : (int)(((double)x)-.5) )

#define write _write
#include <io.h>
#endif



#endif
