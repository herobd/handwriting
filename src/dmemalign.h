#ifndef DMEMALIGN_H
#define DMEMALIGN_H

#ifdef DEBUG
#include <stdio.h>
#include <errno.h>
#endif

#ifdef _WIN32
#include <malloc.h>
///allocate a block of bufSize bytes, aligned on align-byte boundary
/**align must be a power of 2 and a multiple of sizeof(void*)*/
inline void* daligned_malloc(size_t bufSize, size_t align){
  return _aligned_malloc(bufSize, align);
}
///free a block of memory previously allocated with daligned_malloc()
inline void daligned_free(void *pbuf){
  _aligned_free(pbuf);
}
#else
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif
#include <stdlib.h>
///allocate a block of bufSize bytes, aligned on align-byte boundary
/**align must be a power of 2 and a multiple of sizeof(void*)
*/
inline void* daligned_malloc(size_t bufSize, size_t align){
  void *ptmp;
  int retVal;

  retVal = posix_memalign(&ptmp, align, bufSize);
  if(0 == retVal)
    return ptmp;
#ifdef DEBUG
  if(EINVAL == retVal)
    fprintf(stderr, "daligned_malloc() posix_memalign returned EINVAL(%d) "
	    "(align=%d must be power of 2 and multiple of sizeof(void*)=%d\n",
	    EINVAL, align, sizeof(void*));
#endif
  return NULL;
}
///free a block of memory previously allocated with daligned_malloc()
inline void daligned_free(void *pbuf){
  free(pbuf);
}
#endif



#endif
