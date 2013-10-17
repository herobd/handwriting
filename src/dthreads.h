#ifndef DTHREADS_H
#define DTHREADS_H

///dthreads.h provides access to basic pthread functionality
/**pthread functionality is provided by including appropriate headers,
 * and (if compiling on windows) defining some macros so that the
 * pthread api can be used instead of windows api for creating
 * threads.
 */

#ifdef _WIN32
#ifndef WIN32
#define WIN32
#endif
#endif

#ifdef WIN32
#include <windows.h>
#include <process.h>

typedef HANDLE pthread_mutex_t;
typedef HANDLE pthread_t;
typedef HANDLE sem_t;
typedef void pthread_attr_t;

// pthread macros from Kenric B. White
// (some were modified by me (Douglas J. Kennard))
#define pthread_mutex_lock(A) (WAIT_FAILED == \
                               WaitForSingleObject(*(A),INFINITE))//kbw,me
#define pthread_mutex_unlock(A) (!ReleaseMutex( *(A) )) //kbw,me
#define sem_wait(A) ((WaitForSingleObject( *(A), INFINITE)),0)//kbw,me
#define sem_post(A) (!ReleaseSemaphore( *(A), 1, NULL))//kbw,me
#define sem_init(A,B,C) (-1 * (NULL == \
                        ((*(A)) = (CreateSemaphore(NULL,(C),9999999,NULL)))))//kbw,me
#define pthread_mutex_init(A,B) ((*(A) = (CreateMutex(NULL, FALSE, "A"))),0)//kbw,me

// this may not work, I haven't tried it.  It may not even compile
inline int pthread_create(pthread_t *thread, pthread_attr_t *attr,
			  void *(*start_routine)(void *), void *arg){
  (*thread) = _beginthreadex(NULL, 0,
			     (unsigned (__stdcall *)(void *))start_routine,
			     arg, 0, NULL);
  if(0 == (*thread))
    return 1;
  return 0;
}

// this may not work, I haven't tried it.  It may not even compile
inline int pthread_mutex_destroy(pthread_mutex_t mtx){
  return !(closehandle(mtx));
}


inline int getNumCPUs(){
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  return si.dwNumberOfProcessors;
}

#endif


#ifndef WIN32
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <pthread.h>
#include <semaphore.h>

#include <unistd.h> // needed for getNumCPUs()
inline int getNumCPUs(){
  return sysconf(_SC_NPROCESSORS_ONLN);
}

#endif



#endif
