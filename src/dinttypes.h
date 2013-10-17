#ifndef DINTTYPES_H
#define DINTTYPES_H

typedef unsigned char D_uint8;
typedef unsigned short D_uint16;
typedef unsigned int D_uint32;
typedef signed char D_sint8;
typedef signed short D_sint16;
typedef signed int D_sint32;
typedef unsigned long long D_uint64;

#ifdef _WIN32
typedef int ssize_t;
#endif

#endif
