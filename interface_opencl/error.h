#ifndef ERROR_H
#define ERROR_H

#include "common_magma.h"

// overloaded C++ functions to deal with errors
void magma_xerror( cl_int         err, const char* func, const char* file, int line );
void magma_xerror( clblasStatus   err, const char* func, const char* file, int line );
//void magma_xerror( magma_int_t  err, const char* func, const char* file, int line );

#ifdef __cplusplus
extern "C" {
#endif

// In magma.h, we also provide magma_strerror.
const char* magma_clGetErrorString( cl_int error );
const char* magma_clblasGetErrorString( clblasStatus error );

#ifdef __cplusplus
}
#endif

#ifdef NDEBUG
#define check_error( err )                     ((void)0)
#define check_xerror( err, func, file, line )  ((void)0)
#else
#define check_error( err )                     magma_xerror( err, __func__, __FILE__, __LINE__ )
#define check_xerror( err, func, file, line )  magma_xerror( err, func, file, line )
#endif

#endif        //  #ifndef ERROR_H
