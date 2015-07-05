#ifndef MAGMA_ZLASWP_H
#define MAGMA_ZLASWP_H

/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @precisions normal z -> s d c

       auto-converted from zlaswp.cu
       
       @author Stan Tomov
       @author Mathieu Faverge
       @author Ichitaro Yamazaki
       @author Mark Gates
*/
//#include "common_magma.h"

// MAX_PIVOTS is maximum number of pivots to apply in each kernel launch
// NTHREADS is number of threads in a block
// 64 and 256 are better on Kepler; 
//#define MAX_PIVOTS 64
//#define NTHREADS   256
#define MAX_PIVOTS 32
#define NTHREADS   64

typedef struct {
    int npivots;
    int ipiv[MAX_PIVOTS];
} zlaswp_params_t;

#endif // MAGMA_ZLASWP_H
