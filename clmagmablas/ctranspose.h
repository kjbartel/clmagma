#ifndef MAGMA_CTRANSPOSE_H
#define MAGMA_CTRANSPOSE_H

/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from ztranspose.h normal z -> c, Sat Nov 15 00:21:35 2014

       auto-converted from ctranspose.cu

       @author Stan Tomov
       @author Mark Gates
*/
//#include "common_magma.h"

#define PRECISION_c

#if defined(PRECISION_z)
    #define NX 16
#else
    #define NX 32
#endif

#define NB 32
#define NY 8

#endif // MAGMA_CTRANSPOSE_H
