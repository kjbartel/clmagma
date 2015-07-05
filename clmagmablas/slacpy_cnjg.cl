/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zlacpy_cnjg.cl normal z -> s, Sat Nov 15 00:21:35 2014
       
       @author Ichitaro Yamazaki
       @author Stan Tomov
*/
#include "kernels_header.h"

/*********************************************************
 *
 * SWAP BLAS: permute to set of N elements
 *
 ********************************************************/
/*
 *  First version: line per line
 */

__kernel void slacpy_cnjg_kernel(
    int n, 
    __global float *A1, int offset_A1, int lda1,
    __global float *A2, int offset_A2, int lda2 )
{
    int x = get_global_id(0);
    int offset1 = x * lda1;
    int offset2 = x * lda2;

    if( x < n )
    {
        __global float *a1  = A1 + offset_A1 + offset1;
        __global float *a2  = A2 + offset_A2 + offset2;
        *a2 = MAGMA_S_CNJG(*a1);
    }
}
