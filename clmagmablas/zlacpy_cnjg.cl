/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @precisions normal z -> s d c
       
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

__kernel void zlacpy_cnjg_kernel(
    int n, 
    __global magmaDoubleComplex *A1, int offset_A1, int lda1,
    __global magmaDoubleComplex *A2, int offset_A2, int lda2 )
{
    int x = get_global_id(0);
    int offset1 = x * lda1;
    int offset2 = x * lda2;

    if( x < n )
    {
        __global magmaDoubleComplex *a1  = A1 + offset_A1 + offset1;
        __global magmaDoubleComplex *a2  = A2 + offset_A2 + offset2;
        *a2 = MAGMA_Z_CNJG(*a1);
    }
}
