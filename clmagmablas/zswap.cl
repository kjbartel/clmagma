/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014
       
       @author Mark Gates

       @precisions normal z -> s d c

       auto-converted from zswap.cu

*/
#include "kernels_header.h"
#include "zswap.h"


/* Vector is divided into ceil(n/nb) blocks.
   Each thread swaps one element, x[tid] <---> y[tid].
*/
__kernel void zswap_kernel(
    int n,
    __global magmaDoubleComplex *x, unsigned long x_offset, int incx,
    __global magmaDoubleComplex *y, unsigned long y_offset, int incy )
{
    x += x_offset;
    y += y_offset;

    magmaDoubleComplex tmp;
    int ind = get_local_id(0) + get_local_size(0)*get_group_id(0);
    if ( ind < n ) {
        x += ind*incx;
        y += ind*incy;
        tmp = *x;
        *x  = *y;
        *y  = tmp;
    }
}
