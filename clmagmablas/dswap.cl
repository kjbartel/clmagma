/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014
       
       @author Mark Gates

       @generated from zswap.cl normal z -> d, Sat Nov 15 00:21:35 2014

       auto-converted from dswap.cu

*/
#include "kernels_header.h"
#include "dswap.h"


/* Vector is divided into ceil(n/nb) blocks.
   Each thread swaps one element, x[tid] <---> y[tid].
*/
__kernel void dswap_kernel(
    int n,
    __global double *x, unsigned long x_offset, int incx,
    __global double *y, unsigned long y_offset, int incy )
{
    x += x_offset;
    y += y_offset;

    double tmp;
    int ind = get_local_id(0) + get_local_size(0)*get_group_id(0);
    if ( ind < n ) {
        x += ind*incx;
        y += ind*incy;
        tmp = *x;
        *x  = *y;
        *y  = tmp;
    }
}
