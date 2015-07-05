/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zcaxpycp.cl mixed zc -> ds, Sat Nov 15 00:21:35 2014

       auto-converted from dsaxpycp.cu

*/
#include "kernels_header.h"
#include "dsaxpycp.h"

// adds   x += r (including conversion to double)  --and--
// copies w = b
// each thread does one index, x[i] and w[i]
__kernel void
dsaxpycp_kernel(
    int m, __global float *r, unsigned long r_offset, 
    __global double *x, unsigned long x_offset,
    const __global double *b, unsigned long b_offset, 
    __global double *w, unsigned long w_offset )
{
    r += r_offset;
    x += x_offset;
    b += b_offset;
    w += w_offset;

    const int i = get_local_id(0) + get_group_id(0)*NB;
    if ( i < m ) {
        x[i] = MAGMA_D_ADD( x[i], MAGMA_D_MAKE( (double)MAGMA_D_REAL( r[i] ) , 
                                                (double)MAGMA_D_IMAG( r[i] ) ));
        w[i] = b[i];
    }
}


// adds   x += r  --and--
// copies r = b
// each thread does one index, x[i] and r[i]
__kernel void
daxpycp_kernel(
    int m, __global double *r, unsigned long r_offset, 
    __global double *x, unsigned long x_offset,
    const __global double *b, unsigned long b_offset)
{
    r += r_offset;
    x += x_offset;
    b += b_offset;

    const int i = get_local_id(0) + get_group_id(0)*NB;
    if ( i < m ) {
        x[i] = MAGMA_D_ADD( x[i], r[i] );
        r[i] = b[i];
    }
}
