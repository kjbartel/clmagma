/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @precisions mixed zc -> ds

       auto-converted from zcaxpycp.cu

*/
#include "kernels_header.h"
#include "zcaxpycp.h"

// adds   x += r (including conversion to double)  --and--
// copies w = b
// each thread does one index, x[i] and w[i]
__kernel void
zcaxpycp_kernel(
    int m, __global magmaFloatComplex *r, unsigned long r_offset, 
    __global magmaDoubleComplex *x, unsigned long x_offset,
    const __global magmaDoubleComplex *b, unsigned long b_offset, 
    __global magmaDoubleComplex *w, unsigned long w_offset )
{
    r += r_offset;
    x += x_offset;
    b += b_offset;
    w += w_offset;

    const int i = get_local_id(0) + get_group_id(0)*NB;
    if ( i < m ) {
        x[i] = MAGMA_Z_ADD( x[i], MAGMA_Z_MAKE( (double)MAGMA_Z_REAL( r[i] ) , 
                                                (double)MAGMA_Z_IMAG( r[i] ) ));
        w[i] = b[i];
    }
}


// adds   x += r  --and--
// copies r = b
// each thread does one index, x[i] and r[i]
__kernel void
zaxpycp_kernel(
    int m, __global magmaDoubleComplex *r, unsigned long r_offset, 
    __global magmaDoubleComplex *x, unsigned long x_offset,
    const __global magmaDoubleComplex *b, unsigned long b_offset)
{
    r += r_offset;
    x += x_offset;
    b += b_offset;

    const int i = get_local_id(0) + get_group_id(0)*NB;
    if ( i < m ) {
        x[i] = MAGMA_Z_ADD( x[i], r[i] );
        r[i] = b[i];
    }
}
