/*
 *   -- clMAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      @date November 2014
 */
#include "kernels_header.h"

// empty kernel, benchmark in iwocl 2013
__kernel void empty_kernel(
    int i0, int i1, int i2, int i3, int i4, 
    int i5, int i6, int i7, int i8, int i9,
    double d0, double d1, double d2, double d3, double d4, 
    __global double *dA,
    __global double *dB,
    __global double *dC )
{
    int tid = get_local_id(0);

    for( int i=0; i < i0; i++ ) {
        dC[i+tid] += d1*dC[i+tid] + d2*dA[i]*dB[i];
    }
    barrier( CLK_LOCAL_MEM_FENCE );
    for( int i=0; i < i0; i++ ) {
        dC[i+tid] += d1*dC[i+tid] + d2*dA[i]*dB[i];
    }
}
