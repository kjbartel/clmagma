/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zlag2c.cl mixed zc -> ds, Sat Nov 15 00:21:35 2014

       auto-converted from dlag2s.cu
       @author Mark Gates
*/
#include "kernels_header.h"
#include "dlag2s.h"

// TODO get rid of global variable!
// TODO __kernel int flag = 0;


/*
    Divides matrix into ceil( m/BLK_X ) x ceil( n/BLK_Y ) blocks.
    Each block has BLK_X threads.
    Each thread loops across one row, updating BLK_Y entries.
    
    Code similar to dlat2s and zlaset.
*/
__kernel
void dlag2s_kernel(
    int m, int n,
    const __global double *A, unsigned long A_offset, int lda,
    __global float *SA, unsigned long SA_offset,       int ldsa,
    double rmax )
{
    A += A_offset;
    SA += SA_offset;

    double tmp;
    // double neg_rmax = - rmax;
    
    int ind = get_group_id(0)*BLK_X + get_local_id(0);
    int iby = get_group_id(1)*BLK_Y;
    /* check if full block-column */
    bool full = (iby + BLK_Y <= n);
    /* do only rows inside matrix */
    if ( ind < m ) {
        A  += ind + iby*lda;
        SA += ind + iby*ldsa;
        if ( full ) {
            // full block-column
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                tmp = A[j*lda];
                // if (   (MAGMA_D_REAL(tmp) < neg_rmax) || (MAGMA_D_REAL(tmp) > rmax)
                //     #if defined(PRECISION_z) || defined(PRECISION_c)
                //     || (MAGMA_D_IMAG(tmp) < neg_rmax) || (MAGMA_D_IMAG(tmp) > rmax)
                //     #endif
                //     )
                // {
                //     flag = 1;
                // }
                SA[j*ldsa] = MAGMA_S_MAKE( MAGMA_D_REAL(tmp), MAGMA_D_IMAG(tmp) );
            }
        }
        else {
            // partial block-column
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                tmp = A[j*lda];
                // if (   (MAGMA_D_REAL(tmp) < neg_rmax) || (MAGMA_D_REAL(tmp) > rmax)
                //     #if defined(PRECISION_z) || defined(PRECISION_c)
                //     || (MAGMA_D_IMAG(tmp) < neg_rmax) || (MAGMA_D_IMAG(tmp) > rmax)
                //     #endif
                //     )
                // {
                //     flag = 1;
                // }
                SA[j*ldsa] = MAGMA_S_MAKE( MAGMA_D_REAL(tmp), MAGMA_D_IMAG(tmp) );
            }
        }
    }
}
