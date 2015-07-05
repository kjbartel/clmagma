/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @precisions mixed zc -> ds

       auto-converted from clag2z.cu
       @author Mark Gates
*/
#include "kernels_header.h"
#include "clag2z.h"


/*
    Divides matrix into ceil( m/BLK_X ) x ceil( n/BLK_Y ) blocks.
    Each block has BLK_X threads.
    Each thread loops across one row, updating BLK_Y entries.
    
    Code similar to clat2z and zlaset.
*/
__kernel
void clag2z_kernel(
    int m, int n,
    const __global magmaFloatComplex *SA, unsigned long SA_offset, int ldsa,
    __global magmaDoubleComplex       *A, unsigned long A_offset, int lda )
{
    SA += SA_offset;
    A += A_offset;

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
                A[j*lda] = MAGMA_Z_MAKE( MAGMA_C_REAL( SA[j*ldsa] ), MAGMA_C_IMAG( SA[j*ldsa] ));
            }
        }
        else {
            // partial block-column
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                A[j*lda] = MAGMA_Z_MAKE( MAGMA_C_REAL( SA[j*ldsa] ), MAGMA_C_IMAG( SA[j*ldsa] ));
            }
        }
    }
}
