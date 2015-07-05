/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zlascl_2x2.cl normal z -> c, Sat Nov 15 00:21:35 2014

       @author Ichitaro Yamazaki
       @author Stan Tomov
*/
#include "kernels_header.h"

#define A(i,j) (A[offset_A + (i) + (j)*lda])
#define W(i,j) (W[offset_W + (i) + (j)*ldw])

// each thread block does one NB x n block row of A.
// each thread does one row, starting from left edge and moving right.
__kernel void clascl_2x2_full_kernel(
    int m, 
    __global const magmaFloatComplex* W, int offset_W, int ldw, 
    __global magmaFloatComplex* A, int offset_A, int lda)
{
    int ind = get_global_id(0);
 
    magmaFloatComplex D21 = W( 1, 0 );
    magmaFloatComplex D11 = MAGMA_C_DIV( W( 1, 1 ), D21 );
    magmaFloatComplex D22 = MAGMA_C_DIV( W( 0, 0 ), MAGMA_C_CNJG( D21 ) );
    float T = 1.0 / ( MAGMA_C_REAL( MAGMA_C_MUL(D11,D22) ) - 1.0 );
    D21 = MAGMA_C_DIV( MAGMA_C_MAKE(T,0.0), D21 );

    if (ind < m) {
        A( ind, 0 ) = MAGMA_C_MUL( MAGMA_C_CNJG( D21 ), 
                                   MAGMA_C_SUB( MAGMA_C_MUL(D11,W( 2+ind, 0 )), 
                                                W( 2+ind, 1 ) ));
        A( ind, 1 ) = MAGMA_C_MUL( D21, 
                                   MAGMA_C_SUB( MAGMA_C_MUL(D22,W( 2+ind, 1 )),
                                                W( 2+ind, 0 ) ));
    }
}


// each thread block does one NB x n block row of A.
// each thread does one row, starting from left edge and moving right to diagonal.
__kernel void clascl_2x2_lower_kernel(
    int m,
    __global const magmaFloatComplex* W, int offset_W, int ldw,
    __global magmaFloatComplex* A, int offset_A, int lda)
{
    int ind = get_global_id(0);

    magmaFloatComplex D21 = W( 1, 0 );
    magmaFloatComplex D11 = MAGMA_C_DIV( W( 1, 1 ), D21 );
    magmaFloatComplex D22 = MAGMA_C_DIV( W( 0, 0 ), MAGMA_C_CNJG( D21 ) );
    float T = 1.0 / ( MAGMA_C_REAL( MAGMA_C_MUL(D11,D22) ) - 1.0 );
    D21 = MAGMA_C_DIV( MAGMA_C_MAKE(T,0.0), D21 );

    if (ind < m) {
    A( ind, 0 ) = MAGMA_C_MUL( MAGMA_C_CNJG( D21 ),
                                   MAGMA_C_SUB( MAGMA_C_MUL(D11,W( 2+ind, 0 )),
                                                W( 2+ind, 1 ) ));
        A( ind, 1 ) = MAGMA_C_MUL( D21,
                                   MAGMA_C_SUB( MAGMA_C_MUL(D22,W( 2+ind, 1 )),
                                                W( 2+ind, 0 ) ));
    }
}


// each thread block does one NB x n block row of A.
// each thread does one row, starting from right edge and moving left to diagonal.
__kernel void clascl_2x2_upper_kernel(
    int m,
    __global const magmaFloatComplex* W, int offset_W, int ldw,
    __global magmaFloatComplex* A, int offset_A, int lda)
{
    int ind = get_global_id(0);

    magmaFloatComplex D21 = W( m, 1 );
    magmaFloatComplex D11 = MAGMA_C_DIV( W( m+1, 1 ), MAGMA_C_CNJG(D21) );
    magmaFloatComplex D22 = MAGMA_C_DIV( W( m, 0 ), D21 );
    float T = 1.0 / ( MAGMA_C_REAL( MAGMA_C_MUL(D11,D22) ) - 1.0 );
    D21 = MAGMA_C_DIV( MAGMA_C_MAKE(T,0.0), D21 );

    if (ind < m) {
        A( ind, 0 ) = MAGMA_C_MUL( D21,
                                   MAGMA_C_SUB( MAGMA_C_MUL(D11, W(ind,0)),
                                                W( ind, 1 ) ));
        A( ind, 1 ) = MAGMA_C_MUL( MAGMA_C_CNJG(D21),
                                   MAGMA_C_SUB( MAGMA_C_MUL(D22, W(ind,1)),
                                                W(ind,0 ) ));
    }
}
