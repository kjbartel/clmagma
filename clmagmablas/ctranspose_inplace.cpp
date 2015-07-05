/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from ztranspose_inplace.cpp normal z -> c, Sat Nov 15 00:21:35 2014

       auto-converted from ctranspose_inplace.cu

       @author Stan Tomov
       @author Mark Gates
*/
#include "clmagma_runtime.h"
#include "common_magma.h"
#include "ctranspose_inplace.h"


/**
    Purpose
    -------
    ctranspose_inplace_q transposes a square N-by-N matrix in-place.
    
    Same as ctranspose_inplace, but adds queue argument.
    
    Arguments
    ---------
    @param[in]
    n       INTEGER
            The number of rows & columns of the matrix dA.  N >= 0.
    
    @param[in]
    dA      COMPLEX array, dimension (LDDA,N)
            The N-by-N matrix dA.
            On exit, dA(j,i) = dA_original(i,j), for 0 <= i,j < N.
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= N.
    
    @param[in]
    queue   magma_queue_t
            Queue to execute in.
    
    @ingroup magma_caux2
    ********************************************************************/
extern "C" void
magmablas_ctranspose_inplace(
    magma_int_t n,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue )
{
    cl_kernel kernel;
    cl_int err;
    int i;

    magma_int_t info = 0;
    if ( n < 0 )
        info = -1;
    else if ( ldda < n )
        info = -3;
    
    if ( info != 0 ) {
        magma_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    size_t threads[2] = { NB, NB };
    int nblock = (n + NB - 1)/NB;
    
    // need 1/2 * (nblock+1) * nblock to cover lower triangle and diagonal of matrix.
    // block assignment differs depending on whether nblock is odd or even.
    if( nblock % 2 == 1 ) {
        size_t grid[2] = { nblock, (nblock+1)/2 };
        grid[0] *= threads[0];
        grid[1] *= threads[1];
        
        kernel = g_runtime.get_kernel( "ctranspose_inplace_odd" );
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
            err |= clSetKernelArg( kernel, i++, sizeof(dA       ), &dA        );
            err |= clSetKernelArg( kernel, i++, sizeof(dA_offset), &dA_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldda     ), &ldda      );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 2, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    }
    else {
        size_t grid[2] = { nblock+1, nblock/2 };
        grid[0] *= threads[0];
        grid[1] *= threads[1];
        
        kernel = g_runtime.get_kernel( "ctranspose_inplace_even" );
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
            err |= clSetKernelArg( kernel, i++, sizeof(dA       ), &dA        );
            err |= clSetKernelArg( kernel, i++, sizeof(dA_offset), &dA_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldda     ), &ldda      );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 2, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    }
}
