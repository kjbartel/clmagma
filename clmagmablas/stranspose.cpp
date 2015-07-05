/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from ztranspose.cpp normal z -> s, Sat Nov 15 00:21:35 2014

       auto-converted from stranspose.cu

       @author Stan Tomov
       @author Mark Gates
*/
#include "clmagma_runtime.h"
#include "common_magma.h"
#include "stranspose.h"


/**
    Purpose
    -------
    stranspose_q copies and transposes a matrix dA to matrix dAT.
    
    Same as stranspose, but adds queue argument.
        
    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix dA.  M >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix dA.  N >= 0.
    
    @param[in]
    dA      REAL array, dimension (LDDA,N)
            The M-by-N matrix dA.
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= M.
    
    @param[in]
    dAT     REAL array, dimension (LDDA,N)
            The N-by-M matrix dAT.
    
    @param[in]
    lddat   INTEGER
            The leading dimension of the array dAT.  LDDAT >= N.
    
    @param[in]
    queue   magma_queue_t
            Queue to execute in.
    
    @ingroup magma_saux2
    ********************************************************************/
extern "C" void
magmablas_stranspose(
    magma_int_t m, magma_int_t n,
    magmaFloat_const_ptr dA, size_t dA_offset,  magma_int_t ldda,
    magmaFloat_ptr       dAT, size_t dAT_offset, magma_int_t lddat,
    magma_queue_t queue )
{
    cl_kernel kernel;
    cl_int err;
    int i;

    magma_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( ldda < m )
        info = -4;
    else if ( lddat < n )
        info = -6;
    
    if ( info != 0 ) {
        magma_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    size_t threads[2] = { NX, NY };
    size_t grid[2] = { (m+NB-1)/NB, (n+NB-1)/NB };
    grid[0] *= threads[0];
    grid[1] *= threads[1];
    kernel = g_runtime.get_kernel( "stranspose_kernel" );
    if ( kernel != NULL ) {
        err = 0;
        i   = 0;
        err |= clSetKernelArg( kernel, i++, sizeof(m         ), &m          );
        err |= clSetKernelArg( kernel, i++, sizeof(n         ), &n          );
        err |= clSetKernelArg( kernel, i++, sizeof(dA        ), &dA         );
        err |= clSetKernelArg( kernel, i++, sizeof(dA_offset ), &dA_offset  );
        err |= clSetKernelArg( kernel, i++, sizeof(ldda      ), &ldda       );
        err |= clSetKernelArg( kernel, i++, sizeof(dAT       ), &dAT        );
        err |= clSetKernelArg( kernel, i++, sizeof(dAT_offset), &dAT_offset );
        err |= clSetKernelArg( kernel, i++, sizeof(lddat     ), &lddat      );
        check_error( err );

        err = clEnqueueNDRangeKernel( queue, kernel, 2, NULL, grid, threads, 0, NULL, NULL );
        check_error( err );
    }
}
