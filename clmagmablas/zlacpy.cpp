/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @precisions normal z -> s d c

       auto-converted from zlacpy.cu
       @author Mark Gates
       @author Azzam Haidar
*/
#include "clmagma_runtime.h"
#include "common_magma.h"
#include "zlacpy.h"


/**
    Purpose
    -------
    ZLACPY_Q copies all or part of a two-dimensional matrix dA to another
    matrix dB.
    
    This is the same as ZLACPY, but adds queue argument.
    
    Arguments
    ---------
    
    @param[in]
    uplo    magma_uplo_t
            Specifies the part of the matrix dA to be copied to dB.
      -     = MagmaUpper:      Upper triangular part
      -     = MagmaLower:      Lower triangular part
            Otherwise:  All of the matrix dA
    
    @param[in]
    m       INTEGER
            The number of rows of the matrix dA.  M >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix dA.  N >= 0.
    
    @param[in]
    dA      COMPLEX_16 array, dimension (LDDA,N)
            The m by n matrix dA.
            If UPLO = MagmaUpper, only the upper triangle or trapezoid is accessed;
            if UPLO = MagmaLower, only the lower triangle or trapezoid is accessed.
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
    
    @param[out]
    dB      COMPLEX_16 array, dimension (LDDB,N)
            The m by n matrix dB.
            On exit, dB = dA in the locations specified by UPLO.
    
    @param[in]
    lddb    INTEGER
            The leading dimension of the array dB.  LDDB >= max(1,M).
    
    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    @ingroup magma_zaux2
    ********************************************************************/
extern "C" void
magmablas_zlacpy(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr       dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue )
{
    cl_kernel kernel;
    cl_int err;
    int i;

    magma_int_t info = 0;
    if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( ldda < max(1,m))
        info = -5;
    else if ( lddb < max(1,m))
        info = -7;
    
    if ( info != 0 ) {
        magma_xerbla( __func__, -(info) );
        return;
    }
    
    if ( m == 0 || n == 0 )
        return;
    
    size_t threads[2] = { BLK_X, 1 };
    size_t grid[2] = { (m + BLK_X - 1)/BLK_X, (n + BLK_Y - 1)/BLK_Y };
    grid[0] *= threads[0];
    grid[1] *= threads[1];
    
    if ( uplo == MagmaLower ) {
        kernel = g_runtime.get_kernel( "zlacpy_kernel_lower" );
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            err |= clSetKernelArg( kernel, i++, sizeof(m        ), &m         );
            err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
            err |= clSetKernelArg( kernel, i++, sizeof(dA       ), &dA        );
            err |= clSetKernelArg( kernel, i++, sizeof(dA_offset), &dA_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldda     ), &ldda      );
            err |= clSetKernelArg( kernel, i++, sizeof(dB       ), &dB        );
            err |= clSetKernelArg( kernel, i++, sizeof(dB_offset), &dB_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(lddb     ), &lddb      );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 2, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    }
    else if ( uplo == MagmaUpper ) {
        kernel = g_runtime.get_kernel( "zlacpy_kernel_upper" );
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            err |= clSetKernelArg( kernel, i++, sizeof(m        ), &m         );
            err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
            err |= clSetKernelArg( kernel, i++, sizeof(dA       ), &dA        );
            err |= clSetKernelArg( kernel, i++, sizeof(dA_offset), &dA_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldda     ), &ldda      );
            err |= clSetKernelArg( kernel, i++, sizeof(dB       ), &dB        );
            err |= clSetKernelArg( kernel, i++, sizeof(dB_offset), &dB_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(lddb     ), &lddb      );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 2, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    }
    else {
        kernel = g_runtime.get_kernel( "zlacpy_kernel_full" );
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            err |= clSetKernelArg( kernel, i++, sizeof(m        ), &m         );
            err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
            err |= clSetKernelArg( kernel, i++, sizeof(dA       ), &dA        );
            err |= clSetKernelArg( kernel, i++, sizeof(dA_offset), &dA_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldda     ), &ldda      );
            err |= clSetKernelArg( kernel, i++, sizeof(dB       ), &dB        );
            err |= clSetKernelArg( kernel, i++, sizeof(dB_offset), &dB_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(lddb     ), &lddb      );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 2, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    }
}


/**
    Purpose
    -------
    ZLACPY_BATCHED_Q copies all or part of each two-dimensional matrix
    dAarray[i] to matrix dBarray[i], for 0 <= i < batchcount.
    
    This is the same as ZLACPY_BATCHED, but adds queue argument.
    
    Arguments
    ---------
    
    @param[in]
    uplo    magma_uplo_t
            Specifies the part of each matrix dA to be copied to dB.
      -     = MagmaUpper:      Upper triangular part
      -     = MagmaLower:      Lower triangular part
            Otherwise:  All of each matrix dA
    
    @param[in]
    m       INTEGER
            The number of rows of each matrix dA.  M >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of each matrix dA.  N >= 0.
    
    @param[in]
    dAarray COMPLEX_16* array, dimension (batchCount)
            array of pointers to the matrices dA, where each dA is of dimension (LDDA,N)
            The m by n matrix dA.
            If UPLO = MagmaUpper, only the upper triangle or trapezoid is accessed;
            if UPLO = MagmaLower, only the lower triangle or trapezoid is accessed.
    
    @param[in]
    ldda    INTEGER
            The leading dimension of each array dA.  LDDA >= max(1,M).
    
    @param[out]
    dBarray COMPLEX_16* array, dimension (batchCount)
            array of pointers to the matrices dB, where each dB is of dimension (LDDB,N)
            The m by n matrix dB.
            On exit, dB = dA in the locations specified by UPLO.
    
    @param[in]
    lddb    INTEGER
            The leading dimension of each array dB.  LDDB >= max(1,M).
    
    @param[in]
    batchCount  Number of matrices in dAarray and dBarray.
    
    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    @ingroup magma_zaux2
    ********************************************************************/
extern "C" void
magmablas_zlacpy_batched(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr const *dAarray, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr             *dBarray, size_t dB_offset, magma_int_t lddb,
    magma_int_t batchCount, magma_queue_t queue )
{
    cl_kernel kernel;
    cl_int err;
    int i;

    magma_int_t info = 0;
    if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( ldda < max(1,m))
        info = -5;
    else if ( lddb < max(1,m))
        info = -7;
    else if ( batchCount < 0 )
        info = -8;
    
    if ( info != 0 ) {
        magma_xerbla( __func__, -(info) );
        return;
    }
    
    if ( m == 0 || n == 0 || batchCount == 0 )
        return;
    
    size_t threads[3] = { BLK_X, 1, 1 };
    size_t grid[3] = { (m + BLK_X - 1)/BLK_X, (n + BLK_Y - 1)/BLK_Y, batchCount };
    grid[0] *= threads[0];
    grid[1] *= threads[1];
    grid[2] *= threads[2];
    
    if ( uplo == MagmaLower ) {
        kernel = g_runtime.get_kernel( "zlacpy_kernel_batched_lower" );
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            err |= clSetKernelArg( kernel, i++, sizeof(m        ), &m         );
            err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
            err |= clSetKernelArg( kernel, i++, sizeof(dAarray  ), &dAarray   );
            err |= clSetKernelArg( kernel, i++, sizeof(dA_offset), &dA_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldda     ), &ldda      );
            err |= clSetKernelArg( kernel, i++, sizeof(dBarray  ), &dBarray   );
            err |= clSetKernelArg( kernel, i++, sizeof(dB_offset), &dB_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(lddb     ), &lddb      );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 3, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    }
    else if ( uplo == MagmaUpper ) {
        kernel = g_runtime.get_kernel( "zlacpy_kernel_batched_upper" );
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            err |= clSetKernelArg( kernel, i++, sizeof(m        ), &m         );
            err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
            err |= clSetKernelArg( kernel, i++, sizeof(dAarray  ), &dAarray   );
            err |= clSetKernelArg( kernel, i++, sizeof(dA_offset), &dA_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldda     ), &ldda      );
            err |= clSetKernelArg( kernel, i++, sizeof(dBarray  ), &dBarray   );
            err |= clSetKernelArg( kernel, i++, sizeof(dB_offset), &dB_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(lddb     ), &lddb      );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 3, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    }
    else {
        kernel = g_runtime.get_kernel( "zlacpy_kernel_batched_full" );
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            err |= clSetKernelArg( kernel, i++, sizeof(m        ), &m         );
            err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
            err |= clSetKernelArg( kernel, i++, sizeof(dAarray  ), &dAarray   );
            err |= clSetKernelArg( kernel, i++, sizeof(dA_offset), &dA_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldda     ), &ldda      );
            err |= clSetKernelArg( kernel, i++, sizeof(dBarray  ), &dBarray   );
            err |= clSetKernelArg( kernel, i++, sizeof(dB_offset), &dB_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(lddb     ), &lddb      );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 3, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    }
}
