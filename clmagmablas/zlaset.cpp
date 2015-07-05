/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @author Mark Gates
       
       @precisions normal z -> s d c

       auto-converted from zlaset.cu

*/
#include "clmagma_runtime.h"
#include "common_magma.h"
#include "zlaset.h"


/**
    Purpose
    -------
    ZLASET_STREAM initializes a 2-D array A to DIAG on the diagonal and
    OFFDIAG on the off-diagonals.
    
    This is the same as ZLASET, but adds queue argument.
    
    Arguments
    ---------
    @param[in]
    uplo    magma_uplo_t
            Specifies the part of the matrix dA to be set.
      -     = MagmaUpper:      Upper triangular part
      -     = MagmaLower:      Lower triangular part
            Otherwise:         All of the matrix dA is set.
    
    @param[in]
    m       INTEGER
            The number of rows of the matrix dA.  M >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix dA.  N >= 0.
    
    @param[in]
    offdiag COMPLEX_16
            The scalar OFFDIAG. (In LAPACK this is called ALPHA.)
    
    @param[in]
    diag    COMPLEX_16
            The scalar DIAG. (In LAPACK this is called BETA.)
    
    @param[in]
    dA      COMPLEX_16 array, dimension (LDDA,N)
            The M-by-N matrix dA.
            If UPLO = MagmaUpper, only the upper triangle or trapezoid is accessed;
            if UPLO = MagmaLower, only the lower triangle or trapezoid is accessed.
            On exit, A(i,j) = OFFDIAG, 1 <= i <= m, 1 <= j <= n, i != j;
                     A(i,i) = DIAG,    1 <= i <= min(m,n)
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
    
    @param[in]
    queue   magma_queue_t
            Queue to execute in.
    
    @ingroup magma_zaux2
    ********************************************************************/
extern "C"
void magmablas_zlaset(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magmaDoubleComplex offdiag, magmaDoubleComplex diag,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue)
{
    cl_kernel kernel;
    cl_int err;
    int i;

    magma_int_t info = 0;
    if ( uplo != MagmaLower && uplo != MagmaUpper && uplo != MagmaFull )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( ldda < max(1,m) )
        info = -7;
    
    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    if ( m == 0 || n == 0 ) {
        return;
    }
    
    size_t threads[2] = { BLK_X, 1 };
    size_t grid[2] = { (m-1)/BLK_X + 1, (n-1)/BLK_Y + 1 };
    grid[0] *= threads[0];
    grid[1] *= threads[1];
    
    if (uplo == MagmaLower) {
        kernel = g_runtime.get_kernel( "zlaset_lower" );
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            err |= clSetKernelArg( kernel, i++, sizeof(m        ), &m         );
            err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
            err |= clSetKernelArg( kernel, i++, sizeof(offdiag  ), &offdiag   );
            err |= clSetKernelArg( kernel, i++, sizeof(diag     ), &diag      );
            err |= clSetKernelArg( kernel, i++, sizeof(dA       ), &dA        );
            err |= clSetKernelArg( kernel, i++, sizeof(dA_offset), &dA_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldda     ), &ldda      );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 2, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    }
    else if (uplo == MagmaUpper) {
        kernel = g_runtime.get_kernel( "zlaset_upper" );
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            err |= clSetKernelArg( kernel, i++, sizeof(m        ), &m         );
            err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
            err |= clSetKernelArg( kernel, i++, sizeof(offdiag  ), &offdiag   );
            err |= clSetKernelArg( kernel, i++, sizeof(diag     ), &diag      );
            err |= clSetKernelArg( kernel, i++, sizeof(dA       ), &dA        );
            err |= clSetKernelArg( kernel, i++, sizeof(dA_offset), &dA_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldda     ), &ldda      );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 2, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    }
    else {
        kernel = g_runtime.get_kernel( "zlaset_full" );
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            err |= clSetKernelArg( kernel, i++, sizeof(m        ), &m         );
            err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
            err |= clSetKernelArg( kernel, i++, sizeof(offdiag  ), &offdiag   );
            err |= clSetKernelArg( kernel, i++, sizeof(diag     ), &diag      );
            err |= clSetKernelArg( kernel, i++, sizeof(dA       ), &dA        );
            err |= clSetKernelArg( kernel, i++, sizeof(dA_offset), &dA_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldda     ), &ldda      );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 2, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    }
}
