/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @precisions mixed zc -> ds

       auto-converted from zlat2c.cu
       @author Mark Gates
*/
#include "clmagma_runtime.h"
#include "common_magma.h"
#include "zlat2c.h"


/**
    Purpose
    -------
    ZLAT2C converts a double-complex matrix, A,
                 to a single-complex matrix, SA.
    
    RMAX is the overflow for the single-complex arithmetic.
    ZLAT2C checks that all the entries of A are between -RMAX and
    RMAX. If not, the conversion is aborted and a flag is raised.
        
    Arguments
    ---------
    @param[in]
    uplo    magma_uplo_t
            Specifies the part of the matrix A to be converted.
      -     = MagmaUpper:      Upper triangular part
      -     = MagmaLower:      Lower triangular part
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  n >= 0.
    
    @param[in]
    A       COMPLEX_16 array, dimension (LDA,n)
            On entry, the n-by-n coefficient matrix A.
    
    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,n).
    
    @param[out]
    SA      COMPLEX array, dimension (LDSA,n)
            On exit, if INFO=0, the n-by-n coefficient matrix SA;
            if INFO > 0, the content of SA is unspecified.
    
    @param[in]
    ldsa    INTEGER
            The leading dimension of the array SA.  LDSA >= max(1,n).
    
    @param[out]
    info    INTEGER
      -     = 0:  successful exit.
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     = 1:  an entry of the matrix A is greater than the COMPLEX
                  overflow threshold, in this case, the content
                  of SA on exit is unspecified.
    
    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    @ingroup magma_zaux2
    ********************************************************************/
extern "C" void
magmablas_zlat2c(
    magma_uplo_t uplo, magma_int_t n,
    magmaDoubleComplex_const_ptr  A, size_t A_offset, magma_int_t lda,
    magmaFloatComplex_ptr        SA, size_t SA_offset, magma_int_t ldsa,
    magma_queue_t queue,
    magma_int_t *info )
{
    cl_kernel kernel;
    cl_int err;
    int i;

    magma_ptr dflag;
    magma_imalloc(&dflag, 1); 

    *info = 0;
    if ( uplo != MagmaLower && uplo != MagmaUpper )
        *info = -1;
    else if ( n < 0 )
        *info = -2;
    else if ( lda < max(1,n) )
        *info = -4;
    else if ( ldsa < max(1,n) )
        *info = -6;
    
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return; //*info;
    }

    /* quick return */
    if ( n == 0 ) {
        return;
    }
    
    double rmax = (double)lapackf77_slamch("O");

    size_t threads[2] = { BLK_X, 1 };
    size_t    grid[2] = { (n+BLK_X-1)/BLK_X, (n+BLK_Y-1)/BLK_Y };
    grid[0] *= threads[0];
    grid[1] *= threads[1];
    // cudaMemcpyToSymbol( flag, info, sizeof(flag) );    // flag = 0
    magma_setvector(1, sizeof(magma_int_t), info, 1, dflag,0, 1, queue);
    
    if (uplo == MagmaLower)
        kernel = g_runtime.get_kernel( "zlat2c_lower" );
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            err  = clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
            err |= clSetKernelArg( kernel, i++, sizeof(A        ), &A         );
            err |= clSetKernelArg( kernel, i++, sizeof(A_offset ), &A_offset  );
            err |= clSetKernelArg( kernel, i++, sizeof(lda      ), &lda       );
            err |= clSetKernelArg( kernel, i++, sizeof(SA       ), &SA        );
            err |= clSetKernelArg( kernel, i++, sizeof(SA_offset), &SA_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldsa     ), &ldsa      );
            err |= clSetKernelArg( kernel, i++, sizeof(rmax     ), &rmax      );
            err |= clSetKernelArg( kernel, i++, sizeof(dflag    ), &dflag     );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 2, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    else if (uplo == MagmaUpper)
        kernel = g_runtime.get_kernel( "zlat2c_upper" );
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            err  = clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
            err |= clSetKernelArg( kernel, i++, sizeof(A        ), &A         );
            err |= clSetKernelArg( kernel, i++, sizeof(A_offset ), &A_offset  );
            err |= clSetKernelArg( kernel, i++, sizeof(lda      ), &lda       );
            err |= clSetKernelArg( kernel, i++, sizeof(SA       ), &SA        );
            err |= clSetKernelArg( kernel, i++, sizeof(SA_offset), &SA_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldsa     ), &ldsa      );
            err |= clSetKernelArg( kernel, i++, sizeof(rmax     ), &rmax      );
            err |= clSetKernelArg( kernel, i++, sizeof(dflag    ), &dflag     );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 2, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    
        // cudaMemcpyFromSymbol( info, flag, sizeof(flag) );  // info = flag
        magma_getvector(1, sizeof(magma_int_t), dflag,0, 1, info, 1, queue);

        magma_free(dflag);
}
