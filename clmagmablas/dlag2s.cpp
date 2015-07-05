/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zlag2c.cpp mixed zc -> ds, Sat Nov 15 00:21:35 2014

       auto-converted from dlag2s.cu
       @author Mark Gates
*/
#include "clmagma_runtime.h"
#include "common_magma.h"
#include "dlag2s.h"


/**
    Purpose
    -------
    DLAG2S_STREAM converts a double-real matrix, A,
                        to a single-real matrix, SA.
    
    RMAX is the overflow for the single-real arithmetic.
    DLAG2S checks that all the entries of A are between -RMAX and
    RMAX. If not, the conversion is aborted and a flag is raised.
    
    This is the same as DLAG2S, but adds queue argument.
        
    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of lines of the matrix A.  m >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  n >= 0.
    
    @param[in]
    A       DOUBLE PRECISION array, dimension (LDA,n)
            On entry, the m-by-n coefficient matrix A.
    
    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,m).
    
    @param[out]
    SA      SINGLE PRECISION array, dimension (LDSA,n)
            On exit, if INFO=0, the m-by-n coefficient matrix SA;
            if INFO > 0, the content of SA is unspecified.
    
    @param[in]
    ldsa    INTEGER
            The leading dimension of the array SA.  LDSA >= max(1,m).
    
    @param[out]
    info    INTEGER
      -     = 0:  successful exit.
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     = 1:  an entry of the matrix A is greater than the SINGLE PRECISION
                  overflow threshold, in this case, the content
                  of SA on exit is unspecified.
    
    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    @ingroup magma_daux2
    ********************************************************************/
extern "C" void
magmablas_dlag2s(
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr A, size_t A_offset, magma_int_t lda,
    magmaFloat_ptr SA, size_t SA_offset,       magma_int_t ldsa,
    magma_queue_t queue,
    magma_int_t *info )
{
    cl_kernel kernel;
    cl_int err;
    int i;

    *info = 0;
    if ( m < 0 )
        *info = -1;
    else if ( n < 0 )
        *info = -2;
    else if ( lda < max(1,m) )
        *info = -4;
    else if ( ldsa < max(1,m) )
        *info = -6;
    
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return; //*info;
    }

    /* quick return */
    if ( m == 0 || n == 0 ) {
        return;
    }
    
    double rmax = (double)lapackf77_slamch("O");

    size_t threads[2] = { BLK_X, 1 };
    size_t grid[2] = { (m+BLK_X-1)/BLK_X, (n+BLK_Y-1)/BLK_Y };
    grid[0] *= threads[0];
    grid[1] *= threads[1];
    
    // TODO cudaMemcpyToSymbol( flag, info, sizeof(flag) );    // flag = 0
    
    kernel = g_runtime.get_kernel( "dlag2s_kernel" );
    if ( kernel != NULL ) {
        err = 0;
        i   = 0;
        err |= clSetKernelArg( kernel, i++, sizeof(m        ), &m         );
        err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
        err |= clSetKernelArg( kernel, i++, sizeof(A        ), &A         );
        err |= clSetKernelArg( kernel, i++, sizeof(A_offset ), &A_offset  );
        err |= clSetKernelArg( kernel, i++, sizeof(lda      ), &lda       );
        err |= clSetKernelArg( kernel, i++, sizeof(SA       ), &SA        );
        err |= clSetKernelArg( kernel, i++, sizeof(SA_offset), &SA_offset );
        err |= clSetKernelArg( kernel, i++, sizeof(ldsa     ), &ldsa      );
        err |= clSetKernelArg( kernel, i++, sizeof(rmax     ), &rmax      );
        check_error( err );

        err = clEnqueueNDRangeKernel( queue, kernel, 2, NULL, grid, threads, 0, NULL, NULL );
        check_error( err );
    }
    
    // TODO cudaMemcpyFromSymbol( info, flag, sizeof(flag) );  // info = flag
}
