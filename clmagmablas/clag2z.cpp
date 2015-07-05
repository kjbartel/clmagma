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
#include "clmagma_runtime.h"
#include "common_magma.h"
#include "clag2z.h"


/**
    Purpose
    -------
    CLAG2Z_STREAM converts a single-complex matrix, SA,
                        to a double-complex matrix, A.

    Note that while it is possible to overflow while converting
    from double to single, it is not possible to overflow when
    converting from single to double.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of lines of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in]
    SA      REAL array, dimension (LDSA,N)
            On entry, the M-by-N coefficient matrix SA.

    @param[in]
    ldsa    INTEGER
            The leading dimension of the array SA.  LDSA >= max(1,M).

    @param[out]
    A       DOUBLE PRECISION array, dimension (LDA,N)
            On exit, the M-by-N coefficient matrix A.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
    
    @param[in]
    queue   magma_queue_t
            Queue to execute in.
    
    @ingroup magma_caux2
    ********************************************************************/
extern "C" void
magmablas_clag2z(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr SA, size_t SA_offset, magma_int_t ldsa,
    magmaDoubleComplex_ptr       A, size_t A_offset, magma_int_t lda,
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
    else if ( ldsa < max(1,m) )
        *info = -4;
    else if ( lda < max(1,m) )
        *info = -6;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return; //*info;
    }

    /* quick return */
    if ( m == 0 || n == 0 ) {
        return;
    }

    size_t threads[2] = { BLK_X, 1 };
    size_t grid[2] = { (m+BLK_X-1)/BLK_X, (n+BLK_Y-1)/BLK_Y };
    grid[0] *= threads[0];
    grid[1] *= threads[1];
    kernel = g_runtime.get_kernel( "clag2z_kernel" );
    if ( kernel != NULL ) {
        err = 0;
        i   = 0;
        err |= clSetKernelArg( kernel, i++, sizeof(m        ), &m         );
        err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
        err |= clSetKernelArg( kernel, i++, sizeof(SA       ), &SA        );
        err |= clSetKernelArg( kernel, i++, sizeof(SA_offset), &SA_offset );
        err |= clSetKernelArg( kernel, i++, sizeof(ldsa     ), &ldsa      );
        err |= clSetKernelArg( kernel, i++, sizeof(A        ), &A         );
        err |= clSetKernelArg( kernel, i++, sizeof(A_offset ), &A_offset  );
        err |= clSetKernelArg( kernel, i++, sizeof(lda      ), &lda       );
        check_error( err );

        err = clEnqueueNDRangeKernel( queue, kernel, 2, NULL, grid, threads, 0, NULL, NULL );
        check_error( err );
    }
}
