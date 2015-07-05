/*
   -- clMAGMA (version 1.3.0) --
   Univ. of Tennessee, Knoxville
   Univ. of California, Berkeley
   Univ. of Colorado, Denver
   @date November 2014

   @generated from zpotri_gpu.cpp normal z -> c, Sat Nov 15 00:21:37 2014

 */
#include "common_magma.h"

#define dA(i, j)  dA, (offset_dA + (j)*ldda + (i))

extern "C" magma_int_t
magma_cpotri_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex_ptr dA, size_t offset_dA, magma_int_t ldda,
    magma_queue_t queues[2],
    magma_int_t *info)
{
/*  -- MAGMA (version 1.3.0) --
    Univ. of Tennessee, Knoxville
    Univ. of California, Berkeley
    Univ. of Colorado, Denver
    @date November 2014

    Purpose
    =======
    CPOTRI computes the inverse of a real symmetric positive definite
    matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
    computed by CPOTRF.

    Arguments
    =========
    UPLO    (input) CHARACTER*1
            = 'U':  Upper triangle of A is stored;
            = 'L':  Lower triangle of A is stored.

    N       (input) INTEGER
            The order of the matrix A.  N >= 0.

    dA      (input/output) COMPLEX array on the GPU, dimension (LDDA,N)
            On entry, the triangular factor U or L from the Cholesky
            factorization A = U**T*U or A = L*L**T, as computed by
            CPOTRF.
            On exit, the upper or lower triangle of the (symmetric)
            inverse of A, overwriting the input factor U or L.

    LDDA    (input) INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).
            INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
            > 0:  if INFO = i, the (i,i) element of the factor U or L is
            zero, and the inverse could not be computed.

    ===================================================================== */

    /* Local variables */
    *info = 0;
    if ((uplo != MagmaUpper) && (uplo != MagmaLower))
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (ldda < max(1,n))
        *info = -4;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if ( n == 0 )
        return *info;

    /* Invert the triangular Cholesky factor U or L */
    magma_ctrtri_gpu( uplo, MagmaNonUnit, n, dA, offset_dA, ldda, queues, info );
    
    if ( *info == 0 ) {
        /* Form inv(U) * inv(U)**T or inv(L)**T * inv(L) */
        magma_clauum_gpu( uplo, n, dA, offset_dA, ldda, queues[0], info );
    }

    return *info;
} /* magma_cpotri */
