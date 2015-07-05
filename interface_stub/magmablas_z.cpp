/*
 *   -- clMAGMA (version 0.2.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @precisions normal z -> s d c
 */

// defines stub functions, if no GPU implementation is available.
#include <stdio.h>
#include <stdlib.h>

#include "magma.h"

#if !defined(HAVE_CUBLAS) && !defined(HAVE_clAmdBlas) && !defined(HAVE_CBLAS)

// ========================================
// copying sub-matrices (contiguous columns)
magma_err_t
magma_zsetmatrix(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex const* hA_src, size_t hA_offset, magma_int_t ldha,
    magmaDoubleComplex_ptr    dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue )
{
    printf( "%s( %d, %d, ptr, %ld, %d, ptr, %ld, %d, queue )\n", __func__,
        m, n, hA_offset, ldha, dA_offset, ldda );
    return 0;
}

// --------------------
magma_err_t
magma_zgetmatrix(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex*          hA_dst, size_t hA_offset, magma_int_t ldha,
    magma_queue_t queue )
{
    printf( "%s( %d, %d, ptr, %ld, %d, ptr, %ld, %d, queue )\n", __func__,
        m, n, dA_offset, ldda, hA_offset, ldha );
    return 0;
}

// --------------------
magma_err_t
magma_zsetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex const* hA_src, size_t hA_offset, magma_int_t ldha,
    magmaDoubleComplex_ptr    dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue )
{
    printf( "%s( %d, %d, ptr, %ld, %d, ptr, %ld, %d, queue )\n", __func__,
        m, n, hA_offset, ldha, dA_offset, ldda );
    return 0;
}

// --------------------
magma_err_t
magma_zgetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex*          hA_dst, size_t hA_offset, magma_int_t ldha,
    magma_queue_t queue )
{
    printf( "%s( %d, %d, ptr, %ld, %d, ptr, %ld, %d, queue )\n", __func__,
        m, n, dA_offset, ldda, hA_offset, ldha );
    return 0;
}

// ========================================
// BLAS functions
magma_err_t
magma_zgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaDoubleComplex alpha, magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDoubleComplex_const_ptr dB, size_t dB_offset, magma_int_t ldb,
    magmaDoubleComplex beta,  magmaDoubleComplex_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t handle )
{
    printf( "%s( %d, %d, %d, %d, %d, %f, ptr, %ld, %d, ptr, %ld, %d, %f, ptr, %ld, %d, handle )\n",
        __func__,
        transA, transB, m, n, k,
        MAGMA_Z_REAL(alpha), dA_offset, lda,
                             dB_offset, ldb,
        MAGMA_Z_REAL(beta),  dC_offset, ldc );
    return 0;
}

// --------------------
magma_err_t
magma_zgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex alpha, magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDoubleComplex_const_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaDoubleComplex beta,  magmaDoubleComplex_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t handle )
{
    printf( "%s( %d,      %d, %d,    %f, ptr, %ld, %d, ptr, %ld, %d, %f, ptr, %ld, %d, handle )\n",
        __func__,
        transA, m, n,
        MAGMA_Z_REAL(alpha), dA_offset, lda,
                             dx_offset, incx,
        MAGMA_Z_REAL(beta),  dy_offset, incy );
    return 0;
}

// --------------------
magma_err_t
magma_zhemm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex alpha, magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDoubleComplex_const_ptr dB, size_t dB_offset, magma_int_t ldb,
    magmaDoubleComplex beta,  magmaDoubleComplex_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t handle )
{
    printf( "%s( %d, %d, %d, %d,    %f, ptr, %ld, %d, ptr, %ld, %d, %f, ptr, %ld, %d, handle )\n",
        __func__,
        side, uplo, m, n,
        MAGMA_Z_REAL(alpha), dA_offset, lda,
                             dB_offset, ldb,
        MAGMA_Z_REAL(beta),  dC_offset, ldc );
    return 0;
}

// --------------------
magma_err_t
magma_zhemv(
    magma_uplo_t uplo,
    magma_int_t n,
    magmaDoubleComplex alpha, magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDoubleComplex_const_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaDoubleComplex beta,  magmaDoubleComplex_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t handle )
{
    printf( "%s(      %d, %d,       %f, ptr, %ld, %d, ptr, %ld, %d, %f, ptr, %ld, %d, handle )\n",
        __func__,
        uplo, n,
        MAGMA_Z_REAL(alpha), dA_offset, lda,
                             dx_offset, incx,
        MAGMA_Z_REAL(beta),  dy_offset, incy );
    return 0;
}

// --------------------
magma_err_t
magma_zherk(
    magma_uplo_t uplo,
    magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha, magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
    double beta,  magmaDoubleComplex_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t handle )
{
    printf( "%s(   %d, %d, %d, %d,   %f, ptr, %ld, %d, %f, ptr, %ld, %d, handle )\n",
        __func__,
        uplo, trans, n, k,
        alpha, dA_offset, lda,
        beta,  dC_offset, ldc );
    return 0;
}

// --------------------
magma_err_t
magma_ztrsm(
    magma_side_t  side,
    magma_uplo_t  uplo,
    magma_trans_t trans,
    magma_diag_t  diag,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex alpha, magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDoubleComplex_ptr       dB, size_t dB_offset, magma_int_t ldb,
    magma_queue_t handle )
{
    printf( "%s( %d, %d, %d, %d, %d, %d, %f, ptr, %ld, %d, ptr, %ld, %d, handle )\n",
        __func__,
        side, uplo, trans, diag, m, n,
        MAGMA_Z_REAL(alpha), dA_offset, lda,
                             dB_offset, ldb );
    return 0;
}

#endif // not HAVE_CUBLAS and not HAVE_clAmdBlas and not HAVE_CBLAS
