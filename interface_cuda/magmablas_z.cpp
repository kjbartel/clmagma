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

#include <stdlib.h>
#include <stdio.h>

#include "magma.h"

#ifdef HAVE_CUBLAS

// ========================================
// copying sub-matrices (contiguous columns)
magma_err_t
magma_zsetmatrix(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex const* hA_src, size_t hA_offset, magma_int_t ldha,
    magmaDoubleComplex_ptr    dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue )
{
    cublasStatus_t stat;
    stat = cublasSetMatrix(
        m, n, sizeof(magmaDoubleComplex),
        hA_src + hA_offset, ldha,
        dA_dst + dA_offset, ldda );
    return stat;
}

// --------------------
magma_err_t
magma_zgetmatrix(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex*          hA_dst, size_t hA_offset, magma_int_t ldha,
    magma_queue_t queue )
{
    cublasStatus_t stat;
    stat = cublasGetMatrix(
        m, n, sizeof(magmaDoubleComplex),
        dA_src + dA_offset, ldda,
        hA_dst + hA_offset, ldha );
    return stat;
}

// --------------------
magma_err_t
magma_zsetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex const* hA_src, size_t hA_offset, magma_int_t ldha,
    magmaDoubleComplex_ptr    dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue )
{
    cudaStream_t stream;
    cublasStatus_t stat;
    stat = cublasGetStream( queue, &stream );
    stat = cublasSetMatrixAsync(
        m, n, sizeof(magmaDoubleComplex),
        hA_src + hA_offset, ldha,
        dA_dst + dA_offset, ldda, stream );
    return stat;
}

// --------------------
magma_err_t
magma_zgetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex*          hA_dst, size_t hA_offset, magma_int_t ldha,
    magma_queue_t queue )
{
    cudaStream_t stream;
    cublasStatus_t stat;
    stat = cublasGetStream( queue, &stream );
    stat = cublasGetMatrixAsync(
        m, n, sizeof(magmaDoubleComplex),
        dA_src + dA_offset, ldda,
        hA_dst + hA_offset, ldha, stream );
    return stat;
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
    magma_queue_t queue )
{
    cublasStatus_t stat;
    stat = cublasZgemm(
        queue,
        cublas_trans_const( transA ),
        cublas_trans_const( transB ),
        m, n, k,
        &alpha, dA + dA_offset, lda,
                dB + dB_offset, ldb,
        &beta,  dC + dC_offset, ldc );
    return stat;
}

// --------------------
magma_err_t
magma_zgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex alpha, magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDoubleComplex_const_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaDoubleComplex beta,  magmaDoubleComplex_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue )
{
    cublasStatus_t stat;
    stat = cublasZgemv(
        queue,
        cublas_trans_const( transA ),
        m, n,
        &alpha, dA + dA_offset, lda,
                dx + dx_offset, incx,
        &beta,  dy + dy_offset, incy );
    return stat;
}

// --------------------
magma_err_t
magma_zhemm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex alpha, magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDoubleComplex_const_ptr dB, size_t dB_offset, magma_int_t ldb,
    magmaDoubleComplex beta,  magmaDoubleComplex_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t queue )
{
    cublasStatus_t stat;
    stat = cublasZhemm(
        queue,
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        m, n,
        &alpha, dA + dA_offset, lda,
                dB + dB_offset, ldb,
        &beta,  dC + dC_offset, ldc );
    return stat;
}

// --------------------
magma_err_t
magma_zhemv(
    magma_uplo_t uplo,
    magma_int_t n,
    magmaDoubleComplex alpha, magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDoubleComplex_const_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaDoubleComplex beta,  magmaDoubleComplex_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue )
{
    cublasStatus_t stat;
    stat = cublasZhemv(
        queue,
        cublas_uplo_const( uplo ),
        n,
        &alpha, dA + dA_offset, lda,
                dx + dx_offset, incx,
        &beta,  dy + dy_offset, incy );
    return stat;
}

// --------------------
magma_err_t
magma_zherk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha, magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
    double beta,  magmaDoubleComplex_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t queue )
{
    cublasStatus_t stat;
    stat = cublasZherk(
        queue,
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        &alpha, dA + dA_offset, lda,
        &beta,  dC + dC_offset, ldc );
    return stat;
}

// --------------------
magma_err_t
magma_ztrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex alpha, magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDoubleComplex_ptr       dB, size_t dB_offset, magma_int_t ldb,
    magma_queue_t queue )
{
    cublasStatus_t stat;
    stat = cublasZtrsm(
        queue,
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        cublas_diag_const( diag ),
        m, n,
        &alpha, dA + dA_offset, lda,
                dB + dB_offset, ldb );
    return stat;
}

#endif // HAVE_CUBLAS
