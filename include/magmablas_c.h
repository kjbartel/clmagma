/*
 *   -- clMAGMA (version 0.2.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @generated c Thu May 24 17:09:40 2012
 */

#ifndef MAGMA_BLAS_C_H
#define MAGMA_BLAS_C_H

#include "magma_types.h"

#ifdef __cplusplus
extern "C" {
#endif

// ========================================
// copying sub-matrices (contiguous columns)
magma_err_t
magma_csetmatrix(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex const* hA_src, size_t hA_offset, magma_int_t ldha,
    magmaFloatComplex_ptr    dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

magma_err_t
magma_cgetmatrix(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex*          hA_dst, size_t hA_offset, magma_int_t ldha,
    magma_queue_t queue );

magma_err_t
magma_csetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex const* hA_src, size_t hA_offset, magma_int_t ldha,
    magmaFloatComplex_ptr    dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue, magma_event_t *event );

magma_err_t
magma_cgetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex*          hA_dst, size_t hA_offset, magma_int_t ldha,
    magma_queue_t queue, magma_event_t *event );

magma_err_t
magma_ccopymatrix(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr    dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

// ========================================
// matrix transpose and swapping functions
magma_err_t
magma_cinplace_transpose(
    magmaFloatComplex_ptr dA, size_t dA_offset, int lda, int n,
    magma_queue_t queue );

magma_err_t
magma_ctranspose2(
    magmaFloatComplex_ptr odata, size_t odata_offset, int ldo,
    magmaFloatComplex_ptr idata, size_t idata_offset, int ldi,
    int m, int n,
    magma_queue_t queue );

magma_err_t
magma_ctranspose(
    magmaFloatComplex_ptr odata, int offo, int ldo,
    magmaFloatComplex_ptr idata, int offi, int ldi,
    int m, int n,
    magma_queue_t queue );

magma_err_t
magma_cpermute_long2(
    magmaFloatComplex_ptr dAT, size_t dAT_offset, int lda,
    int *ipiv, int nb, int ind,
    magma_queue_t queue );


// ========================================
// BLAS functions
magma_err_t
magma_cgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloatComplex alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaFloatComplex_const_ptr dB, size_t dB_offset, magma_int_t ldb,
    magmaFloatComplex beta,  magmaFloatComplex_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t queue );

magma_err_t
magma_cgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaFloatComplex_const_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaFloatComplex beta,  magmaFloatComplex_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

magma_err_t
magma_chemm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaFloatComplex_const_ptr dB, size_t dB_offset, magma_int_t ldb,
    magmaFloatComplex beta,  magmaFloatComplex_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t queue );

magma_err_t
magma_chemv(
    magma_uplo_t uplo,
    magma_int_t n,
    magmaFloatComplex alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaFloatComplex_const_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaFloatComplex beta,  magmaFloatComplex_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

magma_err_t
magma_cherk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
    float beta,  magmaFloatComplex_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t queue );

magma_err_t
magma_ctrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaFloatComplex_ptr       dB, size_t dB_offset, magma_int_t ldb,
    magma_queue_t queue );

magma_err_t
magma_ctrsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag, 
	magma_int_t n, 
	magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda, 
	magmaFloatComplex_ptr dx, size_t dx_offset, magma_int_t incx,
	magma_queue_t queue );

magma_err_t
magma_ctrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaFloatComplex_ptr       dB, size_t dB_offset, magma_int_t ldb,
    magma_queue_t queue );

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_BLAS_H
