/*
 *   -- clMAGMA (version 0.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @generated d Wed Apr  4 01:12:51 2012
 */

#ifndef MAGMA_BLAS_D_H
#define MAGMA_BLAS_D_H

#include "magma_types.h"

#ifdef __cplusplus
extern "C" {
#endif

// ========================================
// copying sub-matrices (contiguous columns)
magma_err_t
magma_dsetmatrix(
	magma_int_t m, magma_int_t n,
	double const* hA_src, size_t hA_offset, magma_int_t ldha,
	magmaDouble_ptr    dA_dst, size_t dA_offset, magma_int_t ldda,
	magma_queue_t queue );

magma_err_t
magma_dgetmatrix(
	magma_int_t m, magma_int_t n,
	magmaDouble_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
	double*          hA_dst, size_t hA_offset, magma_int_t ldha,
	magma_queue_t queue );

magma_err_t
magma_dsetmatrix_async(
	magma_int_t m, magma_int_t n,
	double const* hA_src, size_t hA_offset, magma_int_t ldha,
	magmaDouble_ptr    dA_dst, size_t dA_offset, magma_int_t ldda,
	magma_queue_t queue, magma_event_t *event );

magma_err_t
magma_dgetmatrix_async(
	magma_int_t m, magma_int_t n,
	magmaDouble_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
	double*          hA_dst, size_t hA_offset, magma_int_t ldha,
	magma_queue_t queue, magma_event_t *event );


// ========================================
// matrix transpose and swapping functions
magma_err_t
magma_dinplace_transpose(
    magmaDouble_ptr dA, size_t dA_offset, int lda, int n,
    magma_queue_t queue );

magma_err_t
magma_dtranspose2(
    magmaDouble_ptr odata, size_t odata_offset, int ldo,
    magmaDouble_ptr idata, size_t idata_offset, int ldi,
    int m, int n,
    magma_queue_t queue );

magma_err_t
magma_dtranspose(
    magmaDouble_ptr odata, int offo, int ldo,
    magmaDouble_ptr idata, int offi, int ldi,
    int m, int n,
    magma_queue_t queue );

magma_err_t
magma_dpermute_long2(
    magmaDouble_ptr dAT, size_t dAT_offset, int lda,
    int *ipiv, int nb, int ind,
    magma_queue_t queue );


// ========================================
// BLAS functions
magma_err_t
magma_dgemm(
	magma_trans_t transA, magma_trans_t transB,
	magma_int_t m, magma_int_t n, magma_int_t k,
	double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
	                          magmaDouble_const_ptr dB, size_t dB_offset, magma_int_t ldb,
	double beta,  magmaDouble_ptr       dC, size_t dC_offset, magma_int_t ldc,
	magma_queue_t queue );

magma_err_t
magma_dgemv(
	magma_trans_t transA,
	magma_int_t m, magma_int_t n,
	double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
	                          magmaDouble_const_ptr dx, size_t dx_offset, magma_int_t incx,
	double beta,  magmaDouble_ptr       dy, size_t dy_offset, magma_int_t incy,
	magma_queue_t queue );

magma_err_t
magma_dsymm(
	magma_side_t side, magma_uplo_t uplo,
	magma_int_t m, magma_int_t n,
	double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
	                          magmaDouble_const_ptr dB, size_t dB_offset, magma_int_t ldb,
	double beta,  magmaDouble_ptr       dC, size_t dC_offset, magma_int_t ldc,
	magma_queue_t queue );

magma_err_t
magma_dsymv(
	magma_uplo_t uplo,
	magma_int_t n,
	double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
	                          magmaDouble_const_ptr dx, size_t dx_offset, magma_int_t incx,
	double beta,  magmaDouble_ptr       dy, size_t dy_offset, magma_int_t incy,
	magma_queue_t queue );

magma_err_t
magma_dsyrk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
    double beta,  magmaDouble_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t queue );

magma_err_t
magma_dtrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDouble_ptr       dB, size_t dB_offset, magma_int_t ldb,
    magma_queue_t queue );

magma_err_t
magma_dtrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDouble_ptr       dB, size_t dB_offset, magma_int_t ldb,
    magma_queue_t queue );

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_BLAS_H
