/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @precisions normal z -> s d c
*/

#ifndef MAGMABLAS_Z_H
#define MAGMABLAS_Z_H

#include "magma_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */

extern "C" magma_int_t
magma_zlarfbx_gpu(
    int m, int k, magmaDoubleComplex_ptr V, size_t v_offset, int ldv,
    magmaDoubleComplex_ptr dT, size_t dT_offset, int ldt, magmaDoubleComplex_ptr c, size_t c_offset, 
    magmaDoubleComplex_ptr dwork, size_t dwork_offset,
    magma_queue_t queue);

extern "C" magma_int_t
magma_zlarfgtx_gpu(
    int n, magmaDoubleComplex_ptr dx0, size_t dx0_offset, magmaDoubleComplex_ptr dx, size_t dx_offset, 
    magmaDoubleComplex_ptr dtau, size_t dtau_offset, magmaDouble_ptr dxnorm, size_t dxnorm_offset, 
    magmaDoubleComplex_ptr dA, size_t dA_offset, int it,
    magmaDoubleComplex_ptr V, size_t V_offset, int ldv, magmaDoubleComplex_ptr T, size_t T_offset, int ldt,
    magmaDoubleComplex_ptr dwork, size_t dwork_offset,
    magma_queue_t queue);

extern "C" magma_int_t
magmablas_dznrm2(
    int m, int num, magmaDoubleComplex_ptr da, size_t da_offset, magma_int_t ldda, magmaDouble_ptr dxnorm, size_t dxnorm_offset, 
    magma_queue_t queue);

extern "C" magma_int_t
magmablas_dznrm2_adjust(
    int k, magmaDouble_ptr xnorm, size_t xnorm_offset, magmaDoubleComplex_ptr c, size_t c_offset,
    magma_queue_t queue);
    
extern "C" magma_int_t
magmablas_zgemm_reduce(
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaDoubleComplex alpha, const magmaDoubleComplex_ptr d_A, size_t d_A_offset, magma_int_t lda,
    const magmaDoubleComplex_ptr d_B, size_t d_B_offset, magma_int_t ldb,
    magmaDoubleComplex beta,        magmaDoubleComplex_ptr d_C, size_t d_C_offset, magma_int_t ldc,
    magma_queue_t queue);

// iwocl 2013 benchmark
void
magmablas_empty(
    magmaDouble_ptr dA,
    magmaDouble_ptr dB,
    magmaDouble_ptr dC,
    magma_queue_t queue );

  /*
   * Transpose functions
   */

void
magmablas_ztranspose_inplace(
    magma_int_t n,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

void
magmablas_ztranspose(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_ptr dA,  size_t dA_offset,  magma_int_t ldda,
    magmaDoubleComplex_ptr dAT, size_t dAT_offset, magma_int_t lddat,
    magma_queue_t queue );


  /*
   * Multi-GPU copy functions
   */

void
magma_zgetmatrix_1D_col_bcyclic(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_ptr dA[], magma_int_t ldda,
    magmaDoubleComplex    *hA,   magma_int_t lda,
    magma_int_t ngpu, magma_int_t nb,
    magma_queue_t queues[] );

void
magma_zsetmatrix_1D_col_bcyclic(
    magma_int_t m, magma_int_t n,
    const magmaDoubleComplex *hA,   magma_int_t lda,
    magmaDoubleComplex_ptr    dA[], magma_int_t ldda,
    magma_int_t ngpu, magma_int_t nb,
    magma_queue_t queues[] );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magmablas_zlacpy(
    magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr       dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

void
magmablas_zlacpy_cnjg(
    magma_int_t n,
    magmaDoubleComplex_ptr dA1T, size_t dA1T_offset, magma_int_t lda1,
    magmaDoubleComplex_ptr dA2T, size_t dA2T_offset, magma_int_t lda2,
    magma_queue_t queue);

double
magmablas_zlange(
    magma_norm_t norm,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dwork, size_t dwork_offset, magma_queue_t queue );

double
magmablas_zlanhe(
    magma_norm_t norm, magma_uplo_t uplo,
    magma_int_t n,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dwork, size_t dwork_offset, magma_queue_t queue );

double
magmablas_zlansy(
    magma_norm_t norm, magma_uplo_t uplo,
    magma_int_t n,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dwork, size_t dwork_offset, magma_queue_t queue );

void
magmablas_zlascl_2x2(
    magma_type_t type, magma_int_t m,
    const magmaDoubleComplex_ptr dW, size_t dW_offset, magma_int_t lddw,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_int_t *info, magma_queue_t queue );

// TODO update to offdiag, diag interface
void
magmablas_zlaset(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magmaDoubleComplex offdiag, magmaDoubleComplex diag,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

void
magmablas_zlaswp(
    magma_int_t n,
    magmaDoubleComplex_ptr dAT, size_t dAT_offset, magma_int_t ldda,
    magma_int_t k1, magma_int_t k2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_queue_t queue );

void
magmablas_zlaswp2(
    magma_int_t n,
    magmaDoubleComplex_ptr dAT, size_t dAT_offset, magma_int_t ldda,
    magma_int_t k1, magma_int_t k2,
    magmaInt_const_ptr d_ipiv, size_t d_ipiv_offset, magma_int_t inci,
    magma_queue_t queue );

void
magmablas_zlaswpx(
    magma_int_t n,
    magmaDoubleComplex_ptr dA, size_t dAT_offset, magma_int_t ldx, magma_int_t ldy,
    magma_int_t k1, magma_int_t k2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_queue_t queue );


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magmablas_zswap(
    magma_int_t n,
    magmaDoubleComplex_ptr dx, size_t offset_dx, magma_int_t incx,
    magmaDoubleComplex_ptr dy, size_t offset_dy, magma_int_t incy,
    magma_queue_t queue );


  /*
   * Level 2 BLAS (alphabetical order)
   */


  /*
   * Level 3 BLAS (alphabetical order)
   */


  /*
   * Wrappers for platform independence.
   * These wrap CUBLAS or AMD OpenCL BLAS functions.
   */

// ========================================
// copying vectors
// set copies host to device
// get copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

void
magma_zsetvector(
    magma_int_t n,
    magmaDoubleComplex const*    hx_src,                   magma_int_t incx,
    magmaDoubleComplex_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_zsetvector_async(
    magma_int_t n,
    magmaDoubleComplex const*    hx_src,                   magma_int_t incx,
    magmaDoubleComplex_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );

void
magma_zgetvector(
    magma_int_t n,
    magmaDoubleComplex_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaDoubleComplex*          hy_dst,                   magma_int_t incy,
    magma_queue_t queue );

void
magma_zgetvector_async(
    magma_int_t n,
    magmaDoubleComplex_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaDoubleComplex*          hy_dst,                   magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );

void
magma_zcopyvector(
    magma_int_t n,
    magmaDoubleComplex_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaDoubleComplex_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_zcopyvector_async(
    magma_int_t n,
    magmaDoubleComplex_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaDoubleComplex_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

void
magma_zsetmatrix(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex const*    hA_src,                   magma_int_t ldha,
    magmaDoubleComplex_ptr       dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

void
magma_zsetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex const*    hA_src,                   magma_int_t ldha,
    magmaDoubleComplex_ptr       dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue, magma_event_t *event );

void
magma_zgetmatrix(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex*          hA_dst,                   magma_int_t ldha,
    magma_queue_t queue );

void
magma_zgetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex*          hA_dst,                   magma_int_t ldha,
    magma_queue_t queue, magma_event_t *event );

void
magma_zcopymatrix(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

void
magma_zcopymatrix_async(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue, magma_event_t *event );


// ========================================
// Level 1 BLAS (alphabetical order)

magma_int_t
magma_izamax(
    magma_int_t n,
    const magmaDoubleComplex_ptr dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

void
magma_zcopy(
    magma_int_t n,
    magmaDoubleComplex_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaDoubleComplex_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_zscal(
    magma_int_t n,
    magmaDoubleComplex alpha,
    magmaDoubleComplex_ptr dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

void
magma_zdscal(
    magma_int_t n,
    double alpha,
    magmaDoubleComplex_ptr dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

void
magma_zswap(
    magma_int_t n,
    magmaDoubleComplex_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaDoubleComplex_ptr dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_zgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex alpha,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_const_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaDoubleComplex beta,
    magmaDoubleComplex_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_zhemv(
    magma_uplo_t uplo,
    magma_int_t n,
    magmaDoubleComplex alpha,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_const_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaDoubleComplex beta,
    magmaDoubleComplex_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_ztrsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t n,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr       dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_zgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaDoubleComplex alpha,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_const_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaDoubleComplex beta,
    magmaDoubleComplex_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_zhemm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex alpha,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_const_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaDoubleComplex beta,
    magmaDoubleComplex_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_zher2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    magmaDoubleComplex alpha,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_const_ptr dB, size_t dB_offset, magma_int_t lddb,
    double beta,
    magmaDoubleComplex_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_zherk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    double beta,
    magmaDoubleComplex_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_ztrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex alpha,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr       dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

void
magma_ztrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex alpha,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr       dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );


#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif  /* MAGMABLAS_Z_H */
