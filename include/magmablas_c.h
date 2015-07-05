/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from magmablas_z.h normal z -> c, Sat Nov 15 00:21:34 2014
*/

#ifndef MAGMABLAS_C_H
#define MAGMABLAS_C_H

#include "magma_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */

extern "C" magma_int_t
magma_clarfbx_gpu(
    int m, int k, magmaFloatComplex_ptr V, size_t v_offset, int ldv,
    magmaFloatComplex_ptr dT, size_t dT_offset, int ldt, magmaFloatComplex_ptr c, size_t c_offset, 
    magmaFloatComplex_ptr dwork, size_t dwork_offset,
    magma_queue_t queue);

extern "C" magma_int_t
magma_clarfgtx_gpu(
    int n, magmaFloatComplex_ptr dx0, size_t dx0_offset, magmaFloatComplex_ptr dx, size_t dx_offset, 
    magmaFloatComplex_ptr dtau, size_t dtau_offset, magmaFloat_ptr dxnorm, size_t dxnorm_offset, 
    magmaFloatComplex_ptr dA, size_t dA_offset, int it,
    magmaFloatComplex_ptr V, size_t V_offset, int ldv, magmaFloatComplex_ptr T, size_t T_offset, int ldt,
    magmaFloatComplex_ptr dwork, size_t dwork_offset,
    magma_queue_t queue);

extern "C" magma_int_t
magmablas_scnrm2(
    int m, int num, magmaFloatComplex_ptr da, size_t da_offset, magma_int_t ldda, magmaFloat_ptr dxnorm, size_t dxnorm_offset, 
    magma_queue_t queue);

extern "C" magma_int_t
magmablas_scnrm2_adjust(
    int k, magmaFloat_ptr xnorm, size_t xnorm_offset, magmaFloatComplex_ptr c, size_t c_offset,
    magma_queue_t queue);
    
extern "C" magma_int_t
magmablas_cgemm_reduce(
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloatComplex alpha, const magmaFloatComplex_ptr d_A, size_t d_A_offset, magma_int_t lda,
    const magmaFloatComplex_ptr d_B, size_t d_B_offset, magma_int_t ldb,
    magmaFloatComplex beta,        magmaFloatComplex_ptr d_C, size_t d_C_offset, magma_int_t ldc,
    magma_queue_t queue);

// iwocl 2013 benchmark
void
magmablas_empty(
    magmaFloat_ptr dA,
    magmaFloat_ptr dB,
    magmaFloat_ptr dC,
    magma_queue_t queue );

  /*
   * Transpose functions
   */

void
magmablas_ctranspose_inplace(
    magma_int_t n,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

void
magmablas_ctranspose(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr dA,  size_t dA_offset,  magma_int_t ldda,
    magmaFloatComplex_ptr dAT, size_t dAT_offset, magma_int_t lddat,
    magma_queue_t queue );


  /*
   * Multi-GPU copy functions
   */

void
magma_cgetmatrix_1D_col_bcyclic(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr dA[], magma_int_t ldda,
    magmaFloatComplex    *hA,   magma_int_t lda,
    magma_int_t ngpu, magma_int_t nb,
    magma_queue_t queues[] );

void
magma_csetmatrix_1D_col_bcyclic(
    magma_int_t m, magma_int_t n,
    const magmaFloatComplex *hA,   magma_int_t lda,
    magmaFloatComplex_ptr    dA[], magma_int_t ldda,
    magma_int_t ngpu, magma_int_t nb,
    magma_queue_t queues[] );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magmablas_clacpy(
    magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr       dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

void
magmablas_clacpy_cnjg(
    magma_int_t n,
    magmaFloatComplex_ptr dA1T, size_t dA1T_offset, magma_int_t lda1,
    magmaFloatComplex_ptr dA2T, size_t dA2T_offset, magma_int_t lda2,
    magma_queue_t queue);

float
magmablas_clange(
    magma_norm_t norm,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr dwork, size_t dwork_offset, magma_queue_t queue );

float
magmablas_clanhe(
    magma_norm_t norm, magma_uplo_t uplo,
    magma_int_t n,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr dwork, size_t dwork_offset, magma_queue_t queue );

float
magmablas_clansy(
    magma_norm_t norm, magma_uplo_t uplo,
    magma_int_t n,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr dwork, size_t dwork_offset, magma_queue_t queue );

void
magmablas_clascl_2x2(
    magma_type_t type, magma_int_t m,
    const magmaFloatComplex_ptr dW, size_t dW_offset, magma_int_t lddw,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_int_t *info, magma_queue_t queue );

// TODO update to offdiag, diag interface
void
magmablas_claset(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magmaFloatComplex offdiag, magmaFloatComplex diag,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

void
magmablas_claswp(
    magma_int_t n,
    magmaFloatComplex_ptr dAT, size_t dAT_offset, magma_int_t ldda,
    magma_int_t k1, magma_int_t k2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_queue_t queue );

void
magmablas_claswp2(
    magma_int_t n,
    magmaFloatComplex_ptr dAT, size_t dAT_offset, magma_int_t ldda,
    magma_int_t k1, magma_int_t k2,
    magmaInt_const_ptr d_ipiv, size_t d_ipiv_offset, magma_int_t inci,
    magma_queue_t queue );

void
magmablas_claswpx(
    magma_int_t n,
    magmaFloatComplex_ptr dA, size_t dAT_offset, magma_int_t ldx, magma_int_t ldy,
    magma_int_t k1, magma_int_t k2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_queue_t queue );


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magmablas_cswap(
    magma_int_t n,
    magmaFloatComplex_ptr dx, size_t offset_dx, magma_int_t incx,
    magmaFloatComplex_ptr dy, size_t offset_dy, magma_int_t incy,
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
magma_csetvector(
    magma_int_t n,
    magmaFloatComplex const*    hx_src,                   magma_int_t incx,
    magmaFloatComplex_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_csetvector_async(
    magma_int_t n,
    magmaFloatComplex const*    hx_src,                   magma_int_t incx,
    magmaFloatComplex_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );

void
magma_cgetvector(
    magma_int_t n,
    magmaFloatComplex_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaFloatComplex*          hy_dst,                   magma_int_t incy,
    magma_queue_t queue );

void
magma_cgetvector_async(
    magma_int_t n,
    magmaFloatComplex_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaFloatComplex*          hy_dst,                   magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );

void
magma_ccopyvector(
    magma_int_t n,
    magmaFloatComplex_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaFloatComplex_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_ccopyvector_async(
    magma_int_t n,
    magmaFloatComplex_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaFloatComplex_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

void
magma_csetmatrix(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex const*    hA_src,                   magma_int_t ldha,
    magmaFloatComplex_ptr       dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

void
magma_csetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex const*    hA_src,                   magma_int_t ldha,
    magmaFloatComplex_ptr       dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue, magma_event_t *event );

void
magma_cgetmatrix(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex*          hA_dst,                   magma_int_t ldha,
    magma_queue_t queue );

void
magma_cgetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex*          hA_dst,                   magma_int_t ldha,
    magma_queue_t queue, magma_event_t *event );

void
magma_ccopymatrix(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

void
magma_ccopymatrix_async(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue, magma_event_t *event );


// ========================================
// Level 1 BLAS (alphabetical order)

magma_int_t
magma_icamax(
    magma_int_t n,
    const magmaFloatComplex_ptr dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

void
magma_ccopy(
    magma_int_t n,
    magmaFloatComplex_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaFloatComplex_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_cscal(
    magma_int_t n,
    magmaFloatComplex alpha,
    magmaFloatComplex_ptr dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

void
magma_csscal(
    magma_int_t n,
    float alpha,
    magmaFloatComplex_ptr dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

void
magma_cswap(
    magma_int_t n,
    magmaFloatComplex_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaFloatComplex_ptr dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_cgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex alpha,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_const_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaFloatComplex beta,
    magmaFloatComplex_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_chemv(
    magma_uplo_t uplo,
    magma_int_t n,
    magmaFloatComplex alpha,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_const_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaFloatComplex beta,
    magmaFloatComplex_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_ctrsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t n,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr       dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_cgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloatComplex alpha,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_const_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaFloatComplex beta,
    magmaFloatComplex_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_chemm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex alpha,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_const_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaFloatComplex beta,
    magmaFloatComplex_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_cher2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    magmaFloatComplex alpha,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_const_ptr dB, size_t dB_offset, magma_int_t lddb,
    float beta,
    magmaFloatComplex_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_cherk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    float beta,
    magmaFloatComplex_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_ctrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex alpha,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr       dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

void
magma_ctrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex alpha,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr       dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );


#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif  /* MAGMABLAS_C_H */
