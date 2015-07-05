/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from magmablas_z.h normal z -> s, Sat Nov 15 00:21:34 2014
*/

#ifndef MAGMABLAS_S_H
#define MAGMABLAS_S_H

#include "magma_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */

extern "C" magma_int_t
magma_slarfbx_gpu(
    int m, int k, magmaFloat_ptr V, size_t v_offset, int ldv,
    magmaFloat_ptr dT, size_t dT_offset, int ldt, magmaFloat_ptr c, size_t c_offset, 
    magmaFloat_ptr dwork, size_t dwork_offset,
    magma_queue_t queue);

extern "C" magma_int_t
magma_slarfgtx_gpu(
    int n, magmaFloat_ptr dx0, size_t dx0_offset, magmaFloat_ptr dx, size_t dx_offset, 
    magmaFloat_ptr dtau, size_t dtau_offset, magmaFloat_ptr dxnorm, size_t dxnorm_offset, 
    magmaFloat_ptr dA, size_t dA_offset, int it,
    magmaFloat_ptr V, size_t V_offset, int ldv, magmaFloat_ptr T, size_t T_offset, int ldt,
    magmaFloat_ptr dwork, size_t dwork_offset,
    magma_queue_t queue);

extern "C" magma_int_t
magmablas_snrm2(
    int m, int num, magmaFloat_ptr da, size_t da_offset, magma_int_t ldda, magmaFloat_ptr dxnorm, size_t dxnorm_offset, 
    magma_queue_t queue);

extern "C" magma_int_t
magmablas_snrm2_adjust(
    int k, magmaFloat_ptr xnorm, size_t xnorm_offset, magmaFloat_ptr c, size_t c_offset,
    magma_queue_t queue);
    
extern "C" magma_int_t
magmablas_sgemm_reduce(
    magma_int_t m, magma_int_t n, magma_int_t k,
    float alpha, const magmaFloat_ptr d_A, size_t d_A_offset, magma_int_t lda,
    const magmaFloat_ptr d_B, size_t d_B_offset, magma_int_t ldb,
    float beta,        magmaFloat_ptr d_C, size_t d_C_offset, magma_int_t ldc,
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
magmablas_stranspose_inplace(
    magma_int_t n,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

void
magmablas_stranspose(
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr dA,  size_t dA_offset,  magma_int_t ldda,
    magmaFloat_ptr dAT, size_t dAT_offset, magma_int_t lddat,
    magma_queue_t queue );


  /*
   * Multi-GPU copy functions
   */

void
magma_sgetmatrix_1D_col_bcyclic(
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr dA[], magma_int_t ldda,
    float    *hA,   magma_int_t lda,
    magma_int_t ngpu, magma_int_t nb,
    magma_queue_t queues[] );

void
magma_ssetmatrix_1D_col_bcyclic(
    magma_int_t m, magma_int_t n,
    const float *hA,   magma_int_t lda,
    magmaFloat_ptr    dA[], magma_int_t ldda,
    magma_int_t ngpu, magma_int_t nb,
    magma_queue_t queues[] );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magmablas_slacpy(
    magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr       dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

void
magmablas_slacpy_cnjg(
    magma_int_t n,
    magmaFloat_ptr dA1T, size_t dA1T_offset, magma_int_t lda1,
    magmaFloat_ptr dA2T, size_t dA2T_offset, magma_int_t lda2,
    magma_queue_t queue);

float
magmablas_slange(
    magma_norm_t norm,
    magma_int_t m, magma_int_t n,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr dwork, size_t dwork_offset, magma_queue_t queue );

float
magmablas_slansy(
    magma_norm_t norm, magma_uplo_t uplo,
    magma_int_t n,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr dwork, size_t dwork_offset, magma_queue_t queue );

float
magmablas_slansy(
    magma_norm_t norm, magma_uplo_t uplo,
    magma_int_t n,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr dwork, size_t dwork_offset, magma_queue_t queue );

void
magmablas_slascl_2x2(
    magma_type_t type, magma_int_t m,
    const magmaFloat_ptr dW, size_t dW_offset, magma_int_t lddw,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_int_t *info, magma_queue_t queue );

// TODO update to offdiag, diag interface
void
magmablas_slaset(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    float offdiag, float diag,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

void
magmablas_slaswp(
    magma_int_t n,
    magmaFloat_ptr dAT, size_t dAT_offset, magma_int_t ldda,
    magma_int_t k1, magma_int_t k2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_queue_t queue );

void
magmablas_slaswp2(
    magma_int_t n,
    magmaFloat_ptr dAT, size_t dAT_offset, magma_int_t ldda,
    magma_int_t k1, magma_int_t k2,
    magmaInt_const_ptr d_ipiv, size_t d_ipiv_offset, magma_int_t inci,
    magma_queue_t queue );

void
magmablas_slaswpx(
    magma_int_t n,
    magmaFloat_ptr dA, size_t dAT_offset, magma_int_t ldx, magma_int_t ldy,
    magma_int_t k1, magma_int_t k2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_queue_t queue );


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magmablas_sswap(
    magma_int_t n,
    magmaFloat_ptr dx, size_t offset_dx, magma_int_t incx,
    magmaFloat_ptr dy, size_t offset_dy, magma_int_t incy,
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
magma_ssetvector(
    magma_int_t n,
    float const*    hx_src,                   magma_int_t incx,
    magmaFloat_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_ssetvector_async(
    magma_int_t n,
    float const*    hx_src,                   magma_int_t incx,
    magmaFloat_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );

void
magma_sgetvector(
    magma_int_t n,
    magmaFloat_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    float*          hy_dst,                   magma_int_t incy,
    magma_queue_t queue );

void
magma_sgetvector_async(
    magma_int_t n,
    magmaFloat_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    float*          hy_dst,                   magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );

void
magma_scopyvector(
    magma_int_t n,
    magmaFloat_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaFloat_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_scopyvector_async(
    magma_int_t n,
    magmaFloat_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaFloat_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

void
magma_ssetmatrix(
    magma_int_t m, magma_int_t n,
    float const*    hA_src,                   magma_int_t ldha,
    magmaFloat_ptr       dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

void
magma_ssetmatrix_async(
    magma_int_t m, magma_int_t n,
    float const*    hA_src,                   magma_int_t ldha,
    magmaFloat_ptr       dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue, magma_event_t *event );

void
magma_sgetmatrix(
    magma_int_t m, magma_int_t n,
    magmaFloat_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    float*          hA_dst,                   magma_int_t ldha,
    magma_queue_t queue );

void
magma_sgetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaFloat_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    float*          hA_dst,                   magma_int_t ldha,
    magma_queue_t queue, magma_event_t *event );

void
magma_scopymatrix(
    magma_int_t m, magma_int_t n,
    magmaFloat_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

void
magma_scopymatrix_async(
    magma_int_t m, magma_int_t n,
    magmaFloat_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue, magma_event_t *event );


// ========================================
// Level 1 BLAS (alphabetical order)

magma_int_t
magma_isamax(
    magma_int_t n,
    const magmaFloat_ptr dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

void
magma_scopy(
    magma_int_t n,
    magmaFloat_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaFloat_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_sscal(
    magma_int_t n,
    float alpha,
    magmaFloat_ptr dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

void
magma_sscal(
    magma_int_t n,
    float alpha,
    magmaFloat_ptr dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

void
magma_sswap(
    magma_int_t n,
    magmaFloat_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaFloat_ptr dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_sgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    float alpha,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_const_ptr dx, size_t dx_offset, magma_int_t incx,
    float beta,
    magmaFloat_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_ssymv(
    magma_uplo_t uplo,
    magma_int_t n,
    float alpha,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_const_ptr dx, size_t dx_offset, magma_int_t incx,
    float beta,
    magmaFloat_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_strsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t n,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr       dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_sgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    float alpha,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_const_ptr dB, size_t dB_offset, magma_int_t lddb,
    float beta,
    magmaFloat_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_ssymm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    float alpha,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_const_ptr dB, size_t dB_offset, magma_int_t lddb,
    float beta,
    magmaFloat_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_ssyr2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_const_ptr dB, size_t dB_offset, magma_int_t lddb,
    float beta,
    magmaFloat_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_ssyrk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    float beta,
    magmaFloat_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_strmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    float alpha,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr       dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

void
magma_strsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    float alpha,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr       dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );


#ifdef __cplusplus
}
#endif

#undef REAL

#endif  /* MAGMABLAS_S_H */
