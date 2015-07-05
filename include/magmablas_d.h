/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from magmablas_z.h normal z -> d, Sat Nov 15 00:21:34 2014
*/

#ifndef MAGMABLAS_D_H
#define MAGMABLAS_D_H

#include "magma_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */

extern "C" magma_int_t
magma_dlarfbx_gpu(
    int m, int k, magmaDouble_ptr V, size_t v_offset, int ldv,
    magmaDouble_ptr dT, size_t dT_offset, int ldt, magmaDouble_ptr c, size_t c_offset, 
    magmaDouble_ptr dwork, size_t dwork_offset,
    magma_queue_t queue);

extern "C" magma_int_t
magma_dlarfgtx_gpu(
    int n, magmaDouble_ptr dx0, size_t dx0_offset, magmaDouble_ptr dx, size_t dx_offset, 
    magmaDouble_ptr dtau, size_t dtau_offset, magmaDouble_ptr dxnorm, size_t dxnorm_offset, 
    magmaDouble_ptr dA, size_t dA_offset, int it,
    magmaDouble_ptr V, size_t V_offset, int ldv, magmaDouble_ptr T, size_t T_offset, int ldt,
    magmaDouble_ptr dwork, size_t dwork_offset,
    magma_queue_t queue);

extern "C" magma_int_t
magmablas_dnrm2(
    int m, int num, magmaDouble_ptr da, size_t da_offset, magma_int_t ldda, magmaDouble_ptr dxnorm, size_t dxnorm_offset, 
    magma_queue_t queue);

extern "C" magma_int_t
magmablas_dnrm2_adjust(
    int k, magmaDouble_ptr xnorm, size_t xnorm_offset, magmaDouble_ptr c, size_t c_offset,
    magma_queue_t queue);
    
extern "C" magma_int_t
magmablas_dgemm_reduce(
    magma_int_t m, magma_int_t n, magma_int_t k,
    double alpha, const magmaDouble_ptr d_A, size_t d_A_offset, magma_int_t lda,
    const magmaDouble_ptr d_B, size_t d_B_offset, magma_int_t ldb,
    double beta,        magmaDouble_ptr d_C, size_t d_C_offset, magma_int_t ldc,
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
magmablas_dtranspose_inplace(
    magma_int_t n,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

void
magmablas_dtranspose(
    magma_int_t m, magma_int_t n,
    magmaDouble_ptr dA,  size_t dA_offset,  magma_int_t ldda,
    magmaDouble_ptr dAT, size_t dAT_offset, magma_int_t lddat,
    magma_queue_t queue );


  /*
   * Multi-GPU copy functions
   */

void
magma_dgetmatrix_1D_col_bcyclic(
    magma_int_t m, magma_int_t n,
    magmaDouble_ptr dA[], magma_int_t ldda,
    double    *hA,   magma_int_t lda,
    magma_int_t ngpu, magma_int_t nb,
    magma_queue_t queues[] );

void
magma_dsetmatrix_1D_col_bcyclic(
    magma_int_t m, magma_int_t n,
    const double *hA,   magma_int_t lda,
    magmaDouble_ptr    dA[], magma_int_t ldda,
    magma_int_t ngpu, magma_int_t nb,
    magma_queue_t queues[] );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magmablas_dlacpy(
    magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr       dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

void
magmablas_dlacpy_cnjg(
    magma_int_t n,
    magmaDouble_ptr dA1T, size_t dA1T_offset, magma_int_t lda1,
    magmaDouble_ptr dA2T, size_t dA2T_offset, magma_int_t lda2,
    magma_queue_t queue);

double
magmablas_dlange(
    magma_norm_t norm,
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dwork, size_t dwork_offset, magma_queue_t queue );

double
magmablas_dlansy(
    magma_norm_t norm, magma_uplo_t uplo,
    magma_int_t n,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dwork, size_t dwork_offset, magma_queue_t queue );

double
magmablas_dlansy(
    magma_norm_t norm, magma_uplo_t uplo,
    magma_int_t n,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dwork, size_t dwork_offset, magma_queue_t queue );

void
magmablas_dlascl_2x2(
    magma_type_t type, magma_int_t m,
    const magmaDouble_ptr dW, size_t dW_offset, magma_int_t lddw,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_int_t *info, magma_queue_t queue );

// TODO update to offdiag, diag interface
void
magmablas_dlaset(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    double offdiag, double diag,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

void
magmablas_dlaswp(
    magma_int_t n,
    magmaDouble_ptr dAT, size_t dAT_offset, magma_int_t ldda,
    magma_int_t k1, magma_int_t k2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_queue_t queue );

void
magmablas_dlaswp2(
    magma_int_t n,
    magmaDouble_ptr dAT, size_t dAT_offset, magma_int_t ldda,
    magma_int_t k1, magma_int_t k2,
    magmaInt_const_ptr d_ipiv, size_t d_ipiv_offset, magma_int_t inci,
    magma_queue_t queue );

void
magmablas_dlaswpx(
    magma_int_t n,
    magmaDouble_ptr dA, size_t dAT_offset, magma_int_t ldx, magma_int_t ldy,
    magma_int_t k1, magma_int_t k2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_queue_t queue );


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magmablas_dswap(
    magma_int_t n,
    magmaDouble_ptr dx, size_t offset_dx, magma_int_t incx,
    magmaDouble_ptr dy, size_t offset_dy, magma_int_t incy,
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
magma_dsetvector(
    magma_int_t n,
    double const*    hx_src,                   magma_int_t incx,
    magmaDouble_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_dsetvector_async(
    magma_int_t n,
    double const*    hx_src,                   magma_int_t incx,
    magmaDouble_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );

void
magma_dgetvector(
    magma_int_t n,
    magmaDouble_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    double*          hy_dst,                   magma_int_t incy,
    magma_queue_t queue );

void
magma_dgetvector_async(
    magma_int_t n,
    magmaDouble_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    double*          hy_dst,                   magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );

void
magma_dcopyvector(
    magma_int_t n,
    magmaDouble_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaDouble_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_dcopyvector_async(
    magma_int_t n,
    magmaDouble_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaDouble_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

void
magma_dsetmatrix(
    magma_int_t m, magma_int_t n,
    double const*    hA_src,                   magma_int_t ldha,
    magmaDouble_ptr       dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

void
magma_dsetmatrix_async(
    magma_int_t m, magma_int_t n,
    double const*    hA_src,                   magma_int_t ldha,
    magmaDouble_ptr       dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue, magma_event_t *event );

void
magma_dgetmatrix(
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    double*          hA_dst,                   magma_int_t ldha,
    magma_queue_t queue );

void
magma_dgetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    double*          hA_dst,                   magma_int_t ldha,
    magma_queue_t queue, magma_event_t *event );

void
magma_dcopymatrix(
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

void
magma_dcopymatrix_async(
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue, magma_event_t *event );


// ========================================
// Level 1 BLAS (alphabetical order)

magma_int_t
magma_idamax(
    magma_int_t n,
    const magmaDouble_ptr dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

void
magma_dcopy(
    magma_int_t n,
    magmaDouble_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaDouble_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_dscal(
    magma_int_t n,
    double alpha,
    magmaDouble_ptr dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

void
magma_dscal(
    magma_int_t n,
    double alpha,
    magmaDouble_ptr dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

void
magma_dswap(
    magma_int_t n,
    magmaDouble_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaDouble_ptr dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_dgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    double alpha,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_const_ptr dx, size_t dx_offset, magma_int_t incx,
    double beta,
    magmaDouble_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_dsymv(
    magma_uplo_t uplo,
    magma_int_t n,
    double alpha,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_const_ptr dx, size_t dx_offset, magma_int_t incx,
    double beta,
    magmaDouble_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_dtrsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t n,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr       dx, size_t dx_offset, magma_int_t incx,
    magma_queue_t queue );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_dgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    double alpha,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_const_ptr dB, size_t dB_offset, magma_int_t lddb,
    double beta,
    magmaDouble_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_dsymm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    double alpha,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_const_ptr dB, size_t dB_offset, magma_int_t lddb,
    double beta,
    magmaDouble_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_dsyr2k(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_const_ptr dB, size_t dB_offset, magma_int_t lddb,
    double beta,
    magmaDouble_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_dsyrk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    double beta,
    magmaDouble_ptr       dC, size_t dC_offset, magma_int_t lddc,
    magma_queue_t queue );

void
magma_dtrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    double alpha,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr       dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

void
magma_dtrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    double alpha,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr       dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );


#ifdef __cplusplus
}
#endif

#undef REAL

#endif  /* MAGMABLAS_D_H */
