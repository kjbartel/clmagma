/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from magma_z.h normal z -> s, Sat Nov 15 00:21:34 2014
*/

#ifndef MAGMA_S_H
#define MAGMA_S_H

#include "magma_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA Auxiliary functions to get the NB used
*/
magma_int_t magma_get_slaex3_m_nb();       // defined in slaex3_m.cpp

magma_int_t magma_get_spotrf_nb( magma_int_t m );
magma_int_t magma_get_sgetrf_nb( magma_int_t m );
magma_int_t magma_get_sgetri_nb( magma_int_t m );
magma_int_t magma_get_sgeqp3_nb( magma_int_t m );
magma_int_t magma_get_sgeqrf_nb( magma_int_t m );
magma_int_t magma_get_sgeqlf_nb( magma_int_t m );
magma_int_t magma_get_sgehrd_nb( magma_int_t m );
magma_int_t magma_get_ssytrd_nb( magma_int_t m );
magma_int_t magma_get_ssytrf_nb( magma_int_t m );
magma_int_t magma_get_sgelqf_nb( magma_int_t m );
magma_int_t magma_get_sgebrd_nb( magma_int_t m );
magma_int_t magma_get_ssygst_nb( magma_int_t m );
magma_int_t magma_get_sgesvd_nb( magma_int_t m );
magma_int_t magma_get_ssygst_nb_m( magma_int_t m );
magma_int_t magma_get_sbulge_nb( magma_int_t m, magma_int_t nbthreads );
magma_int_t magma_get_sbulge_nb_mgpu( magma_int_t m );
magma_int_t magma_sbulge_get_Vblksiz( magma_int_t m, magma_int_t nb, magma_int_t nbthreads );
magma_int_t magma_get_sbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU (alphabetical order)
*/

/* magma_sgebrd */
/* magma_sgebrd */
/* magma_sgebrd */
magma_int_t
magma_sgebrd(
    magma_int_t m, magma_int_t n,
    float *A, magma_int_t lda,
    float *d, float *e,
    float *tauq, float *taup,
    float *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sgeev */
/* magma_sgeev */
/* magma_sgeev */
magma_int_t
magma_sgeev(
    magma_vec_t jobvl, magma_vec_t jobvr, magma_int_t n,
    float *A, magma_int_t lda,
    #ifdef COMPLEX
    float *w,
    #else
    float *wr, float *wi,
    #endif
    float *VL, magma_int_t ldvl,
    float *VR, magma_int_t ldvr,
    float *work, magma_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sgehrd */
/* magma_sgehrd */
/* magma_sgehrd */
magma_int_t
magma_sgehrd(
    magma_int_t n, magma_int_t ilo, magma_int_t ihi,
    float *A, magma_int_t lda,
    float *tau,
    float *work, magma_int_t lwork,
    magmaFloat_ptr dT, size_t dT_offset,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sgeqrf */
/* magma_sgeqrf */
/* magma_sgeqrf */
magma_int_t
magma_sgeqrf(
    magma_int_t m, magma_int_t n,
    float *A, magma_int_t lda,
    float *tau,
    float *work, magma_int_t lwork,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_sgesv */
/* magma_sgesv */
/* magma_sgesv */
magma_int_t
magma_sgesv(
    magma_int_t n, magma_int_t nrhs,
    float *A, magma_int_t lda, magma_int_t *ipiv,
    float *B, magma_int_t ldb,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_sgesvd */
/* magma_sgesvd */
/* magma_sgesvd */
magma_int_t
magma_sgesvd(
    magma_vec_t jobu, magma_vec_t jobvt, magma_int_t m, magma_int_t n,
    float *A,    magma_int_t lda, float *s,
    float *U,    magma_int_t ldu,
    float *VT,   magma_int_t ldvt,
    float *work, magma_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sgetrf */
/* magma_sgetrf */
/* magma_sgetrf */
magma_int_t
magma_sgetrf(
    magma_int_t m, magma_int_t n,
    float *A, magma_int_t lda, magma_int_t *ipiv,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_ssyevd */
/* magma_ssyevd */
/* magma_ssyevd */
magma_int_t
magma_ssyevd(
    magma_vec_t jobz, magma_uplo_t uplo, magma_int_t n,
    float *A, magma_int_t lda,
    float *w,
    float *work, magma_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_int_t lrwork,
    #endif
    magma_int_t *iwork, magma_int_t liwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_ssytrd */
/* magma_ssytrd */
/* magma_ssytrd */
magma_int_t
magma_ssytrd(
    magma_uplo_t uplo, magma_int_t n,
    float *A, magma_int_t lda,
    float *d, float *e, float *tau,
    float *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_ssytrf */
/* magma_ssytrf */
/* magma_ssytrf */
magma_int_t
magma_ssytrf(
    magma_uplo_t uplo, magma_int_t n,
    float *A, magma_int_t lda,
    magma_int_t *ipiv,
    magma_queue_t queue,
    magma_int_t *info);

#ifdef REAL
// only applicable to real [sd] precisions
/* magma_slaex0 */
/* magma_slaex0 */
/* magma_slaex0 */
magma_int_t
magma_slaex0(
    magma_int_t n, float *d, float *e,
    float *Q, magma_int_t ldq,
    float *work, magma_int_t *iwork,
    magmaFloat_ptr dwork,
    magma_range_t range, float vl, float vu, magma_int_t il, magma_int_t iu,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_slaex1 */
/* magma_slaex1 */
/* magma_slaex1 */
magma_int_t
magma_slaex1(
    magma_int_t n, float *d,
    float *Q, magma_int_t ldq,
    magma_int_t *indxq, float rho, magma_int_t cutpnt,
    float *work, magma_int_t *iwork,
    magmaFloat_ptr dwork,
    magma_range_t range, float vl, float vu, magma_int_t il, magma_int_t iu,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_slaex3 */
/* magma_slaex3 */
/* magma_slaex3 */
magma_int_t
magma_slaex3(
    magma_int_t k, magma_int_t n, magma_int_t n1, float *d,
    float *Q, magma_int_t ldq,
    float rho,
    float *dlamda, float *Q2, magma_int_t *indx,
    magma_int_t *ctot, float *w, float *s, magma_int_t *indxq,
    magmaFloat_ptr dwork,
    magma_range_t range, float vl, float vu, magma_int_t il, magma_int_t iu,
    magma_queue_t queue,
    magma_int_t *info);
#endif  // REAL

/* magma_slasyf_gpu */
/* magma_slasyf_gpu */
/* magma_slasyf_gpu */
magma_int_t
magma_slasyf_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nb, magma_int_t *kb,
    float *hA, magma_int_t lda,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_int_t *ipiv,
    magmaFloat_ptr dW, size_t dW_offset, magma_int_t lddw,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_slahr2 */
/* magma_slahr2 */
/* magma_slahr2 */
magma_int_t
magma_slahr2(
    magma_int_t n, magma_int_t k, magma_int_t nb,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr dV, size_t dV_offset, magma_int_t lddv,
    float *A,  magma_int_t lda,
    float *tau,
    float *T,  magma_int_t ldt,
    float *Y,  magma_int_t ldy,
    magma_queue_t queue);

/* magma_slahru */
/* magma_slahru */
/* magma_slahru */
magma_int_t
magma_slahru(
    magma_int_t n, magma_int_t ihi, magma_int_t k, magma_int_t nb,
    float     *A, magma_int_t lda,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr dY, size_t dY_offset, magma_int_t lddy,
    magmaFloat_ptr dV, size_t dV_offset, magma_int_t lddv,
    magmaFloat_ptr dT, size_t dT_offset,
    magmaFloat_ptr dwork, size_t dwork_offset,
    magma_queue_t queue);

/* magma_slatrd */
/* magma_slatrd */
/* magma_slatrd */
magma_int_t
magma_slatrd(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nb,
    float *A,  magma_int_t lda,
    float *e, float *tau,
    float *W,  magma_int_t ldw,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr dW, size_t dW_offset, magma_int_t lddw,
    magma_queue_t queue);

/* magma_sposv */
/* magma_sposv */
/* magma_sposv */
magma_int_t
magma_sposv(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    float *A, magma_int_t lda,
    float *B, magma_int_t ldb,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_spotrf */
/* magma_spotrf */
/* magma_spotrf */
magma_int_t
magma_spotrf(
    magma_uplo_t uplo, magma_int_t n,
    float *A, magma_int_t lda,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_sstedx */
/* magma_sstedx */
/* magma_sstedx */
magma_int_t
magma_sstedx(
    magma_range_t range, magma_int_t n, float vl, float vu,
    magma_int_t il, magma_int_t iu, float *d, float *e,
    float *Z, magma_int_t ldz,
    float *rwork, magma_int_t lrwork,
    magma_int_t *iwork, magma_int_t liwork,
    magmaFloat_ptr dwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sorghr */
/* magma_sorghr */
/* magma_sorghr */
magma_int_t
magma_sorghr(
    magma_int_t n, magma_int_t ilo, magma_int_t ihi,
    float *A, magma_int_t lda,
    float *tau,
    magmaFloat_ptr dT, size_t dT_offset, magma_int_t nb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sorgqr */
/* magma_sorgqr */
/* magma_sorgqr */
magma_int_t
magma_sorgqr(
    magma_int_t m, magma_int_t n, magma_int_t k,
    float *A, magma_int_t lda,
    float *tau,
    magmaFloat_ptr dT, size_t dT_offset, magma_int_t nb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sormql */
/* magma_sormql */
/* magma_sormql */
magma_int_t
magma_sormql(
    magma_side_t side, magma_trans_t trans,
    magma_int_t m, magma_int_t n, magma_int_t k,
    float *A, magma_int_t lda,
    float *tau,
    float *C, magma_int_t ldc,
    float *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sormqr */
/* magma_sormqr */
/* magma_sormqr */
magma_int_t
magma_sormqr(
    magma_side_t side, magma_trans_t trans,
    magma_int_t m, magma_int_t n, magma_int_t k,
    float *A, magma_int_t lda,
    float *tau,
    float *C, magma_int_t ldc,
    float *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sormtr */
/* magma_sormtr */
/* magma_sormtr */
magma_int_t
magma_sormtr(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t m, magma_int_t n,
    float *A,    magma_int_t lda,
    float *tau,
    float *C,    magma_int_t ldc,
    float *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on GPU (alphabetical order)
*/

/* magma_sgels_gpu */
/* magma_sgels_gpu */
/* magma_sgels_gpu */
magma_int_t
magma_sgels_gpu(
    magma_trans_t trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr dB, size_t dB_offset, magma_int_t lddb,
    float *hwork, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sgeqr2x3_gpu */
/* magma_sgeqr2x3_gpu */
/* magma_sgeqr2x3_gpu */
magma_int_t
magma_sgeqr2x3_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr dA,    size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr dtau,  size_t dtau_offset,
    magmaFloat_ptr dT,    size_t dT_offset,
    magmaFloat_ptr ddA,   size_t ddA_offset,
    magmaFloat_ptr        dwork, size_t dwork_offset,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sgeqrf_gpu */
/* magma_sgeqrf_gpu */
/* magma_sgeqrf_gpu */
magma_int_t
magma_sgeqrf_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr dA, size_t dA_offset,  magma_int_t ldda,
    float *tau,
    magmaFloat_ptr dT, size_t dT_offset,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sgeqrf_msub */
/* magma_sgeqrf_msub */
/* magma_sgeqrf_msub */
magma_int_t
magma_sgeqrf_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr dlA[], magma_int_t ldda,
    float *tau,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_sgeqrf2_2q_gpu */
/* magma_sgeqrf2_2q_gpu */
/* magma_sgeqrf2_2q_gpu */
magma_int_t
magma_sgeqrf2_2q_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    float *tau,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_sgeqrf2_gpu */
/* magma_sgeqrf2_gpu */
/* magma_sgeqrf2_gpu */
magma_int_t
magma_sgeqrf2_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    float *tau,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_sgeqrf2_mgpu */
/* magma_sgeqrf2_mgpu */
/* magma_sgeqrf2_mgpu */
magma_int_t
magma_sgeqrf2_mgpu(
    magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr dlA[], magma_int_t ldda,
    float *tau,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_sgeqrs_gpu */
/* magma_sgeqrs_gpu */
/* magma_sgeqrs_gpu */
magma_int_t
magma_sgeqrs_gpu(
    magma_int_t m, magma_int_t n, magma_int_t nrhs,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    float *tau,
    magmaFloat_ptr dT, size_t dT_offset,
    magmaFloat_ptr dB, size_t dB_offset, magma_int_t lddb,
    float *hwork, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sgesv_gpu */
/* magma_sgesv_gpu */
/* magma_sgesv_gpu */
magma_int_t
magma_sgesv_gpu(
    magma_int_t n, magma_int_t nrhs,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magmaFloat_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sgetrf_gpu */
/* magma_sgetrf_gpu */
/* magma_sgetrf_gpu */
magma_int_t
magma_sgetrf_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sgetrf_mgpu */
/* magma_sgetrf_mgpu */
/* magma_sgetrf_mgpu */
magma_int_t
magma_sgetrf_mgpu(
    magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr d_lA[], size_t dlA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_sgetrf_msub */
/* magma_sgetrf_msub */
/* magma_sgetrf_msub */
magma_int_t
magma_sgetrf_msub(
    magma_trans_t trans, magma_int_t num_subs, magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr d_lA[], size_t dlA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_sgetrf2_gpu */
/* magma_sgetrf2_gpu */
/* magma_sgetrf2_gpu */
magma_int_t
magma_sgetrf2_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_sgetrf2_mgpu */
/* magma_sgetrf2_mgpu */
/* magma_sgetrf2_mgpu */
magma_int_t
magma_sgetrf2_mgpu(
    magma_int_t ngpu,
    magma_int_t m, magma_int_t n, magma_int_t nb, magma_int_t offset,
    magmaFloat_ptr d_lAT[], size_t dlAT_offset, magma_int_t lddat, magma_int_t *ipiv,
    magmaFloat_ptr d_lAP[], size_t dlAP_offset,
    float *W, magma_int_t ldw,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_sgetrf2_msub */
/* magma_sgetrf2_msub */
/* magma_sgetrf2_msub */
magma_int_t
magma_sgetrf2_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_int_t m, magma_int_t n, magma_int_t nb, magma_int_t offset,
    magmaFloat_ptr d_lAT[], size_t dlAT_offset, magma_int_t lddat, magma_int_t *ipiv,
    magmaFloat_ptr d_panel[],
    magmaFloat_ptr d_lAP[], size_t dlAP_offset,
    float *W, magma_int_t ldw,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_sgetri_gpu */
/* magma_sgetri_gpu */
/* magma_sgetri_gpu */
magma_int_t
magma_sgetri_gpu(
    magma_int_t n,
    magmaFloat_ptr dA,    size_t dA_offset,    magma_int_t ldda, magma_int_t *ipiv,
    magmaFloat_ptr dwork, size_t dwork_offset, magma_int_t lwork,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_sgetrs_gpu */
/* magma_sgetrs_gpu */
/* magma_sgetrs_gpu */
magma_int_t
magma_sgetrs_gpu(
    magma_trans_t trans, magma_int_t n, magma_int_t nrhs,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magmaFloat_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_slabrd_gpu */
/* magma_slabrd_gpu */
/* magma_slabrd_gpu */
magma_int_t
magma_slabrd_gpu(
    magma_int_t m, magma_int_t n, magma_int_t nb,
    float     *A,                   magma_int_t lda,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    float *d, float *e, float *tauq, float *taup,
    float     *X,                   magma_int_t ldx,
    magmaFloat_ptr dX, size_t dX_offset, magma_int_t lddx,
    float     *Y,                   magma_int_t ldy,
    magmaFloat_ptr dY, size_t dY_offset, magma_int_t lddy,
    magma_queue_t queue);

/* magma_slarfb_gpu */
/* magma_slarfb_gpu */
/* magma_slarfb_gpu */
magma_int_t
magma_slarfb_gpu(
    magma_side_t side, magma_trans_t trans, magma_direct_t direct, magma_storev_t storev,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloat_const_ptr dV, size_t dV_offset,    magma_int_t lddv,
    magmaFloat_const_ptr dT, size_t dT_offset,    magma_int_t lddt,
    magmaFloat_ptr dC,       size_t dC_offset,    magma_int_t lddc,
    magmaFloat_ptr dwork,    size_t dwork_offset, magma_int_t ldwork,
    magma_queue_t queue);

/* magma_slauum_gpu */
/* magma_slauum_gpu */
/* magma_slauum_gpu */
magma_int_t
magma_slauum_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_sposv_gpu */
/* magma_sposv_gpu */
/* magma_sposv_gpu */
magma_int_t
magma_sposv_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_spotrf_gpu */
/* magma_spotrf_gpu */
/* magma_spotrf_gpu */
magma_int_t
magma_spotrf_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_spotrf_mgpu */
/* magma_spotrf_mgpu */
/* magma_spotrf_mgpu */
magma_int_t
magma_spotrf_mgpu(
    magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t n,
    magmaFloat_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_spotrf_msub */
/* magma_spotrf_msub */
/* magma_spotrf_msub */
magma_int_t
magma_spotrf_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t n,
    magmaFloat_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_spotrf2_gpu */
/* magma_spotrf2_gpu */
/* magma_spotrf2_gpu */
magma_int_t
magma_spotrf2_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_spotrf2_mgpu */
/* magma_spotrf2_mgpu */
/* magma_spotrf2_mgpu */
magma_int_t
magma_spotrf2_mgpu(
    magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
    magmaFloat_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr d_lP[],                   magma_int_t lddp,
    float *A, magma_int_t lda, magma_int_t h,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_spotrf2_msub */
/* magma_spotrf2_msub */
/* magma_spotrf2_msub */
magma_int_t
magma_spotrf2_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
    magmaFloat_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr d_lP[],                   magma_int_t lddp,
    float *A, magma_int_t lda, magma_int_t h,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_spotri_gpu */
/* magma_spotri_gpu */
/* magma_spotri_gpu */
magma_int_t
magma_spotri_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_spotrs_gpu */
/* magma_spotrs_gpu */
/* magma_spotrs_gpu */
magma_int_t
magma_spotrs_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_strtri_gpu */
/* magma_strtri_gpu */
/* magma_strtri_gpu */
magma_int_t
magma_strtri_gpu(
    magma_uplo_t uplo, magma_diag_t diag, magma_int_t n,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_sormqr_gpu */
/* magma_sormqr_gpu */
/* magma_sormqr_gpu */
magma_int_t
magma_sormqr_gpu(
    magma_side_t side, magma_trans_t trans,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
    float *tau,
    magmaFloat_ptr dC, size_t dC_offset, magma_int_t lddc,
    float *hwork, magma_int_t lwork,
    magmaFloat_ptr dT, size_t dT_offset, magma_int_t nb,
    magma_queue_t queue,
    magma_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA utility function definitions
*/

extern const float MAGMA_S_NAN;
extern const float MAGMA_S_INF;

/* magma_snan_inf */
/* magma_snan_inf */
/* magma_snan_inf */
magma_int_t
magma_snan_inf(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    const float *A, magma_int_t lda,
    magma_int_t *cnt_nan,
    magma_int_t *cnt_inf);

/* magma_snan_inf_gpu */
/* magma_snan_inf_gpu */
/* magma_snan_inf_gpu */
magma_int_t
magma_snan_inf_gpu(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magmaFloat_const_ptr dA, magma_int_t dA_offset, magma_int_t ldda,
    magma_int_t *cnt_nan,
    magma_int_t *cnt_inf,
    magma_queue_t queue);

void magma_sprint(
    magma_int_t m, magma_int_t n,
    const float *A, magma_int_t lda);

void magma_sprint_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloat_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue);

void spanel_to_q(
    magma_uplo_t uplo, magma_int_t ib,
    float *A, magma_int_t lda,
    float *work);

void sq_to_panel(
    magma_uplo_t uplo, magma_int_t ib,
    float *A, magma_int_t lda,
    float *work);

#ifdef __cplusplus
}
#endif

#undef REAL

#endif /* MAGMA_S_H */
