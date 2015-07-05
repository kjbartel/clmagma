/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from magma_z.h normal z -> d, Sat Nov 15 00:21:34 2014
*/

#ifndef MAGMA_D_H
#define MAGMA_D_H

#include "magma_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA Auxiliary functions to get the NB used
*/
magma_int_t magma_get_dlaex3_m_nb();       // defined in dlaex3_m.cpp

magma_int_t magma_get_dpotrf_nb( magma_int_t m );
magma_int_t magma_get_dgetrf_nb( magma_int_t m );
magma_int_t magma_get_dgetri_nb( magma_int_t m );
magma_int_t magma_get_dgeqp3_nb( magma_int_t m );
magma_int_t magma_get_dgeqrf_nb( magma_int_t m );
magma_int_t magma_get_dgeqlf_nb( magma_int_t m );
magma_int_t magma_get_dgehrd_nb( magma_int_t m );
magma_int_t magma_get_dsytrd_nb( magma_int_t m );
magma_int_t magma_get_dsytrf_nb( magma_int_t m );
magma_int_t magma_get_dgelqf_nb( magma_int_t m );
magma_int_t magma_get_dgebrd_nb( magma_int_t m );
magma_int_t magma_get_dsygst_nb( magma_int_t m );
magma_int_t magma_get_dgesvd_nb( magma_int_t m );
magma_int_t magma_get_dsygst_nb_m( magma_int_t m );
magma_int_t magma_get_dbulge_nb( magma_int_t m, magma_int_t nbthreads );
magma_int_t magma_get_dbulge_nb_mgpu( magma_int_t m );
magma_int_t magma_dbulge_get_Vblksiz( magma_int_t m, magma_int_t nb, magma_int_t nbthreads );
magma_int_t magma_get_dbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU (alphabetical order)
*/

/* magma_dgebrd */
/* magma_dgebrd */
/* magma_dgebrd */
magma_int_t
magma_dgebrd(
    magma_int_t m, magma_int_t n,
    double *A, magma_int_t lda,
    double *d, double *e,
    double *tauq, double *taup,
    double *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dgeev */
/* magma_dgeev */
/* magma_dgeev */
magma_int_t
magma_dgeev(
    magma_vec_t jobvl, magma_vec_t jobvr, magma_int_t n,
    double *A, magma_int_t lda,
    #ifdef COMPLEX
    double *w,
    #else
    double *wr, double *wi,
    #endif
    double *VL, magma_int_t ldvl,
    double *VR, magma_int_t ldvr,
    double *work, magma_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dgehrd */
/* magma_dgehrd */
/* magma_dgehrd */
magma_int_t
magma_dgehrd(
    magma_int_t n, magma_int_t ilo, magma_int_t ihi,
    double *A, magma_int_t lda,
    double *tau,
    double *work, magma_int_t lwork,
    magmaDouble_ptr dT, size_t dT_offset,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dgeqrf */
/* magma_dgeqrf */
/* magma_dgeqrf */
magma_int_t
magma_dgeqrf(
    magma_int_t m, magma_int_t n,
    double *A, magma_int_t lda,
    double *tau,
    double *work, magma_int_t lwork,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_dgesv */
/* magma_dgesv */
/* magma_dgesv */
magma_int_t
magma_dgesv(
    magma_int_t n, magma_int_t nrhs,
    double *A, magma_int_t lda, magma_int_t *ipiv,
    double *B, magma_int_t ldb,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_dgesvd */
/* magma_dgesvd */
/* magma_dgesvd */
magma_int_t
magma_dgesvd(
    magma_vec_t jobu, magma_vec_t jobvt, magma_int_t m, magma_int_t n,
    double *A,    magma_int_t lda, double *s,
    double *U,    magma_int_t ldu,
    double *VT,   magma_int_t ldvt,
    double *work, magma_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dgetrf */
/* magma_dgetrf */
/* magma_dgetrf */
magma_int_t
magma_dgetrf(
    magma_int_t m, magma_int_t n,
    double *A, magma_int_t lda, magma_int_t *ipiv,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_dsyevd */
/* magma_dsyevd */
/* magma_dsyevd */
magma_int_t
magma_dsyevd(
    magma_vec_t jobz, magma_uplo_t uplo, magma_int_t n,
    double *A, magma_int_t lda,
    double *w,
    double *work, magma_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_int_t lrwork,
    #endif
    magma_int_t *iwork, magma_int_t liwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dsytrd */
/* magma_dsytrd */
/* magma_dsytrd */
magma_int_t
magma_dsytrd(
    magma_uplo_t uplo, magma_int_t n,
    double *A, magma_int_t lda,
    double *d, double *e, double *tau,
    double *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dsytrf */
/* magma_dsytrf */
/* magma_dsytrf */
magma_int_t
magma_dsytrf(
    magma_uplo_t uplo, magma_int_t n,
    double *A, magma_int_t lda,
    magma_int_t *ipiv,
    magma_queue_t queue,
    magma_int_t *info);

#ifdef REAL
// only applicable to real [sd] precisions
/* magma_dlaex0 */
/* magma_dlaex0 */
/* magma_dlaex0 */
magma_int_t
magma_dlaex0(
    magma_int_t n, double *d, double *e,
    double *Q, magma_int_t ldq,
    double *work, magma_int_t *iwork,
    magmaDouble_ptr dwork,
    magma_range_t range, double vl, double vu, magma_int_t il, magma_int_t iu,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dlaex1 */
/* magma_dlaex1 */
/* magma_dlaex1 */
magma_int_t
magma_dlaex1(
    magma_int_t n, double *d,
    double *Q, magma_int_t ldq,
    magma_int_t *indxq, double rho, magma_int_t cutpnt,
    double *work, magma_int_t *iwork,
    magmaDouble_ptr dwork,
    magma_range_t range, double vl, double vu, magma_int_t il, magma_int_t iu,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dlaex3 */
/* magma_dlaex3 */
/* magma_dlaex3 */
magma_int_t
magma_dlaex3(
    magma_int_t k, magma_int_t n, magma_int_t n1, double *d,
    double *Q, magma_int_t ldq,
    double rho,
    double *dlamda, double *Q2, magma_int_t *indx,
    magma_int_t *ctot, double *w, double *s, magma_int_t *indxq,
    magmaDouble_ptr dwork,
    magma_range_t range, double vl, double vu, magma_int_t il, magma_int_t iu,
    magma_queue_t queue,
    magma_int_t *info);
#endif  // REAL

/* magma_dlasyf_gpu */
/* magma_dlasyf_gpu */
/* magma_dlasyf_gpu */
magma_int_t
magma_dlasyf_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nb, magma_int_t *kb,
    double *hA, magma_int_t lda,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_int_t *ipiv,
    magmaDouble_ptr dW, size_t dW_offset, magma_int_t lddw,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dlahr2 */
/* magma_dlahr2 */
/* magma_dlahr2 */
magma_int_t
magma_dlahr2(
    magma_int_t n, magma_int_t k, magma_int_t nb,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dV, size_t dV_offset, magma_int_t lddv,
    double *A,  magma_int_t lda,
    double *tau,
    double *T,  magma_int_t ldt,
    double *Y,  magma_int_t ldy,
    magma_queue_t queue);

/* magma_dlahru */
/* magma_dlahru */
/* magma_dlahru */
magma_int_t
magma_dlahru(
    magma_int_t n, magma_int_t ihi, magma_int_t k, magma_int_t nb,
    double     *A, magma_int_t lda,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dY, size_t dY_offset, magma_int_t lddy,
    magmaDouble_ptr dV, size_t dV_offset, magma_int_t lddv,
    magmaDouble_ptr dT, size_t dT_offset,
    magmaDouble_ptr dwork, size_t dwork_offset,
    magma_queue_t queue);

/* magma_dlatrd */
/* magma_dlatrd */
/* magma_dlatrd */
magma_int_t
magma_dlatrd(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nb,
    double *A,  magma_int_t lda,
    double *e, double *tau,
    double *W,  magma_int_t ldw,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dW, size_t dW_offset, magma_int_t lddw,
    magma_queue_t queue);

/* magma_dposv */
/* magma_dposv */
/* magma_dposv */
magma_int_t
magma_dposv(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    double *A, magma_int_t lda,
    double *B, magma_int_t ldb,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_dpotrf */
/* magma_dpotrf */
/* magma_dpotrf */
magma_int_t
magma_dpotrf(
    magma_uplo_t uplo, magma_int_t n,
    double *A, magma_int_t lda,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_dstedx */
/* magma_dstedx */
/* magma_dstedx */
magma_int_t
magma_dstedx(
    magma_range_t range, magma_int_t n, double vl, double vu,
    magma_int_t il, magma_int_t iu, double *d, double *e,
    double *Z, magma_int_t ldz,
    double *rwork, magma_int_t lrwork,
    magma_int_t *iwork, magma_int_t liwork,
    magmaDouble_ptr dwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dorghr */
/* magma_dorghr */
/* magma_dorghr */
magma_int_t
magma_dorghr(
    magma_int_t n, magma_int_t ilo, magma_int_t ihi,
    double *A, magma_int_t lda,
    double *tau,
    magmaDouble_ptr dT, size_t dT_offset, magma_int_t nb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dorgqr */
/* magma_dorgqr */
/* magma_dorgqr */
magma_int_t
magma_dorgqr(
    magma_int_t m, magma_int_t n, magma_int_t k,
    double *A, magma_int_t lda,
    double *tau,
    magmaDouble_ptr dT, size_t dT_offset, magma_int_t nb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dormql */
/* magma_dormql */
/* magma_dormql */
magma_int_t
magma_dormql(
    magma_side_t side, magma_trans_t trans,
    magma_int_t m, magma_int_t n, magma_int_t k,
    double *A, magma_int_t lda,
    double *tau,
    double *C, magma_int_t ldc,
    double *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dormqr */
/* magma_dormqr */
/* magma_dormqr */
magma_int_t
magma_dormqr(
    magma_side_t side, magma_trans_t trans,
    magma_int_t m, magma_int_t n, magma_int_t k,
    double *A, magma_int_t lda,
    double *tau,
    double *C, magma_int_t ldc,
    double *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dormtr */
/* magma_dormtr */
/* magma_dormtr */
magma_int_t
magma_dormtr(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t m, magma_int_t n,
    double *A,    magma_int_t lda,
    double *tau,
    double *C,    magma_int_t ldc,
    double *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on GPU (alphabetical order)
*/

/* magma_dgels_gpu */
/* magma_dgels_gpu */
/* magma_dgels_gpu */
magma_int_t
magma_dgels_gpu(
    magma_trans_t trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dB, size_t dB_offset, magma_int_t lddb,
    double *hwork, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dgeqr2x3_gpu */
/* magma_dgeqr2x3_gpu */
/* magma_dgeqr2x3_gpu */
magma_int_t
magma_dgeqr2x3_gpu(
    magma_int_t m, magma_int_t n,
    magmaDouble_ptr dA,    size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dtau,  size_t dtau_offset,
    magmaDouble_ptr dT,    size_t dT_offset,
    magmaDouble_ptr ddA,   size_t ddA_offset,
    magmaDouble_ptr        dwork, size_t dwork_offset,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dgeqrf_gpu */
/* magma_dgeqrf_gpu */
/* magma_dgeqrf_gpu */
magma_int_t
magma_dgeqrf_gpu(
    magma_int_t m, magma_int_t n,
    magmaDouble_ptr dA, size_t dA_offset,  magma_int_t ldda,
    double *tau,
    magmaDouble_ptr dT, size_t dT_offset,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dgeqrf_msub */
/* magma_dgeqrf_msub */
/* magma_dgeqrf_msub */
magma_int_t
magma_dgeqrf_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaDouble_ptr dlA[], magma_int_t ldda,
    double *tau,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_dgeqrf2_2q_gpu */
/* magma_dgeqrf2_2q_gpu */
/* magma_dgeqrf2_2q_gpu */
magma_int_t
magma_dgeqrf2_2q_gpu(
    magma_int_t m, magma_int_t n,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    double *tau,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_dgeqrf2_gpu */
/* magma_dgeqrf2_gpu */
/* magma_dgeqrf2_gpu */
magma_int_t
magma_dgeqrf2_gpu(
    magma_int_t m, magma_int_t n,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    double *tau,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_dgeqrf2_mgpu */
/* magma_dgeqrf2_mgpu */
/* magma_dgeqrf2_mgpu */
magma_int_t
magma_dgeqrf2_mgpu(
    magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaDouble_ptr dlA[], magma_int_t ldda,
    double *tau,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_dgeqrs_gpu */
/* magma_dgeqrs_gpu */
/* magma_dgeqrs_gpu */
magma_int_t
magma_dgeqrs_gpu(
    magma_int_t m, magma_int_t n, magma_int_t nrhs,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    double *tau,
    magmaDouble_ptr dT, size_t dT_offset,
    magmaDouble_ptr dB, size_t dB_offset, magma_int_t lddb,
    double *hwork, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dgesv_gpu */
/* magma_dgesv_gpu */
/* magma_dgesv_gpu */
magma_int_t
magma_dgesv_gpu(
    magma_int_t n, magma_int_t nrhs,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magmaDouble_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dgetrf_gpu */
/* magma_dgetrf_gpu */
/* magma_dgetrf_gpu */
magma_int_t
magma_dgetrf_gpu(
    magma_int_t m, magma_int_t n,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dgetrf_mgpu */
/* magma_dgetrf_mgpu */
/* magma_dgetrf_mgpu */
magma_int_t
magma_dgetrf_mgpu(
    magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaDouble_ptr d_lA[], size_t dlA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_dgetrf_msub */
/* magma_dgetrf_msub */
/* magma_dgetrf_msub */
magma_int_t
magma_dgetrf_msub(
    magma_trans_t trans, magma_int_t num_subs, magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaDouble_ptr d_lA[], size_t dlA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_dgetrf2_gpu */
/* magma_dgetrf2_gpu */
/* magma_dgetrf2_gpu */
magma_int_t
magma_dgetrf2_gpu(
    magma_int_t m, magma_int_t n,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_dgetrf2_mgpu */
/* magma_dgetrf2_mgpu */
/* magma_dgetrf2_mgpu */
magma_int_t
magma_dgetrf2_mgpu(
    magma_int_t ngpu,
    magma_int_t m, magma_int_t n, magma_int_t nb, magma_int_t offset,
    magmaDouble_ptr d_lAT[], size_t dlAT_offset, magma_int_t lddat, magma_int_t *ipiv,
    magmaDouble_ptr d_lAP[], size_t dlAP_offset,
    double *W, magma_int_t ldw,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_dgetrf2_msub */
/* magma_dgetrf2_msub */
/* magma_dgetrf2_msub */
magma_int_t
magma_dgetrf2_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_int_t m, magma_int_t n, magma_int_t nb, magma_int_t offset,
    magmaDouble_ptr d_lAT[], size_t dlAT_offset, magma_int_t lddat, magma_int_t *ipiv,
    magmaDouble_ptr d_panel[],
    magmaDouble_ptr d_lAP[], size_t dlAP_offset,
    double *W, magma_int_t ldw,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_dgetri_gpu */
/* magma_dgetri_gpu */
/* magma_dgetri_gpu */
magma_int_t
magma_dgetri_gpu(
    magma_int_t n,
    magmaDouble_ptr dA,    size_t dA_offset,    magma_int_t ldda, magma_int_t *ipiv,
    magmaDouble_ptr dwork, size_t dwork_offset, magma_int_t lwork,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_dgetrs_gpu */
/* magma_dgetrs_gpu */
/* magma_dgetrs_gpu */
magma_int_t
magma_dgetrs_gpu(
    magma_trans_t trans, magma_int_t n, magma_int_t nrhs,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magmaDouble_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dlabrd_gpu */
/* magma_dlabrd_gpu */
/* magma_dlabrd_gpu */
magma_int_t
magma_dlabrd_gpu(
    magma_int_t m, magma_int_t n, magma_int_t nb,
    double     *A,                   magma_int_t lda,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    double *d, double *e, double *tauq, double *taup,
    double     *X,                   magma_int_t ldx,
    magmaDouble_ptr dX, size_t dX_offset, magma_int_t lddx,
    double     *Y,                   magma_int_t ldy,
    magmaDouble_ptr dY, size_t dY_offset, magma_int_t lddy,
    magma_queue_t queue);

/* magma_dlarfb_gpu */
/* magma_dlarfb_gpu */
/* magma_dlarfb_gpu */
magma_int_t
magma_dlarfb_gpu(
    magma_side_t side, magma_trans_t trans, magma_direct_t direct, magma_storev_t storev,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaDouble_const_ptr dV, size_t dV_offset,    magma_int_t lddv,
    magmaDouble_const_ptr dT, size_t dT_offset,    magma_int_t lddt,
    magmaDouble_ptr dC,       size_t dC_offset,    magma_int_t lddc,
    magmaDouble_ptr dwork,    size_t dwork_offset, magma_int_t ldwork,
    magma_queue_t queue);

/* magma_dlauum_gpu */
/* magma_dlauum_gpu */
/* magma_dlauum_gpu */
magma_int_t
magma_dlauum_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dposv_gpu */
/* magma_dposv_gpu */
/* magma_dposv_gpu */
magma_int_t
magma_dposv_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dpotrf_gpu */
/* magma_dpotrf_gpu */
/* magma_dpotrf_gpu */
magma_int_t
magma_dpotrf_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dpotrf_mgpu */
/* magma_dpotrf_mgpu */
/* magma_dpotrf_mgpu */
magma_int_t
magma_dpotrf_mgpu(
    magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t n,
    magmaDouble_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_dpotrf_msub */
/* magma_dpotrf_msub */
/* magma_dpotrf_msub */
magma_int_t
magma_dpotrf_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t n,
    magmaDouble_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_dpotrf2_gpu */
/* magma_dpotrf2_gpu */
/* magma_dpotrf2_gpu */
magma_int_t
magma_dpotrf2_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_dpotrf2_mgpu */
/* magma_dpotrf2_mgpu */
/* magma_dpotrf2_mgpu */
magma_int_t
magma_dpotrf2_mgpu(
    magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
    magmaDouble_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr d_lP[],                   magma_int_t lddp,
    double *A, magma_int_t lda, magma_int_t h,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_dpotrf2_msub */
/* magma_dpotrf2_msub */
/* magma_dpotrf2_msub */
magma_int_t
magma_dpotrf2_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
    magmaDouble_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr d_lP[],                   magma_int_t lddp,
    double *A, magma_int_t lda, magma_int_t h,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_dpotri_gpu */
/* magma_dpotri_gpu */
/* magma_dpotri_gpu */
magma_int_t
magma_dpotri_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_dpotrs_gpu */
/* magma_dpotrs_gpu */
/* magma_dpotrs_gpu */
magma_int_t
magma_dpotrs_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_dtrtri_gpu */
/* magma_dtrtri_gpu */
/* magma_dtrtri_gpu */
magma_int_t
magma_dtrtri_gpu(
    magma_uplo_t uplo, magma_diag_t diag, magma_int_t n,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_dormqr_gpu */
/* magma_dormqr_gpu */
/* magma_dormqr_gpu */
magma_int_t
magma_dormqr_gpu(
    magma_side_t side, magma_trans_t trans,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    double *tau,
    magmaDouble_ptr dC, size_t dC_offset, magma_int_t lddc,
    double *hwork, magma_int_t lwork,
    magmaDouble_ptr dT, size_t dT_offset, magma_int_t nb,
    magma_queue_t queue,
    magma_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA utility function definitions
*/

extern const double MAGMA_D_NAN;
extern const double MAGMA_D_INF;

/* magma_dnan_inf */
/* magma_dnan_inf */
/* magma_dnan_inf */
magma_int_t
magma_dnan_inf(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    const double *A, magma_int_t lda,
    magma_int_t *cnt_nan,
    magma_int_t *cnt_inf);

/* magma_dnan_inf_gpu */
/* magma_dnan_inf_gpu */
/* magma_dnan_inf_gpu */
magma_int_t
magma_dnan_inf_gpu(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA, magma_int_t dA_offset, magma_int_t ldda,
    magma_int_t *cnt_nan,
    magma_int_t *cnt_inf,
    magma_queue_t queue);

void magma_dprint(
    magma_int_t m, magma_int_t n,
    const double *A, magma_int_t lda);

void magma_dprint_gpu(
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue);

void dpanel_to_q(
    magma_uplo_t uplo, magma_int_t ib,
    double *A, magma_int_t lda,
    double *work);

void dq_to_panel(
    magma_uplo_t uplo, magma_int_t ib,
    double *A, magma_int_t lda,
    double *work);

#ifdef __cplusplus
}
#endif

#undef REAL

#endif /* MAGMA_D_H */
