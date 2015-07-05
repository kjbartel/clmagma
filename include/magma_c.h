/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from magma_z.h normal z -> c, Sat Nov 15 00:21:34 2014
*/

#ifndef MAGMA_C_H
#define MAGMA_C_H

#include "magma_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA Auxiliary functions to get the NB used
*/
magma_int_t magma_get_slaex3_m_nb();       // defined in slaex3_m.cpp

magma_int_t magma_get_cpotrf_nb( magma_int_t m );
magma_int_t magma_get_cgetrf_nb( magma_int_t m );
magma_int_t magma_get_cgetri_nb( magma_int_t m );
magma_int_t magma_get_cgeqp3_nb( magma_int_t m );
magma_int_t magma_get_cgeqrf_nb( magma_int_t m );
magma_int_t magma_get_cgeqlf_nb( magma_int_t m );
magma_int_t magma_get_cgehrd_nb( magma_int_t m );
magma_int_t magma_get_chetrd_nb( magma_int_t m );
magma_int_t magma_get_chetrf_nb( magma_int_t m );
magma_int_t magma_get_cgelqf_nb( magma_int_t m );
magma_int_t magma_get_cgebrd_nb( magma_int_t m );
magma_int_t magma_get_chegst_nb( magma_int_t m );
magma_int_t magma_get_cgesvd_nb( magma_int_t m );
magma_int_t magma_get_chegst_nb_m( magma_int_t m );
magma_int_t magma_get_cbulge_nb( magma_int_t m, magma_int_t nbthreads );
magma_int_t magma_get_cbulge_nb_mgpu( magma_int_t m );
magma_int_t magma_cbulge_get_Vblksiz( magma_int_t m, magma_int_t nb, magma_int_t nbthreads );
magma_int_t magma_get_cbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU (alphabetical order)
*/

/* magma_cgebrd */
/* magma_cgebrd */
/* magma_cgebrd */
magma_int_t
magma_cgebrd(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex *A, magma_int_t lda,
    float *d, float *e,
    magmaFloatComplex *tauq, magmaFloatComplex *taup,
    magmaFloatComplex *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cgeev */
/* magma_cgeev */
/* magma_cgeev */
magma_int_t
magma_cgeev(
    magma_vec_t jobvl, magma_vec_t jobvr, magma_int_t n,
    magmaFloatComplex *A, magma_int_t lda,
    #ifdef COMPLEX
    magmaFloatComplex *w,
    #else
    float *wr, float *wi,
    #endif
    magmaFloatComplex *VL, magma_int_t ldvl,
    magmaFloatComplex *VR, magma_int_t ldvr,
    magmaFloatComplex *work, magma_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cgehrd */
/* magma_cgehrd */
/* magma_cgehrd */
magma_int_t
magma_cgehrd(
    magma_int_t n, magma_int_t ilo, magma_int_t ihi,
    magmaFloatComplex *A, magma_int_t lda,
    magmaFloatComplex *tau,
    magmaFloatComplex *work, magma_int_t lwork,
    magmaFloatComplex_ptr dT, size_t dT_offset,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cgeqrf */
/* magma_cgeqrf */
/* magma_cgeqrf */
magma_int_t
magma_cgeqrf(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex *A, magma_int_t lda,
    magmaFloatComplex *tau,
    magmaFloatComplex *work, magma_int_t lwork,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_cgesv */
/* magma_cgesv */
/* magma_cgesv */
magma_int_t
magma_cgesv(
    magma_int_t n, magma_int_t nrhs,
    magmaFloatComplex *A, magma_int_t lda, magma_int_t *ipiv,
    magmaFloatComplex *B, magma_int_t ldb,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_cgesvd */
/* magma_cgesvd */
/* magma_cgesvd */
magma_int_t
magma_cgesvd(
    magma_vec_t jobu, magma_vec_t jobvt, magma_int_t m, magma_int_t n,
    magmaFloatComplex *A,    magma_int_t lda, float *s,
    magmaFloatComplex *U,    magma_int_t ldu,
    magmaFloatComplex *VT,   magma_int_t ldvt,
    magmaFloatComplex *work, magma_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cgetrf */
/* magma_cgetrf */
/* magma_cgetrf */
magma_int_t
magma_cgetrf(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex *A, magma_int_t lda, magma_int_t *ipiv,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_cheevd */
/* magma_cheevd */
/* magma_cheevd */
magma_int_t
magma_cheevd(
    magma_vec_t jobz, magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex *A, magma_int_t lda,
    float *w,
    magmaFloatComplex *work, magma_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_int_t lrwork,
    #endif
    magma_int_t *iwork, magma_int_t liwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_chetrd */
/* magma_chetrd */
/* magma_chetrd */
magma_int_t
magma_chetrd(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex *A, magma_int_t lda,
    float *d, float *e, magmaFloatComplex *tau,
    magmaFloatComplex *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_chetrf */
/* magma_chetrf */
/* magma_chetrf */
magma_int_t
magma_chetrf(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex *A, magma_int_t lda,
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

/* magma_clahef_gpu */
/* magma_clahef_gpu */
/* magma_clahef_gpu */
magma_int_t
magma_clahef_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nb, magma_int_t *kb,
    magmaFloatComplex *hA, magma_int_t lda,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_int_t *ipiv,
    magmaFloatComplex_ptr dW, size_t dW_offset, magma_int_t lddw,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_clahr2 */
/* magma_clahr2 */
/* magma_clahr2 */
magma_int_t
magma_clahr2(
    magma_int_t n, magma_int_t k, magma_int_t nb,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr dV, size_t dV_offset, magma_int_t lddv,
    magmaFloatComplex *A,  magma_int_t lda,
    magmaFloatComplex *tau,
    magmaFloatComplex *T,  magma_int_t ldt,
    magmaFloatComplex *Y,  magma_int_t ldy,
    magma_queue_t queue);

/* magma_clahru */
/* magma_clahru */
/* magma_clahru */
magma_int_t
magma_clahru(
    magma_int_t n, magma_int_t ihi, magma_int_t k, magma_int_t nb,
    magmaFloatComplex     *A, magma_int_t lda,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr dY, size_t dY_offset, magma_int_t lddy,
    magmaFloatComplex_ptr dV, size_t dV_offset, magma_int_t lddv,
    magmaFloatComplex_ptr dT, size_t dT_offset,
    magmaFloatComplex_ptr dwork, size_t dwork_offset,
    magma_queue_t queue);

/* magma_clatrd */
/* magma_clatrd */
/* magma_clatrd */
magma_int_t
magma_clatrd(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nb,
    magmaFloatComplex *A,  magma_int_t lda,
    float *e, magmaFloatComplex *tau,
    magmaFloatComplex *W,  magma_int_t ldw,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr dW, size_t dW_offset, magma_int_t lddw,
    magma_queue_t queue);

/* magma_cposv */
/* magma_cposv */
/* magma_cposv */
magma_int_t
magma_cposv(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    magmaFloatComplex *A, magma_int_t lda,
    magmaFloatComplex *B, magma_int_t ldb,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_cpotrf */
/* magma_cpotrf */
/* magma_cpotrf */
magma_int_t
magma_cpotrf(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex *A, magma_int_t lda,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_cstedx */
/* magma_cstedx */
/* magma_cstedx */
magma_int_t
magma_cstedx(
    magma_range_t range, magma_int_t n, float vl, float vu,
    magma_int_t il, magma_int_t iu, float *d, float *e,
    magmaFloatComplex *Z, magma_int_t ldz,
    float *rwork, magma_int_t lrwork,
    magma_int_t *iwork, magma_int_t liwork,
    magmaFloat_ptr dwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cunghr */
/* magma_cunghr */
/* magma_cunghr */
magma_int_t
magma_cunghr(
    magma_int_t n, magma_int_t ilo, magma_int_t ihi,
    magmaFloatComplex *A, magma_int_t lda,
    magmaFloatComplex *tau,
    magmaFloatComplex_ptr dT, size_t dT_offset, magma_int_t nb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cungqr */
/* magma_cungqr */
/* magma_cungqr */
magma_int_t
magma_cungqr(
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloatComplex *A, magma_int_t lda,
    magmaFloatComplex *tau,
    magmaFloatComplex_ptr dT, size_t dT_offset, magma_int_t nb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cunmql */
/* magma_cunmql */
/* magma_cunmql */
magma_int_t
magma_cunmql(
    magma_side_t side, magma_trans_t trans,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloatComplex *A, magma_int_t lda,
    magmaFloatComplex *tau,
    magmaFloatComplex *C, magma_int_t ldc,
    magmaFloatComplex *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cunmqr */
/* magma_cunmqr */
/* magma_cunmqr */
magma_int_t
magma_cunmqr(
    magma_side_t side, magma_trans_t trans,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloatComplex *A, magma_int_t lda,
    magmaFloatComplex *tau,
    magmaFloatComplex *C, magma_int_t ldc,
    magmaFloatComplex *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cunmtr */
/* magma_cunmtr */
/* magma_cunmtr */
magma_int_t
magma_cunmtr(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex *A,    magma_int_t lda,
    magmaFloatComplex *tau,
    magmaFloatComplex *C,    magma_int_t ldc,
    magmaFloatComplex *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on GPU (alphabetical order)
*/

/* magma_cgels_gpu */
/* magma_cgels_gpu */
/* magma_cgels_gpu */
magma_int_t
magma_cgels_gpu(
    magma_trans_t trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaFloatComplex *hwork, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cgeqr2x3_gpu */
/* magma_cgeqr2x3_gpu */
/* magma_cgeqr2x3_gpu */
magma_int_t
magma_cgeqr2x3_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr dA,    size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr dtau,  size_t dtau_offset,
    magmaFloatComplex_ptr dT,    size_t dT_offset,
    magmaFloatComplex_ptr ddA,   size_t ddA_offset,
    magmaFloat_ptr        dwork, size_t dwork_offset,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cgeqrf_gpu */
/* magma_cgeqrf_gpu */
/* magma_cgeqrf_gpu */
magma_int_t
magma_cgeqrf_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr dA, size_t dA_offset,  magma_int_t ldda,
    magmaFloatComplex *tau,
    magmaFloatComplex_ptr dT, size_t dT_offset,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cgeqrf_msub */
/* magma_cgeqrf_msub */
/* magma_cgeqrf_msub */
magma_int_t
magma_cgeqrf_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr dlA[], magma_int_t ldda,
    magmaFloatComplex *tau,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_cgeqrf2_2q_gpu */
/* magma_cgeqrf2_2q_gpu */
/* magma_cgeqrf2_2q_gpu */
magma_int_t
magma_cgeqrf2_2q_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex *tau,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_cgeqrf2_gpu */
/* magma_cgeqrf2_gpu */
/* magma_cgeqrf2_gpu */
magma_int_t
magma_cgeqrf2_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex *tau,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_cgeqrf2_mgpu */
/* magma_cgeqrf2_mgpu */
/* magma_cgeqrf2_mgpu */
magma_int_t
magma_cgeqrf2_mgpu(
    magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr dlA[], magma_int_t ldda,
    magmaFloatComplex *tau,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_cgeqrs_gpu */
/* magma_cgeqrs_gpu */
/* magma_cgeqrs_gpu */
magma_int_t
magma_cgeqrs_gpu(
    magma_int_t m, magma_int_t n, magma_int_t nrhs,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex *tau,
    magmaFloatComplex_ptr dT, size_t dT_offset,
    magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaFloatComplex *hwork, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cgesv_gpu */
/* magma_cgesv_gpu */
/* magma_cgesv_gpu */
magma_int_t
magma_cgesv_gpu(
    magma_int_t n, magma_int_t nrhs,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cgetrf_gpu */
/* magma_cgetrf_gpu */
/* magma_cgetrf_gpu */
magma_int_t
magma_cgetrf_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cgetrf_mgpu */
/* magma_cgetrf_mgpu */
/* magma_cgetrf_mgpu */
magma_int_t
magma_cgetrf_mgpu(
    magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr d_lA[], size_t dlA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_cgetrf_msub */
/* magma_cgetrf_msub */
/* magma_cgetrf_msub */
magma_int_t
magma_cgetrf_msub(
    magma_trans_t trans, magma_int_t num_subs, magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr d_lA[], size_t dlA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_cgetrf2_gpu */
/* magma_cgetrf2_gpu */
/* magma_cgetrf2_gpu */
magma_int_t
magma_cgetrf2_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_cgetrf2_mgpu */
/* magma_cgetrf2_mgpu */
/* magma_cgetrf2_mgpu */
magma_int_t
magma_cgetrf2_mgpu(
    magma_int_t ngpu,
    magma_int_t m, magma_int_t n, magma_int_t nb, magma_int_t offset,
    magmaFloatComplex_ptr d_lAT[], size_t dlAT_offset, magma_int_t lddat, magma_int_t *ipiv,
    magmaFloatComplex_ptr d_lAP[], size_t dlAP_offset,
    magmaFloatComplex *W, magma_int_t ldw,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_cgetrf2_msub */
/* magma_cgetrf2_msub */
/* magma_cgetrf2_msub */
magma_int_t
magma_cgetrf2_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_int_t m, magma_int_t n, magma_int_t nb, magma_int_t offset,
    magmaFloatComplex_ptr d_lAT[], size_t dlAT_offset, magma_int_t lddat, magma_int_t *ipiv,
    magmaFloatComplex_ptr d_panel[],
    magmaFloatComplex_ptr d_lAP[], size_t dlAP_offset,
    magmaFloatComplex *W, magma_int_t ldw,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_cgetri_gpu */
/* magma_cgetri_gpu */
/* magma_cgetri_gpu */
magma_int_t
magma_cgetri_gpu(
    magma_int_t n,
    magmaFloatComplex_ptr dA,    size_t dA_offset,    magma_int_t ldda, magma_int_t *ipiv,
    magmaFloatComplex_ptr dwork, size_t dwork_offset, magma_int_t lwork,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_cgetrs_gpu */
/* magma_cgetrs_gpu */
/* magma_cgetrs_gpu */
magma_int_t
magma_cgetrs_gpu(
    magma_trans_t trans, magma_int_t n, magma_int_t nrhs,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_clabrd_gpu */
/* magma_clabrd_gpu */
/* magma_clabrd_gpu */
magma_int_t
magma_clabrd_gpu(
    magma_int_t m, magma_int_t n, magma_int_t nb,
    magmaFloatComplex     *A,                   magma_int_t lda,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    float *d, float *e, magmaFloatComplex *tauq, magmaFloatComplex *taup,
    magmaFloatComplex     *X,                   magma_int_t ldx,
    magmaFloatComplex_ptr dX, size_t dX_offset, magma_int_t lddx,
    magmaFloatComplex     *Y,                   magma_int_t ldy,
    magmaFloatComplex_ptr dY, size_t dY_offset, magma_int_t lddy,
    magma_queue_t queue);

/* magma_clarfb_gpu */
/* magma_clarfb_gpu */
/* magma_clarfb_gpu */
magma_int_t
magma_clarfb_gpu(
    magma_side_t side, magma_trans_t trans, magma_direct_t direct, magma_storev_t storev,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloatComplex_const_ptr dV, size_t dV_offset,    magma_int_t lddv,
    magmaFloatComplex_const_ptr dT, size_t dT_offset,    magma_int_t lddt,
    magmaFloatComplex_ptr dC,       size_t dC_offset,    magma_int_t lddc,
    magmaFloatComplex_ptr dwork,    size_t dwork_offset, magma_int_t ldwork,
    magma_queue_t queue);

/* magma_clauum_gpu */
/* magma_clauum_gpu */
/* magma_clauum_gpu */
magma_int_t
magma_clauum_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cposv_gpu */
/* magma_cposv_gpu */
/* magma_cposv_gpu */
magma_int_t
magma_cposv_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cpotrf_gpu */
/* magma_cpotrf_gpu */
/* magma_cpotrf_gpu */
magma_int_t
magma_cpotrf_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_cpotrf_mgpu */
/* magma_cpotrf_mgpu */
/* magma_cpotrf_mgpu */
magma_int_t
magma_cpotrf_mgpu(
    magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_cpotrf_msub */
/* magma_cpotrf_msub */
/* magma_cpotrf_msub */
magma_int_t
magma_cpotrf_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_cpotrf2_gpu */
/* magma_cpotrf2_gpu */
/* magma_cpotrf2_gpu */
magma_int_t
magma_cpotrf2_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_cpotrf2_mgpu */
/* magma_cpotrf2_mgpu */
/* magma_cpotrf2_mgpu */
magma_int_t
magma_cpotrf2_mgpu(
    magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
    magmaFloatComplex_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr d_lP[],                   magma_int_t lddp,
    magmaFloatComplex *A, magma_int_t lda, magma_int_t h,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_cpotrf2_msub */
/* magma_cpotrf2_msub */
/* magma_cpotrf2_msub */
magma_int_t
magma_cpotrf2_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
    magmaFloatComplex_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr d_lP[],                   magma_int_t lddp,
    magmaFloatComplex *A, magma_int_t lda, magma_int_t h,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_cpotri_gpu */
/* magma_cpotri_gpu */
/* magma_cpotri_gpu */
magma_int_t
magma_cpotri_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_cpotrs_gpu */
/* magma_cpotrs_gpu */
/* magma_cpotrs_gpu */
magma_int_t
magma_cpotrs_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_ctrtri_gpu */
/* magma_ctrtri_gpu */
/* magma_ctrtri_gpu */
magma_int_t
magma_ctrtri_gpu(
    magma_uplo_t uplo, magma_diag_t diag, magma_int_t n,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_cunmqr_gpu */
/* magma_cunmqr_gpu */
/* magma_cunmqr_gpu */
magma_int_t
magma_cunmqr_gpu(
    magma_side_t side, magma_trans_t trans,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex *tau,
    magmaFloatComplex_ptr dC, size_t dC_offset, magma_int_t lddc,
    magmaFloatComplex *hwork, magma_int_t lwork,
    magmaFloatComplex_ptr dT, size_t dT_offset, magma_int_t nb,
    magma_queue_t queue,
    magma_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA utility function definitions
*/

extern const magmaFloatComplex MAGMA_C_NAN;
extern const magmaFloatComplex MAGMA_C_INF;

/* magma_cnan_inf */
/* magma_cnan_inf */
/* magma_cnan_inf */
magma_int_t
magma_cnan_inf(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    const magmaFloatComplex *A, magma_int_t lda,
    magma_int_t *cnt_nan,
    magma_int_t *cnt_inf);

/* magma_cnan_inf_gpu */
/* magma_cnan_inf_gpu */
/* magma_cnan_inf_gpu */
magma_int_t
magma_cnan_inf_gpu(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr dA, magma_int_t dA_offset, magma_int_t ldda,
    magma_int_t *cnt_nan,
    magma_int_t *cnt_inf,
    magma_queue_t queue);

void magma_cprint(
    magma_int_t m, magma_int_t n,
    const magmaFloatComplex *A, magma_int_t lda);

void magma_cprint_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue);

void cpanel_to_q(
    magma_uplo_t uplo, magma_int_t ib,
    magmaFloatComplex *A, magma_int_t lda,
    magmaFloatComplex *work);

void cq_to_panel(
    magma_uplo_t uplo, magma_int_t ib,
    magmaFloatComplex *A, magma_int_t lda,
    magmaFloatComplex *work);

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif /* MAGMA_C_H */
