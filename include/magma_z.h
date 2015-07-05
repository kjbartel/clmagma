/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @precisions normal z -> s d c
*/

#ifndef MAGMA_Z_H
#define MAGMA_Z_H

#include "magma_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA Auxiliary functions to get the NB used
*/
magma_int_t magma_get_dlaex3_m_nb();       // defined in dlaex3_m.cpp

magma_int_t magma_get_zpotrf_nb( magma_int_t m );
magma_int_t magma_get_zgetrf_nb( magma_int_t m );
magma_int_t magma_get_zgetri_nb( magma_int_t m );
magma_int_t magma_get_zgeqp3_nb( magma_int_t m );
magma_int_t magma_get_zgeqrf_nb( magma_int_t m );
magma_int_t magma_get_zgeqlf_nb( magma_int_t m );
magma_int_t magma_get_zgehrd_nb( magma_int_t m );
magma_int_t magma_get_zhetrd_nb( magma_int_t m );
magma_int_t magma_get_zhetrf_nb( magma_int_t m );
magma_int_t magma_get_zgelqf_nb( magma_int_t m );
magma_int_t magma_get_zgebrd_nb( magma_int_t m );
magma_int_t magma_get_zhegst_nb( magma_int_t m );
magma_int_t magma_get_zgesvd_nb( magma_int_t m );
magma_int_t magma_get_zhegst_nb_m( magma_int_t m );
magma_int_t magma_get_zbulge_nb( magma_int_t m, magma_int_t nbthreads );
magma_int_t magma_get_zbulge_nb_mgpu( magma_int_t m );
magma_int_t magma_zbulge_get_Vblksiz( magma_int_t m, magma_int_t nb, magma_int_t nbthreads );
magma_int_t magma_get_zbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU (alphabetical order)
*/

/* magma_zgebrd */
/* magma_zgebrd */
/* magma_zgebrd */
magma_int_t
magma_zgebrd(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex *A, magma_int_t lda,
    double *d, double *e,
    magmaDoubleComplex *tauq, magmaDoubleComplex *taup,
    magmaDoubleComplex *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zgeev */
/* magma_zgeev */
/* magma_zgeev */
magma_int_t
magma_zgeev(
    magma_vec_t jobvl, magma_vec_t jobvr, magma_int_t n,
    magmaDoubleComplex *A, magma_int_t lda,
    #ifdef COMPLEX
    magmaDoubleComplex *w,
    #else
    double *wr, double *wi,
    #endif
    magmaDoubleComplex *VL, magma_int_t ldvl,
    magmaDoubleComplex *VR, magma_int_t ldvr,
    magmaDoubleComplex *work, magma_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zgehrd */
/* magma_zgehrd */
/* magma_zgehrd */
magma_int_t
magma_zgehrd(
    magma_int_t n, magma_int_t ilo, magma_int_t ihi,
    magmaDoubleComplex *A, magma_int_t lda,
    magmaDoubleComplex *tau,
    magmaDoubleComplex *work, magma_int_t lwork,
    magmaDoubleComplex_ptr dT, size_t dT_offset,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zgeqrf */
/* magma_zgeqrf */
/* magma_zgeqrf */
magma_int_t
magma_zgeqrf(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex *A, magma_int_t lda,
    magmaDoubleComplex *tau,
    magmaDoubleComplex *work, magma_int_t lwork,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_zgesv */
/* magma_zgesv */
/* magma_zgesv */
magma_int_t
magma_zgesv(
    magma_int_t n, magma_int_t nrhs,
    magmaDoubleComplex *A, magma_int_t lda, magma_int_t *ipiv,
    magmaDoubleComplex *B, magma_int_t ldb,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_zgesvd */
/* magma_zgesvd */
/* magma_zgesvd */
magma_int_t
magma_zgesvd(
    magma_vec_t jobu, magma_vec_t jobvt, magma_int_t m, magma_int_t n,
    magmaDoubleComplex *A,    magma_int_t lda, double *s,
    magmaDoubleComplex *U,    magma_int_t ldu,
    magmaDoubleComplex *VT,   magma_int_t ldvt,
    magmaDoubleComplex *work, magma_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zgetrf */
/* magma_zgetrf */
/* magma_zgetrf */
magma_int_t
magma_zgetrf(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex *A, magma_int_t lda, magma_int_t *ipiv,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_zheevd */
/* magma_zheevd */
/* magma_zheevd */
magma_int_t
magma_zheevd(
    magma_vec_t jobz, magma_uplo_t uplo, magma_int_t n,
    magmaDoubleComplex *A, magma_int_t lda,
    double *w,
    magmaDoubleComplex *work, magma_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_int_t lrwork,
    #endif
    magma_int_t *iwork, magma_int_t liwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zhetrd */
/* magma_zhetrd */
/* magma_zhetrd */
magma_int_t
magma_zhetrd(
    magma_uplo_t uplo, magma_int_t n,
    magmaDoubleComplex *A, magma_int_t lda,
    double *d, double *e, magmaDoubleComplex *tau,
    magmaDoubleComplex *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zhetrf */
/* magma_zhetrf */
/* magma_zhetrf */
magma_int_t
magma_zhetrf(
    magma_uplo_t uplo, magma_int_t n,
    magmaDoubleComplex *A, magma_int_t lda,
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

/* magma_zlahef_gpu */
/* magma_zlahef_gpu */
/* magma_zlahef_gpu */
magma_int_t
magma_zlahef_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nb, magma_int_t *kb,
    magmaDoubleComplex *hA, magma_int_t lda,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_int_t *ipiv,
    magmaDoubleComplex_ptr dW, size_t dW_offset, magma_int_t lddw,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zlahr2 */
/* magma_zlahr2 */
/* magma_zlahr2 */
magma_int_t
magma_zlahr2(
    magma_int_t n, magma_int_t k, magma_int_t nb,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr dV, size_t dV_offset, magma_int_t lddv,
    magmaDoubleComplex *A,  magma_int_t lda,
    magmaDoubleComplex *tau,
    magmaDoubleComplex *T,  magma_int_t ldt,
    magmaDoubleComplex *Y,  magma_int_t ldy,
    magma_queue_t queue);

/* magma_zlahru */
/* magma_zlahru */
/* magma_zlahru */
magma_int_t
magma_zlahru(
    magma_int_t n, magma_int_t ihi, magma_int_t k, magma_int_t nb,
    magmaDoubleComplex     *A, magma_int_t lda,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr dY, size_t dY_offset, magma_int_t lddy,
    magmaDoubleComplex_ptr dV, size_t dV_offset, magma_int_t lddv,
    magmaDoubleComplex_ptr dT, size_t dT_offset,
    magmaDoubleComplex_ptr dwork, size_t dwork_offset,
    magma_queue_t queue);

/* magma_zlatrd */
/* magma_zlatrd */
/* magma_zlatrd */
magma_int_t
magma_zlatrd(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nb,
    magmaDoubleComplex *A,  magma_int_t lda,
    double *e, magmaDoubleComplex *tau,
    magmaDoubleComplex *W,  magma_int_t ldw,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr dW, size_t dW_offset, magma_int_t lddw,
    magma_queue_t queue);

/* magma_zposv */
/* magma_zposv */
/* magma_zposv */
magma_int_t
magma_zposv(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    magmaDoubleComplex *A, magma_int_t lda,
    magmaDoubleComplex *B, magma_int_t ldb,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_zpotrf */
/* magma_zpotrf */
/* magma_zpotrf */
magma_int_t
magma_zpotrf(
    magma_uplo_t uplo, magma_int_t n,
    magmaDoubleComplex *A, magma_int_t lda,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_zstedx */
/* magma_zstedx */
/* magma_zstedx */
magma_int_t
magma_zstedx(
    magma_range_t range, magma_int_t n, double vl, double vu,
    magma_int_t il, magma_int_t iu, double *d, double *e,
    magmaDoubleComplex *Z, magma_int_t ldz,
    double *rwork, magma_int_t lrwork,
    magma_int_t *iwork, magma_int_t liwork,
    magmaDouble_ptr dwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zunghr */
/* magma_zunghr */
/* magma_zunghr */
magma_int_t
magma_zunghr(
    magma_int_t n, magma_int_t ilo, magma_int_t ihi,
    magmaDoubleComplex *A, magma_int_t lda,
    magmaDoubleComplex *tau,
    magmaDoubleComplex_ptr dT, size_t dT_offset, magma_int_t nb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zungqr */
/* magma_zungqr */
/* magma_zungqr */
magma_int_t
magma_zungqr(
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaDoubleComplex *A, magma_int_t lda,
    magmaDoubleComplex *tau,
    magmaDoubleComplex_ptr dT, size_t dT_offset, magma_int_t nb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zunmql */
/* magma_zunmql */
/* magma_zunmql */
magma_int_t
magma_zunmql(
    magma_side_t side, magma_trans_t trans,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaDoubleComplex *A, magma_int_t lda,
    magmaDoubleComplex *tau,
    magmaDoubleComplex *C, magma_int_t ldc,
    magmaDoubleComplex *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zunmqr */
/* magma_zunmqr */
/* magma_zunmqr */
magma_int_t
magma_zunmqr(
    magma_side_t side, magma_trans_t trans,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaDoubleComplex *A, magma_int_t lda,
    magmaDoubleComplex *tau,
    magmaDoubleComplex *C, magma_int_t ldc,
    magmaDoubleComplex *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zunmtr */
/* magma_zunmtr */
/* magma_zunmtr */
magma_int_t
magma_zunmtr(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex *A,    magma_int_t lda,
    magmaDoubleComplex *tau,
    magmaDoubleComplex *C,    magma_int_t ldc,
    magmaDoubleComplex *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on GPU (alphabetical order)
*/

/* magma_zgels_gpu */
/* magma_zgels_gpu */
/* magma_zgels_gpu */
magma_int_t
magma_zgels_gpu(
    magma_trans_t trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaDoubleComplex *hwork, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zgeqr2x3_gpu */
/* magma_zgeqr2x3_gpu */
/* magma_zgeqr2x3_gpu */
magma_int_t
magma_zgeqr2x3_gpu(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_ptr dA,    size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr dtau,  size_t dtau_offset,
    magmaDoubleComplex_ptr dT,    size_t dT_offset,
    magmaDoubleComplex_ptr ddA,   size_t ddA_offset,
    magmaDouble_ptr        dwork, size_t dwork_offset,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zgeqrf_gpu */
/* magma_zgeqrf_gpu */
/* magma_zgeqrf_gpu */
magma_int_t
magma_zgeqrf_gpu(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_ptr dA, size_t dA_offset,  magma_int_t ldda,
    magmaDoubleComplex *tau,
    magmaDoubleComplex_ptr dT, size_t dT_offset,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zgeqrf_msub */
/* magma_zgeqrf_msub */
/* magma_zgeqrf_msub */
magma_int_t
magma_zgeqrf_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_ptr dlA[], magma_int_t ldda,
    magmaDoubleComplex *tau,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_zgeqrf2_2q_gpu */
/* magma_zgeqrf2_2q_gpu */
/* magma_zgeqrf2_2q_gpu */
magma_int_t
magma_zgeqrf2_2q_gpu(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex *tau,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_zgeqrf2_gpu */
/* magma_zgeqrf2_gpu */
/* magma_zgeqrf2_gpu */
magma_int_t
magma_zgeqrf2_gpu(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex *tau,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_zgeqrf2_mgpu */
/* magma_zgeqrf2_mgpu */
/* magma_zgeqrf2_mgpu */
magma_int_t
magma_zgeqrf2_mgpu(
    magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_ptr dlA[], magma_int_t ldda,
    magmaDoubleComplex *tau,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_zgeqrs_gpu */
/* magma_zgeqrs_gpu */
/* magma_zgeqrs_gpu */
magma_int_t
magma_zgeqrs_gpu(
    magma_int_t m, magma_int_t n, magma_int_t nrhs,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex *tau,
    magmaDoubleComplex_ptr dT, size_t dT_offset,
    magmaDoubleComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaDoubleComplex *hwork, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zgesv_gpu */
/* magma_zgesv_gpu */
/* magma_zgesv_gpu */
magma_int_t
magma_zgesv_gpu(
    magma_int_t n, magma_int_t nrhs,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magmaDoubleComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zgetrf_gpu */
/* magma_zgetrf_gpu */
/* magma_zgetrf_gpu */
magma_int_t
magma_zgetrf_gpu(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zgetrf_mgpu */
/* magma_zgetrf_mgpu */
/* magma_zgetrf_mgpu */
magma_int_t
magma_zgetrf_mgpu(
    magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_ptr d_lA[], size_t dlA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_zgetrf_msub */
/* magma_zgetrf_msub */
/* magma_zgetrf_msub */
magma_int_t
magma_zgetrf_msub(
    magma_trans_t trans, magma_int_t num_subs, magma_int_t ngpu,
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_ptr d_lA[], size_t dlA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_zgetrf2_gpu */
/* magma_zgetrf2_gpu */
/* magma_zgetrf2_gpu */
magma_int_t
magma_zgetrf2_gpu(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_zgetrf2_mgpu */
/* magma_zgetrf2_mgpu */
/* magma_zgetrf2_mgpu */
magma_int_t
magma_zgetrf2_mgpu(
    magma_int_t ngpu,
    magma_int_t m, magma_int_t n, magma_int_t nb, magma_int_t offset,
    magmaDoubleComplex_ptr d_lAT[], size_t dlAT_offset, magma_int_t lddat, magma_int_t *ipiv,
    magmaDoubleComplex_ptr d_lAP[], size_t dlAP_offset,
    magmaDoubleComplex *W, magma_int_t ldw,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_zgetrf2_msub */
/* magma_zgetrf2_msub */
/* magma_zgetrf2_msub */
magma_int_t
magma_zgetrf2_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_int_t m, magma_int_t n, magma_int_t nb, magma_int_t offset,
    magmaDoubleComplex_ptr d_lAT[], size_t dlAT_offset, magma_int_t lddat, magma_int_t *ipiv,
    magmaDoubleComplex_ptr d_panel[],
    magmaDoubleComplex_ptr d_lAP[], size_t dlAP_offset,
    magmaDoubleComplex *W, magma_int_t ldw,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_zgetri_gpu */
/* magma_zgetri_gpu */
/* magma_zgetri_gpu */
magma_int_t
magma_zgetri_gpu(
    magma_int_t n,
    magmaDoubleComplex_ptr dA,    size_t dA_offset,    magma_int_t ldda, magma_int_t *ipiv,
    magmaDoubleComplex_ptr dwork, size_t dwork_offset, magma_int_t lwork,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_zgetrs_gpu */
/* magma_zgetrs_gpu */
/* magma_zgetrs_gpu */
magma_int_t
magma_zgetrs_gpu(
    magma_trans_t trans, magma_int_t n, magma_int_t nrhs,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda, magma_int_t *ipiv,
    magmaDoubleComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zlabrd_gpu */
/* magma_zlabrd_gpu */
/* magma_zlabrd_gpu */
magma_int_t
magma_zlabrd_gpu(
    magma_int_t m, magma_int_t n, magma_int_t nb,
    magmaDoubleComplex     *A,                   magma_int_t lda,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    double *d, double *e, magmaDoubleComplex *tauq, magmaDoubleComplex *taup,
    magmaDoubleComplex     *X,                   magma_int_t ldx,
    magmaDoubleComplex_ptr dX, size_t dX_offset, magma_int_t lddx,
    magmaDoubleComplex     *Y,                   magma_int_t ldy,
    magmaDoubleComplex_ptr dY, size_t dY_offset, magma_int_t lddy,
    magma_queue_t queue);

/* magma_zlarfb_gpu */
/* magma_zlarfb_gpu */
/* magma_zlarfb_gpu */
magma_int_t
magma_zlarfb_gpu(
    magma_side_t side, magma_trans_t trans, magma_direct_t direct, magma_storev_t storev,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaDoubleComplex_const_ptr dV, size_t dV_offset,    magma_int_t lddv,
    magmaDoubleComplex_const_ptr dT, size_t dT_offset,    magma_int_t lddt,
    magmaDoubleComplex_ptr dC,       size_t dC_offset,    magma_int_t lddc,
    magmaDoubleComplex_ptr dwork,    size_t dwork_offset, magma_int_t ldwork,
    magma_queue_t queue);

/* magma_zlauum_gpu */
/* magma_zlauum_gpu */
/* magma_zlauum_gpu */
magma_int_t
magma_zlauum_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zposv_gpu */
/* magma_zposv_gpu */
/* magma_zposv_gpu */
magma_int_t
magma_zposv_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zpotrf_gpu */
/* magma_zpotrf_gpu */
/* magma_zpotrf_gpu */
magma_int_t
magma_zpotrf_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_zpotrf_mgpu */
/* magma_zpotrf_mgpu */
/* magma_zpotrf_mgpu */
magma_int_t
magma_zpotrf_mgpu(
    magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t n,
    magmaDoubleComplex_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_zpotrf_msub */
/* magma_zpotrf_msub */
/* magma_zpotrf_msub */
magma_int_t
magma_zpotrf_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t n,
    magmaDoubleComplex_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_zpotrf2_gpu */
/* magma_zpotrf2_gpu */
/* magma_zpotrf2_gpu */
magma_int_t
magma_zpotrf2_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_zpotrf2_mgpu */
/* magma_zpotrf2_mgpu */
/* magma_zpotrf2_mgpu */
magma_int_t
magma_zpotrf2_mgpu(
    magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
    magmaDoubleComplex_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr d_lP[],                   magma_int_t lddp,
    magmaDoubleComplex *A, magma_int_t lda, magma_int_t h,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_zpotrf2_msub */
/* magma_zpotrf2_msub */
/* magma_zpotrf2_msub */
magma_int_t
magma_zpotrf2_msub(
    magma_int_t num_subs, magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magma_int_t off_i, magma_int_t off_j, magma_int_t nb,
    magmaDoubleComplex_ptr d_lA[], size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr d_lP[],                   magma_int_t lddp,
    magmaDoubleComplex *A, magma_int_t lda, magma_int_t h,
    magma_queue_t queues[],
    magma_int_t *info);

/* magma_zpotri_gpu */
/* magma_zpotri_gpu */
/* magma_zpotri_gpu */
magma_int_t
magma_zpotri_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_zpotrs_gpu */
/* magma_zpotrs_gpu */
/* magma_zpotrs_gpu */
magma_int_t
magma_zpotrs_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue,
    magma_int_t *info);

/* magma_ztrtri_gpu */
/* magma_ztrtri_gpu */
/* magma_ztrtri_gpu */
magma_int_t
magma_ztrtri_gpu(
    magma_uplo_t uplo, magma_diag_t diag, magma_int_t n,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queues[2],
    magma_int_t *info);

/* magma_zunmqr_gpu */
/* magma_zunmqr_gpu */
/* magma_zunmqr_gpu */
magma_int_t
magma_zunmqr_gpu(
    magma_side_t side, magma_trans_t trans,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex *tau,
    magmaDoubleComplex_ptr dC, size_t dC_offset, magma_int_t lddc,
    magmaDoubleComplex *hwork, magma_int_t lwork,
    magmaDoubleComplex_ptr dT, size_t dT_offset, magma_int_t nb,
    magma_queue_t queue,
    magma_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA utility function definitions
*/

extern const magmaDoubleComplex MAGMA_Z_NAN;
extern const magmaDoubleComplex MAGMA_Z_INF;

/* magma_znan_inf */
/* magma_znan_inf */
/* magma_znan_inf */
magma_int_t
magma_znan_inf(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    const magmaDoubleComplex *A, magma_int_t lda,
    magma_int_t *cnt_nan,
    magma_int_t *cnt_inf);

/* magma_znan_inf_gpu */
/* magma_znan_inf_gpu */
/* magma_znan_inf_gpu */
magma_int_t
magma_znan_inf_gpu(
    magma_uplo_t uplo, magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr dA, magma_int_t dA_offset, magma_int_t ldda,
    magma_int_t *cnt_nan,
    magma_int_t *cnt_inf,
    magma_queue_t queue);

void magma_zprint(
    magma_int_t m, magma_int_t n,
    const magmaDoubleComplex *A, magma_int_t lda);

void magma_zprint_gpu(
    magma_int_t m, magma_int_t n,
    magmaDoubleComplex_const_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue);

void zpanel_to_q(
    magma_uplo_t uplo, magma_int_t ib,
    magmaDoubleComplex *A, magma_int_t lda,
    magmaDoubleComplex *work);

void zq_to_panel(
    magma_uplo_t uplo, magma_int_t ib,
    magmaDoubleComplex *A, magma_int_t lda,
    magmaDoubleComplex *work);

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif /* MAGMA_Z_H */
