/*
 *   -- clMAGMA (version 0.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @generated c Wed Jun 27 23:49:46 2012
 */

#ifndef MAGMA_C_H
#define MAGMA_C_H

#include "magma_types.h"

#define PRECISION_c

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/


/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/

magma_err_t
magma_cgebrd(magma_int_t m, magma_int_t n,
             magmaFloatComplex *a, magma_int_t lda, float *d, float *e,
             magmaFloatComplex *tauq, magmaFloatComplex *taup,
             magmaFloatComplex *work, magma_int_t lwork,
             magma_int_t *info, magma_queue_t queue);

magma_err_t
magma_cgeqrf2_gpu(
        magma_int_t m, magma_int_t n, 
        magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda, 
        magmaFloatComplex *tau, magma_int_t *info, magma_queue_t queue);

magma_err_t
magma_cgetrf_gpu(
        magma_int_t m, magma_int_t n, 
        magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda, 
        magma_int_t *ipiv, magma_int_t *info, magma_queue_t queue);

magma_err_t
magma_clarfb_gpu( 
        int side, int trans, int direct, int storev, 
        magma_int_t m, magma_int_t n, magma_int_t k,
        magmaFloatComplex_ptr dV, size_t dV_offset, magma_int_t ldv,
        magmaFloatComplex_ptr dT, size_t dT_offset, magma_int_t ldt, 
        magmaFloatComplex_ptr dC, size_t dC_offset, magma_int_t ldc,
        magmaFloatComplex_ptr dwork, size_t dwork_offset, magma_int_t ldwork,
        magma_queue_t queue);

magma_err_t
magma_cpotrf_gpu(
        int uplo,
        magma_int_t n, 
        magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda, 
        magma_int_t *info, magma_queue_t queue);

magma_err_t
magma_cpotrs_gpu(
		magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
        magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
		magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
		magma_err_t *info, magma_queue_t queue );

magma_err_t
magma_cposv_gpu( 
		magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
        magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
		magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
		magma_err_t *info, magma_queue_t queue );

magma_err_t
magma_cgetrs_gpu(magma_trans_t trans, magma_int_t n, magma_int_t nrhs, 
		magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
		magma_int_t *ipiv, 
		magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb, 
		magma_int_t *info, magma_queue_t queue);

magma_err_t
magma_cgesv_gpu( magma_int_t n, magma_int_t nrhs,
                 magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
				 magma_int_t *ipiv,
				 magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
				 magma_err_t *info, magma_queue_t queue );

magma_int_t
magma_cunmqr_gpu(magma_side_t side, magma_trans_t trans,
                 magma_int_t m, magma_int_t n, magma_int_t k,
                 magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda, 
                 magmaFloatComplex *tau,
                 magmaFloatComplex_ptr dC, size_t dC_offset, magma_int_t lddc,
                 magmaFloatComplex *hwork, magma_int_t lwork,
                 magmaFloatComplex_ptr dT, size_t dT_offset, magma_int_t nb, 
                 magma_int_t *info, magma_queue_t queue);

magma_err_t
magma_cgeqrs_gpu(magma_int_t m, magma_int_t n, magma_int_t nrhs,
                 magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda, 
                 magmaFloatComplex *tau,   magmaFloatComplex_ptr dT, size_t dT_offset, 
				 magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb, 
                 magmaFloatComplex *hwork, magma_int_t lwork, 
                 magma_int_t *info, magma_queue_t queue);

magma_err_t
magma_cgeqrf_gpu( magma_int_t m, magma_int_t n, 
                  magmaFloatComplex_ptr dA, size_t dA_offset,  magma_int_t ldda,
                  magmaFloatComplex *tau, magmaFloatComplex_ptr dT, size_t dT_offset, 
                  magma_int_t *info, magma_queue_t queue);

magma_int_t
magma_cgels_gpu( magma_trans_t trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
                 magmaFloatComplex_ptr dA, size_t dA_offset,  magma_int_t ldda, 
                 magmaFloatComplex_ptr dB, size_t dB_offset,  magma_int_t lddb, 
                 magmaFloatComplex *hwork, magma_int_t lwork, 
                 magma_int_t *info, magma_queue_t queue );

magma_int_t
magma_cgehrd(	magma_int_t n, magma_int_t ilo, magma_int_t ihi, 
		magmaFloatComplex *a, magma_int_t lda, 
		magmaFloatComplex *tau, 
		magmaFloatComplex *work, magma_int_t lwork, 
		magmaFloatComplex_ptr dT, size_t dT_offset, 
		magma_int_t *info, magma_queue_t queue );

magma_int_t
magma_clabrd_gpu( magma_int_t m, magma_int_t n, magma_int_t nb,
                  magmaFloatComplex *a, magma_int_t lda,
                  magmaFloatComplex_ptr da, size_t da_offset, magma_int_t ldda,
                  float *d, float *e, magmaFloatComplex *tauq, magmaFloatComplex *taup,
                  magmaFloatComplex *x, magma_int_t ldx,
                  magmaFloatComplex_ptr dx, size_t dx_offset, magma_int_t lddx,
                  magmaFloatComplex *y, magma_int_t ldy,
                  magmaFloatComplex_ptr dy, size_t dy_offset, magma_int_t lddy,
                  magma_queue_t queue );

magma_err_t
magma_clahr2(	magma_int_t n, magma_int_t k, magma_int_t nb,
		magmaFloatComplex_ptr da, size_t da_offset, magmaFloatComplex_ptr dv, size_t dv_offset, 
		magmaFloatComplex *a, magma_int_t lda, 
		magmaFloatComplex *tau, magmaFloatComplex *t, magma_int_t ldt, 
		magmaFloatComplex *y, magma_int_t ldy, 
		magma_queue_t queue	);

magma_err_t
magma_clahru(	magma_int_t n, magma_int_t ihi, magma_int_t k, magma_int_t nb, 
		magmaFloatComplex *a, magma_int_t lda, 
		magmaFloatComplex_ptr d_a, size_t d_a_offset, magmaFloatComplex_ptr y, size_t y_offset, 
		magmaFloatComplex_ptr v, size_t v_offset, magmaFloatComplex_ptr d_t, size_t dt_offset, 
		magmaFloatComplex_ptr d_work, size_t d_work_offset, magma_queue_t queue );

magma_err_t
magma_cunghr(	magma_int_t n, magma_int_t ilo, magma_int_t ihi, 
		magmaFloatComplex *a, magma_int_t lda, 
		magmaFloatComplex *tau, 
		magmaFloatComplex_ptr dT, size_t dT_offset, magma_int_t nb, 
		magma_int_t *info, magma_queue_t queue );

magma_err_t
magma_cungqr(	magma_int_t m, magma_int_t n, magma_int_t k,
		magmaFloatComplex *a, magma_int_t lda,
		magmaFloatComplex *tau, magmaFloatComplex_ptr dT, size_t dT_offset,
		magma_int_t nb, magma_int_t *info, magma_queue_t queue );

magma_err_t 
magma_clatrd(	char uplo, magma_int_t n, magma_int_t nb, 
				magmaFloatComplex *a,  magma_int_t lda, 
				float *e, magmaFloatComplex *tau, 
				magmaFloatComplex *w,  magma_int_t ldw, 
				magmaFloatComplex_ptr da, size_t da_offset, magma_int_t ldda, 
				magmaFloatComplex_ptr dw, size_t dw_offset, magma_int_t lddw, magma_queue_t queue);

 magma_err_t
 magma_chetrd(	char uplo, magma_int_t n,
				magmaFloatComplex *a, magma_int_t lda,
			    float *d, float *e, magmaFloatComplex *tau,
				magmaFloatComplex *work, magma_int_t lwork,
				magma_int_t *info, magma_queue_t queue);

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA utility function definitions
*/

void magma_cprint    ( magma_int_t m, magma_int_t n, magmaFloatComplex     *A, magma_int_t lda  );
void magma_cprint_gpu( magma_int_t m, magma_int_t n, magmaFloatComplex_ptr dA, magma_int_t ldda );

#ifdef __cplusplus
}
#endif

#undef PRECISION_c
#endif /* MAGMA_C_H */
