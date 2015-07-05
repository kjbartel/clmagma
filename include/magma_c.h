/*
 *   -- clMAGMA (version 0.2.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @generated c Thu May 24 17:09:40 2012
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

extern "C" magma_err_t
magma_cpotrs_gpu(
		magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
        magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
		magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
		magma_err_t *info, magma_queue_t queue );

extern "C" magma_err_t
magma_cposv_gpu( 
		magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
        magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
		magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
		magma_err_t *info, magma_queue_t queue );

extern "C" magma_err_t
magma_cgetrs_gpu(magma_trans_t trans, magma_int_t n, magma_int_t nrhs, 
		magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
		magma_int_t *ipiv, 
		magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb, 
		magma_int_t *info, magma_queue_t queue);

extern "C" magma_err_t
magma_cgesv_gpu( magma_int_t n, magma_int_t nrhs,
                 magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
				 magma_int_t *ipiv,
				 magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
				 magma_err_t *info, magma_queue_t queue );

extern "C" magma_int_t
magma_cunmqr_gpu(magma_side_t side, magma_trans_t trans,
                 magma_int_t m, magma_int_t n, magma_int_t k,
                 magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda, 
                 magmaFloatComplex *tau,
                 magmaFloatComplex_ptr dC, size_t dC_offset, magma_int_t lddc,
                 magmaFloatComplex *hwork, magma_int_t lwork,
                 magmaFloatComplex_ptr dT, size_t dT_offset, magma_int_t nb, 
                 magma_int_t *info, magma_queue_t queue);

extern "C" magma_err_t
magma_cgeqrs_gpu(magma_int_t m, magma_int_t n, magma_int_t nrhs,
                 magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda, 
                 magmaFloatComplex *tau,   magmaFloatComplex_ptr dT, size_t dT_offset, 
				 magmaFloatComplex_ptr dB, size_t dB_offset, magma_int_t lddb, 
                 magmaFloatComplex *hwork, magma_int_t lwork, 
                 magma_int_t *info, magma_queue_t queue);

extern "C" magma_err_t
magma_cgeqrf_gpu( magma_int_t m, magma_int_t n, 
                  magmaFloatComplex_ptr dA, size_t dA_offset,  magma_int_t ldda,
                  magmaFloatComplex *tau, magmaFloatComplex_ptr dT, size_t dT_offset, 
                  magma_int_t *info, magma_queue_t queue);

extern "C" magma_int_t
magma_cgels_gpu( magma_trans_t trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
                 magmaFloatComplex_ptr dA, size_t dA_offset,  magma_int_t ldda, 
                 magmaFloatComplex_ptr dB, size_t dB_offset,  magma_int_t lddb, 
                 magmaFloatComplex *hwork, magma_int_t lwork, 
                 magma_int_t *info, magma_queue_t queue );
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
