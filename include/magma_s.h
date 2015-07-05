/*
 *   -- clMAGMA (version 0.2.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @generated s Thu May 24 17:09:39 2012
 */

#ifndef MAGMA_S_H
#define MAGMA_S_H

#include "magma_types.h"

#define PRECISION_s

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
magma_sgeqrf2_gpu(
        magma_int_t m, magma_int_t n, 
        magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda, 
        float *tau, magma_int_t *info, magma_queue_t queue);

magma_err_t
magma_sgetrf_gpu(
        magma_int_t m, magma_int_t n, 
        magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda, 
        magma_int_t *ipiv, magma_int_t *info, magma_queue_t queue);

magma_err_t
magma_slarfb_gpu( 
        int side, int trans, int direct, int storev, 
        magma_int_t m, magma_int_t n, magma_int_t k,
        magmaFloat_ptr dV, size_t dV_offset, magma_int_t ldv,
        magmaFloat_ptr dT, size_t dT_offset, magma_int_t ldt, 
        magmaFloat_ptr dC, size_t dC_offset, magma_int_t ldc,
        magmaFloat_ptr dwork, size_t dwork_offset, magma_int_t ldwork,
        magma_queue_t queue);

magma_err_t
magma_spotrf_gpu(
        int uplo,
        magma_int_t n, 
        magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda, 
        magma_int_t *info, magma_queue_t queue);

extern "C" magma_err_t
magma_spotrs_gpu(
		magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
        magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
		magmaFloat_ptr dB, size_t dB_offset, magma_int_t lddb,
		magma_err_t *info, magma_queue_t queue );

extern "C" magma_err_t
magma_sposv_gpu( 
		magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
        magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
		magmaFloat_ptr dB, size_t dB_offset, magma_int_t lddb,
		magma_err_t *info, magma_queue_t queue );

extern "C" magma_err_t
magma_sgetrs_gpu(magma_trans_t trans, magma_int_t n, magma_int_t nrhs, 
		magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
		magma_int_t *ipiv, 
		magmaFloat_ptr dB, size_t dB_offset, magma_int_t lddb, 
		magma_int_t *info, magma_queue_t queue);

extern "C" magma_err_t
magma_sgesv_gpu( magma_int_t n, magma_int_t nrhs,
                 magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda,
				 magma_int_t *ipiv,
				 magmaFloat_ptr dB, size_t dB_offset, magma_int_t lddb,
				 magma_err_t *info, magma_queue_t queue );

extern "C" magma_int_t
magma_sormqr_gpu(magma_side_t side, magma_trans_t trans,
                 magma_int_t m, magma_int_t n, magma_int_t k,
                 magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda, 
                 float *tau,
                 magmaFloat_ptr dC, size_t dC_offset, magma_int_t lddc,
                 float *hwork, magma_int_t lwork,
                 magmaFloat_ptr dT, size_t dT_offset, magma_int_t nb, 
                 magma_int_t *info, magma_queue_t queue);

extern "C" magma_err_t
magma_sgeqrs_gpu(magma_int_t m, magma_int_t n, magma_int_t nrhs,
                 magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda, 
                 float *tau,   magmaFloat_ptr dT, size_t dT_offset, 
				 magmaFloat_ptr dB, size_t dB_offset, magma_int_t lddb, 
                 float *hwork, magma_int_t lwork, 
                 magma_int_t *info, magma_queue_t queue);

extern "C" magma_err_t
magma_sgeqrf_gpu( magma_int_t m, magma_int_t n, 
                  magmaFloat_ptr dA, size_t dA_offset,  magma_int_t ldda,
                  float *tau, magmaFloat_ptr dT, size_t dT_offset, 
                  magma_int_t *info, magma_queue_t queue);

extern "C" magma_int_t
magma_sgels_gpu( magma_trans_t trans, magma_int_t m, magma_int_t n, magma_int_t nrhs,
                 magmaFloat_ptr dA, size_t dA_offset,  magma_int_t ldda, 
                 magmaFloat_ptr dB, size_t dB_offset,  magma_int_t lddb, 
                 float *hwork, magma_int_t lwork, 
                 magma_int_t *info, magma_queue_t queue );
/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA utility function definitions
*/

void magma_sprint    ( magma_int_t m, magma_int_t n, float     *A, magma_int_t lda  );
void magma_sprint_gpu( magma_int_t m, magma_int_t n, magmaFloat_ptr dA, magma_int_t ldda );

#ifdef __cplusplus
}
#endif

#undef PRECISION_s
#endif /* MAGMA_S_H */
