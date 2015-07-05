/*
 *   -- clMAGMA (version 0.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @generated s Wed Apr  4 01:12:51 2012
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

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA utility function definitions
*/

void magma_sprint    ( magma_int_t m, magma_int_t n, float  *A, magma_int_t lda  );
void magma_sprint_gpu( magma_int_t m, magma_int_t n, float *dA, magma_int_t ldda );

#ifdef __cplusplus
}
#endif

#undef PRECISION_s
#endif /* MAGMA_S_H */
