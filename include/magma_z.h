/*
 *   -- clMAGMA (version 0.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @precisions normal z -> s d c
 */

#ifndef MAGMA_Z_H
#define MAGMA_Z_H

#include "magma_types.h"

#define PRECISION_z

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
magma_zgeqrf2_gpu(
        magma_int_t m, magma_int_t n, 
        magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda, 
        magmaDoubleComplex *tau, magma_int_t *info, magma_queue_t queue);

magma_err_t
magma_zgetrf_gpu(
        magma_int_t m, magma_int_t n, 
        magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda, 
        magma_int_t *ipiv, magma_int_t *info, magma_queue_t queue);

magma_err_t
magma_zlarfb_gpu( 
        int side, int trans, int direct, int storev, 
        magma_int_t m, magma_int_t n, magma_int_t k,
        magmaDoubleComplex_ptr dV, size_t dV_offset, magma_int_t ldv,
        magmaDoubleComplex_ptr dT, size_t dT_offset, magma_int_t ldt, 
        magmaDoubleComplex_ptr dC, size_t dC_offset, magma_int_t ldc,
        magmaDoubleComplex_ptr dwork, size_t dwork_offset, magma_int_t ldwork,
	magma_queue_t queue);

magma_err_t
magma_zpotrf_gpu(
        int uplo,
        magma_int_t n, 
        magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda, 
	magma_int_t *info, magma_queue_t queue);

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA utility function definitions
*/

void magma_zprint    ( magma_int_t m, magma_int_t n, magmaDoubleComplex  *A, magma_int_t lda  );
void magma_zprint_gpu( magma_int_t m, magma_int_t n, magmaDoubleComplex *dA, magma_int_t ldda );

#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* MAGMA_Z_H */
