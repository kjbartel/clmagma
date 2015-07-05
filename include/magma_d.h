/*
 *   -- clMAGMA (version 0.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @generated d Wed Apr  4 01:12:51 2012
 */

#ifndef MAGMA_D_H
#define MAGMA_D_H

#include "magma_types.h"

#define PRECISION_d

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
magma_dgeqrf2_gpu(
        magma_int_t m, magma_int_t n, 
        magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda, 
        double *tau, magma_int_t *info, magma_queue_t queue);

magma_err_t
magma_dgetrf_gpu(
        magma_int_t m, magma_int_t n, 
        magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda, 
        magma_int_t *ipiv, magma_int_t *info, magma_queue_t queue);

magma_err_t
magma_dlarfb_gpu( 
        int side, int trans, int direct, int storev, 
        magma_int_t m, magma_int_t n, magma_int_t k,
        magmaDouble_ptr dV, size_t dV_offset, magma_int_t ldv,
        magmaDouble_ptr dT, size_t dT_offset, magma_int_t ldt, 
        magmaDouble_ptr dC, size_t dC_offset, magma_int_t ldc,
        magmaDouble_ptr dwork, size_t dwork_offset, magma_int_t ldwork,
	magma_queue_t queue);

magma_err_t
magma_dpotrf_gpu(
        int uplo,
        magma_int_t n, 
        magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda, 
	magma_int_t *info, magma_queue_t queue);

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA utility function definitions
*/

void magma_dprint    ( magma_int_t m, magma_int_t n, double  *A, magma_int_t lda  );
void magma_dprint_gpu( magma_int_t m, magma_int_t n, double *dA, magma_int_t ldda );

#ifdef __cplusplus
}
#endif

#undef PRECISION_d
#endif /* MAGMA_D_H */
