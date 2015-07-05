/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @precisions mixed zc -> ds
*/

#ifndef MAGMA_ZC_H
#define MAGMA_ZC_H

#include "magma_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Mixed precision */
magma_int_t
magma_zcgesv_gpu(
    magma_trans_t trans, magma_int_t n, magma_int_t nrhs,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_int_t *ipiv,
    magmaInt_ptr dipiv,
    magmaDoubleComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaDoubleComplex_ptr dX, size_t dX_offset, magma_int_t lddx,
    magmaDoubleComplex_ptr dworkd, size_t dworkd_offset,
    magmaFloatComplex_ptr dworks, size_t dworks_offset,
    magma_int_t *iter,
    magma_queue_t queue,
    magma_int_t *info );

magma_int_t
magma_zcgetrs_gpu(
    magma_trans_t trans, magma_int_t n, magma_int_t nrhs,
    magmaFloatComplex_ptr  dA, size_t dA_offset, magma_int_t ldda,
    magmaInt_ptr        dipiv, size_t dipiv_offset,
    magmaDoubleComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaDoubleComplex_ptr dX, size_t dX_offset, magma_int_t lddx,
    magmaFloatComplex_ptr dSX, size_t dSX_offset,
    magma_queue_t queue,
    magma_int_t *info );

magma_int_t
magma_zcposv_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaDoubleComplex_ptr dX, size_t dX_offset, magma_int_t lddx,
    magmaDoubleComplex_ptr dworkd, size_t dworkd_offset, 
    magmaFloatComplex_ptr  dworks, size_t dworks_offset,
    magma_int_t *iter,
    magma_queue_t queue,
    magma_int_t *info );

magma_int_t
magma_zcgeqrsv_gpu(
    magma_int_t m, magma_int_t n, magma_int_t nrhs,
    magmaDoubleComplex_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDoubleComplex_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaDoubleComplex_ptr dX, size_t dX_offset, magma_int_t lddx,
    magma_int_t *iter,
    magma_queue_t queue,
    magma_int_t *info );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_ZC_H */
