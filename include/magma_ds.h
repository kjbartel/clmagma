/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from magma_zc.h mixed zc -> ds, Sat Nov 15 00:21:34 2014
*/

#ifndef MAGMA_DS_H
#define MAGMA_DS_H

#include "magma_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Mixed precision */
magma_int_t
magma_dsgesv_gpu(
    magma_trans_t trans, magma_int_t n, magma_int_t nrhs,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magma_int_t *ipiv,
    magmaInt_ptr dipiv,
    magmaDouble_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaDouble_ptr dX, size_t dX_offset, magma_int_t lddx,
    magmaDouble_ptr dworkd, size_t dworkd_offset,
    magmaFloat_ptr dworks, size_t dworks_offset,
    magma_int_t *iter,
    magma_queue_t queue,
    magma_int_t *info );

magma_int_t
magma_dsgetrs_gpu(
    magma_trans_t trans, magma_int_t n, magma_int_t nrhs,
    magmaFloat_ptr  dA, size_t dA_offset, magma_int_t ldda,
    magmaInt_ptr        dipiv, size_t dipiv_offset,
    magmaDouble_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaDouble_ptr dX, size_t dX_offset, magma_int_t lddx,
    magmaFloat_ptr dSX, size_t dSX_offset,
    magma_queue_t queue,
    magma_int_t *info );

magma_int_t
magma_dsposv_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nrhs,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaDouble_ptr dX, size_t dX_offset, magma_int_t lddx,
    magmaDouble_ptr dworkd, size_t dworkd_offset, 
    magmaFloat_ptr  dworks, size_t dworks_offset,
    magma_int_t *iter,
    magma_queue_t queue,
    magma_int_t *info );

magma_int_t
magma_dsgeqrsv_gpu(
    magma_int_t m, magma_int_t n, magma_int_t nrhs,
    magmaDouble_ptr dA, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr dB, size_t dB_offset, magma_int_t lddb,
    magmaDouble_ptr dX, size_t dX_offset, magma_int_t lddx,
    magma_int_t *iter,
    magma_queue_t queue,
    magma_int_t *info );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_DS_H */
