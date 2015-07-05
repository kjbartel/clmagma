/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014
*/

#ifndef MAGMABLAS_H
#define MAGMABLAS_H

#include "magmablas_z.h"
#include "magmablas_c.h"
#include "magmablas_d.h"
#include "magmablas_s.h"
#include "magmablas_zc.h"
#include "magmablas_ds.h"

#ifdef __cplusplus
extern "C" {
#endif


// ========================================
// copying vectors
// set copies host to device
// get copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)

void
magma_setvector(
    magma_int_t n, magma_int_t elemSize,
    void const*     hx_src,                   magma_int_t incx,
    magma_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_setvector_async(
    magma_int_t n, magma_int_t elemSize,
    void const*     hx_src,                   magma_int_t incx,
    magma_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );

void
magma_getvector(
    magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    void*           hy_dst,                   magma_int_t incy,
    magma_queue_t queue );

void
magma_getvector_async(
    magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    void*           hy_dst,                   magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );

void
magma_copyvector(
    magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magma_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue );

void
magma_copyvector_async(
    magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magma_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)

void
magma_setmatrix(
    magma_int_t m, magma_int_t n, magma_int_t elemSize,
    void const*     hA_src,                   magma_int_t ldha,
    magma_ptr       dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue );

void
magma_setmatrix_async(
    magma_int_t m, magma_int_t n, magma_int_t elemSize,
    void const*     hA_src,                   magma_int_t ldha,
    magma_ptr       dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue, magma_event_t *event );

void
magma_getmatrix(
    magma_int_t m, magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    void*           hA_dst,                   magma_int_t ldha,
    magma_queue_t queue );

void
magma_getmatrix_async(
    magma_int_t m, magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    void*           hA_dst,                   magma_int_t ldha,
    magma_queue_t queue, magma_event_t *event );

void
magma_copymatrix(
    magma_int_t m, magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magma_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue );

void
magma_copymatrix_async(
    magma_int_t m, magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magma_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue, magma_event_t *event );


#ifdef __cplusplus
}
#endif

#endif // MAGMABLAS_H
