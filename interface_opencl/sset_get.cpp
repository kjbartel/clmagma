/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014
 
       @author Mark Gates
       @generated from zset_get.cpp normal z -> s, Sat Nov 15 00:21:38 2014
*/

#include <stdlib.h>
#include <stdio.h>

#include "magma.h"
#include "error.h"

#if defined(HAVE_clBLAS)


// ========================================
// globals, defined in interface.c
extern magma_event_t* g_event;

// ========================================
// copying vectors
extern "C" void
magma_ssetvector(
    magma_int_t n,
    float const* hx_src,                   magma_int_t incx,
    magmaFloat_ptr    dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueWriteBuffer(
            queue, dy_dst, CL_TRUE,
            dy_offset*sizeof(float), n*sizeof(float),
            hx_src, 0, NULL, g_event);
        check_error( err );
    }
    else {
        magma_int_t ldha = incx;
        magma_int_t lddb = incy;
        magma_ssetmatrix( 1, n,
            hx_src,            ldha,
            dy_dst, dy_offset, lddb,
            queue);
    }
}

// --------------------
extern "C" void
magma_ssetvector_async(
    magma_int_t n,
    float const* hx_src,                   magma_int_t incx,
    magmaFloat_ptr    dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueWriteBuffer(
            queue, dy_dst, CL_FALSE,
            dy_offset*sizeof(float), n*sizeof(float),
            hx_src, 0, NULL, event);
        check_error( err );
    }
    else {
        magma_int_t ldha = incx;
        magma_int_t lddb = incy;
        magma_ssetmatrix_async( 1, n,
            hx_src,            ldha,
            dy_dst, dy_offset, lddb,
            queue, event);
    }
}

// --------------------
extern "C" void
magma_sgetvector(
    magma_int_t n,
    magmaFloat_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    float*          hy_dst,                   magma_int_t incy,
    magma_queue_t queue )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueReadBuffer(
            queue, dx_src, CL_TRUE,
            dx_offset*sizeof(float), n*sizeof(float),
            hy_dst, 0, NULL, g_event);
        check_error( err );
    }
    else {
        magma_int_t ldda = incx;
        magma_int_t ldhb = incy;
        magma_sgetmatrix( 1, n,
            dx_src, dx_offset, ldda,
            hy_dst,            ldhb,
            queue);
    }
}

// --------------------
extern "C" void
magma_sgetvector_async(
    magma_int_t n,
    magmaFloat_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    float*          hy_dst,                   magma_int_t incy,
    magma_queue_t queue, magma_event_t *event )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueReadBuffer(
            queue, dx_src, CL_FALSE,
            dx_offset*sizeof(float), n*sizeof(float),
            hy_dst, 0, NULL, event);
        check_error( err );
    }
    else {
        magma_int_t ldda = incx;
        magma_int_t ldhb = incy;
        magma_sgetmatrix_async( 1, n,
            dx_src, dx_offset, ldda,
            hy_dst,            ldhb,
            queue, event);
    }
}

// --------------------
extern "C" void
magma_scopyvector(
    magma_int_t n,
    magmaFloat_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaFloat_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueReadBuffer(
            queue, dx_src, CL_TRUE,
            dx_offset*sizeof(float), n*sizeof(float),
            dy_dst, dy_offset*sizeof(float), NULL, g_event);
        check_error( err );
    }
    else {
        magma_int_t ldda = incx;
        magma_int_t lddb = incy;
        magma_scopymatrix( 1, n,
            dx_src, dx_offset, ldda,
            dy_dst, dy_offset, lddb,
            queue);
    }
}

// --------------------
extern "C" void
magma_scopyvector_async(
    magma_int_t n,
    magmaFloat_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaFloat_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueReadBuffer(
            queue, dx_src, CL_FALSE,
            dx_offset*sizeof(float), n*sizeof(float),
            dy_dst, dy_offset*sizeof(float), NULL, event);
        check_error( err );
    }
    else {
        magma_int_t ldda = incx;
        magma_int_t lddb = incy;
        magma_scopymatrix_async( 1, n,
            dx_src, dx_offset, ldda,
            dy_dst, dy_offset, lddb,
            queue, event);
    }
}


// ========================================
// copying sub-matrices (contiguous columns)
// OpenCL takes queue even for blocking transfers, oddly.
extern "C" void
magma_ssetmatrix(
    magma_int_t m, magma_int_t n,
    float const* hA_src,                   magma_int_t ldha,
    magmaFloat_ptr    dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue )
{
    if (m <= 0 || n <= 0)
        return;

    size_t buffer_origin[3] = { dB_offset*sizeof(float), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(float), n, 1 };
    cl_int err = clEnqueueWriteBufferRect(
        queue, dB_dst, CL_TRUE,  // blocking
        buffer_origin, host_orig, region,
        lddb*sizeof(float), 0,
        ldha*sizeof(float), 0,
        hA_src, 0, NULL, g_event );
    check_error( err );
}

// --------------------
extern "C" void
magma_ssetmatrix_async(
    magma_int_t m, magma_int_t n,
    float const* hA_src,                   magma_int_t ldha,
    magmaFloat_ptr    dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue, magma_event_t *event )
{
    if (m <= 0 || n <= 0)
        return;

    size_t buffer_origin[3] = { dB_offset*sizeof(float), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(float), n, 1 };
    cl_int err = clEnqueueWriteBufferRect(
        queue, dB_dst, CL_FALSE,  // non-blocking
        buffer_origin, host_orig, region,
        lddb*sizeof(float), 0,
        ldha*sizeof(float), 0,
        hA_src, 0, NULL, event );
    clFlush(queue);
    check_error( err );
}

// --------------------
extern "C" void
magma_sgetmatrix(
    magma_int_t m, magma_int_t n,
    magmaFloat_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    float*          hB_dst,                   magma_int_t ldhb,
    magma_queue_t queue )
{
    if (m <= 0 || n <= 0)
       return;

    size_t buffer_origin[3] = { dA_offset*sizeof(float), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(float), n, 1 };
    cl_int err = clEnqueueReadBufferRect(
        queue, dA_src, CL_TRUE,  // blocking
        buffer_origin, host_orig, region,
        ldda*sizeof(float), 0,
        ldhb*sizeof(float), 0,
        hB_dst, 0, NULL, g_event );
    check_error( err );
}

// --------------------
extern "C" void
magma_sgetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaFloat_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    float*          hB_dst,                   magma_int_t ldhb,
    magma_queue_t queue, magma_event_t *event )
{
    if (m <= 0 || n <= 0)
        return;

    size_t buffer_origin[3] = { dA_offset*sizeof(float), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(float), n, 1 };
    cl_int err = clEnqueueReadBufferRect(
        queue, dA_src, CL_FALSE,  // non-blocking
        buffer_origin, host_orig, region,
        ldda*sizeof(float), 0,
        ldhb*sizeof(float), 0,
        hB_dst, 0, NULL, event );
    clFlush(queue);
    check_error( err );
}

// --------------------
extern "C" void
magma_scopymatrix(
    magma_int_t m, magma_int_t n,
    magmaFloat_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue )
{
    if (m <= 0 || n <= 0)
        return;

    size_t src_origin[3] = { dA_offset*sizeof(float), 0, 0 };
    size_t dst_orig[3]   = { dB_offset*sizeof(float), 0, 0 };
    size_t region[3]     = { m*sizeof(float), n, 1 };
    cl_int err = clEnqueueCopyBufferRect(
        queue, dA_src, dB_dst,
        src_origin, dst_orig, region,
        ldda*sizeof(float), 0,
        lddb*sizeof(float), 0,
        0, NULL, g_event );
    check_error( err );
}

// --------------------
extern "C" void
magma_scopymatrix_async(
    magma_int_t m, magma_int_t n,
    magmaFloat_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaFloat_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue, magma_event_t *event )
{
    if (m <= 0 || n <= 0)
        return;

    // TODO how to make non-blocking?
    size_t src_origin[3] = { dA_offset*sizeof(float), 0, 0 };
    size_t dst_orig[3]   = { dB_offset*sizeof(float), 0, 0 };
    size_t region[3]     = { m*sizeof(float), n, 1 };
    cl_int err = clEnqueueCopyBufferRect(
        queue, dA_src, dB_dst,
        src_origin, dst_orig, region,
        ldda*sizeof(float), 0,
        lddb*sizeof(float), 0,
        0, NULL, event );
    check_error( err );
}

#endif // HAVE_clBLAS
