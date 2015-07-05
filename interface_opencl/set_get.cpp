/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014
 
       @author Mark Gates
       @precisions normal z -> s d c
*/

#include <stdlib.h>
#include <stdio.h>

#include "magma.h"
#include "error.h"

#if defined(HAVE_clBLAS)

// generic, type-independent routines to copy data.
// type-safe versions which avoid the user needing sizeof(...) are in [sdcz]set_get.cpp

// ========================================
// globals, defined in interface.c
extern magma_event_t* g_event;

// ========================================
// copying vectors
extern "C" void
magma_setvector(
    magma_int_t n, magma_int_t elemSize,
    void const* hx_src,                   magma_int_t incx,
    magma_ptr   dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueWriteBuffer(
            queue, dy_dst, CL_TRUE,
            dy_offset*elemSize, n*elemSize,
            hx_src, 0, NULL, g_event);
        check_error( err );
    }
    else {
        magma_int_t ldha = incx;
        magma_int_t lddb = incy;
        magma_setmatrix( 1, n, elemSize,
            hx_src,            ldha,
            dy_dst, dy_offset, lddb,
            queue);
    }
}

// --------------------
extern "C" void
magma_setvector_async(
    magma_int_t n, magma_int_t elemSize,
    void const* hx_src,                   magma_int_t incx,
    magma_ptr   dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueWriteBuffer(
            queue, dy_dst, CL_FALSE,
            dy_offset*elemSize, n*elemSize,
            hx_src, 0, NULL, event);
        check_error( err );
    }
    else {
        magma_int_t ldha = incx;
        magma_int_t lddb = incy;
        magma_setmatrix_async( 1, n, elemSize,
            hx_src,            ldha,
            dy_dst, dy_offset, lddb,
            queue, event);
    }
}

// --------------------
extern "C" void
magma_getvector(
    magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    void*           hy_dst,                   magma_int_t incy,
    magma_queue_t queue )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueReadBuffer(
            queue, dx_src, CL_TRUE,
            dx_offset*elemSize, n*elemSize,
            hy_dst, 0, NULL, g_event);
        check_error( err );
    }
    else {
        magma_int_t ldda = incx;
        magma_int_t ldhb = incy;
        magma_getmatrix( 1, n, elemSize,
            dx_src, dx_offset, ldda,
            hy_dst,            ldhb,
            queue);
    }
}

// --------------------
extern "C" void
magma_getvector_async(
    magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    void*           hy_dst,                   magma_int_t incy,
    magma_queue_t queue, magma_event_t *event )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueReadBuffer(
            queue, dx_src, CL_FALSE,
            dx_offset*elemSize, n*elemSize,
            hy_dst, 0, NULL, event);
        check_error( err );
    }
    else {
        magma_int_t ldda = incx;
        magma_int_t ldhb = incy;
        magma_getmatrix_async( 1, n, elemSize,
            dx_src, dx_offset, ldda,
            hy_dst,            ldhb,
            queue, event);
    }
}

// --------------------
extern "C" void
magma_copyvector(
    magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magma_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueReadBuffer(
            queue, dx_src, CL_TRUE,
            dx_offset*elemSize, n*elemSize,
            dy_dst, dy_offset*elemSize, NULL, g_event);
        check_error( err );
    }
    else {
        magma_int_t ldda = incx;
        magma_int_t lddb = incy;
        magma_copymatrix( 1, n, elemSize,
            dx_src, dx_offset, ldda,
            dy_dst, dy_offset, lddb,
            queue);
    }
}

// --------------------
extern "C" void
magma_copyvector_async(
    magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magma_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueReadBuffer(
            queue, dx_src, CL_FALSE,
            dx_offset*elemSize, n*elemSize,
            dy_dst, dy_offset*elemSize, NULL, event);
        check_error( err );
    }
    else {
        magma_int_t ldda = incx;
        magma_int_t lddb = incy;
        magma_copymatrix_async( 1, n, elemSize,
            dx_src, dx_offset, ldda,
            dy_dst, dy_offset, lddb,
            queue, event);
    }
}


// ========================================
// copying sub-matrices (contiguous columns)
// OpenCL takes queue even for blocking transfers, oddly.
extern "C" void
magma_setmatrix(
    magma_int_t m, magma_int_t n, magma_int_t elemSize,
    void const* hA_src,                   magma_int_t ldha,
    magma_ptr   dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue )
{
    if (m <= 0 || n <= 0)
        return;

    size_t buffer_origin[3] = { dB_offset*elemSize, 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*elemSize, n, 1 };
    cl_int err = clEnqueueWriteBufferRect(
        queue, dB_dst, CL_TRUE,  // blocking
        buffer_origin, host_orig, region,
        lddb*elemSize, 0,
        ldha*elemSize, 0,
        hA_src, 0, NULL, g_event );
    check_error( err );
}

// --------------------
extern "C" void
magma_setmatrix_async(
    magma_int_t m, magma_int_t n, magma_int_t elemSize,
    void const* hA_src,                   magma_int_t ldha,
    magma_ptr   dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue, magma_event_t *event )
{
    if (m <= 0 || n <= 0)
        return;

    size_t buffer_origin[3] = { dB_offset*elemSize, 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*elemSize, n, 1 };
    cl_int err = clEnqueueWriteBufferRect(
        queue, dB_dst, CL_FALSE,  // non-blocking
        buffer_origin, host_orig, region,
        lddb*elemSize, 0,
        ldha*elemSize, 0,
        hA_src, 0, NULL, event );
    clFlush( queue );
    check_error( err );
}

// --------------------
extern "C" void
magma_getmatrix(
    magma_int_t m, magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    void*           hB_dst,                   magma_int_t ldhb,
    magma_queue_t queue )
{
    if (m <= 0 || n <= 0)
       return;

    size_t buffer_origin[3] = { dA_offset*elemSize, 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*elemSize, n, 1 };
    cl_int err = clEnqueueReadBufferRect(
        queue, dA_src, CL_TRUE,  // blocking
        buffer_origin, host_orig, region,
        ldda*elemSize, 0,
        ldhb*elemSize, 0,
        hB_dst, 0, NULL, g_event );
    check_error( err );
}

// --------------------
extern "C" void
magma_getmatrix_async(
    magma_int_t m, magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    void*           hB_dst,                   magma_int_t ldhb,
    magma_queue_t queue, magma_event_t *event )
{
    if (m <= 0 || n <= 0)
        return;

    size_t buffer_origin[3] = { dA_offset*elemSize, 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*elemSize, n, 1 };
    cl_int err = clEnqueueReadBufferRect(
        queue, dA_src, CL_FALSE,  // non-blocking
        buffer_origin, host_orig, region,
        ldda*elemSize, 0,
        ldhb*elemSize, 0,
        hB_dst, 0, NULL, event );
    clFlush( queue );
    check_error( err );
}

// --------------------
extern "C" void
magma_copymatrix(
    magma_int_t m, magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magma_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue )
{
    if (m <= 0 || n <= 0)
        return;

    size_t src_origin[3] = { dA_offset*elemSize, 0, 0 };
    size_t dst_orig[3]   = { dB_offset*elemSize, 0, 0 };
    size_t region[3]     = { m*elemSize, n, 1 };
    cl_int err = clEnqueueCopyBufferRect(
        queue, dA_src, dB_dst,
        src_origin, dst_orig, region,
        ldda*elemSize, 0,
        lddb*elemSize, 0,
        0, NULL, g_event );
    check_error( err );
}

// --------------------
extern "C" void
magma_copymatrix_async(
    magma_int_t m, magma_int_t n, magma_int_t elemSize,
    magma_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magma_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue, magma_event_t *event )
{
    if (m <= 0 || n <= 0)
        return;

    // TODO how to make non-blocking?
    size_t src_origin[3] = { dA_offset*elemSize, 0, 0 };
    size_t dst_orig[3]   = { dB_offset*elemSize, 0, 0 };
    size_t region[3]     = { m*elemSize, n, 1 };
    cl_int err = clEnqueueCopyBufferRect(
        queue, dA_src, dB_dst,
        src_origin, dst_orig, region,
        ldda*elemSize, 0,
        lddb*elemSize, 0,
        0, NULL, event );
    check_error( err );
}

#endif // HAVE_clBLAS
