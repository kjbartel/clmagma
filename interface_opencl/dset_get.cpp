/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014
 
       @author Mark Gates
       @generated from zset_get.cpp normal z -> d, Sat Nov 15 00:21:38 2014
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
magma_dsetvector(
    magma_int_t n,
    double const* hx_src,                   magma_int_t incx,
    magmaDouble_ptr    dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueWriteBuffer(
            queue, dy_dst, CL_TRUE,
            dy_offset*sizeof(double), n*sizeof(double),
            hx_src, 0, NULL, g_event);
        check_error( err );
    }
    else {
        magma_int_t ldha = incx;
        magma_int_t lddb = incy;
        magma_dsetmatrix( 1, n,
            hx_src,            ldha,
            dy_dst, dy_offset, lddb,
            queue);
    }
}

// --------------------
extern "C" void
magma_dsetvector_async(
    magma_int_t n,
    double const* hx_src,                   magma_int_t incx,
    magmaDouble_ptr    dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueWriteBuffer(
            queue, dy_dst, CL_FALSE,
            dy_offset*sizeof(double), n*sizeof(double),
            hx_src, 0, NULL, event);
        check_error( err );
    }
    else {
        magma_int_t ldha = incx;
        magma_int_t lddb = incy;
        magma_dsetmatrix_async( 1, n,
            hx_src,            ldha,
            dy_dst, dy_offset, lddb,
            queue, event);
    }
}

// --------------------
extern "C" void
magma_dgetvector(
    magma_int_t n,
    magmaDouble_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    double*          hy_dst,                   magma_int_t incy,
    magma_queue_t queue )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueReadBuffer(
            queue, dx_src, CL_TRUE,
            dx_offset*sizeof(double), n*sizeof(double),
            hy_dst, 0, NULL, g_event);
        check_error( err );
    }
    else {
        magma_int_t ldda = incx;
        magma_int_t ldhb = incy;
        magma_dgetmatrix( 1, n,
            dx_src, dx_offset, ldda,
            hy_dst,            ldhb,
            queue);
    }
}

// --------------------
extern "C" void
magma_dgetvector_async(
    magma_int_t n,
    magmaDouble_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    double*          hy_dst,                   magma_int_t incy,
    magma_queue_t queue, magma_event_t *event )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueReadBuffer(
            queue, dx_src, CL_FALSE,
            dx_offset*sizeof(double), n*sizeof(double),
            hy_dst, 0, NULL, event);
        check_error( err );
    }
    else {
        magma_int_t ldda = incx;
        magma_int_t ldhb = incy;
        magma_dgetmatrix_async( 1, n,
            dx_src, dx_offset, ldda,
            hy_dst,            ldhb,
            queue, event);
    }
}

// --------------------
extern "C" void
magma_dcopyvector(
    magma_int_t n,
    magmaDouble_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaDouble_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueReadBuffer(
            queue, dx_src, CL_TRUE,
            dx_offset*sizeof(double), n*sizeof(double),
            dy_dst, dy_offset*sizeof(double), NULL, g_event);
        check_error( err );
    }
    else {
        magma_int_t ldda = incx;
        magma_int_t lddb = incy;
        magma_dcopymatrix( 1, n,
            dx_src, dx_offset, ldda,
            dy_dst, dy_offset, lddb,
            queue);
    }
}

// --------------------
extern "C" void
magma_dcopyvector_async(
    magma_int_t n,
    magmaDouble_const_ptr dx_src, size_t dx_offset, magma_int_t incx,
    magmaDouble_ptr       dy_dst, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue, magma_event_t *event )
{
    if (n <= 0)
        return;

    if (incx == 1 && incy == 1) {
        cl_int err = clEnqueueReadBuffer(
            queue, dx_src, CL_FALSE,
            dx_offset*sizeof(double), n*sizeof(double),
            dy_dst, dy_offset*sizeof(double), NULL, event);
        check_error( err );
    }
    else {
        magma_int_t ldda = incx;
        magma_int_t lddb = incy;
        magma_dcopymatrix_async( 1, n,
            dx_src, dx_offset, ldda,
            dy_dst, dy_offset, lddb,
            queue, event);
    }
}


// ========================================
// copying sub-matrices (contiguous columns)
// OpenCL takes queue even for blocking transfers, oddly.
extern "C" void
magma_dsetmatrix(
    magma_int_t m, magma_int_t n,
    double const* hA_src,                   magma_int_t ldha,
    magmaDouble_ptr    dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue )
{
    if (m <= 0 || n <= 0)
        return;

    size_t buffer_origin[3] = { dB_offset*sizeof(double), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(double), n, 1 };
    cl_int err = clEnqueueWriteBufferRect(
        queue, dB_dst, CL_TRUE,  // blocking
        buffer_origin, host_orig, region,
        lddb*sizeof(double), 0,
        ldha*sizeof(double), 0,
        hA_src, 0, NULL, g_event );
    check_error( err );
}

// --------------------
extern "C" void
magma_dsetmatrix_async(
    magma_int_t m, magma_int_t n,
    double const* hA_src,                   magma_int_t ldha,
    magmaDouble_ptr    dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue, magma_event_t *event )
{
    if (m <= 0 || n <= 0)
        return;

    size_t buffer_origin[3] = { dB_offset*sizeof(double), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(double), n, 1 };
    cl_int err = clEnqueueWriteBufferRect(
        queue, dB_dst, CL_FALSE,  // non-blocking
        buffer_origin, host_orig, region,
        lddb*sizeof(double), 0,
        ldha*sizeof(double), 0,
        hA_src, 0, NULL, event );
    clFlush(queue);
    check_error( err );
}

// --------------------
extern "C" void
magma_dgetmatrix(
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    double*          hB_dst,                   magma_int_t ldhb,
    magma_queue_t queue )
{
    if (m <= 0 || n <= 0)
       return;

    size_t buffer_origin[3] = { dA_offset*sizeof(double), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(double), n, 1 };
    cl_int err = clEnqueueReadBufferRect(
        queue, dA_src, CL_TRUE,  // blocking
        buffer_origin, host_orig, region,
        ldda*sizeof(double), 0,
        ldhb*sizeof(double), 0,
        hB_dst, 0, NULL, g_event );
    check_error( err );
}

// --------------------
extern "C" void
magma_dgetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    double*          hB_dst,                   magma_int_t ldhb,
    magma_queue_t queue, magma_event_t *event )
{
    if (m <= 0 || n <= 0)
        return;

    size_t buffer_origin[3] = { dA_offset*sizeof(double), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(double), n, 1 };
    cl_int err = clEnqueueReadBufferRect(
        queue, dA_src, CL_FALSE,  // non-blocking
        buffer_origin, host_orig, region,
        ldda*sizeof(double), 0,
        ldhb*sizeof(double), 0,
        hB_dst, 0, NULL, event );
    clFlush(queue);
    check_error( err );
}

// --------------------
extern "C" void
magma_dcopymatrix(
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue )
{
    if (m <= 0 || n <= 0)
        return;

    size_t src_origin[3] = { dA_offset*sizeof(double), 0, 0 };
    size_t dst_orig[3]   = { dB_offset*sizeof(double), 0, 0 };
    size_t region[3]     = { m*sizeof(double), n, 1 };
    cl_int err = clEnqueueCopyBufferRect(
        queue, dA_src, dB_dst,
        src_origin, dst_orig, region,
        ldda*sizeof(double), 0,
        lddb*sizeof(double), 0,
        0, NULL, g_event );
    check_error( err );
}

// --------------------
extern "C" void
magma_dcopymatrix_async(
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaDouble_ptr       dB_dst, size_t dB_offset, magma_int_t lddb,
    magma_queue_t queue, magma_event_t *event )
{
    if (m <= 0 || n <= 0)
        return;

    // TODO how to make non-blocking?
    size_t src_origin[3] = { dA_offset*sizeof(double), 0, 0 };
    size_t dst_orig[3]   = { dB_offset*sizeof(double), 0, 0 };
    size_t region[3]     = { m*sizeof(double), n, 1 };
    cl_int err = clEnqueueCopyBufferRect(
        queue, dA_src, dB_dst,
        src_origin, dst_orig, region,
        ldda*sizeof(double), 0,
        lddb*sizeof(double), 0,
        0, NULL, event );
    check_error( err );
}

#endif // HAVE_clBLAS
