/*
 *   -- clMAGMA (version 0.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @generated d Wed Apr  4 01:12:54 2012
 */

#include <stdlib.h>
#include <stdio.h>

#include "magma.h"

#ifdef HAVE_AMDBLAS

// ========================================
// globals, defined in interface.c
extern cl_platform_id gPlatform;
extern cl_context     gContext;


// ========================================
// copying sub-matrices (contiguous columns)
// OpenCL takes queue even for blocking transfers, oddly.
magma_err_t
magma_dsetmatrix(
    magma_int_t m, magma_int_t n,
    double const* hA_src, size_t hA_offset, magma_int_t ldha,
    magmaDouble_ptr    dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue )
{
    size_t buffer_origin[3] = { dA_offset*sizeof(double), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(double), n, 1 };
    cl_int err = clEnqueueWriteBufferRect(
        queue, dA_dst, CL_TRUE,  // blocking
        buffer_origin, host_orig, region,
        ldda*sizeof(double), 0,
        ldha*sizeof(double), 0,
        hA_src, 0, NULL, NULL );
    return err;
}

// --------------------
magma_err_t
magma_dgetmatrix(
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    double*          hA_dst, size_t hA_offset, magma_int_t ldha,
    magma_queue_t queue )
{
    size_t buffer_origin[3] = { dA_offset*sizeof(double), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(double), n, 1 };
    cl_int err = clEnqueueReadBufferRect(
        queue, dA_src, CL_TRUE,  // blocking
        buffer_origin, host_orig, region,
        ldda*sizeof(double), 0,
        ldha*sizeof(double), 0,
        hA_dst, 0, NULL, NULL );
    return err;
}

// --------------------
magma_err_t
magma_dsetmatrix_async(
    magma_int_t m, magma_int_t n,
    double const* hA_src, size_t hA_offset, magma_int_t ldha,
    magmaDouble_ptr    dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue, magma_event_t *event )
{
    size_t buffer_origin[3] = { dA_offset*sizeof(double), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(double), n, 1 };
    cl_int err = clEnqueueWriteBufferRect(
        queue, dA_dst, CL_FALSE,  // non-blocking
        buffer_origin, host_orig, region,
        ldda*sizeof(double), 0,
        ldha*sizeof(double), 0,
        hA_src, 0, NULL, event );
    return err;
}

// --------------------
magma_err_t
magma_dgetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaDouble_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    double*          hA_dst, size_t hA_offset, magma_int_t ldha,
    magma_queue_t queue, magma_event_t *event )
{
    size_t buffer_origin[3] = { dA_offset*sizeof(double), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(double), n, 1 };
    cl_int err = clEnqueueReadBufferRect(
        queue, dA_src, CL_FALSE,  // non-blocking
        buffer_origin, host_orig, region,
        ldda*sizeof(double), 0,
        ldha*sizeof(double), 0,
        hA_dst, 0, NULL, event );
    return err;
}

// ========================================
// BLAS functions
magma_err_t
magma_dgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDouble_const_ptr dB, size_t dB_offset, magma_int_t ldb,
    double beta,  magmaDouble_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasDgemmEx(
        clAmdBlasColumnMajor,
        amdblas_trans_const( transA ),
        amdblas_trans_const( transB ),
        m, n, k,
        alpha, dA, dA_offset, lda,
               dB, dB_offset, ldb,
        beta,  dC, dC_offset, ldc,
        1, &queue, 0, NULL, NULL );
    return err;
}

// --------------------
magma_err_t
magma_dgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDouble_const_ptr dx, size_t dx_offset, magma_int_t incx,
    double beta,  magmaDouble_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasDgemvEx(
        clAmdBlasColumnMajor,
        amdblas_trans_const( transA ),
        m, n,
        alpha, dA, dA_offset, lda,
               dx, dx_offset, incx,
        beta,  dy, dy_offset, incy,
        1, &queue, 0, NULL, NULL );
    return err;
}

// --------------------
magma_err_t
magma_dsymm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDouble_const_ptr dB, size_t dB_offset, magma_int_t ldb,
    double beta,  magmaDouble_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasDsymm(
        clAmdBlasColumnMajor,
        amdblas_side_const( side ),
        amdblas_uplo_const( uplo ),
        m, n,
        alpha, dA, dA_offset, lda,
               dB, dB_offset, ldb,
        beta,  dC, dC_offset, ldc,
        1, &queue, 0, NULL, NULL );
    return err;
}

// --------------------
magma_err_t
magma_dsymv(
    magma_uplo_t uplo,
    magma_int_t n,
    double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDouble_const_ptr dx, size_t dx_offset, magma_int_t incx,
    double beta,  magmaDouble_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasDsymvEx(
        clAmdBlasColumnMajor,
        amdblas_uplo_const( uplo ),
        n,
        alpha, dA, dA_offset, lda,
               dx, dx_offset, incx,
        beta,  dy, dy_offset, incy,
        1, &queue, 0, NULL, NULL );
    return err;
}

// --------------------
magma_err_t
magma_dsyrk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
    double beta,  magmaDouble_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasDsyrkEx(
        clAmdBlasColumnMajor,
        amdblas_uplo_const( uplo ),
        amdblas_trans_const( trans ),
        n, k,
        alpha, dA, dA_offset, lda,
        beta,  dC, dC_offset, ldc,
        1, &queue, 0, NULL, NULL );
    return err;
}

// --------------------
magma_err_t
magma_dtrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDouble_ptr       dB, size_t dB_offset, magma_int_t ldb,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasDtrsmEx(
        clAmdBlasColumnMajor,
        amdblas_side_const( side ),
        amdblas_uplo_const( uplo ),
        amdblas_trans_const( trans ),
        amdblas_diag_const( diag ),
        m, n,
        alpha, dA, dA_offset, lda,
               dB, dB_offset, ldb,
        1, &queue, 0, NULL, NULL );
    return err;
}

// --------------------
magma_err_t
magma_dtrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    double alpha, magmaDouble_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaDouble_ptr       dB, size_t dB_offset, magma_int_t ldb,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasDtrmmEx(
        clAmdBlasColumnMajor,
        amdblas_side_const( side ),
        amdblas_uplo_const( uplo ),
        amdblas_trans_const( trans ),
        amdblas_diag_const( diag ),
        m, n,
        alpha, dA, dA_offset, lda,
               dB, dB_offset, ldb,
        1, &queue, 0, NULL, NULL );
    return err;
}

#endif // HAVE_AMDBLAS
