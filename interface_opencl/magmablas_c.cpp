/*
 *   -- clMAGMA (version 0.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 * @generated c Wed Apr  4 01:12:54 2012
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
magma_csetmatrix(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex const* hA_src, size_t hA_offset, magma_int_t ldha,
    magmaFloatComplex_ptr    dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue )
{
    size_t buffer_origin[3] = { dA_offset*sizeof(magmaFloatComplex), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(magmaFloatComplex), n, 1 };
    cl_int err = clEnqueueWriteBufferRect(
        queue, dA_dst, CL_TRUE,  // blocking
        buffer_origin, host_orig, region,
        ldda*sizeof(magmaFloatComplex), 0,
        ldha*sizeof(magmaFloatComplex), 0,
        hA_src, 0, NULL, NULL );
    return err;
}

// --------------------
magma_err_t
magma_cgetmatrix(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex*          hA_dst, size_t hA_offset, magma_int_t ldha,
    magma_queue_t queue )
{
    size_t buffer_origin[3] = { dA_offset*sizeof(magmaFloatComplex), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(magmaFloatComplex), n, 1 };
    cl_int err = clEnqueueReadBufferRect(
        queue, dA_src, CL_TRUE,  // blocking
        buffer_origin, host_orig, region,
        ldda*sizeof(magmaFloatComplex), 0,
        ldha*sizeof(magmaFloatComplex), 0,
        hA_dst, 0, NULL, NULL );
    return err;
}

// --------------------
magma_err_t
magma_csetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex const* hA_src, size_t hA_offset, magma_int_t ldha,
    magmaFloatComplex_ptr    dA_dst, size_t dA_offset, magma_int_t ldda,
    magma_queue_t queue, magma_event_t *event )
{
    size_t buffer_origin[3] = { dA_offset*sizeof(magmaFloatComplex), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(magmaFloatComplex), n, 1 };
    cl_int err = clEnqueueWriteBufferRect(
        queue, dA_dst, CL_FALSE,  // non-blocking
        buffer_origin, host_orig, region,
        ldda*sizeof(magmaFloatComplex), 0,
        ldha*sizeof(magmaFloatComplex), 0,
        hA_src, 0, NULL, event );
    return err;
}

// --------------------
magma_err_t
magma_cgetmatrix_async(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_const_ptr dA_src, size_t dA_offset, magma_int_t ldda,
    magmaFloatComplex*          hA_dst, size_t hA_offset, magma_int_t ldha,
    magma_queue_t queue, magma_event_t *event )
{
    size_t buffer_origin[3] = { dA_offset*sizeof(magmaFloatComplex), 0, 0 };
    size_t host_orig[3]     = { 0, 0, 0 };
    size_t region[3]        = { m*sizeof(magmaFloatComplex), n, 1 };
    cl_int err = clEnqueueReadBufferRect(
        queue, dA_src, CL_FALSE,  // non-blocking
        buffer_origin, host_orig, region,
        ldda*sizeof(magmaFloatComplex), 0,
        ldha*sizeof(magmaFloatComplex), 0,
        hA_dst, 0, NULL, event );
    return err;
}

// ========================================
// BLAS functions
magma_err_t
magma_cgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloatComplex alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaFloatComplex_const_ptr dB, size_t dB_offset, magma_int_t ldb,
    magmaFloatComplex beta,  magmaFloatComplex_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasCgemmEx(
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
magma_cgemv(
    magma_trans_t transA,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaFloatComplex_const_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaFloatComplex beta,  magmaFloatComplex_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasCgemvEx(
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
magma_chemm(
    magma_side_t side, magma_uplo_t uplo,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaFloatComplex_const_ptr dB, size_t dB_offset, magma_int_t ldb,
    magmaFloatComplex beta,  magmaFloatComplex_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasChemm(
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
magma_chemv(
    magma_uplo_t uplo,
    magma_int_t n,
    magmaFloatComplex alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaFloatComplex_const_ptr dx, size_t dx_offset, magma_int_t incx,
    magmaFloatComplex beta,  magmaFloatComplex_ptr       dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasChemv(
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
magma_cherk(
    magma_uplo_t uplo, magma_trans_t trans,
    magma_int_t n, magma_int_t k,
    float alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
    float beta,  magmaFloatComplex_ptr       dC, size_t dC_offset, magma_int_t ldc,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasCherk(
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
magma_ctrsm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaFloatComplex_ptr       dB, size_t dB_offset, magma_int_t ldb,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasCtrsmEx(
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
magma_ctrmm(
    magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t m, magma_int_t n,
    magmaFloatComplex alpha, magmaFloatComplex_const_ptr dA, size_t dA_offset, magma_int_t lda,
                              magmaFloatComplex_ptr       dB, size_t dB_offset, magma_int_t ldb,
    magma_queue_t queue )
{
    cl_int err = clAmdBlasCtrmmEx(
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
