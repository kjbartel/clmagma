/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014
       
       @author Mark Gates

       @generated from zswap.cpp normal z -> d, Sat Nov 15 00:21:35 2014

       auto-converted from dswap.cu

*/
#include "clmagma_runtime.h"
#include "common_magma.h"
#include "dswap.h"


/**
    Purpose:
    =============
    Swap vector x and y; \f$ x <-> y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in,out]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      DOUBLE_PRECISION array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_dblas1
    ********************************************************************/
extern "C" void 
magmablas_dswap(
    magma_int_t n,
    magmaDouble_ptr dx, size_t dx_offset, magma_int_t incx, 
    magmaDouble_ptr dy, size_t dy_offset, magma_int_t incy,
    magma_queue_t queue )
{
    cl_kernel kernel;
    cl_int err;
    int i;

    size_t grid[1] = { (n+NB-1) / NB };
    size_t threads[1] = { NB };
    grid[0] *= threads[0];
    kernel = g_runtime.get_kernel( "dswap_kernel" );
    if ( kernel != NULL ) {
        err = 0;
        i   = 0;
        err |= clSetKernelArg( kernel, i++, sizeof(n        ), &n         );
        err |= clSetKernelArg( kernel, i++, sizeof(dx       ), &dx        );
        err |= clSetKernelArg( kernel, i++, sizeof(dx_offset), &dx_offset );
        err |= clSetKernelArg( kernel, i++, sizeof(incx     ), &incx      );
        err |= clSetKernelArg( kernel, i++, sizeof(dy       ), &dy        );
        err |= clSetKernelArg( kernel, i++, sizeof(dy_offset), &dy_offset );
        err |= clSetKernelArg( kernel, i++, sizeof(incy     ), &incy      );
        check_error( err );

        err = clEnqueueNDRangeKernel( queue, kernel, 1, NULL, grid, threads, 0, NULL, NULL );
        check_error( err );
    }
}
