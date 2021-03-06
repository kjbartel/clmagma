/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @precisions mixed zc -> ds

       auto-converted from zcaxpycp.cu

*/
#include "clmagma_runtime.h"
#include "common_magma.h"
#include "zcaxpycp.h"


// ----------------------------------------------------------------------
// adds   x += r (including conversion to double)  --and--
// copies w = b
extern "C" void
magmablas_zcaxpycp(
    magma_int_t m,
    magmaFloatComplex_ptr r, size_t r_offset,
    magmaDoubleComplex_ptr x, size_t x_offset,
    magmaDoubleComplex_const_ptr b, size_t b_offset,
    magmaDoubleComplex_ptr w, size_t w_offset,
    magma_queue_t queue )
{
    cl_kernel kernel;
    cl_int err;
    int i;

    size_t threads[1] = { NB };
    size_t grid[1] = { (m + NB - 1)/NB };
    grid[0] *= threads[0];
    kernel = g_runtime.get_kernel( "zcaxpycp_kernel" );
    if ( kernel != NULL ) {
        err = 0;
        i   = 0;
        err  = clSetKernelArg( kernel, i++, sizeof(m       ), &m        );
        err |= clSetKernelArg( kernel, i++, sizeof(r       ), &r        );
        err |= clSetKernelArg( kernel, i++, sizeof(r_offset), &r_offset );
        err |= clSetKernelArg( kernel, i++, sizeof(x       ), &x        );
        err |= clSetKernelArg( kernel, i++, sizeof(x_offset), &x_offset );
        err |= clSetKernelArg( kernel, i++, sizeof(b       ), &b        );
        err |= clSetKernelArg( kernel, i++, sizeof(b_offset), &b_offset );
        err |= clSetKernelArg( kernel, i++, sizeof(w       ), &w        );
        err |= clSetKernelArg( kernel, i++, sizeof(w_offset), &w_offset );
        check_error( err );

        err = clEnqueueNDRangeKernel( queue, kernel, 1, NULL, grid, threads, 0, NULL, NULL );
        check_error( err );
    }
}


// ----------------------------------------------------------------------
// adds   x += r  --and--
// copies r = b
extern "C" void
magmablas_zaxpycp(
    magma_int_t m,
    magmaDoubleComplex_ptr r, size_t r_offset,
    magmaDoubleComplex_ptr x, size_t x_offset,
    magmaDoubleComplex_const_ptr b, size_t b_offset,
    magma_queue_t queue )
{
    cl_kernel kernel;
    cl_int err;
    int i;

    size_t threads[1] = { NB };
    size_t grid[1] = { (m + NB - 1)/NB };
    grid[0] *= threads[0];
    kernel = g_runtime.get_kernel( "zaxpycp_kernel" );
    if ( kernel != NULL ) {
        err = 0;
        i   = 0;
        err  = clSetKernelArg( kernel, i++, sizeof(m       ), &m        );
        err |= clSetKernelArg( kernel, i++, sizeof(r       ), &r        );
        err |= clSetKernelArg( kernel, i++, sizeof(r_offset), &r_offset );
        err |= clSetKernelArg( kernel, i++, sizeof(x       ), &x        );
        err |= clSetKernelArg( kernel, i++, sizeof(x_offset), &x_offset );
        err |= clSetKernelArg( kernel, i++, sizeof(b       ), &b        );
        err |= clSetKernelArg( kernel, i++, sizeof(b_offset), &b_offset );
        check_error( err );

        err = clEnqueueNDRangeKernel( queue, kernel, 1, NULL, grid, threads, 0, NULL, NULL );
        check_error( err );
    }
}
