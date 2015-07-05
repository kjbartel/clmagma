/*
 *   -- clMAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      @date November 2014
 */
#include "clmagma_runtime.h"
#include "common_magma.h"


#define BLOCK_SIZE 64


// empty kernel calling, benchmarkde for overhead for iwocl 2013
// (updated to current formatting standards, no precision generation)
extern "C" void
magmablas_empty(
    magmaDouble_ptr dA,
    magmaDouble_ptr dB,
    magmaDouble_ptr dC,
    magma_queue_t queue )
{
    cl_kernel kernel;
    cl_int err;
    int i;
    
    int m=1, n=1;
    int i0=0, i1=1, i2=1, i3=1, i4=1, i5=1, i6=1, i7=1, i8=1, i9=1;
    double d0=1, d1=1, d2=1, d3=1, d4=1;

    size_t threads[1] = { BLOCK_SIZE };
    size_t grid[1] = { (n+BLOCK_SIZE-1)/BLOCK_SIZE };
    grid[0] *= threads[0];
    
    kernel = g_runtime.get_kernel( "empty_kernel" );
    if ( kernel != NULL ) {
        err = 0;
        int i = 0;
        err  = clSetKernelArg( kernel, i++, sizeof(i0), &i0 );
        err |= clSetKernelArg( kernel, i++, sizeof(i1), &i1 );
        err |= clSetKernelArg( kernel, i++, sizeof(i2), &i2 );
        err |= clSetKernelArg( kernel, i++, sizeof(i3), &i3 );
        err |= clSetKernelArg( kernel, i++, sizeof(i4), &i4 );
        err |= clSetKernelArg( kernel, i++, sizeof(i5), &i5 );
        err |= clSetKernelArg( kernel, i++, sizeof(i6), &i6 );
        err |= clSetKernelArg( kernel, i++, sizeof(i7), &i7 );
        err |= clSetKernelArg( kernel, i++, sizeof(i8), &i8 );
        err |= clSetKernelArg( kernel, i++, sizeof(i9), &i9 );
        
        err |= clSetKernelArg( kernel, i++, sizeof(d0), &d0 );
        err |= clSetKernelArg( kernel, i++, sizeof(d1), &d1 );
        err |= clSetKernelArg( kernel, i++, sizeof(d2), &d2 );
        err |= clSetKernelArg( kernel, i++, sizeof(d3), &d3 );
        err |= clSetKernelArg( kernel, i++, sizeof(d4), &d4 );
        
        err |= clSetKernelArg( kernel, i++, sizeof(dA), &dA );
        err |= clSetKernelArg( kernel, i++, sizeof(dB), &dB );
        err |= clSetKernelArg( kernel, i++, sizeof(dC), &dC );
        check_error( err );
        
        err = clEnqueueNDRangeKernel( queue, kernel, 1, NULL, grid, threads, 0, NULL, NULL );
        check_error( err );
    }
}
