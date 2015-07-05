/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zlacpy_cnjg.cpp normal z -> s, Sat Nov 15 00:21:35 2014

       @author Ichitaro Yamazaki
       @author Stan Tomov
*/
#include "clmagma_runtime.h"
#include "common_magma.h"
#include "magmablas.h"

/*********************************************************
 *
 * SWAP BLAS: permute to set of N elements
 *
 ********************************************************/
/*
 *  First version: line per line
 */

extern "C" void 
magmablas_slacpy_cnjg(
    magma_int_t n, 
    magmaFloat_ptr dA1T, size_t dA1T_offset, magma_int_t lda1, 
    magmaFloat_ptr dA2T, size_t dA2T_offset, magma_int_t lda2,
    magma_queue_t queue)
{
    int blocksize = 64;
    size_t LocalWorkSize[1]  = {blocksize};
    size_t GlobalWorkSize[1] = { ( (n+blocksize-1) / blocksize) * LocalWorkSize[0] };

    if ( n <=0 ) {
      return;
    }

    cl_int ciErrNum;
    cl_kernel ckKernel = NULL;
    ckKernel = g_runtime.get_kernel( "slacpy_cnjg_kernel" );
    if(!ckKernel){
        printf ("Error: cannot locate kernel in line %d, file %s\n", __LINE__, __FILE__);
        return;
    }

    int offset_A = (int)dA1T_offset;
    int offset_B = (int)dA2T_offset;
    int nn = 0;
    ciErrNum  = clSetKernelArg( ckKernel, nn++, sizeof(cl_int), (void*)&n        );
    ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_mem), (void*)&dA1T     );
    ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_int), (void*)&offset_A );
    ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_int), (void*)&lda1     );
    ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_mem), (void*)&dA2T     );
    ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_int), (void*)&offset_B );
    ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_int), (void*)&lda2     );
    if (ciErrNum != CL_SUCCESS){
      printf("Error: clSetKernelArg at %d in file %s, %s\n", __LINE__, __FILE__, 
             magma_strerror(ciErrNum));
      return;
    }

    // launch kernel                                                                         
    ciErrNum = clEnqueueNDRangeKernel(
                  queue, ckKernel, 1, NULL, GlobalWorkSize, LocalWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        printf("Error: clEnqueueNDRangeKernel at %d in file %s \"%s\"\n",
               __LINE__, __FILE__, magma_strerror(ciErrNum));
        return;
    }
}
