/*
 *   -- clMAGMA (version 1.0.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      August 2012
 *
 * @generated s Wed Oct 24 00:32:57 2012
 */

#include <stdio.h>

#include "magmablas.h"
#include "CL_MAGMA_RT.h"

#define BLOCK_SIZE 64

/*********************************************************
 *
 * SWAP BLAS: permute to set of N elements
 *
 ********************************************************/
/*
 *  First version: line per line
 */
typedef struct {
/*
    magmaFloat_ptr A1;
    magmaFloat_ptr A2;
*/
    int n, offset_dA1, lda1, offset_dA2, lda2;
} magmagpu_sswap_params_t;


extern "C" void 
magmablas_sswap( magma_int_t n, magmaFloat_ptr dA1T, size_t offset_dA1T, magma_int_t lda1, 
                 magmaFloat_ptr dA2T, size_t offset_dA2T, magma_int_t lda2, magma_queue_t queue)
{
    int  blocksize = BLOCK_SIZE;
	size_t LocalWorkSize[1] = {blocksize};
	size_t GlobalWorkSize[1] = {((n+blocksize-1)/blocksize)*LocalWorkSize[0]};

    magmagpu_sswap_params_t params = { n, (int)offset_dA1T, lda1, (int)offset_dA2T, lda2 };
	
	cl_int ciErrNum;                // Error code var
	cl_kernel ckKernel=NULL;
	
	ckKernel = rt->KernelPool["magmagpu_sswap"];
	if (!ckKernel)
	{
		printf ("Error: cannot locate kernel in line %d, file %s\n", __LINE__, __FILE__);
		return;
	}
	
	int nn = 0;
	ciErrNum  = clSetKernelArg( ckKernel, nn++, sizeof(cl_mem), (void*)&dA1T   );
	ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_mem), (void*)&dA2T   );
	ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(magmagpu_sswap_params_t),           (void*)&params	);
	if (ciErrNum != CL_SUCCESS)
	{
		printf("Error: clSetKernelArg at %d in file %s, %s\n", __LINE__, __FILE__, rt->GetErrorCode(ciErrNum));
		return;
	}
	// launch kernel
	ciErrNum = clEnqueueNDRangeKernel(queue, ckKernel, 1, NULL, GlobalWorkSize, LocalWorkSize, 0, NULL, NULL);
	if (ciErrNum != CL_SUCCESS)
	{
		printf("Error: clEnqueueNDRangeKernel at %d in file %s \"%s\"\n",
			__LINE__, __FILE__, rt->GetErrorCode(ciErrNum));
		return;
	}
}

