#include <stdio.h>

#include "magmablas.h"
#include "CL_MAGMA_RT.h"

#define BLOCK_SIZE 64

typedef struct {
        int n;
		int lda;
		int j0;
		int npivots;
        short ipiv[BLOCK_SIZE];
} slaswp_params_t2;

typedef struct 
{
	float a;
	float b;
	float c;
	float d;
} float4;


// ----------------------------------------
void slaswp3(
	cl_mem dA, size_t offset, slaswp_params_t2 params, magma_queue_t queue )
{
	/*
	int blocksize = 64;
	dim3 blocks = (params.n+blocksize-1) / blocksize;
	myslaswp2<<< blocks, blocksize, 0, magma_stream >>>( params );
	*/
	
	cl_int ciErrNum;                // Error code var
	cl_kernel ckKernel=NULL;
	
	ckKernel = rt.KernelPool["myslaswp2"];
	if (!ckKernel)
	{
		printf ("Error: cannot locate kernel in line %d, file %s\n", __LINE__, __FILE__);
		return;
	}
	
	int nn = 0;
	ciErrNum  = clSetKernelArg( ckKernel, nn++, sizeof(cl_mem),           (void*)&dA     );
	ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_int),           (void*)&offset );
	ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(slaswp_params_t2), (void*)&params );
	if (ciErrNum != CL_SUCCESS)
	{
		printf("Error: clSetKernelArg at %d in file %s, %s\n", __LINE__, __FILE__, rt.GetErrorCode(ciErrNum));
		return;
	}
	
	size_t GlobalWorkSize[2]={0,0}, LocalWorkSize[2]={0,0};
	
	LocalWorkSize[0] = BLOCK_SIZE;
	LocalWorkSize[1] = 1;
	
	GlobalWorkSize[0] = ((params.n+BLOCK_SIZE-1) / BLOCK_SIZE)*LocalWorkSize[0];
	GlobalWorkSize[1] = 1;
	
	// launch kernel
	ciErrNum = clEnqueueNDRangeKernel(
		queue, ckKernel, 2, NULL, GlobalWorkSize, LocalWorkSize, 0, NULL, NULL);
	if (ciErrNum != CL_SUCCESS)
	{
		printf("Error: clEnqueueNDRangeKernel at %d in file %s \"%s\"\n",
			__LINE__, __FILE__, rt.GetErrorCode(ciErrNum));
		return;
	}
}


// ----------------------------------------
magma_err_t
magma_spermute_long2(
	cl_mem dAT, size_t dAT_offset, int lda,
	int *ipiv, int nb, int ind,
	magma_queue_t queue )
{
	int k;
	
	for( k = 0; k < nb-BLOCK_SIZE; k += BLOCK_SIZE )
	{
		slaswp_params_t2 params = { lda, lda, ind + k, BLOCK_SIZE };
		for( int j = 0; j < BLOCK_SIZE; j++ )
		{
			params.ipiv[j] = ipiv[ind + k + j] - k - 1;
			ipiv[ind + k + j] += ind;
		}
		slaswp3( dAT, dAT_offset, params, queue );
	}
	
	int num_pivots = nb - k;
	
	slaswp_params_t2 params = { lda, lda, ind + k, num_pivots };
	for( int j = 0; j < num_pivots; j++ )
	{
		params.ipiv[j] = ipiv[ind + k + j] - k - 1;
		ipiv[ind + k + j] += ind;
	}
	slaswp3( dAT, dAT_offset, params, queue );
	
	return MAGMA_SUCCESS;
}
