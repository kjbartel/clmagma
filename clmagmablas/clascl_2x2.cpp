/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zlascl_2x2.cpp normal z -> c, Sat Nov 15 00:21:35 2014

       @author Ichitaro Yamazaki
       @author Stan Tomov
*/
#include "clmagma_runtime.h"
#include "common_magma.h"
#include "magmablas.h"

#define NB 64


/**
    Purpose
    -------
    CLASCL2 scales the M by N complex matrix A by the real diagonal matrix dD.
    TYPE specifies that A may be full, upper triangular, lower triangular.

    Arguments
    ---------
    \param[in]
    type    magma_type_t
            TYPE indices the storage type of the input matrix A.
            = MagmaFull:   full matrix.
            = MagmaLower:  lower triangular matrix.
            = MagmaUpper:  upper triangular matrix.
            Other formats that LAPACK supports, MAGMA does not currently support.

    \param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    \param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    \param[in]
    dD      REAL vector, dimension (M)
            The diagonal matrix containing the scalar factors. Stored as a vector.

    \param[in,out]
    dA      COMPLEX array, dimension (LDDA,N)
            The matrix to be scaled by dD.  See TYPE for the
            storage type.

    \param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    \param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value.

    @ingroup magma_caux2
    ********************************************************************/
extern "C" void
magmablas_clascl_2x2(
    magma_type_t type, magma_int_t m, 
    const magmaFloatComplex_ptr dW, size_t dW_offset, magma_int_t lddw, 
    magmaFloatComplex_ptr dA, size_t dA_offset, magma_int_t ldda, 
    magma_int_t *info, magma_queue_t queue )
{
    *info = 0;
    if ( type != MagmaLower && type != MagmaUpper )
        *info = -1;
    else if ( m < 0 )
        *info = -2;
    else if ( ldda < max(1,m) )
        *info = -4;
    
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return;  //info;
    }
    
    size_t  LocalWorkSize[1] = { NB };
    size_t GlobalWorkSize[1] = { ((m + NB - 1)/NB)*NB };

    cl_int ciErrNum;
    cl_kernel ckKernel = NULL;

    if (type == MagmaLower) {
      ckKernel = g_runtime.get_kernel( "clascl_2x2_lower_kernel" );
    }
    else {
      ckKernel = g_runtime.get_kernel( "clascl_2x2_upper_kernel" );
    }

    if(!ckKernel){
      printf ("Error: cannot locate kernel in line %d, file %s\n", __LINE__, __FILE__);
      return;
    }

    int offset_A = (int)dW_offset;
    int offset_B = (int)dA_offset;
    int nn = 0;
    ciErrNum  = clSetKernelArg( ckKernel, nn++, sizeof(cl_int), (void*)&m        );
    ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_mem), (void*)&dW       );
    ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_int), (void*)&offset_A );
    ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_int), (void*)&lddw     );
    ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_mem), (void*)&dA       );
    ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_int), (void*)&offset_B );
    ciErrNum |= clSetKernelArg( ckKernel, nn++, sizeof(cl_int), (void*)&ldda     );
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
