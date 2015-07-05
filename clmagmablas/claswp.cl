/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zlaswp.cl normal z -> c, Sat Nov 15 00:21:35 2014

       auto-converted from claswp.cu
       
       @author Stan Tomov
       @author Mathieu Faverge
       @author Ichitaro Yamazaki
       @author Mark Gates
*/
#include "kernels_header.h"
#include "claswp.h"


// Matrix A is stored row-wise in dAT.
// Divide matrix A into block-columns of NTHREADS columns each.
// Each GPU block processes one block-column of A.
// Each thread goes down a column of A,
// swapping rows according to pivots stored in params.
__kernel void claswp_kernel(
    int n, __global magmaFloatComplex *dAT, unsigned long dAT_offset, int ldda, claswp_params_t params )
{
    dAT += dAT_offset;

    int tid = get_local_id(0) + get_local_size(0)*get_group_id(0);
    if ( tid < n ) {
        dAT += tid;
        __global magmaFloatComplex *A1  = dAT;
        
        for( int i1 = 0; i1 < params.npivots; ++i1 ) {
            int i2 = params.ipiv[i1];
            __global magmaFloatComplex *A2 = dAT + i2*ldda;
            magmaFloatComplex temp = *A1;
            *A1 = *A2;
            *A2 = temp;
            A1 += ldda;  // A1 = dA + i1*ldx
        }
    }
}






// ------------------------------------------------------------
// Extended version has stride in both directions (ldx, ldy)
// to handle both row-wise and column-wise storage.

// Matrix A is stored row or column-wise in dA.
// Divide matrix A into block-columns of NTHREADS columns each.
// Each GPU block processes one block-column of A.
// Each thread goes down a column of A,
// swapping rows according to pivots stored in params.
__kernel void claswpx_kernel(
    int n, __global magmaFloatComplex *dA, unsigned long dA_offset, int ldx, int ldy, claswp_params_t params )
{
    dA += dA_offset;

    int tid = get_local_id(0) + get_local_size(0)*get_group_id(0);
    if ( tid < n ) {
        dA += tid*ldy;
        __global magmaFloatComplex *A1  = dA;
        
        for( int i1 = 0; i1 < params.npivots; ++i1 ) {
            int i2 = params.ipiv[i1];
            __global magmaFloatComplex *A2 = dA + i2*ldx;
            magmaFloatComplex temp = *A1;
            *A1 = *A2;
            *A2 = temp;
            A1 += ldx;  // A1 = dA + i1*ldx
        }
    }
}






// ------------------------------------------------------------
// This version takes d_ipiv on the GPU. Thus it does not pass pivots
// as an argument using a structure, avoiding all the argument size
// limitations of CUDA and OpenCL. It also needs just one kernel launch
// with all the pivots, instead of multiple kernel launches with small
// batches of pivots. On Fermi, it is faster than magmablas_claswp
// (including copying pivots to the GPU).

__kernel void claswp2_kernel(
    int n, __global magmaFloatComplex *dAT, unsigned long dAT_offset, int ldda, int npivots,
    const __global magma_int_t *d_ipiv, unsigned long d_ipiv_offset, int inci )
{
    dAT += dAT_offset;
    d_ipiv += d_ipiv_offset;

    int tid = get_local_id(0) + get_local_size(0)*get_group_id(0);
    if ( tid < n ) {
        dAT += tid;
        __global magmaFloatComplex *A1  = dAT;
        
        for( int i1 = 0; i1 < npivots; ++i1 ) {
            int i2 = d_ipiv[i1*inci] - 1;  // Fortran index
            __global magmaFloatComplex *A2 = dAT + i2*ldda;
            magmaFloatComplex temp = *A1;
            *A1 = *A2;
            *A2 = temp;
            A1 += ldda;  // A1 = dA + i1*ldx
        }
    }
}
