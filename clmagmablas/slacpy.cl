/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zlacpy.cl normal z -> s, Sat Nov 15 00:21:35 2014

       auto-converted from slacpy.cu
       @author Mark Gates
       @author Azzam Haidar
*/
#include "kernels_header.h"
#include "slacpy.h"

/*
    Divides matrix into ceil( m/BLK_X ) x ceil( n/BLK_Y ) blocks.
    Each block has BLK_X threads.
    Each thread loops across one row, updating BLK_Y entries.

    Code similar to slaset.
*/
 __kernel
void slacpy_device_full(
    int m, int n,
    const __global float *dA, unsigned long dA_offset, int ldda,
    __global float       *dB, unsigned long dB_offset, int lddb )
{
    dA += dA_offset;
    dB += dB_offset;

    int ind = get_group_id(0)*BLK_X + get_local_id(0);
    int iby = get_group_id(1)*BLK_Y;
    /* check if full block-column */
    bool full = (iby + BLK_Y <= n);
    /* do only rows inside matrix */
    if ( ind < m ) {
        dA += ind + iby*ldda;
        dB += ind + iby*lddb;
        if ( full ) {
            // full block-column
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
        else {
            // partial block-column
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
    }
}


/*
    Similar to slacpy_full, but updates only the diagonal and below.
    Blocks that are fully above the diagonal exit immediately.

    Code similar to slaset.
*/
 __kernel
void slacpy_device_lower(
    int m, int n,
    const __global float *dA, unsigned long dA_offset, int ldda,
    __global float       *dB, unsigned long dB_offset, int lddb )
{
    dA += dA_offset;
    dB += dB_offset;

    int ind = get_group_id(0)*BLK_X + get_local_id(0);
    int iby = get_group_id(1)*BLK_Y;
    /* check if full block-column && (below diag) */
    bool full = (iby + BLK_Y <= n && (ind >= iby + BLK_Y));
    /* do only rows inside matrix, and blocks not above diag */
    if ( ind < m && ind + BLK_X > iby ) {
        dA += ind + iby*ldda;
        dB += ind + iby*lddb;
        if ( full ) {
            // full block-column, off-diagonal block
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
        else {
            // either partial block-column or diagonal block
            for( int j=0; j < BLK_Y && iby+j < n && ind >= iby+j; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
    }
}


/*
    Similar to slacpy_full, but updates only the diagonal and above.
    Blocks that are fully below the diagonal exit immediately.

    Code similar to slaset.
*/
 __kernel
void slacpy_device_upper(
    int m, int n,
    const __global float *dA, unsigned long dA_offset, int ldda,
    __global float       *dB, unsigned long dB_offset, int lddb )
{
    dA += dA_offset;
    dB += dB_offset;

    int ind = get_group_id(0)*BLK_X + get_local_id(0);
    int iby = get_group_id(1)*BLK_Y;
    /* check if full block-column && (above diag) */
    bool full = (iby + BLK_Y <= n && (ind + BLK_X <= iby));
    /* do only rows inside matrix, and blocks not below diag */
    if ( ind < m && ind < iby + BLK_Y ) {
        dA += ind + iby*ldda;
        dB += ind + iby*lddb;
        if ( full ) {
            // full block-column, off-diagonal block
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
        else {
            // either partial block-column or diagonal block
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                if ( ind <= iby+j ) {
                    dB[j*lddb] = dA[j*ldda];
                }
            }
        }
    }
}

/*
    kernel wrapper to call the device function.
*/
__kernel
void slacpy_kernel_full(
    int m, int n,
    const __global float *dA, unsigned long dA_offset, int ldda,
    __global float       *dB, unsigned long dB_offset, int lddb )
{
    slacpy_device_full(m, n, dA, dA_offset, ldda, dB, dB_offset, lddb);
}

__kernel
void slacpy_kernel_lower(
    int m, int n,
    const __global float *dA, unsigned long dA_offset, int ldda,
    __global float       *dB, unsigned long dB_offset, int lddb )
{
    slacpy_device_lower(m, n, dA, dA_offset, ldda, dB, dB_offset, lddb);
}

__kernel
void slacpy_kernel_upper(
    int m, int n,
    const __global float *dA, unsigned long dA_offset, int ldda,
    __global float       *dB, unsigned long dB_offset, int lddb )
{
    slacpy_device_upper(m, n, dA, dA_offset, ldda, dB, dB_offset, lddb);
}


/*
    kernel wrapper to call the device function for the batched routine.
*/
/*
__kernel
void slacpy_kernel_batched_full(
    int m, int n,
    __global float const * __global  const *dAarray, unsigned long dA_offset, int ldda,
    __global float       * __global        *dBarray, unsigned long dB_offset, int lddb )
{
    int batchid = get_group_id(2);
    slacpy_device_full(m, n, dAarray[batchid], dA_offset, ldda, dBarray[batchid], dB_offset, lddb);
}

__kernel
void slacpy_kernel_batched_lower(
    int m, int n,
    __global float const * __global  const *dAarray, unsigned long dA_offset, int ldda,
    __global float       * __global        *dBarray, unsigned long dB_offset, int lddb )
{
    int batchid = get_group_id(2);
    slacpy_device_lower(m, n, dAarray[batchid], dA_offset, ldda, dBarray[batchid], dB_offset, lddb);
}

__kernel
void slacpy_kernel_batched_upper(
    int m, int n,
    __global float const * __global  const *dAarray, unsigned long dA_offset, int ldda,
    __global float       * __global        *dBarray, unsigned long dB_offset, int lddb )
{
    int batchid = get_group_id(2);
    slacpy_device_upper(m, n, dAarray[batchid], dA_offset, ldda, dBarray[batchid], dB_offset, lddb);
}
*/
