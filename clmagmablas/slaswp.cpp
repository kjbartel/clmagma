/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zlaswp.cpp normal z -> s, Sat Nov 15 00:21:35 2014

       auto-converted from slaswp.cu
       
       @author Stan Tomov
       @author Mathieu Faverge
       @author Ichitaro Yamazaki
       @author Mark Gates
*/
#include "clmagma_runtime.h"
#include "common_magma.h"
#include "slaswp.h"


/**
    Purpose:
    =============
    SLASWP performs a series of row interchanges on the matrix A.
    One row interchange is initiated for each of rows K1 through K2 of A.
    
    ** Unlike LAPACK, here A is stored row-wise (hence dAT). **
    Otherwise, this is identical to LAPACK's interface.
    
    Arguments:
    ==========
    \param[in]
    n        INTEGER
             The number of columns of the matrix A.
    
    \param[in,out]
    dAT      REAL array on GPU, stored row-wise, dimension (LDDA,N)
             On entry, the matrix of column dimension N to which the row
             interchanges will be applied.
             On exit, the permuted matrix.
    
    \param[in]
    ldda     INTEGER
             The leading dimension of the array A. ldda >= n.
    
    \param[in]
    k1       INTEGER
             The first element of IPIV for which a row interchange will
             be done. (Fortran one-based index: 1 <= k1 <= n.)
    
    \param[in]
    k2       INTEGER
             The last element of IPIV for which a row interchange will
             be done. (Fortran one-based index: 1 <= k2 <= n.)
    
    \param[in]
    ipiv     INTEGER array, on CPU, dimension (K2*abs(INCI))
             The vector of pivot indices.  Only the elements in positions
             K1 through K2 of IPIV are accessed.
             IPIV(K) = L implies rows K and L are to be interchanged.
    
    \param[in]
    inci     INTEGER
             The increment between successive values of IPIV.
             Currently, INCI > 0.
             TODO: If INCI is negative, the pivots are applied in reverse order.

    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    @ingroup magma_saux2
    ********************************************************************/
// It is used in sgessm, sgetrf_incpiv.
extern "C" void
magmablas_slaswp(
    magma_int_t n,
    magmaFloat_ptr dAT, size_t dAT_offset, magma_int_t ldda,
    magma_int_t k1, magma_int_t k2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_queue_t queue )
{
    cl_kernel kernel;
    cl_int err;
    int i;
    
    magma_int_t info = 0;
    if ( n < 0 )
        info = -1;
    else if ( k1 < 1 || k1 > n )
        info = -4;
    else if ( k2 < 1 || k2 > n )
        info = -5;
    else if ( inci <= 0 )
        info = -7;

    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    size_t grid[1] = { (n + NTHREADS - 1) / NTHREADS };
    size_t threads[1] = { NTHREADS };
    grid[0] *= threads[0];
    slaswp_params_t params;
    
    kernel = g_runtime.get_kernel( "slaswp_kernel" );
    
    for( int k = k1-1; k < k2; k += MAX_PIVOTS ) {
        int npivots = min( MAX_PIVOTS, k2-k );
        params.npivots = npivots;
        for( int j = 0; j < npivots; ++j ) {
            params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
        }
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            size_t k_offset = dAT_offset + k*ldda;
            err |= clSetKernelArg( kernel, i++, sizeof(n       ), &n        );
            err |= clSetKernelArg( kernel, i++, sizeof(dAT     ), &dAT      );
            err |= clSetKernelArg( kernel, i++, sizeof(k_offset), &k_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldda    ), &ldda     );
            err |= clSetKernelArg( kernel, i++, sizeof(params  ), &params   );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 1, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    }
}


/**
    Purpose:
    =============
    SLASWPX performs a series of row interchanges on the matrix A.
    One row interchange is initiated for each of rows K1 through K2 of A.
    
    ** Unlike LAPACK, here A is stored either row-wise or column-wise,
       depending on ldx and ldy. **
    Otherwise, this is identical to LAPACK's interface.
    
    Arguments:
    ==========
    \param[in]
    n        INTEGER
             The number of columns of the matrix A.
    
    \param[in,out]
    dA       REAL array on GPU, dimension (*,*)
             On entry, the matrix of column dimension N to which the row
             interchanges will be applied.
             On exit, the permuted matrix.
    
    \param[in]
    ldx      INTEGER
             Stride between elements in same column.
    
    \param[in]
    ldy      INTEGER
             Stride between elements in same row.
             For A stored row-wise,    set ldx=ldda and ldy=1.
             For A stored column-wise, set ldx=1    and ldy=ldda.
    
    \param[in]
    k1       INTEGER
             The first element of IPIV for which a row interchange will
             be done. (One based index.)
    
    \param[in]
    k2       INTEGER
             The last element of IPIV for which a row interchange will
             be done. (One based index.)
    
    \param[in]
    ipiv     INTEGER array, on CPU, dimension (K2*abs(INCI))
             The vector of pivot indices.  Only the elements in positions
             K1 through K2 of IPIV are accessed.
             IPIV(K) = L implies rows K and L are to be interchanged.
    
    \param[in]
    inci     INTEGER
             The increment between successive values of IPIV.
             Currently, IPIV > 0.
             TODO: If IPIV is negative, the pivots are applied in reverse order.

    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    @ingroup magma_saux2
    ********************************************************************/
extern "C" void
magmablas_slaswpx(
    magma_int_t n,
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldx, magma_int_t ldy,
    magma_int_t k1, magma_int_t k2,
    const magma_int_t *ipiv, magma_int_t inci,
    magma_queue_t queue )
{
    cl_kernel kernel;
    cl_int err;
    int i;

    magma_int_t info = 0;
    if ( n < 0 )
        info = -1;
    else if ( k1 < 0 )
        info = -4;  
    else if ( k2 < 0 || k2 < k1 )
        info = -5;
    else if ( inci <= 0 )
        info = -7;

    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    size_t grid[1] = { (n + NTHREADS - 1) / NTHREADS };
    size_t threads[1] = { NTHREADS };
    grid[0] *= threads[0];
    slaswp_params_t params;
    
    kernel = g_runtime.get_kernel( "slaswpx_kernel" );
    
    for( int k = k1-1; k < k2; k += MAX_PIVOTS ) {
        int npivots = min( MAX_PIVOTS, k2-k );
        params.npivots = npivots;
        for( int j = 0; j < npivots; ++j ) {
            params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
        }
        if ( kernel != NULL ) {
            err = 0;
            i   = 0;
            size_t k_offset = dA_offset + k*ldx;
            err |= clSetKernelArg( kernel, i++, sizeof(n       ), &n        );
            err |= clSetKernelArg( kernel, i++, sizeof(dA      ), &dA       );
            err |= clSetKernelArg( kernel, i++, sizeof(k_offset), &k_offset );
            err |= clSetKernelArg( kernel, i++, sizeof(ldx     ), &ldx      );
            err |= clSetKernelArg( kernel, i++, sizeof(ldy     ), &ldy      );
            err |= clSetKernelArg( kernel, i++, sizeof(params  ), &params   );
            check_error( err );

            err = clEnqueueNDRangeKernel( queue, kernel, 1, NULL, grid, threads, 0, NULL, NULL );
            check_error( err );
        }
    }
}


/**
    Purpose:
    =============
    SLASWP2 performs a series of row interchanges on the matrix A.
    One row interchange is initiated for each of rows K1 through K2 of A.
    
    ** Unlike LAPACK, here A is stored row-wise (hence dAT). **
    Otherwise, this is identical to LAPACK's interface.
    
    Here, d_ipiv is passed in GPU memory.
    
    Arguments:
    ==========
    \param[in]
    n        INTEGER
             The number of columns of the matrix A.
    
    \param[in,out]
    dAT      REAL array on GPU, stored row-wise, dimension (LDDA,*)
             On entry, the matrix of column dimension N to which the row
             interchanges will be applied.
             On exit, the permuted matrix.
    
    \param[in]
    ldda     INTEGER
             The leading dimension of the array A.
             (I.e., stride between elements in a column.)
    
    \param[in]
    k1       INTEGER
             The first element of IPIV for which a row interchange will
             be done. (One based index.)
    
    \param[in]
    k2       INTEGER
             The last element of IPIV for which a row interchange will
             be done. (One based index.)
    
    \param[in]
    d_ipiv   INTEGER array, on GPU, dimension (K2*abs(INCI))
             The vector of pivot indices.  Only the elements in positions
             K1 through K2 of IPIV are accessed.
             IPIV(K) = L implies rows K and L are to be interchanged.
    
    \param[in]
    inci     INTEGER
             The increment between successive values of IPIV.
             Currently, IPIV > 0.
             TODO: If IPIV is negative, the pivots are applied in reverse order.

    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    @ingroup magma_saux2
    ********************************************************************/
extern "C" void
magmablas_slaswp2(
    magma_int_t n,
    magmaFloat_ptr dAT, size_t dAT_offset, magma_int_t ldda,
    magma_int_t k1, magma_int_t k2,
    magmaInt_const_ptr d_ipiv, size_t d_ipiv_offset, magma_int_t inci,
    magma_queue_t queue )
{
    cl_kernel kernel;
    cl_int err;
    int i;

    magma_int_t info = 0;
    if ( n < 0 )
        info = -1;
    else if ( k1 < 0 )
        info = -4;  
    else if ( k2 < 0 || k2 < k1 )
        info = -5;
    else if ( inci <= 0 )
        info = -7;

    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_int_t nb = k2-(k1-1);
    
    size_t grid[1] = { (n + NTHREADS - 1) / NTHREADS };
    size_t threads[1] = { NTHREADS };
    grid[0] *= threads[0];
    kernel = g_runtime.get_kernel( "slaswp2_kernel" );
    if ( kernel != NULL ) {
        err = 0;
        i   = 0;
        err |= clSetKernelArg( kernel, i++, sizeof(n            ), &n             );
        err |= clSetKernelArg( kernel, i++, sizeof(dAT          ), &dAT           );
        err |= clSetKernelArg( kernel, i++, sizeof(dAT_offset   ), &dAT_offset    );
        err |= clSetKernelArg( kernel, i++, sizeof(ldda         ), &ldda          );
        err |= clSetKernelArg( kernel, i++, sizeof(nb           ), &nb            );
        err |= clSetKernelArg( kernel, i++, sizeof(d_ipiv       ), &d_ipiv        );
        err |= clSetKernelArg( kernel, i++, sizeof(d_ipiv_offset), &d_ipiv_offset );
        err |= clSetKernelArg( kernel, i++, sizeof(inci         ), &inci          );
        check_error( err );

        err = clEnqueueNDRangeKernel( queue, kernel, 1, NULL, grid, threads, 0, NULL, NULL );
        check_error( err );
    }
}
