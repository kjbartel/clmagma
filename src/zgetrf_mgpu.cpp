/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @precisions normal z -> s d c

*/
#include <math.h>
#include "common_magma.h"


extern "C" magma_int_t
magma_zgetrf_mgpu(
    magma_int_t ngpu, 
    magma_int_t m, magma_int_t n, 
    magmaDoubleComplex_ptr *d_lA, size_t dlA_offset, magma_int_t ldda,
    magma_int_t *ipiv,
    magma_queue_t *queues,
    magma_int_t *info)
{
/*  -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

    Purpose
    =======
    ZGETRF computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.

    The factorization has the form
        A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 3 BLAS version of the algorithm.

    Arguments
    =========
    NUM_GPUS 
            (input) INTEGER
            The number of GPUS to be used for the factorization.

    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) COMPLEX_16 array on the GPU, dimension (LDDA,N).
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    LDDA     (input) INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    IPIV    (output) INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.
    =====================================================================    */

    magmaDoubleComplex c_one     = MAGMA_Z_ONE;
    magmaDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;

    magma_int_t nb, n_local[MagmaMaxGPUs];
    magma_int_t maxm, mindim;
    magma_int_t d, d2, lddat, ldwork;
    magmaDoubleComplex_ptr d_lAT[MagmaMaxGPUs];
    magmaDoubleComplex_ptr d_panel[MagmaMaxGPUs];
    magmaDoubleComplex *work;

    /* Check arguments */
    *info = 0;
    if (m < 0)
        *info = -2;
    else if (n < 0)
        *info = -3;
    else if (ldda < max(1,m))
        *info = -5;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return *info;

    /* Function Body */
    mindim = min(m, n);
    nb     = magma_get_zgetrf_nb(m);

    if (nb <= 1 || nb >= n) {
        /* Use CPU code. */
        magma_zmalloc_cpu( &work, m * n );
        if ( work == NULL ) {
            *info = MAGMA_ERR_HOST_ALLOC;
            return *info;
        }
        magma_zgetmatrix( m, n, d_lA[0], 0, ldda, work, m, queues[0] );
        lapackf77_zgetrf(&m, &n, work, &m, ipiv, info);
        magma_zsetmatrix( m, n, work, m, d_lA[0], 0, ldda, queues[0] );
        magma_free_cpu(work);
    } else {
        /* Use hybrid blocked code. */
        maxm = ((m + 31)/32)*32;
        if ( ngpu > ceil((double)n/nb) ) {
            printf( " * too many GPUs for the matrix size, using %d GPUs\n", (int) ngpu );
            *info = -1;
            return *info;
        }

        /* allocate workspace for each GPU */
        lddat = (n+nb-1)/nb;                 /* number of block columns         */
        lddat = (lddat+ngpu-1)/ngpu;         /* number of block columns per GPU */
        lddat = nb*lddat;                    /* number of columns per GPU       */
        lddat = ((lddat+31)/32)*32;          /* make it a multiple of 32        */
        for( d=0; d < ngpu; d++ ) {
            /* local-n and local-ld */
            n_local[d] = ((n/nb)/ngpu)*nb;
            if (d < (n/nb)%ngpu)
                n_local[d] += nb;
            else if (d == (n/nb)%ngpu)
                n_local[d] += n%nb;

            /* workspaces */
            if (MAGMA_SUCCESS != magma_zmalloc( &d_panel[d], 3*nb*maxm )) {
                for( d2=0; d2 < d; d2++ ) {
                    magma_free( d_panel[d2] );
                    magma_free( d_lAT[d2]   );
                }
                *info = MAGMA_ERR_DEVICE_ALLOC;
                return *info;
            }

            /* local-matrix storage */
            if (MAGMA_SUCCESS != magma_zmalloc( &d_lAT[d], lddat*maxm )) {
                for( d2=0; d2 <= d; d2++ ) {
                    magma_free( d_panel[d2] );
                }
                for( d2=0; d2 < d; d2++ ) {
                    magma_free( d_lAT[d2] );
                }
                *info = MAGMA_ERR_DEVICE_ALLOC;
                return *info;
            }

            magmablas_ztranspose( m, n_local[d], d_lA[d], 0, ldda, d_lAT[d], 0, lddat, queues[2*d+1] );
        }
        for( d=0; d < ngpu; d++ ) {
            magma_queue_sync(queues[2*d+1]);
        }

        /* cpu workspace */
        ldwork = maxm;
        if (MAGMA_SUCCESS != magma_zmalloc_cpu( &work, ldwork*nb*ngpu )) {
            for( d=0; d < ngpu; d++ ) {
                magma_free( d_panel[d] );
                magma_free( d_lAT[d]   );
            }
            *info = MAGMA_ERR_HOST_ALLOC;
            return *info;
        }

        /* calling multi-gpu interface with allocated workspaces and queues */
        magma_zgetrf2_mgpu(ngpu, m, n, nb, 0, d_lAT, 0, lddat, ipiv, d_panel, 0, work, maxm,
                           queues, info);

        /* clean up */
        for( d=0; d < ngpu; d++ ) {
            /* save on output */
            magmablas_ztranspose( n_local[d], m, d_lAT[d], 0, lddat, d_lA[d], 0, ldda, queues[2*d+1] );
            magma_queue_sync(queues[2*d+1]);
            magma_free( d_lAT[d]   );
            magma_free( d_panel[d] );
        } /* end of for d=1,..,ngpu */
        magma_free_cpu( work );
    }
      
    return *info;       
}
