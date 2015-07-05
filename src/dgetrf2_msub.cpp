/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zgetrf2_msub.cpp normal z -> d, Sat Nov 15 00:21:37 2014

*/
#include <math.h>
#include "common_magma.h"


extern "C" magma_int_t
magma_dgetrf2_msub(
    magma_int_t num_subs, magma_int_t ngpu, 
    magma_int_t m, magma_int_t n, magma_int_t nb, magma_int_t offset,
    magmaDouble_ptr *d_lAT, size_t dlAT_offset, magma_int_t lddat, 
    magma_int_t *ipiv,
    magmaDouble_ptr *d_panel, 
    magmaDouble_ptr *d_lAP, size_t dlAP_offset, 
    double *w, magma_int_t ldw,
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
    DGETRF computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.

    The factorization has the form
       A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 3 BLAS version of the algorithm.
    Use two buffer to send panels..

    Arguments
    =========
    NUM_GPUS 
            (input) INTEGER
            The number of GPUS to be used for the factorization.

    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) DOUBLE_PRECISION array on the GPU, dimension (LDDA,N).
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

#define d_lAT(id,i,j)      d_lAT[(id)], (((offset)+(i)*nb)*lddat + (j)*nb)
#define d_lAT_offset(i, j)              (((offset)+(i)*nb)*lddat + (j)*nb)
#define W(j)     (w +((j)%(1+ngpu))*nb*ldw)

    double c_one     = MAGMA_D_ONE;
    double c_neg_one = MAGMA_D_NEG_ONE;

    magma_int_t tot_subs = num_subs * ngpu;
    magma_int_t block_size = 32;
    magma_int_t iinfo, maxm, mindim;
    magma_int_t i, j, d, dd, rows, cols, s;
    magma_int_t id, j_local, j_local2, nb0, nb1;

    /* local submatrix info */
    magma_int_t ldpan[MagmaMaxSubs * MagmaMaxGPUs],
                n_local[MagmaMaxSubs * MagmaMaxGPUs]; 
    size_t dpanel_local_offset[MagmaMaxSubs * MagmaMaxGPUs];
    magmaDouble_ptr dpanel_local[MagmaMaxSubs * MagmaMaxGPUs];

    /* Check arguments */
    *info = 0;
    if (m < 0)
        *info = -2;
    else if (n < 0)
        *info = -3;
    else if (tot_subs*lddat < max(1,n))
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
    if (tot_subs > ceil((double)n/nb)) {
        *info = -1;
        return *info;
    }
    
    else {
        /* Use hybrid blocked code. */
        maxm  = ((m + block_size-1)/block_size)*block_size;
        
        /* some initializations */
        for (i=0; i < tot_subs; i++) {
            n_local[i] = ((n/nb)/tot_subs)*nb;
            if (i < (n/nb)%tot_subs)
                n_local[i] += nb;
            else if (i == (n/nb)%tot_subs)
                n_local[i] += n%nb;
        }
        
        /* start sending the first panel to cpu */
        nb0 = min(mindim, nb);
        magmablas_dtranspose(  nb0, maxm, d_lAT(0,0,0), lddat, d_lAP[0], dlAP_offset, maxm, queues[2*0+1] );
        magma_dgetmatrix_async( m, nb0,
                                d_lAP[0], dlAP_offset, maxm,
                                W(0), ldw, queues[2*0+1], NULL );
        clFlush(queues[2*0+1]);
        /* ------------------------------------------------------------------------------------- */
        
        s = mindim / nb;
        for (j=0; j < s; j++) {
            /* Set the submatrix ID that holds the current panel */
            id = j%tot_subs;
        
            /* Set the local index where the current panel is */
            j_local = j/tot_subs;
            // cols for gpu panel
            cols  = maxm - j*nb;
            // rows for cpu panel
            rows  = m - j*nb;
        
            /* synchrnoize j-th panel from id-th gpu into work */
            magma_queue_sync( queues[2*(id%ngpu)+1] );
        
            /* j-th panel factorization */
            lapackf77_dgetrf( &rows, &nb, W(j), &ldw, ipiv+j*nb, &iinfo);
            if ((*info == 0) && (iinfo > 0)) {
                *info = iinfo + j*nb;
                //break;
            }
        
            /* start sending the panel to all the gpus */
            d = (j+1)%ngpu;
            for (dd=0; dd < ngpu; dd++) {
                magma_dsetmatrix_async( rows, nb,
                                        W(j), ldw,
                                        d_lAP[d], dlAP_offset+(j%(2+ngpu))*nb*maxm, maxm, 
                                        queues[2*d+1], NULL );
                d = (d+1)%ngpu;
            }
            /* apply the pivoting */
            for( i=j*nb; i < j*nb + nb; ++i ) {
                ipiv[i] += j*nb;
            }
            d = (j+1)%tot_subs;
            for (dd=0; dd < tot_subs; dd++) {
                magmablas_dlaswp( lddat, d_lAT(d,0,0), lddat, j*nb + 1, j*nb + nb, ipiv, 1, queues[2*(d%ngpu)] );
                d = (d+1)%tot_subs;
            }
        
            /* update the trailing-matrix/look-ahead */
            d = (j+1)%tot_subs;
            for (dd=0; dd < tot_subs; dd++) {
                /* storage for panel */
                if (d%ngpu == id%ngpu) {
                    /* the panel belond to this gpu */
                    dpanel_local[d] = d_lAT[id];
                    dpanel_local_offset[d] = d_lAT_offset(j, j_local);
                    ldpan[d] = lddat;
                    /* next column */
                    j_local2 = j_local;
                    if ( d <= id )
                        j_local2++;
                } else {
                    /* the panel belong to another gpu */
                    dpanel_local[d] = d_panel[d%ngpu];  
                    dpanel_local_offset[d] = (j%(2+ngpu))*nb*maxm;
                    ldpan[d] = nb;
                    /* next column */
                    j_local2 = j_local;
                    if ( d < id )
                        j_local2++;
                }
                /* the size of the next column */
                if (s > (j+1)) {
                    nb0 = nb;
                } else {
                    nb0 = n_local[d]-nb*(s/tot_subs);
                    if (d < s%tot_subs)
                        nb0 -= nb;
                }
                if (d == (j+1)%tot_subs) {
                    /* owns the next column, look-ahead the column */
                    nb1 = nb0;
                } else {
                    /* update the entire trailing matrix */
                    nb1 = n_local[d] - j_local2*nb;
                }
                
                /* gpu updating the trailing matrix */
                if (d == (j+1)%tot_subs) { /* look-ahead, this is executed first (j.e., dd=0)  */
                    magma_queue_sync(queues[2*(d%ngpu)]);   /* pivoting done? (overwrite with panel) */
                    magmablas_dtranspose( cols, nb,
                                          d_lAP[d%ngpu], dlAP_offset+(j%(2+ngpu))*nb*maxm, maxm,
                                          dpanel_local[d], dpanel_local_offset[d], ldpan[d], 
                                          queues[2*(d%ngpu)+1] );
                    magma_queue_sync(queues[2*(d%ngpu)+1]); /* panel arrived and transposed for remaining update ? */
        
                    magma_dtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 nb1, nb, c_one,
                                 dpanel_local[d], dpanel_local_offset[d], ldpan[d],
                                 d_lAT(d, j, j_local2), lddat, queues[2*(d%ngpu)+1]);
        
                    magma_dgemm( MagmaNoTrans, MagmaNoTrans, 
                                 nb1, m-(j+1)*nb, nb, 
                                 c_neg_one, d_lAT(d, j,   j_local2),         lddat,
                                            dpanel_local[d], dpanel_local_offset[d]+nb*ldpan[d], ldpan[d], 
                                 c_one,     d_lAT(d, j+1, j_local2),         lddat,
                                 queues[2*(d%ngpu)+1]);
                } else { /* no look-ahead */
                    if (dd < ngpu) {
                        /* synch and transpose only the first time */
                        magma_queue_sync(queues[2*(d%ngpu)+1]); /* panel arrived? */
                        magmablas_dtranspose( cols, nb,
                                              d_lAP[d%ngpu], dlAP_offset+(j%(2+ngpu))*nb*maxm, maxm,
                                              dpanel_local[d], dpanel_local_offset[d], ldpan[d], 
                                              queues[2*(d%ngpu)] );
                    }
        
                    magma_dtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 nb1, nb, c_one,
                                 dpanel_local[d], dpanel_local_offset[d], ldpan[d],
                                 d_lAT(d, j, j_local2), lddat, queues[2*(d%ngpu)]);
                
                    magma_dgemm( MagmaNoTrans, MagmaNoTrans, 
                                 nb1, m-(j+1)*nb, nb, 
                                 c_neg_one, d_lAT(d, j,   j_local2),         lddat,
                                            dpanel_local[d], dpanel_local_offset[d]+nb*ldpan[d], ldpan[d], 
                                 c_one,     d_lAT(d, j+1, j_local2),         lddat,
                                 queues[2*(d%ngpu)]);    
                }
                if (d == (j+1)%tot_subs) {
                    /* Set the local index where the current panel is */
                    int loff    = j+1;
                    int j_local = (j+1)/tot_subs;
                    int ldda    = maxm - (j+1)*nb;
                    int cols    = m - (j+1)*nb;
                    nb0 = min(nb, mindim - (j+1)*nb); /* size of the diagonal block */
                    
                    if (nb0 > 0) {
                        /* transpose the panel for sending it to cpu */
                        magmablas_dtranspose( nb0, ldda,
                                              d_lAT(d,loff,j_local), lddat,
                                              d_lAP[d%ngpu], dlAP_offset + ((j+1)%(2+ngpu))*nb*maxm, ldda, 
                                              queues[2*(d%ngpu)+1] );
                  
                        /* send the panel to cpu */
                        magma_dgetmatrix_async( cols, nb0, 
                                                d_lAP[d%ngpu], dlAP_offset + ((j+1)%(2+ngpu))*nb*maxm, ldda, 
                                                W(j+1), ldw, queues[2*(d%ngpu)+1], NULL );
                    }
                } else {
                    //trace_gpu_end( d, 0 );
                }
                d = (d+1)%tot_subs;
            }
        
            /* update the remaining matrix by gpu owning the next panel */
            if ((j+1) < s) {
                d = (j+1)%tot_subs;
                int j_local = (j+1)/tot_subs;
                int rows  = m - (j+1)*nb;
                
                magma_dtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             n_local[d] - (j_local+1)*nb, nb, 
                             c_one, dpanel_local[d], dpanel_local_offset[d], ldpan[d], 
                                    d_lAT(d,j,j_local+1), lddat, queues[2*(d%ngpu)] );
                    
                magma_dgemm( MagmaNoTrans, MagmaNoTrans, 
                             n_local[d]-(j_local+1)*nb, rows, nb, 
                             c_neg_one, d_lAT(d,j,j_local+1), lddat, 
                                        dpanel_local[d], dpanel_local_offset[d]+nb*ldpan[d], ldpan[d], 
                             c_one,     d_lAT(d,j+1,  j_local+1), lddat, queues[2*(d%ngpu)] );
            }
        } /* end of for j=1..s */
        /* ------------------------------------------------------------------------------ */
        
        /* Set the GPU number that holds the last panel */
        id = s%tot_subs;
        
        /* Set the local index where the last panel is */
        j_local = s/tot_subs;
        
        /* size of the last diagonal-block */
        nb0 = min(m - s*nb, n - s*nb);
        rows = m    - s*nb;
        cols = maxm - s*nb;
        
        if (nb0 > 0) {
        
            /* wait for the last panel on cpu */
            magma_queue_sync( queues[2*(id%ngpu)+1] );
            
            /* factor on cpu */
            lapackf77_dgetrf( &rows, &nb0, W(s), &ldw, ipiv+s*nb, &iinfo );
            if ( (*info == 0) && (iinfo > 0) )
                *info = iinfo + s*nb;
        
            /* send the factor to gpus */
            for (d=0; d < ngpu; d++) {
                magma_dsetmatrix_async( rows, nb0, W(s), ldw,
                                        d_lAP[d], dlAP_offset+(s%(2+ngpu))*nb*maxm, cols, 
                                        queues[2*d+1], NULL );
            }
        
            for( i=s*nb; i < s*nb + nb0; ++i ) {
                ipiv[i] += s*nb;
            }
            for (d=0; d < tot_subs; d++) {
                magmablas_dlaswp( lddat, d_lAT(d,0,0), lddat, s*nb + 1, s*nb + nb0, ipiv, 1, queues[2*(d%ngpu)] );
            }
        
            d = id;
            for (dd=0; dd < tot_subs; dd++) {
                /* wait for the pivoting to be done */
                if (dd < ngpu) {
                    /* synch only the first time */
                    magma_queue_sync( queues[2*(d%ngpu)] );
                }
        
                j_local2 = j_local;
                if (d%ngpu == id%ngpu) {
                    /* the panel belond to this gpu */
                    dpanel_local[d] = d_lAT[id];
                    dpanel_local_offset[d] = d_lAT_offset(s, j_local);
                    if (dd < ngpu) {
                        magmablas_dtranspose( rows, nb0,
                                              d_lAP[d%ngpu], dlAP_offset+(s%(2+ngpu))*nb*maxm, cols, 
                                              dpanel_local[d], dpanel_local_offset[d], lddat, 
                                              queues[2*(d%ngpu)+1] );
                    }
                    /* size of the "extra" block */
                    if (d == id) { /* the last diagonal block belongs to this submatrix */
                        nb1 = nb0;
                    } else if (d < id) {
                        nb1 = nb;
                    } else {
                        nb1 = 0;
                    }
                    if (n_local[d] > j_local*nb+nb1) {
                        magma_dtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                     n_local[d] - (j_local*nb+nb1), nb0, c_one,
                                     dpanel_local[d], dpanel_local_offset[d], lddat, 
                                     d_lAT(d, s, j_local)+nb1, lddat, queues[2*(d%ngpu)+1]);
                    }
                } else if (n_local[d] > j_local2*nb) {
                    /* the panel belong to another gpu */
                    dpanel_local[d] = d_panel[d%ngpu];
                    dpanel_local_offset[d] = (s%(2+ngpu))*nb*maxm;
        
                    /* next column */
                    if (d < ngpu) {
                        /* transpose only the first time */
                        magmablas_dtranspose( rows, nb0,
                                              d_lAP[d%ngpu], dlAP_offset+(s%(2+ngpu))*nb*maxm, cols, 
                                              dpanel_local[d], dpanel_local_offset[d], nb, 
                                              queues[2*(d%ngpu)+1] );
                    }
                    if (d < id)
                        j_local2++;
                    nb1 = n_local[d] - j_local2*nb;
                    magma_dtrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 nb1, nb0, c_one,
                                 dpanel_local[d], dpanel_local_offset[d], nb, 
                                 d_lAT(d,s,j_local2), lddat, queues[2*(d%ngpu)+1]);
                }
                d = (d+1)%tot_subs;
            }
        } /* if( nb0 > 0 ) */

        /* clean up */
        for (d=0; d < ngpu; d++) {
            magma_queue_sync( queues[2*d] );
            magma_queue_sync( queues[2*d+1] );
        } 
    }
    return *info;
    /* End of MAGMA_DGETRF2_MSUB */
}

#undef d_lAT
