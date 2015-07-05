/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from zgetrf2_mgpu.cpp normal z -> c, Sat Nov 15 00:21:37 2014

*/
#include <math.h>
#include "common_magma.h"

extern "C" magma_int_t
magma_cgetrf2_mgpu(
    magma_int_t num_gpus, 
    magma_int_t m, magma_int_t n, magma_int_t nb, magma_int_t offset,
    magmaFloatComplex_ptr *d_lAT, size_t dlAT_offset, magma_int_t lddat, 
    magma_int_t *ipiv,
    magmaFloatComplex_ptr *d_lAP, size_t dlAP_offset, 
    magmaFloatComplex *w, magma_int_t ldw,
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
    CGETRF computes an LU factorization of a general M-by-N matrix A
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

    A       (input/output) COMPLEX array on the GPU, dimension (LDDA,N).
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
    #define W(j) (w+((j)%num_gpus)*nb*ldw)

    magmaFloatComplex c_one     = MAGMA_C_ONE;
    magmaFloatComplex c_neg_one = MAGMA_C_NEG_ONE;

    magma_int_t block_size = 32;
    magma_int_t iinfo, n_local[4]; 
    magma_int_t maxm, mindim;
    magma_int_t i, j, d, dd, rows, cols, s, ldpan[4];
    magma_int_t id, i_local, i_local2, nb0, nb1;
    magmaFloatComplex_ptr d_panel[4], panel_local[4];
    size_t d_panel_offset[4];
    size_t panel_local_offset[4];
    //cudaStream_t streaml[4][2];

    /* Check arguments */
    *info = 0;
    if (m < 0)
    *info = -2;
    else if (n < 0)
    *info = -3;
    else if (num_gpus*lddat < max(1,n))
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
    //nb     = magma_get_cgetrf_nb(m);
    if( num_gpus > ceil((float)n/nb) ) {
      *info = -1;
      return *info;
    }
    
    else{
        /* Use hybrid blocked code. */
        maxm  = ((m + block_size-1)/block_size)*block_size;
        
        /* some initializations */
        for( d=0; d < num_gpus; d++){
          n_local[d] = ((n/nb)/num_gpus)*nb;
          if (d < (n/nb)%num_gpus)
             n_local[d] += nb;
          else if (d == (n/nb)%num_gpus)
             n_local[d] += n%nb;
        
          /* workspaces */
          //d_panel[d] = &(d_lAP[d][nb*maxm]);   /* temporary panel storage */
          d_panel[d] = d_lAP[d];
          d_panel_offset[d] = nb*maxm;
        
        }
        
        /* start sending the panel to cpu */
        nb0 = min(mindim, nb);
        magmablas_ctranspose( nb0, maxm, d_lAT(0,0,0), lddat, d_lAP[0], 0, maxm, queues[2*0+1] );
        magma_cgetmatrix_async( m, nb0,
                                d_lAP[0], 0, maxm,
                                W(0), ldw, queues[2*0+1], NULL );
        clFlush(queues[2*0+1]);
        /* ------------------------------------------------------------------------------------- */
        s = mindim / nb;
        for( j=0; j < s; j++ )
        {
            /* Set the GPU number that holds the current panel */
            id = j%num_gpus;
        
            /* Set the local index where the current panel is */
            i_local = j/num_gpus;
            // cols for gpu panel
            cols  = maxm - j*nb;
            // rows for cpu panel
            rows  = m - j*nb;
        
            /* synchrnoize j-th panel from id-th gpu into work */
            magma_queue_sync( queues[2*id+1] );
        
            /* j-th panel factorization */
            lapackf77_cgetrf( &rows, &nb, W(j), &ldw, ipiv+j*nb, &iinfo);
            
            if ( (*info == 0) && (iinfo > 0) ) {
                *info = iinfo + j*nb;
                //break;
            }
        
            /* start sending the panel to all the gpus */
            d = (j+1)%num_gpus;
            for( dd=0; dd<num_gpus; dd++ ) {
                magma_csetmatrix_async( rows, nb,
                                        W(j), ldw,
                                        d_lAP[d], dlAP_offset, cols, queues[2*d+1], NULL );
                d = (d+1)%num_gpus;
            }
        
            /* apply the pivoting */
            for( i=j*nb; i < j*nb + nb; ++i ) {
                ipiv[i] += j*nb;
            }
            d = (j+1)%num_gpus;
            for( dd=0; dd<num_gpus; dd++ ) {
                magmablas_claswp( lddat, d_lAT(d,0,0), lddat, j*nb + 1, j*nb + nb, ipiv, 1, queues[2*d] );
                d = (d+1)%num_gpus;
            }
        
        
            /* update the trailing-matrix/look-ahead */
            d = (j+1)%num_gpus;
            for( dd=0; dd<num_gpus; dd++ ) {
                /* storage for panel */
                if( d == id ) {
                    /* the panel belond to this gpu */
                    //panel_local[d] = d_lAT(d,j,i_local);
                    panel_local[d] = d_lAT[d];
                    panel_local_offset[d] = d_lAT_offset(j, i_local);
                    ldpan[d] = lddat;
                    /* next column */
                    i_local2 = i_local+1;
                } else {
                    /* the panel belong to another gpu */
                    //panel_local[d] = &d_panel[d][(j%2)*nb*maxm];
                    panel_local[d] = d_panel[d];  
                    panel_local_offset[d] = d_panel_offset[d] + (j%2)*nb*maxm;
                    //panel_local[d] = d_panel[d];
                    ldpan[d] = nb;
                    /* next column */
                    i_local2 = i_local;
                    if( d < id ) i_local2 ++;
                }
                /* the size of the next column */
                if ( s > (j+1) ) {
                    nb0 = nb;
                } else {
                    nb0 = n_local[d]-nb*(s/num_gpus);
                    if( d < s%num_gpus ) nb0 -= nb;
                }
                if( d == (j+1)%num_gpus) {
                    /* owns the next column, look-ahead the column */
                    nb1 = nb0;
                    /* make sure all the pivoting has been applied */
                    //magma_queue_sync(queues[2*d]);
                } else {
                    /* update the entire trailing matrix */
                    nb1 = n_local[d] - i_local2*nb;
        
                    /* synchronization to make sure panel arrived on gpu */
                    //magma_queue_sync(queues[2*d+1]);
                }
                
               /* 
                magma_queue_sync(queues[2*d]);
                magma_queue_sync(queues[2*d+1]);
               */
        
                //magmablas_ctranspose( cols, nb, d_lAP[d], 0, cols, panel_local[d], panel_local_offset[d], ldpan[d], queues[2*d] );
                
                /* gpu updating the trailing matrix */
                if(d == (j+1)%num_gpus){
                magma_queue_sync(queues[2*d]);
                magmablas_ctranspose( cols, nb, d_lAP[d], 0, cols, panel_local[d], panel_local_offset[d], ldpan[d], queues[2*d+1] );
                magma_queue_sync(queues[2*d+1]);
                magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             nb1, nb, c_one,
                             panel_local[d], panel_local_offset[d], ldpan[d],
                             d_lAT(d, j, i_local2), lddat, queues[2*d+1]);
                
                magma_cgemm( MagmaNoTrans, MagmaNoTrans, 
                             nb1, m-(j+1)*nb, nb, 
                             c_neg_one, d_lAT(d, j,   i_local2),         lddat,
                              panel_local[d], panel_local_offset[d]+nb*ldpan[d], ldpan[d], 
                              c_one,     d_lAT(d, j+1, i_local2),         lddat,
                              queues[2*d+1]);
                }else{
                magma_queue_sync(queues[2*d+1]);
                magmablas_ctranspose( cols, nb, d_lAP[d], 0, cols, panel_local[d], panel_local_offset[d], ldpan[d], queues[2*d] );
                magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                             nb1, nb, c_one,
                             panel_local[d], panel_local_offset[d], ldpan[d],
                             d_lAT(d, j, i_local2), lddat, queues[2*d]);
                
                magma_cgemm( MagmaNoTrans, MagmaNoTrans, 
                             nb1, m-(j+1)*nb, nb, 
                             c_neg_one, d_lAT(d, j,   i_local2),         lddat,
                              panel_local[d], panel_local_offset[d]+nb*ldpan[d], ldpan[d], 
                              c_one,     d_lAT(d, j+1, i_local2),         lddat,
                              queues[2*d]);    
                }
                
                if( d == (j+1)%num_gpus ) 
                {
                    /* Set the local index where the current panel is */
                    int loff    = j+1;
                    int i_local = (j+1)/num_gpus;
                    int ldda    = maxm - (j+1)*nb;
                    int cols    = m - (j+1)*nb;
                    nb0 = min(nb, mindim - (j+1)*nb); /* size of the diagonal block */
                    
                    if( nb0 > 0 ) {
                        /* transpose the panel for sending it to cpu */
                        magmablas_ctranspose(   nb0, ldda, d_lAT(d,loff,i_local), lddat, d_lAP[d], 0, ldda, queues[2*d+1] );
                  
                        //clFinish(queues[2*d+1]);
                        /* send the panel to cpu */
                        magma_cgetmatrix_async( cols, nb0, 
                                                d_lAP[d], 0, ldda, 
                                                W(j+1), ldw, queues[2*d+1], NULL );
                    }
                } else {
                      //trace_gpu_end( d, 0 );
                }
                d = (d+1)%num_gpus;
            }
        
            /* update the remaining matrix by gpu owning the next panel */
            if( (j+1) < s ) {
                int i_local = (j+1)/num_gpus;
                int rows  = m - (j+1)*nb;
                
                d = (j+1)%num_gpus;
                
                magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 n_local[d] - (i_local+1)*nb, nb, 
                                 c_one, panel_local[d], panel_local_offset[d], ldpan[d], 
                                 d_lAT(d,j,i_local+1),  lddat, queues[2*d] );
                    
                magma_cgemm( MagmaNoTrans, MagmaNoTrans, 
                              n_local[d]-(i_local+1)*nb, rows, nb, 
                              c_neg_one, d_lAT(d,j,i_local+1),            lddat, 
                              panel_local[d], panel_local_offset[d]+nb*ldpan[d], ldpan[d], 
                              c_one,     d_lAT(d,j+1,  i_local+1),        lddat, queues[2*d] );
            }
        } /* end of for j=1..s */
        /* ------------------------------------------------------------------------------ */
        
        /* Set the GPU number that holds the last panel */
        id = s%num_gpus;
        
        /* Set the local index where the last panel is */
        i_local = s/num_gpus;
        
        /* size of the last diagonal-block */
        nb0 = min(m - s*nb, n - s*nb);
        rows = m    - s*nb;
        cols = maxm - s*nb;
        
        if( nb0 > 0 ) {
        
            /* wait for the last panel on cpu */
            magma_queue_sync( queues[2*id+1]);
            
            /* factor on cpu */
            lapackf77_cgetrf( &rows, &nb0, W(s), &ldw, ipiv+s*nb, &iinfo);
            if ( (*info == 0) && (iinfo > 0) )
                *info = iinfo + s*nb;
        
            /* send the factor to gpus */
            for( d=0; d<num_gpus; d++ ) {
                i_local2 = i_local;
                if( d < id ) i_local2 ++;
        
                if( d == id || n_local[d] > i_local2*nb ) {
                    magma_csetmatrix_async( rows, nb0,
                                            W(s), ldw,
                                            d_lAP[d], 0, cols, queues[2*d+1], NULL );
                }
            }
        
            /* apply the pivoting */
            for( i=s*nb; i < s*nb + nb0; ++i ) {
                ipiv[i] += s*nb;
            }
            for( d=0; d<num_gpus; d++ ) {
                magmablas_claswp( lddat, d_lAT(d,0,0), lddat, s*nb + 1, s*nb + nb0, ipiv, 1, queues[2*d] );
            }
        
            for( d=0; d<num_gpus; d++ ) {
                //magma_queue_sync( queues[2*d+1] );
                /* wait for the pivoting to be done */
                magma_queue_sync( queues[2*d] );
        
                i_local2 = i_local;
                if( d < id ) i_local2++;
                if( d == id ) {
                    /* the panel belond to this gpu */
                    //panel_local[d] = d_lAT(d,s,i_local);
                    panel_local[d] = d_lAT[d];
                    panel_local_offset[d] = d_lAT_offset(s, i_local);
        
                    /* next column */
                    nb1 = n_local[d] - i_local*nb-nb0;
        
                    magmablas_ctranspose( rows, nb0, d_lAP[d], 0, cols, panel_local[d], panel_local_offset[d], lddat, queues[2*d+1] );
        
                    if( nb1 > 0 ){
                        magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                     nb1, nb0, c_one,
                                     panel_local[d], panel_local_offset[d], lddat, 
                                     d_lAT[d], d_lAT_offset(s,i_local)+nb0, lddat, queues[2*d+1]);
                    }
                  } else if( n_local[d] > i_local2*nb ) {
                      /* the panel belong to another gpu */
                      //panel_local[d] = &d_panel[d][(s%2)*nb*maxm];
                      panel_local[d] = d_panel[d];
                      panel_local_offset[d] = d_panel_offset[d] + (s%2)*nb*maxm;
                      //panel_local[d] = d_panel[d];
        
                    /* next column */
                    nb1 = n_local[d] - i_local2*nb;
        
                    magmablas_ctranspose( rows, nb0, d_lAP[d], 0, cols, panel_local[d], panel_local_offset[d], nb, queues[2*d+1] );
                    //cublasCtrsm
                    magma_ctrsm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaUnit, 
                                 nb1, nb0, c_one,
                                 panel_local[d], panel_local_offset[d],    nb, 
                                 d_lAT(d,s,i_local2), lddat, queues[2*d+1]);
                  }
            }
        } /* if( nb0 > 0 ) */
        
        /* clean up */
        for( d=0; d<num_gpus; d++ ) {
                magma_queue_sync( queues[2*d] );
                magma_queue_sync( queues[2*d+1] );
                //magma_queue_destroy(streaml[d][0]);
                //magma_queue_destroy(streaml[d][1]);
        } 
    }
    return *info;
} /* magma_cgetrf2_mgpu */

#undef d_lAT
