/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014
 
       @author Stan Tomov

       @precisions normal z -> s d c

*/
#include "common_magma.h"

extern "C" magma_int_t
magma_zunmqr(
    magma_side_t side, magma_trans_t trans,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaDoubleComplex *a,    magma_int_t lda,
    magmaDoubleComplex *tau,
    magmaDoubleComplex *c,    magma_int_t ldc,
    magmaDoubleComplex *work, magma_int_t lwork,
    magma_queue_t queue,
    magma_int_t *info)
{
/*  -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

    Purpose
    =======
    ZUNMQR overwrites the general complex M-by-N matrix C with

                    SIDE = 'L'     SIDE = 'R'
    TRANS = 'N':      Q * C          C * Q
    TRANS = 'C':      Q**H * C       C * Q**H

    where Q is a complex orthogonal matrix defined as the product of k
    elementary reflectors

          Q = H(1) H(2) . . . H(k)

    as returned by ZGEQRF. Q is of order M if SIDE = 'L' and of order N
    if SIDE = 'R'.

    Arguments
    =========
    SIDE    (input) CHARACTER*1
            = 'L': apply Q or Q**H from the Left;
            = 'R': apply Q or Q**H from the Right.

    TRANS   (input) CHARACTER*1
            = 'N':  No transpose, apply Q;
            = 'C':  Conjugate transpose, apply Q**H.

    M       (input) INTEGER
            The number of rows of the matrix C. M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix C. N >= 0.

    K       (input) INTEGER
            The number of elementary reflectors whose product defines
            the matrix Q.
            If SIDE = 'L', M >= K >= 0;
            if SIDE = 'R', N >= K >= 0.

    A       (input) COMPLEX_16 array, dimension (LDA,K)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,k, as returned by
            ZGEQRF in the first k columns of its array argument A.
            A is modified by the routine but restored on exit.

    LDA     (input) INTEGER
            The leading dimension of the array A.
            If SIDE = 'L', LDA >= max(1,M);
            if SIDE = 'R', LDA >= max(1,N).

    TAU     (input) COMPLEX_16 array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by ZGEQRF.

    C       (input/output) COMPLEX_16 array, dimension (LDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by Q*C or Q**H * C or C * Q**H or C*Q.

    LDC     (input) INTEGER
            The leading dimension of the array C. LDC >= max(1,M).

    WORK    (workspace/output) COMPLEX_16 array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK(0) returns the optimal LWORK.

    LWORK   (input) INTEGER
            The dimension of the array WORK.
            If SIDE = 'L', LWORK >= max(1,N);
            if SIDE = 'R', LWORK >= max(1,M).
            For optimum performance LWORK >= N*NB if SIDE = 'L', and
            LWORK >= M*NB if SIDE = 'R', where NB is the optimal
            blocksize.

            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
    =====================================================================   */
    
    magmaDoubleComplex c_one = MAGMA_Z_ONE;

    /* Allocate work space on the GPU */
    magmaDoubleComplex_ptr dwork, dc;
    magma_zmalloc( &dc, (m)*(n) );
    magma_zmalloc( &dwork, (m + n + 64)*64 );
    
    /* Copy matrix C from the CPU to the GPU */
    magma_zsetmatrix( m, n, c, ldc, dc, 0, m, queue );
    //dc -= (1 + m);
    size_t dc_offset = -(1+m);

    magma_int_t a_offset, c_offset, i__4, lddwork;
    magma_int_t i__;
    magmaDoubleComplex t[2*4160]        /* was [65][64] */;
    magma_int_t i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw;
    int left, notran, lquery;
    magma_int_t iinfo, lwkopt;

    a_offset = 1 + lda;
    a -= a_offset;
    --tau;
    c_offset = 1 + ldc;
    c -= c_offset;

    *info = 0;
    left = (side == MagmaLeft);
    notran = (trans == MagmaNoTrans);
    lquery = (lwork == -1);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }
    if (! left && side != MagmaRight) {
        *info = -1;
    } else if (! notran && trans != MagmaConjTrans) {
        *info = -2;
    } else if (m < 0) {
        *info = -3;
    } else if (n < 0) {
        *info = -4;
    } else if (k < 0 || k > nq) {
        *info = -5;
    } else if (lda < max(1,nq)) {
        *info = -7;
    } else if (ldc < max(1,m)) {
        *info = -10;
    } else if (lwork < max(1,nw) && ! lquery) {
        *info = -12;
    }

    if (*info == 0) {
        /* Determine the block size.  NB may be at most NBMAX, where NBMAX
           is used to define the local array T.    */
        nb = 64;
        lwkopt = max(1,nw) * nb;
        work[0] = MAGMA_Z_MAKE( lwkopt, 0 );
    }

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
      return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0 || k == 0) {
        work[0] = c_one;
        return *info;
    }

    if (nb >= k) {
        /* Use CPU code */
        lapackf77_zunmqr(lapack_const(side), lapack_const(trans), &m, &n, &k, &a[a_offset], &lda, &tau[1],
                         &c[c_offset], &ldc, work, &lwork, &iinfo);
    }
    else {
        /* Use hybrid CPU-GPU code */
        if ( ( left && (! notran) ) ||  ( (! left) && notran ) ) {
            i1 = 1;
            i2 = k;
            i3 = nb;
        } else {
            i1 = (k - 1) / nb * nb + 1;
            i2 = 1;
            i3 = -nb;
        }

        if (left) {
            ni = n;
            jc = 1;
        } else {
            mi = m;
            ic = 1;
        }
        
        for (i__ = i1; i3 < 0 ? i__ >= i2 : i__ <= i2; i__ += i3) {
            ib = min(nb, k - i__ + 1);

            /* Form the triangular factor of the block reflector
               H = H(i) H(i+1) . . . H(i+ib-1) */
            i__4 = nq - i__ + 1;
            lapackf77_zlarft("F", "C", &i__4, &ib, &a[i__ + i__ * lda], &lda,
                             &tau[i__], t, &ib);

            /* 1) Put 0s in the upper triangular part of A;
               2) copy the panel from A to the GPU, and
               3) restore A                                      */
            zpanel_to_q(MagmaUpper, ib, &a[i__ + i__ * lda], lda, t+ib*ib);
            magma_zsetmatrix( i__4, ib, &a[i__ + i__ * lda], lda, dwork, 0, i__4, queue );
            zq_to_panel(MagmaUpper, ib, &a[i__ + i__ * lda], lda, t+ib*ib);

            if (left) {
                /* H or H' is applied to C(i:m,1:n) */
                mi = m - i__ + 1;
                ic = i__;
            }
            else {
                /* H or H' is applied to C(1:m,i:n) */
                ni = n - i__ + 1;
                jc = i__;
            }
            
            if (left)
              lddwork = ni;
            else
              lddwork = mi;

            /* Apply H or H'; First copy T to the GPU */
            magma_zsetmatrix( ib, ib, t, ib, dwork, i__4*ib, ib, queue );
            magma_zlarfb_gpu( side, trans, MagmaForward, MagmaColumnwise,
                              mi, ni, ib,
                              dwork, 0, i__4, dwork, i__4*ib, ib,
                              dc, dc_offset+(ic + jc * m), m,
                              dwork, (i__4*ib + ib*ib), lddwork, queue);
        }

        magma_zgetmatrix( m, n, dc, dc_offset+(1+m), m, &c[c_offset], ldc, queue );
    }
    work[0] = MAGMA_Z_MAKE( lwkopt, 0 );

    //dc += (1 + m);
    magma_free( dc );
    magma_free( dwork );

    return *info;
} /* magma_zunmqr */


