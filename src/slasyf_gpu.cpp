/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @author Ichitaro Yamazaki
       @author Stan Tomov      

       @generated from zlahef_gpu.cpp normal z -> s, Sat Nov 15 00:21:37 2014
*/

#include "common_magma.h"

#define PRECISION_s

/**
    Purpose
    =======

    SLASYF computes a partial factorization of a real symmetric
    matrix A using the Bunch-Kaufman diagonal pivoting method. The
    partial factorization has the form:

    A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  if UPLO = 'U', or:
          ( 0  U22 ) (  0   D  ) ( U12' U22' )

    A  =  ( L11  0 ) (  D   0  ) ( L11' L21' )  if UPLO = 'L'
          ( L21  I ) (  0  A22 ) (  0    I   )

    where the order of D is at most NB. The actual order is returned in
    the argument KB, and is either NB or NB-1, or N if N <= NB.
    Note that U' denotes the conjugate transpose of U.

    SLASYF is an auxiliary routine called by SSYTRF. It uses blocked code
    (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
    A22 (if UPLO = 'L').

    Arguments
    ---------
    @param[in]
    UPLO    CHARACTER
            Specifies whether the upper or lower triangular part of the
            symmetric matrix A is stored:
      -     = 'U':  Upper triangular
      -     = 'L':  Lower triangular

    @param[in]
    N       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    NB      INTEGER
            The maximum number of columns of the matrix A that should be
            factored.  NB should be at least 2 to allow for 2-by-2 pivot
            blocks.

    @param[out]
    KB      INTEGER
            The number of columns of A that were actually factored.
            KB is either NB-1 or NB, or N if N <= NB.

    @param[in,out]
    A       REAL array, dimension (LDA,N)
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading
            n-by-n upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = 'L', the
            leading n-by-n lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
            On exit, A contains details of the partial factorization.

    @param[in]
    LDA     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    ipiv    INTEGER array, dimension (N)
            Details of the interchanges and the block structure of D.
            If UPLO = 'U', only the last KB elements of ipiv are set;
            if UPLO = 'L', only the first KB elements are set.
    \n
            If ipiv(k) > 0, then rows and columns k and ipiv(k) were
            interchanged and D(k,k) is a 1-by-1 diagonal block.
            If UPLO = 'U' and ipiv(k) = ipiv(k-1) < 0, then rows and
            columns k-1 and -ipiv(k) were interchanged and D(k-1:k,k-1:k)
            is a 2-by-2 diagonal block.  If UPLO = 'L' and ipiv(k) =
            ipiv(k+1) < 0, then rows and columns k+1 and -ipiv(k) were
            interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
 
    @param[out]
    W       (workspace) REAL array, dimension (LDW,NB)
 
    @param[in]
    LDW     INTEGER
            The leading dimension of the array W.  LDW >= max(1,N).

    @param[out]
    INFO    INTEGER
      -     = 0: successful exit
      -     > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
                 has been completed, but the block diagonal matrix D is
                 exactly singular.
  
    @ingroup magma_ssytrf_comp
    ********************************************************************/
extern "C" magma_int_t
magma_slasyf_gpu(
    magma_uplo_t uplo, magma_int_t n, magma_int_t nb, magma_int_t *kb, 
    float *hA, magma_int_t lda,  
    magmaFloat_ptr dA, size_t dA_offset, magma_int_t ldda, 
    magma_int_t *ipiv, 
    magmaFloat_ptr dW, size_t dW_offset, magma_int_t lddw, 
    magma_queue_t queue,
    magma_int_t *info) 
{
    /* .. Parameters .. */
    float d_one   = 1.0;
    float d_zero  = 0.0;
    float d_eight = 8.0;
    float d_seven = 7.0;
    #if defined(PRECISION_c)
    float  f_zero =  0.0;
    #endif
    float c_one  =  MAGMA_S_ONE;
    float c_mone = -MAGMA_S_ONE;
    magma_int_t upper = (uplo == MagmaUpper);
    magma_int_t ione = 1;
  
    /* .. Local Scalars .. */
    magma_int_t imax = 0, jmax = 0, kk, kkW, kp, kstep, iinfo;
    float   abs_akk, alpha, colmax, R1, rowmax;
    float Zimax, Z;

    #define dA(i, j)  dA, dA_offset + (j)*ldda  + (i)
    #define dW(i, j)  dW, dW_offset + (j)*lddw  + (i)
    #define  A(i, j) (hA + (j)*lda   + (i))

    /* .. Executable Statements .. */
    *info = 0;

    /* Initialize alpha for use in choosing pivot block size. */
    alpha = ( d_one+sqrt( d_seven ) ) / d_eight;

    magma_event_t event = NULL;
    if( upper ) {
       /* Factorize the trailing columns of A using the upper triangle
          of A and working backwards, and compute the matrix W = U12*D
          for use in updating A11 (note that conjg(W) is actually stored)
 
          K is the main loop index, decreasing from N in steps of 1 or 2
 
          KW is the column of W which corresponds to column K of A   */
       int k, kw = 0;
       for (k = n-1; k+1 > max(n-nb+1, nb); k -= kstep) {
           kw = nb - (n-k);
           /* Copy column K of A to column KW of W and update it */

           magma_scopy( k+1, dA( 0, k ), 1, dW( 0, kw ), 1, queue );

           // set imaginary part of diagonal to be zero
           #if defined(PRECISION_z)
           magma_dsetvector_async( 1, &d_zero, 1, 
                                   dW, 2*(k+ kw*lddw+dW_offset)+1, 1, queue, &event);
           magma_queue_sync( queue );
           #elif defined(PRECISION_c)
           magma_ssetvector_async( 1, &f_zero, 1, 
                                   dW, 2*(k+ kw*lddw+dW_offset)+1, 1, queue, &event);
           magma_queue_sync( queue );
           #endif

           if (k+1 < n) {
                magma_sgemv( MagmaNoTrans, k+1, n-(k+1), c_mone, dA( 0, k+1 ), ldda,
                             dW( k, kw+1 ), lddw, c_one, dW( 0, kw ), ione, queue );

                // set imaginary part of diagonal to be zero
                #if defined(PRECISION_z)
                magma_dsetvector_async( 1, &d_zero, 1,
                                        dW, 2*(k+ kw*lddw+dW_offset)+1, 1, queue, &event );
                magma_queue_sync( queue );
                #elif defined(PRECISION_c)
                magma_ssetvector_async( 1, &f_zero, 1, 
                                        dW, 2*(k+ kw*lddw+dW_offset)+1, 1, queue, &event );
                magma_queue_sync( queue );
                #endif
           }

           kstep = 1;

           /* Determine rows and columns to be interchanged and whether
              a 1-by-1 or 2-by-2 pivot block will be used */
           magma_sgetvector_async( 1, dW( k, kw ), 1, &Z, 1, queue, &event );
           magma_queue_sync( queue );
           abs_akk = fabs( MAGMA_S_REAL( Z ) );

           /* imax is the row-index of the largest off-diagonal element in
              column K, and colmax is its absolute value */
           if( k > 0 ) {
                // magma is one-base
                imax = magma_isamax( k, dW( 0, kw ), 1, queue ) - 1;
                magma_sgetvector( 1, dW( imax, kw ), 1, &Z, 1, queue );
                colmax = MAGMA_S_ABS1( Z );
           } else {
                colmax = d_zero;
           }
           if( max( abs_akk, colmax ) == 0.0 ) {
            
                /* Column K is zero: set INFO and continue */
                if ( *info == 0 ) *info = k;

                kp = k;

                #if defined(PRECISION_z)
                magma_dsetvector_async( 1, &d_zero, 1, 
                                        dA, 2*(k+ k*ldda+dA_offset)+1, 1, queue, &event );
                magma_queue_sync( queue );
                #elif defined(PRECISION_c)
                magma_ssetvector_async( 1, &f_zero, 1, 
                                        dA, 2*(k+ k*ldda+dA_offset)+1, 1, queue, &event );
                magma_queue_sync( queue );
                #endif
           } else {
            if( abs_akk >= alpha*colmax ) {

              /* no interchange, use 1-by-1 pivot block */
              kp = k;
            } else {

              /* Copy column imax to column KW-1 of W and update it */
              magma_scopy( imax+1, dA( 0, imax ), 1, dW( 0, kw-1 ), 1, queue );
              #if defined(PRECISION_z)
              magma_dsetvector_async( 1, &d_zero, 1, 
                                      dW, 2*(imax+ (kw-1)*lddw+dW_offset)+1, 1, queue, &event );
              #elif defined(PRECISION_c)
              magma_ssetvector_async( 1, &f_zero, 1, 
                                      dW, 2*(imax+ (kw-1)*lddw+dW_offset)+1, 1, queue, &event );
              #endif

              #if defined(PRECISION_z) || defined(PRECISION_c)
              magmablas_slacpy_cnjg( k-imax, dA(imax,imax+1), ldda, dW(imax+1,kw-1), 1, queue );
              #else
              magma_scopy( k-imax, dA(imax,imax+1), ldda, dW(imax+1,kw-1), 1, queue );
              #endif
              if( k+1 < n ) {
                 magma_sgemv( MagmaNoTrans, k+1, n-(k+1), c_mone,
                              dA( 0, k+1 ), ldda, dW( imax, kw+1 ), lddw,
                              c_one, dW( 0, kw-1 ), ione, queue );

                 #if defined(PRECISION_z)
                 magma_dsetvector_async( 1, &d_zero, 1, 
                                         dW, 2*(imax+ (kw-1)*lddw+dW_offset)+1, 1, queue, &event );
                 #elif defined(PRECISION_c)
                 magma_ssetvector_async( 1, &f_zero, 1, 
                                         dW, 2*(imax+ (kw-1)*lddw+dW_offset)+1, 1, queue, &event );
                 #endif
              }
              magma_sgetvector_async( 1, dW( imax, kw-1 ), 1, &Zimax, 1, queue, &event );
              magma_queue_sync( queue );

              /* jmax is the column-index of the largest off-diagonal
                element in row imax, and rowmax is its absolute value */
              jmax = imax + magma_isamax( k-imax, dW( imax+1, kw-1 ), 1, queue );
              magma_sgetvector( 1, dW( jmax, kw-1 ), 1, &Z, 1, queue );
              rowmax = MAGMA_S_ABS1( Z );
              if ( imax > 0 ) {
               // magma is one-base                                                                           
               jmax = magma_isamax( imax, dW( 0, kw-1 ), 1, queue ) - 1;
               magma_sgetvector( 1, dW( jmax, kw-1 ), 1, &Z, 1, queue );
               rowmax = max( rowmax, MAGMA_S_ABS1( Z  ) );
              }

              if( abs_akk >= alpha*colmax*( colmax / rowmax ) ) {
 
                     /* no interchange, use 1-by-1 pivot block */
                     kp = k;
              } else if ( fabs( MAGMA_S_REAL( Zimax ) ) >= alpha*rowmax ) {

                     /* interchange rows and columns K and imax, use 1-by-1
                        pivot block */
                     kp = imax;

                     /* copy column KW-1 of W to column KW */
                     magma_scopy( k+1, dW( 0, kw-1 ), 1, dW( 0, kw ), 1, queue );
                } else {

                     /* interchange rows and columns K-1 and imax, use 2-by-2
                        pivot block */
                     kp = imax;
                     kstep = 2;
                }
           }
           kk = k - kstep + 1;
           kkW = nb - (n - kk);

           /* Updated column kp is already stored in column kkW of W */
           if( kp != kk ) {

             /* Interchange rows kk and kp in last kk columns of A and W */
             // note: row-swap A(:,kk)
             magmablas_sswap( n-kk, dA( kk, kk ), ldda, dA( kp, kk ), ldda, queue );
             magmablas_sswap( n-kk, dW( kk, kkW), lddw, dW( kp, kkW), lddw, queue );

             /* Copy non-updated column kk to column kp */
             #if defined(PRECISION_z) || defined(PRECISION_c)
             magmablas_slacpy_cnjg( kk-kp-1, dA( kp+1, kk ), 1, dA( kp, kp+1 ), ldda, queue );
             #else
             magma_scopy( kk-kp-1, dA( kp+1, kk ), 1, dA( kp, kp+1 ), ldda, queue );
             #endif

             // now A(kp,kk) should be A(kk,kk), and copy to A(kp,kp)
             magma_scopy( kp+1, dA( 0, kk ), 1, dA( 0, kp ), 1, queue );
             #if defined(PRECISION_z)
             magma_dsetvector_async( 1, &d_zero, 1, 
                                     dA, 2*(kp+ kp*ldda+dA_offset)+1, 1, queue, &event );
             magma_queue_sync( queue );
             #elif defined(PRECISION_c)
             magma_ssetvector_async( 1, &f_zero, 1, 
                                     dA, 2*(kp+ kp*ldda+dA_offset)+1, 1, queue, &event );
             #endif
           }
           if( kstep == 1 ) {

                /* 1-by-1 pivot block D(k): column KW of W now holds
                      W(k) = U(k)*D(k)
                      where U(k) is the k-th column of U 
                      Store U(k) in column k of A */
                magma_scopy( k+1, dW( 0, kw ), 1, dA( 0, k ), 1, queue );
                if ( k > 0 ) {
                   magma_sgetvector_async( 1, dA( k, k ), 1, &Z, 1, queue, &event );
                   magma_queue_sync( queue );
                   R1 = d_one / MAGMA_S_REAL( Z );
                   magma_sscal( k, R1, dA( 0, k ), 1, queue );

                   /* Conjugate W(k) */
                   #if defined(PRECISION_z) || defined(PRECISION_c)
                   magmablas_slacpy_cnjg( k, dW( 0, kw ), 1, dW( 0, kw ), 1, queue );
                   #endif
                }
           } else {

            /* 2-by-2 pivot block D(k): columns KW and KW-1 of W now hold
              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
              where U(k) and U(k-1) are the k-th and (k-1)-th columns of U */
            if( k > 1 ) {
                /* Store U(k) and U(k-1) in columns k and k-1 of A */
                magmablas_slascl_2x2( MagmaUpper, 
                                      k-1, dW(0, kw-1), lddw, dA(0,k-1), ldda, &iinfo, queue );
             }

             /* Copy D(k) to A */
             magma_scopymatrix( 2, 2, dW( k-1, kw-1 ), lddw, dA( k-1, k-1 ), ldda, queue );

             /* Conjugate W(k) and W(k-1) */
             #if defined(PRECISION_z) || defined(PRECISION_c)
             magmablas_slacpy_cnjg( k,   dW( 0, kw ),   1, dW( 0, kw ),   1, queue );
             magmablas_slacpy_cnjg( k-1, dW( 0, kw-1 ), 1, dW( 0, kw-1 ), 1, queue );
             #endif
           }
          }

          /* Store details of the interchanges in ipiv */
          if( kstep == 1 ) {
            ipiv[ k ] = 1+kp;
          } else {
            ipiv[ k ] = -(1+kp);
            ipiv[ k-1 ] = -(1+kp);
          }
       }
       /* Update the upper triangle of A11 (= A(1:k,1:k)) as
           A11 := A11 - U12*D*U12' = A11 - U12*W'
          computing blocks of NB columns at a time (note that conjg(W) is
          actually stored) */
       kw = nb - (n-k);
       for (int j = ( k / nb )*nb; j >= 0; j -= nb ) {
         int jb = min( nb, k-j+1 );

         #ifdef SYMMETRIC_UPDATE
         /* Update the upper triangle of the diagonal block */
         for (int jj = j; jj < j + jb; jj++) {
            #if defined(PRECISION_z)
            magma_dsetvector_async( 1, &d_zero, 1, 
                                    dA, 2*(jj+ jj*ldda+dA_offset)+1, 1, queue, &event );
            #elif defined(PRECISION_c)
            magma_ssetvector_async( 1, &f_zero, 1,
                                    dA, 2*(jj+ jj*ldda+dA_offset)+1, 1, queue, &event );
            #endif
            magma_sgemv( MagmaNoTrans, jj-j+1, n-(k+1), c_mone,
                         dA( j, k+1 ), ldda, dW( jj, kw+1 ), lddw, c_one,
                         dA( j, jj ), 1, queue );
            #if defined(PRECISION_z)
            magma_dsetvector_async( 1, &d_zero, 1, 
                                    dA, 2*(jj+ jj*ldda+dA_offset)+1, 1, queue, &event );
            #elif defined(PRECISION_c)
            magma_ssetvector_async( 1, &f_zero, 1, 
                                    dA, 2*(jj+ jj*ldda+dA_offset)+1, 1, queue, &event );
            #endif
         }
         /* Update the rectangular superdiagonal block */
         magma_sgemm( MagmaNoTrans, MagmaTrans, j, jb, n-(k+1),
                      c_mone, dA( 0, k+1 ), ldda, dW( j, kw+1 ), lddw,
                      c_one, dA( 0, j ), ldda, queue );
         #else
         #if defined(PRECISION_z)
         magmablas_dlaset(MagmaUpperLower, 1, jb, 
                          0, 0, dA, 2*(j+ j*ldda+dA_offset)+1, 2*(1+ldda), queue );
         #elif defined(PRECISION_c)
         magmablas_slaset(MagmaUpperLower, 1, jb, 
                          0, 0, dA, 2*(j+ j*ldda+dA_offset)+1, 2*(1+ldda), queue );
         #endif
         magma_sgemm( MagmaNoTrans, MagmaTrans, j+jb, jb, n-(k+1),
                      c_mone, dA( 0, k+1 ),  ldda,
                      dW( j, kw+1 ), lddw,
                      c_one,  dA( 0, j ),    ldda, queue );
         #if defined(PRECISION_z)
         magmablas_dlaset(MagmaUpperLower, 1, jb, 
                          0, 0, dA, 2*(j+ j*ldda+dA_offset)+1, 2*(1+ldda), queue );
         #elif defined(PRECISION_c)
         magmablas_slaset(MagmaUpperLower, 1, jb, 
                          0, 0, dA, 2*(j+ j*ldda+dA_offset)+1, 2*(1+ldda), queue );
         #endif
         #endif
       }

       /* Put U12 in standard form by partially undoing the interchanges in columns k+1:n */
       for (int j = k+1; j < n;)
        {
          int jj = j;
          int jp = ipiv[ j ];
          if( jp < 0 ) {
            jp = -jp;
            j = j + 1;
          }
          j = j + 1;
          jp = jp - 1;
          if( jp != jj && j < n )
            magmablas_sswap( n-j, dA( jp, j ), ldda, dA( jj, j ), ldda, queue );
        }

       // copying the panel back to CPU
       magma_sgetmatrix_async( n, n-(k+1), dA(0,k+1), ldda, A(0,k+1), lda, queue, &event );
       magma_queue_sync( queue );

       /* Set KB to the number of columns factorized */
       *kb = n - (k+1);
 
    } else {
       /* Factorize the leading columns of A using the lower triangle
          of A and working forwards, and compute the matrix W = L21*D
          for use in updating A22 (note that conjg(W) is actually stored)

          K is the main loop index, increasing from 1 in steps of 1 or 2 */

       int k;
       for (k = 0; k < min(nb-1,n); k += kstep) {

           /* Copy column K of A to column K of W and update it */
           /* -------------------------------------------------------------- */
           magma_scopy( n-k, dA( k, k ), 1, dW( k, k ), 1, queue );

           // set imaginary part of diagonal to be zero
           #if defined(PRECISION_z)
           magma_dsetvector_async( 1, &d_zero, 1, 
                                   dW, 2*(k*lddw+k+dW_offset)+1, 1, queue, &event);
           magma_queue_sync( queue );
           #elif defined(PRECISION_c)
           magma_ssetvector_async( 1, &f_zero, 1,
                                   dW, 2*(k*lddw+k+dW_offset)+1, 1, queue, &event); 
           magma_queue_sync( queue );
           #endif
           /* -------------------------------------------------------------- */

           magma_sgemv( MagmaNoTrans, n-k, k, c_mone, dA( k, 0 ), ldda, 
                        dW( k, 0 ), lddw, c_one, dW( k, k ), ione, queue );
           // re-set imaginary part of diagonal to be zero
           #if defined(PRECISION_z)
           magma_dsetvector_async( 1, &d_zero, 1, 
                                   dW, 2*(k*lddw+k+dW_offset)+1, 1, queue, &event );
           magma_queue_sync( queue );
           #elif defined(PRECISION_c)
           magma_ssetvector_async( 1, &f_zero, 1,
                                   dW, 2*(k*lddw+k+dW_offset)+1, 1, queue, &event ); 
           magma_queue_sync( queue );
           #endif

           kstep = 1;

           /* Determine rows and columns to be interchanged and whether
              a 1-by-1 or 2-by-2 pivot block will be used */
           magma_sgetvector_async( 1, dW( k, k ), 1, &Z, 1, queue, &event );
           magma_queue_sync( queue );
           abs_akk = fabs( MAGMA_S_REAL( Z ) );

           /* imax is the row-index of the largest off-diagonal element in
              column K, and colmax is its absolute value */
           if( k < n-1 ) {
               // magmablas is one-base
               imax = k + magma_isamax( n-k-1, dW(k+1,k), 1, queue );

               magma_sgetvector( 1, dW( imax,k ), 1, &Z, 1, queue );
               colmax = MAGMA_S_ABS1( Z );
 
           } else {
               colmax = d_zero;
           }

           if ( max( abs_akk, colmax ) == 0.0 ) {

               /* Column K is zero: set INFO and continue */
               if( *info == 0 ) *info = k;
               kp = k;

               // make sure the imaginary part of diagonal is zero
               #if defined(PRECISION_z)
               magma_dsetvector_async( 1, &d_zero, 1, 
                                       dA, 2*(k*ldda+k+dA_offset)+1, 1, queue, &event );
               magma_queue_sync( queue );
               #elif defined(PRECISION_c)
               magma_ssetvector_async( 1, &f_zero, 1,
                                       dA, 2*(k*ldda+k+dA_offset)+1, 1, queue, &event );
               magma_queue_sync( queue );
               #endif
           } else {
               if ( abs_akk >= alpha*colmax ) {

                   /* no interchange, use 1-by-1 pivot block */

                   kp = k;
               } else {
                   /* Copy column imax to column K+1 of W and update it */
                   #if defined(PRECISION_z) || defined(PRECISION_c)
                   magmablas_slacpy_cnjg( imax-k, dA(imax,k), ldda, dW(k,k+1), 1, queue );
                   #else 
                   magma_scopy( imax-k, dA( imax, k ), ldda, dW( k, k+1 ), 1, queue );
                   #endif 

                   magma_scopy( n-imax, dA( imax, imax ), 1, dW( imax, k+1 ), 1, queue );
                   #if defined(PRECISION_z)
                   magma_dsetvector_async( 1, &d_zero, 1, 
                                           dW, 2*((k+1)*lddw+imax+dW_offset)+1, 1, queue, &event);
                   magma_queue_sync( queue );
                   #elif defined(PRECISION_c) 
                   magma_ssetvector_async( 1, &f_zero, 1,
                                           dW, 2*((k+1)*lddw+imax+dW_offset)+1, 1, queue, &event);
                   magma_queue_sync( queue );
                   #endif

                   magma_sgemv( MagmaNoTrans, n-k, k, c_mone, dA( k, 0 ), ldda, 
                                dW( imax, 0 ), lddw, c_one, dW( k, k+1 ), ione, queue );
                   #if defined(PRECISION_z)
                   magma_dsetvector_async( 1, &d_zero, 1, 
                                           dW, 2*((k+1)*lddw+imax+dW_offset)+1, 1, queue, &event);
                   magma_queue_sync( queue );
                   #elif defined(PRECISION_c)
                   magma_ssetvector_async( 1, &f_zero, 1,
                                           dW, 2*((k+1)*lddw+imax+dW_offset)+1, 1, queue, &event);
                   magma_queue_sync( queue );
                   #endif

                   magma_sgetvector_async( 1, dW(imax,k+1), 1, &Zimax, 1, queue, &event);
                   magma_queue_sync( queue );

                   /* jmax is the column-index of the largest off-diagonal
                      element in row imax, and rowmax is its absolute value */
 
                   // magmablas is one-base
                   jmax = k-1 + magma_isamax( imax-k, dW(k, k+1), 1, queue );

                   magma_sgetvector( 1, dW(jmax,k+1), 1, &Z, 1, queue );
                   rowmax = MAGMA_S_ABS1( Z );
                   if( imax < n-1 ) {
                       // magmablas is one-base
                       jmax = imax + magma_isamax( (n-1)-imax, dW(imax+1,k+1), 1, queue);
                       magma_sgetvector( 1, dW(jmax,k+1), 1, &Z, 1, queue );
                       rowmax = max( rowmax, MAGMA_S_ABS1( Z ) );
                   }

                   if( abs_akk >= alpha*colmax*( colmax / rowmax ) ) {

                       /* no interchange, use 1-by-1 pivot block */
                       kp = k;
                   } else if( fabs( MAGMA_S_REAL( Zimax ) ) >= alpha*rowmax ) {

                       /* interchange rows and columns K and imax, use 1-by-1
                          pivot block */
                       kp = imax;

                       /* copy column K+1 of W to column K */
                       magma_scopy( n-k, dW( k, k+1 ), 1, dW( k, k ), 1, queue );
                   } else {

                       /* interchange rows and columns K+1 and imax, use 2-by-2
                          pivot block */
                       kp = imax;
                       kstep = 2;
                   }
               }

               kk = k + kstep - 1;

               /* Updated column kp is already stored in column kk of W */
               if( kp != kk ) {

                   /* Copy non-updated column kk to column kp */
                   /* ------------------------------------------------------------------ */
                   #if defined(PRECISION_z) || defined(PRECISION_c)
                   magmablas_slacpy_cnjg( kp-kk, dA( kk, kk ), 1, dA( kp, kk ), ldda, queue );
                   #else
                   magma_scopy( kp-kk, dA( kk, kk ), 1, dA( kp, kk ), ldda, queue );
                   #endif
                   if ( kp < n ) {
                   magma_scopy( n-kp, dA( kp, kk), 1, dA( kp, kp ), 1, queue );
                   }
                   /* ------------------------------------------------------------------ */

                   /* Interchange rows kk and kp in first kk columns of A and W */
                   magmablas_sswap( kk+1, dA( kk, 0 ), ldda, dA( kp, 0 ), ldda, queue );
                   magmablas_sswap( kk+1, dW( kk, 0 ), lddw, dW( kp, 0 ), lddw, queue );
               }

               if ( kstep == 1 ) {
 
                   /* 1-by-1 pivot block D(k): column k of W now holds

                      W(k) = L(k)*D(k)

                      where L(k) is the k-th column of L
 
                      Store L(k) in column k of A */
                   magma_scopy( n-k, dW( k, k ), 1, dA( k, k ), 1, queue );

                   if ( k < n-1 ) {
                       magma_sgetvector_async( 1, dA(k,k), 1, &Z, 1, queue, &event );
                       magma_queue_sync( queue );
                       R1 = d_one / MAGMA_S_REAL( Z );
                       magma_sscal((n-1)-k, R1, dA( k+1,k ), 1, queue);

                       /* Conjugate W(k) */
                       #if defined(PRECISION_z) || defined(PRECISION_c)
                       magmablas_slacpy_cnjg( (n-1)-k, dW( k+1,k ), 1, dW( k+1,k ), 1, queue );
                       #endif
                   }
               } else {

                   /* 2-by-2 pivot block D(k): columns k and k+1 of W now hold

                   ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

                   where L(k) and L(k+1) are the k-th and (k+1)-th columns
                   of L */
                   magmablas_slascl_2x2( MagmaLower,
                                         n-(k+2), dW(k,k), lddw, dA(k+2,k), ldda, &iinfo, 
                                         queue );

                   /* Copy D(k) to A */
                   magma_scopymatrix( 2, 2, dW( k, k ), lddw, dA( k, k ), ldda, queue );

                   /* Conjugate W(k) and W(k+1) */
                   #if defined(PRECISION_z) || defined(PRECISION_c)
                   magmablas_slacpy_cnjg( (n-1)-k,   dW( k+1,k ),  1, dW( k+1,k ),   1, queue );
                   magmablas_slacpy_cnjg( (n-1)-k-1, dW( k+2,k+1), 1, dW( k+2,k+1 ), 1, queue );
                   #endif
               }
           }

           /* Store details of the interchanges in ipiv */
           if ( kstep == 1 ) {
               ipiv[k] = kp+1;
           } else {
               ipiv[k] = -kp-1;
               ipiv[k+1] = -kp-1;
           }
       } 

       /* Update the lower triangle of A22 (= A(k:n,k:n)) as

          A22 := A22 - L21*D*L21' = A22 - L21*W'

          computing blocks of NB columns at a time (note that conjg(W) is
          actually stored) */
       for( int j = k; j < n; j += nb ) {
           int jb = min( nb, n-j );

           /* Update the lower triangle of the diagonal block */

           #ifdef SYMMETRIC_UPDATE
           for (int jj = j; jj < j + jb; jj++) {
               int jnb = j + jb - jj;

               /* -------------------------------------------------------- */
               magma_sgemv( MagmaNoTrans, jnb, k, c_mone, dA( jj, 0 ), ldda, 
                            dW( jj, 0 ), lddw, c_one, dA( jj, jj ), ione, queue );
               /* -------------------------------------------------------- */
           }

           /* Update the rectangular subdiagonal block */

           if( j+jb < n ) {
               int nk = n - (j+jb);

               /* -------------------------------------------- */
               magma_sgemm( MagmaNoTrans, MagmaTrans, nk, jb, k, 
                            c_mone, dA( j+jb, 0 ), ldda, 
                                    dW( j, 0 ),    lddw,
                            c_one,  dA( j+jb, j ), ldda, queue );
               /* ------------------------------------------- */
           }
           #else

           #if defined(PRECISION_z)
           magmablas_dlaset(MagmaUpperLower, 1, jb,  
                            0, 0, dA, 2*(j*ldda+j+dA_offset)+1, 2*(1+ldda), queue );
           #elif defined(PRECISION_c)
           magmablas_slaset(MagmaUpperLower, 1, jb,
                            0, 0, dA, 2*(j*ldda+j+dA_offset)+1, 2*(1+ldda), queue );  
           #endif
           magma_sgemm( MagmaNoTrans, MagmaTrans, n-j, jb, k, 
                        c_mone, dA( j, 0 ), ldda, 
                                dW( j, 0 ), lddw,
                        c_one,  dA( j, j ), ldda, queue );
           #if defined(PRECISION_z)
           magmablas_dlaset(MagmaUpperLower, 1, jb,  
                            0, 0, dA, 2*(j*ldda+j+dA_offset)+1, 2*(1+ldda), queue );
           #elif defined(PRECISION_c)
           magmablas_slaset(MagmaUpperLower, 1, jb,
                            0, 0, dA, 2*(j*ldda+j+dA_offset)+1, 2*(1+ldda), queue );
           #endif
           #endif
       }

       /* Put L21 in standard form by partially undoing the interchanges
          in columns 1:k-1 */
       for (int j = k; j > 0;) {
           int jj = j;
           int jp = ipiv[j-1];
           if( jp < 0 ) {
               jp = -jp;
               j--;
           }
           j--;
           if ( jp != jj && j >= 1 ) {
               magmablas_sswap( j, dA( jp-1,0 ), ldda, dA( jj-1,0 ), ldda, queue );
           }
       }
       // copying the panel back to CPU
       magma_sgetmatrix_async( n, k, dA(0,0), ldda, A(0,0), lda, queue, &event );
       magma_queue_sync( queue );

       /* Set KB to the number of columns factorized */
       *kb = k;
    }

    return *info;
    /* End of SLASYF */
}

