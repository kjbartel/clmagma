/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from testing_zgetrf_msub.cpp normal z -> c, Sat Nov 15 00:21:40 2014
       @author Mark Gates
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"


// Initialize matrix to random.
// Having this in separate function ensures the same ISEED is always used,
// so we can re-generate the identical matrix.
void init_matrix( int m, int n, magmaFloatComplex *h_A, magma_int_t lda )
{
    magma_int_t ione = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t n2 = lda*n;
    lapackf77_clarnv( &ione, ISEED, &n2, h_A );
}


// On input, A and ipiv is LU factorization of A. On output, A is overwritten.
// Requires m == n.
// Uses init_matrix() to re-generate original A as needed.
// Generates random RHS b and solves Ax=b.
// Returns residual, |Ax - b| / (n |A| |x|).
float get_residual(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex *A, magma_int_t lda,
    magma_int_t *ipiv )
{
    if ( m != n ) {
        printf( "\nERROR: residual check defined only for square matrices\n" );
        return -1;
    }

    const magmaFloatComplex c_one     = MAGMA_C_ONE;
    const magmaFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
    const magma_int_t ione = 1;

    // this seed should be DIFFERENT than used in init_matrix
    // (else x is column of A, so residual can be exactly zero)
    magma_int_t ISEED[4] = {0,0,0,2};
    magma_int_t info = 0;
    magmaFloatComplex *x, *b;

    // initialize RHS
    TESTING_MALLOC_CPU( x, magmaFloatComplex, n );
    TESTING_MALLOC_CPU( b, magmaFloatComplex, n );
    lapackf77_clarnv( &ione, ISEED, &n, b );
    blasf77_ccopy( &n, b, &ione, x, &ione );

    // solve Ax = b
    lapackf77_cgetrs( "Notrans", &n, &ione, A, &lda, ipiv, x, &n, &info );
    if ( info != 0 )
        printf("lapackf77_cgetrs returned error %d: %s.\n",
               (int) info, magma_strerror( info ));

    // reset to original A
    init_matrix( m, n, A, lda );

    // compute r = Ax - b, saved in b
    blasf77_cgemv( "Notrans", &m, &n, &c_one, A, &lda, x, &ione, &c_neg_one, b, &ione );

    // compute residual |Ax - b| / (n*|A|*|x|)
    float norm_x, norm_A, norm_r, work[1];
    norm_A = lapackf77_clange( "F", &m, &n, A, &lda, work );
    norm_r = lapackf77_clange( "F", &n, &ione, b, &n, work );
    norm_x = lapackf77_clange( "F", &n, &ione, x, &n, work );

    //printf( "r=\n" ); magma_cprint( 1, n, b, 1 );
    
    TESTING_FREE_CPU( x );
    TESTING_FREE_CPU( b );

    //printf( "r=%.2e, A=%.2e, x=%.2e, n=%d\n", norm_r, norm_A, norm_x, n );
    return norm_r / (n * norm_A * norm_x);
}


// On input, LU and ipiv is LU factorization of A. On output, LU is overwritten.
// Works for any m, n.
// Uses init_matrix() to re-generate original A as needed.
// Returns error in factorization, |PA - LU| / (n |A|)
// This allocates 3 more matrices to store A, L, and U.
float get_LU_error(magma_int_t M, magma_int_t N,
                    magmaFloatComplex *LU, magma_int_t lda,
                    magma_int_t *ipiv)
{
    magma_int_t min_mn = min(M,N);
    magma_int_t ione   = 1;
    magma_int_t i, j;
    magmaFloatComplex alpha = MAGMA_C_ONE;
    magmaFloatComplex beta  = MAGMA_C_ZERO;
    magmaFloatComplex *A, *L, *U;
    float work[1], matnorm, residual;
                       
    TESTING_MALLOC_CPU( A, magmaFloatComplex, lda*N    );
    TESTING_MALLOC_CPU( L, magmaFloatComplex, M*min_mn );
    TESTING_MALLOC_CPU( U, magmaFloatComplex, min_mn*N );
    memset( L, 0, M*min_mn*sizeof(magmaFloatComplex) );
    memset( U, 0, min_mn*N*sizeof(magmaFloatComplex) );

    // set to original A
    init_matrix( M, N, A, lda );
    lapackf77_claswp( &N, A, &lda, &ione, &min_mn, ipiv, &ione );
    
    // copy LU to L and U, and set diagonal to 1
    lapackf77_clacpy( MagmaLowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_clacpy( MagmaUpperStr, &min_mn, &N, LU, &lda, U, &min_mn );
    for( j=0; j < min_mn; j++ )
        L[j+j*M] = MAGMA_C_MAKE( 1., 0. );
    
    matnorm = lapackf77_clange( "f", &M, &N, A, &lda, work );

    blasf77_cgemm( "N", "N", &M, &N, &min_mn,
                   &alpha, L, &M, U, &min_mn, &beta, LU, &lda );

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_C_SUB( LU[i+j*lda], A[i+j*lda] );
        }
    }
    residual = lapackf77_clange( "f", &M, &N, LU, &lda, work );

    TESTING_FREE_CPU( A );
    TESTING_FREE_CPU( L );
    TESTING_FREE_CPU( U );

    return residual / (matnorm * N);
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgetrf_mgpu
*/
int main( int argc, char** argv)
{
    TESTING_INIT();
    
    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf=0, cpu_time=0;
    float error;
    magmaFloatComplex *h_A, *h_P;
    magmaFloatComplex_ptr d_lA[ MagmaMaxSubs * MagmaMaxGPUs ];
    magma_int_t     *ipiv;
    magma_int_t M, N, n2, lda, ldda, info, min_mn;
    magma_int_t dev, j, k, ngpu, nsub, n_local, nb, nk, ldn_local, maxm;
    magma_int_t status   = 0;

    magma_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    /* Initialize queues */
    magma_queue_t  queues[MagmaMaxGPUs * 2];
    magma_device_t devices[MagmaMaxGPUs];
    magma_int_t num = 0;
    magma_int_t err;
    err = magma_getdevices( devices, MagmaMaxGPUs, &num );
    if ( err != 0 || num < 1 ) {
        fprintf( stderr, "magma_getdevices failed: %d\n", (int) err );
        exit(-1);
    }
    for( dev=0; dev < opts.ngpu; dev++ ) {
        err = magma_queue_create( devices[dev], &queues[2*dev] );
        if ( err != 0 ) {
            fprintf( stderr, "magma_queue_create failed: %d (device %d)\n", (int) err, dev );
            exit(-1);
        }
        err = magma_queue_create( devices[dev], &queues[2*dev+1] );
        if ( err != 0 ) {
            fprintf( stderr, "magma_queue_create failed: %d (device %d)\n", (int) err, dev );
            exit(-1);
        }
    }
    
    printf("trans %s, ngpu %d, nsub %d\n",
           lapack_trans_const(opts.transA), (int) opts.ngpu, (int) opts.nsub );
    if ( opts.check == 2 ) {
        printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |Ax-b|/(N*|A|*|x|)\n");
    }
    else {
        printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |PA-LU|/(N*|A|)\n");
    }
    printf("=========================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            min_mn = min(M, N);
            maxm   = 32*((M+31)/32);
            lda    = M;
            n2     = lda*N;
            nb     = magma_get_cgetrf_nb(M);
            gflops = FLOPS_CGETRF( M, N ) / 1e9;
            
            // nsubs * ngpu must be at least the number of blocks
            ngpu = opts.ngpu;
            nsub = opts.nsub;
            if ( nsub*ngpu > N/nb ) {
                nsub = 1;
                ngpu = 1;
                printf( " * too many GPUs for the matrix size, using %d GPUs and %d submatrices\n", (int) ngpu, (int) nsub );
            }
    
            /* Allocate host memory for the matrix */
            TESTING_MALLOC_CPU( ipiv, magma_int_t, min_mn );
            TESTING_MALLOC_CPU( h_A, magmaFloatComplex, n2 );
            TESTING_MALLOC_CPU( h_P, magmaFloatComplex, lda*nb );
            
            /* Allocate device memory */
            if ( opts.transA == MagmaNoTrans ) {
                ldda = N/nb;                     /* number of block columns         */
                ldda = ldda/(ngpu*nsub);     /* number of block columns per GPU */
                ldda = nb*ldda;                  /* number of columns per GPU       */
                if ( ldda * ngpu*nsub < N ) {
                    /* left over */
                    if ( N-ldda*ngpu*nsub >= nb ) {
                        ldda += nb;
                    } else {
                        ldda += (N-ldda*ngpu*nsub)%nb;
                    }
                }
                ldda = ((ldda+31)/32)*32; /* make it a multiple of 32 */
                for( j=0; j < nsub * ngpu; j++ ) {
                    TESTING_MALLOC_DEV( d_lA[j], magmaFloatComplex, ldda*maxm );
                }
            } else {
                ldda = ((M+31)/32)*32;
                for( j=0; j < nsub * ngpu; j++ ) {
                    n_local = ((N/nb)/(nsub*ngpu))*nb;
                    if ( j < (N/nb)%(nsub*ngpu) ) {
                        n_local += nb;
                    } else if ( j == (N/nb)%(nsub*ngpu) ) {
                        n_local += N%nb;
                    }
                    TESTING_MALLOC_DEV( d_lA[j], magmaFloatComplex, ldda*n_local );
                }
            }

            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                init_matrix( M, N, h_A, lda );
                
                cpu_time = magma_wtime();
                lapackf77_cgetrf( &M, &N, h_A, &lda, ipiv, &info );
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if ( info != 0 )
                    printf("lapackf77_cgetrf returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
            }
    
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            init_matrix( M, N, h_A, lda );
            if ( opts.transA == MagmaNoTrans ) {
                for( j=0; j < N; j += nb ) {
                    k = (j/nb)%(nsub*ngpu);
                    nk = min(nb, N-j);
    
                    /* transpose on CPU, then copy to GPU */
                    int ii,jj;
                    for( ii=0; ii < M; ii++ ) {
                        for( jj=0; jj < nk; jj++ ) {
                            h_P[jj+ii*nk] = h_A[j*lda + ii+jj*lda];
                        }
                    }
                    magma_csetmatrix( nk, M,
                                      h_P, nk,
                                      d_lA[k], j/(nb*nsub*ngpu)*nb, ldda,
                                      queues[2*(k%ngpu)] );
                }
            } else {
                ldda = ((M+31)/32)*32;
                for( j=0; j < N; j += nb ) {
                    k = (j/nb)%(nsub*ngpu);
                    nk = min(nb, N-j);
                    magma_csetmatrix( M, nk,
                                      h_A + j*lda, lda,
                                      d_lA[k], j/(nb*nsub*ngpu)*nb*ldda, ldda,
                                      queues[2*(k%ngpu)] );
                }
            }
            
            gpu_time = magma_wtime();
            magma_cgetrf_msub( opts.transA, nsub, ngpu, M, N, d_lA, 0, ldda, ipiv, queues, &info );
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_cgetrf_mgpu returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));
    
            /* get the matrix from GPUs */
            if ( opts.transA == MagmaNoTrans ) {
                for (j=0; j < N; j+=nb) {
                    k = (j/nb)%(nsub*ngpu);
                    nk = min(nb, N-j);
    
                    /* copy to CPU and then transpose */
                    magma_cgetmatrix( nk, M,
                                      d_lA[k], j/(nb*nsub*ngpu)*nb, ldda,
                                      h_P, nk, queues[2*(k%ngpu)] );
                    int ii, jj;
                    for( ii=0; ii < M; ii++ ) {
                        for( jj=0; jj < nk; jj++ ) {
                            h_A[j*lda + ii+jj*lda] = h_P[jj+ii*nk];
                        }
                    }
                }
            } else {
                for (j=0; j < N; j+=nb) {
                    k = (j/nb)%(nsub*ngpu);
                    nk = min(nb, N-j);
                    magma_cgetmatrix( M, nk,
                                      d_lA[k], j/(nb*nsub*ngpu)*nb*ldda, ldda,
                                      h_A + j*lda, lda, queues[2*(k%ngpu)] );
                }
            }
    
            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (int) M, (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            }
            else {
                printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)",
                       (int) M, (int) N, gpu_perf, gpu_time );
            }
            if ( opts.check == 2 ) {
                error = get_residual( M, N, h_A, lda, ipiv );
                printf("   %8.2e   %s\n", error, (error < tol ? "ok" : "failed"));
                status += ! (error < tol);
            }
            else if ( opts.check ) {
                error = get_LU_error( M, N, h_A, lda, ipiv );
                printf("   %8.2e   %s\n", error, (error < tol ? "ok" : "failed"));
                status += ! (error < tol);
            }
            else {
                printf("     ---  \n");
            }
            
            TESTING_FREE_CPU( ipiv );
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_P );
            for( dev=0; dev < ngpu; dev++ ) {
                for( k=0; k < nsub; k++ ) {
                    TESTING_FREE_DEV( d_lA[dev*nsub + k] );
                }
            }
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    /* Free queues */
    for( dev=0; dev < opts.ngpu; dev++ ) {
        magma_queue_destroy( queues[2*dev] );
        magma_queue_destroy( queues[2*dev+1] );
    }

    TESTING_FINALIZE();
    return status;
}
