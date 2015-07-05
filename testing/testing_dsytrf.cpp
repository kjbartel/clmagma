/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @generated from testing_zhetrf.cpp normal z -> d, Sat Nov 15 00:21:40 2014
       @author Ichitaro Yamazaki
       @author Stan Tomov
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

#include "common_magma.h"

/* ================================================================================================== */
// Initialize matrix to random.
// Having this in separate function ensures the same ISEED is always used,
// so we can re-generate the identical matrix.
void init_matrix( int nopiv, int m, int n, double *h_A, magma_int_t lda )
{
    magma_int_t ione = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t n2 = lda*n;
    //double *A = (double*)malloc(n2*sizeof(double));
    //lapackf77_dlarnv( &ione, ISEED, &n2, A );
    //for (int i=0; i<n; i++) for (int j=0; j<=i; j++) h_A[j+i*lda] = MAGMA_D_MAKE(A[j+i*lda],0.0);
    //free(A);
    //
    lapackf77_dlarnv( &ione, ISEED, &n2, h_A );
    // symmetrize
    for (int i=0; i<n; i++) for (int j=0; j<i; j++) h_A[i+j*lda] = MAGMA_D_CNJG(h_A[j+i*lda]);
    if (nopiv) for (int i=0; i<n; i++) h_A[i+i*lda] = MAGMA_D_MAKE(MAGMA_D_REAL(h_A[i+i*lda]) + n, 0.0);
    else       for (int i=0; i<n; i++) h_A[i+i*lda] = MAGMA_D_MAKE(MAGMA_D_REAL(h_A[i+i*lda]), 0.0);
}


// On input, A and ipiv is LU factorization of A. On output, A is overwritten.
// Requires m == n.
// Uses init_matrix() to re-generate original A as needed.
// Generates random RHS b and solves Ax=b.
// Returns residual, |Ax - b| / (n |A| |x|).
double get_residual(
    int nopiv, magma_uplo_t uplo, magma_int_t n,
    double *A, magma_int_t lda,
    magma_int_t *ipiv )
{
    const double c_one     = MAGMA_D_ONE;
    const double c_neg_one = MAGMA_D_NEG_ONE;
    const magma_int_t ione = 1;
    magma_int_t upper = (uplo == MagmaUpper);
    
    // this seed should be DIFFERENT than used in init_matrix
    // (else x is column of A, so residual can be exactly zero)
    magma_int_t ISEED[4] = {0,0,0,2};
    magma_int_t info = 0;
    double *x, *b;
    
    // initialize RHS
    TESTING_MALLOC_CPU( x, double, n );
    TESTING_MALLOC_CPU( b, double, n );
    lapackf77_dlarnv( &ione, ISEED, &n, b );
    blasf77_dcopy( &n, b, &ione, x, &ione );
    
    // solve Ax = b
    if (nopiv) {
        blasf77_dtrsm( "L", "L", "N", "U", &n, &ione, &c_one,
                       A, &lda, x, &n );
        //for (int i=0; i<n; i++) x[i] /= A[i+i*lda];
        for (int i=0; i<n; i++) x[i] = x[i]/A[i+i*lda]; //MAGMA_D_DIV( x[i], A[i+i*lda] );
        blasf77_dtrsm( "L", "L", "C", "U", &n, &ione, &c_one,
                       A, &lda, x, &n );
    }else {
        dsytrs_( (upper ? MagmaUpperStr: MagmaLowerStr), &n, &ione, A, &lda, ipiv, x, &n, &info );
    }
    if (info != 0)
        printf("lapackf77_dsytrs returned error %d: %s.\n",
               (int) info, magma_strerror( info ));
    // reset to original A
    init_matrix( nopiv, n, n, A, lda );
    
    // compute r = Ax - b, saved in b
    blasf77_dgemv( "Notrans", &n, &n, &c_one, A, &lda, x, &ione, &c_neg_one, b, &ione );
    
    // compute residual |Ax - b| / (n*|A|*|x|)
    double norm_x, norm_A, norm_r, work[1];
    norm_A = lapackf77_dlange( "F", &n, &n, A, &lda, work );
    norm_r = lapackf77_dlange( "F", &n, &ione, b, &n, work );
    norm_x = lapackf77_dlange( "F", &n, &ione, x, &n, work );
    
    TESTING_FREE_CPU( x );
    TESTING_FREE_CPU( b );
    
    return norm_r / (n * norm_A * norm_x);
}


// On input, LU and ipiv is LU factorization of A. On output, LU is overwritten.
// Works for any m, n.
// Uses init_matrix() to re-generate original A as needed.
// Returns error in factorization, |PA - LU| / (n |A|)
// This allocates 3 more matrices to store A, L, and U.
double get_LU_error(magma_int_t M, magma_int_t N,
                    double *LU, magma_int_t lda,
                    magma_int_t *ipiv)
{
    magma_int_t min_mn = min(M,N);
    magma_int_t ione   = 1;
    magma_int_t i, j;
    double alpha = MAGMA_D_ONE;
    double beta  = MAGMA_D_ZERO;
    double *A, *L, *U;
    double work[1], matnorm, residual;
    
    TESTING_MALLOC_CPU( A, double, lda*N    );
    TESTING_MALLOC_CPU( L, double, M*min_mn );
    TESTING_MALLOC_CPU( U, double, min_mn*N );
    memset( L, 0, M*min_mn*sizeof(double) );
    memset( U, 0, min_mn*N*sizeof(double) );

    // set to original A
    init_matrix( 0, M, N, A, lda );
    lapackf77_dlaswp( &N, A, &lda, &ione, &min_mn, ipiv, &ione);
    
    // copy LU to L and U, and set diagonal to 1
    lapackf77_dlacpy( MagmaLowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_dlacpy( MagmaUpperStr, &min_mn, &N, LU, &lda, U, &min_mn );
    for(j=0; j<min_mn; j++)
        L[j+j*M] = MAGMA_D_MAKE( 1., 0. );
    
    matnorm = lapackf77_dlange("f", &M, &N, A, &lda, work);

    blasf77_dgemm("N", "N", &M, &N, &min_mn,
                  &alpha, L, &M, U, &min_mn, &beta, LU, &lda);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_D_SUB( LU[i+j*lda], A[i+j*lda] );
        }
    }
    residual = lapackf77_dlange("f", &M, &N, LU, &lda, work);

    TESTING_FREE_CPU( A );
    TESTING_FREE_CPU( L );
    TESTING_FREE_CPU( U );

    return residual / (matnorm * N);
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dgetrf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf=0, cpu_time=0;
    double          error, error_lapack = 0.0;
    double *h_A, *work, temp;
    magma_int_t     *ipiv;
    magma_int_t     N, n2, lda, lwork, info;
    magma_int_t     status = 0;
    magma_int_t     cpu = 0, gpu = 1, nopiv = 0, row = 0;
    
    magma_opts opts;
    for(int i = 1; i < argc; ++i ) {
        if ( strcmp("--cpu", argv[i]) == 0 ) {
            cpu = 1;
        }
        if ( strcmp("--gpu", argv[i]) == 0 ) {
            gpu = 1;
        }
        if ( strcmp("--row", argv[i]) == 0 ) {
            row = 1;
        }
        if ( strcmp("--nopiv", argv[i]) == 0 ) {
            nopiv = 1;
        }
    }
    parse_opts( argc, argv, &opts );
    magma_uplo_t uplo = opts.uplo;

    if (nopiv)
        printf( "\n No-piv version (A is SPD)" );
    else if (cpu)
        printf( "\n CPU-only version" );
    else if (gpu)
        printf( "\n GPU-only version" );
    else if (row)
        printf( "\n GPU-only version (row-major)" );
    else
        printf( " hybrid CPU-GPU version" );
    printf( " (%s)\n\n",(uplo == MagmaUpper ? "upper" : "lower") );
    
    magma_int_t upper = (uplo == MagmaUpper);
    double tol = opts.tolerance * lapackf77_dlamch("E");

    printf("ngpu %d\n", (int) opts.ngpu );
    if ( opts.check == 2 ) {
        printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |Ax-b|/(N*|A|*|x|)\n");
    }
    else {
        printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |PA-LU|/(N*|A|)\n");
    }
    printf("=========================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            n2     = lda*N;
            gflops = FLOPS_DGETRF( N, N ) / 2e9;
            
            TESTING_MALLOC_CPU( ipiv, magma_int_t, N );
            TESTING_MALLOC_PIN( h_A,  double, n2 );
            
            lwork = -1;
            lapackf77_dsytrf((upper ? MagmaUpperStr: MagmaLowerStr), &N, h_A, &lda, ipiv, &temp, &lwork, &info);
            lwork = max(N*(1+magma_get_dsytrf_nb(N)), (int)MAGMA_D_REAL(temp));
            TESTING_MALLOC_PIN( work, double, lwork );

            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                init_matrix( nopiv, N, N, h_A, lda );
                cpu_time = magma_wtime();
                lapackf77_dsytrf((upper ? MagmaUpperStr: MagmaLowerStr), &N, h_A, &lda, 
                                 ipiv, work, &lwork, &info);
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_dgetrf returned error %d: %s.\n",
                           (int) info, magma_strerror( info ));
                error_lapack = get_residual( nopiv, uplo, N, h_A, lda, ipiv );
            }

            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            init_matrix( nopiv, N, N, h_A, lda );

            gpu_time = magma_wtime();
            if (nopiv) {
                //magma_dsytrf_nopiv( uplo, N, h_A, lda, &info);
            } else if (cpu) {
                //magma_dsytrf_cpu( uplo, N, h_A, lda, ipiv, work, lwork, &info);
            } else if (gpu) {
                magma_dsytrf( uplo, N, h_A, lda, ipiv, opts.queue, &info);
            } else if (row) {
                //magma_dsytrf_gpu_row( uplo, N, h_A, lda, ipiv, work, lwork, &info);
            } else {
                //magma_dsytrf( uplo, N, h_A, lda, ipiv, work, lwork, &info);
            }
            gpu_time = magma_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_dsytrf returned error %d: %s.\n",
                       (int) info, magma_strerror( info ));

            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (int) N, (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            }
            else {
                printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)",
                       (int) N, (int) N, gpu_perf, gpu_time );
            }
            if ( opts.check == 2 ) {
                error = get_residual( nopiv, uplo, N, h_A, lda, ipiv );
                printf("   %8.2e   %s", error, (error < tol ? "ok" : "failed"));
                if (opts.lapack)
                    printf(" (lapack rel.res. = %8.2e)", error_lapack);
                printf("\n");
                status += ! (error < tol);
            }
            else if ( opts.check ) {
                printf( " not yet..\n" ); exit(0);
                //error = get_LU_error( M, N, h_A, lda, ipiv );
                //printf("   %8.2e   %s\n", error, (error < tol ? "ok" : "failed"));
                //status += ! (error < tol);
            }
            else {
                printf("     ---   \n");
            }
            
            TESTING_FREE_CPU( ipiv );
            TESTING_FREE_PIN( work );
            TESTING_FREE_PIN( h_A  );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
