#ifndef TESTINGS_H
#define TESTINGS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "magma.h"

#ifndef min
#define min(a,b)  (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b)  (((a)<(b))?(b):(a))
#endif

#define TESTING_MALLOC( ptr, type, size )                        \
    ptr = (type*) malloc((size) * sizeof(type));                 \
    if ( ptr == 0 ) {                                            \
        fprintf( stderr, "!!!! malloc failed for: %s\n", #ptr ); \
        exit(-1);                                                \
    }

#define TESTING_MALLOC_HOST( ptr, type, size )                                       \
    if ( MAGMA_SUCCESS != magma_malloc_cpu( (void**)&ptr, (size)*sizeof(type) ) ) { \
        fprintf( stderr, "!!!! magma_malloc_cpu failed for: %s\n", #ptr );          \
        exit(-1);                                                                    \
    }

#define TESTING_MALLOC_DEV( ptr, type, size )                              \
    if ( MAGMA_SUCCESS != magma_malloc( &ptr, (size)*sizeof(type) ) ) {    \
        fprintf( stderr, "!!!! magma_malloc_dev failed for: %s\n", #ptr ); \
        exit(-1);                                                          \
    }

#define TESTING_FREE( ptr ) free( ptr )

#define TESTING_FREE_HOST( ptr ) magma_free_cpu( ptr )

#define TESTING_FREE_DEV( ptr ) magma_free( ptr )

// need double precision type that does NOT get converted to float
// by the codegen.py script, to store times without overflow.
typedef double real_Double_t;

#ifdef __cplusplus
extern "C" {
#endif

void magma_zmake_hermitian( magma_int_t N, magmaDoubleComplex* A, magma_int_t lda );
void magma_cmake_hermitian( magma_int_t N, magmaFloatComplex*  A, magma_int_t lda );
void magma_dmake_symmetric( magma_int_t N, double*             A, magma_int_t lda );
void magma_smake_symmetric( magma_int_t N, float*              A, magma_int_t lda );

void magma_zmake_hpd( magma_int_t N, magmaDoubleComplex* A, magma_int_t lda );
void magma_cmake_hpd( magma_int_t N, magmaFloatComplex*  A, magma_int_t lda );
void magma_dmake_hpd( magma_int_t N, double*             A, magma_int_t lda );
void magma_smake_hpd( magma_int_t N, float*              A, magma_int_t lda );

void magma_assert( bool condition, const char* msg, ... );

#define MAX_NTEST 1000

typedef struct magma_opts
{
    // matrix size
    magma_int_t ntest;
    magma_int_t msize[ MAX_NTEST ];
    magma_int_t nsize[ MAX_NTEST ];
    magma_int_t ksize[ MAX_NTEST ];
    magma_int_t mmax;
    magma_int_t nmax;
    magma_int_t kmax;

    // scalars
    magma_int_t device;
    magma_int_t nb;
    magma_int_t nrhs;
    magma_int_t nstream;
    magma_int_t ngpu;
    magma_int_t niter;
    magma_int_t nthread;
    magma_int_t itype;     // hegvd: problem type
    magma_int_t svd_work;  // gesvd
    magma_int_t version;   // hemm_mgpu, hetrd
    double      fraction;  // hegvdx
    double      tolerance;
    magma_int_t panel_nthread; //first dimension for a 2D big panel

    // boolean arguments
    int check;
    int lapack;
    int warmup;
    int all;

    // lapack flags
    magma_uplo_t    uplo;
    magma_trans_t   transA;
    magma_trans_t   transB;
    magma_side_t    side;
    magma_diag_t    diag;
    magma_vec_t     jobu;    // gesvd:  no left  singular vectors
    magma_vec_t     jobvt;   // gesvd:  no right singular vectors
    magma_vec_t     jobz;    // heev:   no eigen vectors
    magma_vec_t     jobvr;   // geev:   no right eigen vectors
    magma_vec_t     jobvl;   // geev:   no left  eigen vectors
} magma_opts;

void parse_opts( int argc, char** argv, magma_opts *opts );

#ifdef __cplusplus
}
#endif

#endif /* TESTINGS_H */
