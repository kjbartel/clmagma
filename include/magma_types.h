/*
 *   -- clMAGMA (version 0.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 */

#ifndef MAGMA_TYPES_H
#define MAGMA_TYPES_H

#include <sys/types.h>
#include <assert.h>


// ========================================
// types
typedef int magma_err_t;
typedef int magma_int_t;


// ========================================
#if HAVE_CUBLAS
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>
    
    typedef cublasHandle_t magma_queue_t;
    typedef cudaEvent_t    magma_event_t;
    typedef int            magma_device_t;
    
    typedef cuDoubleComplex magmaDoubleComplex;
    typedef cuFloatComplex  magmaFloatComplex;
    #define MAGMA_Z_MAKE(r,i)     make_cuDoubleComplex(r, i)
    #define MAGMA_Z_REAL(a)       (a).x
    #define MAGMA_Z_IMAG(a)       (a).y
    #define MAGMA_Z_SET2REAL(a,r) (a).x = (r);   (a).y =  0.0
    #define MAGMA_Z_ADD(a, b)     cuCadd(a, b)
    #define MAGMA_Z_SUB(a, b)     cuCsub(a, b)
    #define MAGMA_Z_MUL(a, b)     cuCmul(a, b)
    #define MAGMA_Z_DIV(a, b)     cuCdiv(a, b)
    #define MAGMA_Z_ABS(a)        cuCabs(a)
    #define MAGMA_Z_CNJG(a)       cuConj(a)
    
    #define MAGMA_C_MAKE(r,i)     make_cuFloatComplex(r, i)
    #define MAGMA_C_REAL(a)       (a).x
    #define MAGMA_C_IMAG(a)       (a).y
    #define MAGMA_C_SET2REAL(a,r) (a).x = (r);   (a).y =  0.0
    #define MAGMA_C_ADD(a, b)     cuCaddf(a, b)
    #define MAGMA_C_SUB(a, b)     cuCsubf(a, b)
    #define MAGMA_C_MUL(a, b)     cuCmulf(a, b)
    #define MAGMA_C_DIV(a, b)     cuCdivf(a, b)
    #define MAGMA_C_ABS(a)        cuCabsf(a)
    #define MAGMA_C_CNJG(a)       cuConjf(a)
    
#elif HAVE_AMDBLAS
    #if defined(__APPLE__) || defined(__MACOSX)
    #include "my_amdblas.h"
    #else
    #include <clAmdBlas.h>
    #endif
    
    typedef cl_command_queue  magma_queue_t;
    typedef cl_event          magma_event_t;
    typedef cl_device_id      magma_device_t;
    
    typedef DoubleComplex magmaDoubleComplex;
    typedef FloatComplex  magmaFloatComplex;
    #define MAGMA_Z_MAKE(r,i)     doubleComplex(r,i)
    #define MAGMA_Z_REAL(a)       (a).x
    #define MAGMA_Z_IMAG(a)       (a).y
    #define MAGMA_Z_SET2REAL(a,r) (a).x = (r);   (a).y =  0.0
    #define MAGMA_Z_ADD(a, b)     MAGMA_Z_MAKE((a).x+(b).x, (a).y+(b).y)
    #define MAGMA_Z_SUB(a, b)     MAGMA_Z_MAKE((a).x-(b).x, (a).y-(b).y)
    #define MAGMA_Z_CNJG(a)       MAGMA_Z_MAKE((a).x, -(a).y)
    
    #define MAGMA_C_MAKE(r,i)     floatComplex(r,i)
    #define MAGMA_C_REAL(a)       (a).x
    #define MAGMA_C_IMAG(a)       (a).y
    #define MAGMA_C_SET2REAL(a,r) (a).x = (r);   (a).y =  0.0
    #define MAGMA_C_ADD(a, b)     MAGMA_C_MAKE((a).x+(b).x, (a).y+(b).y)
    #define MAGMA_C_SUB(a, b)     MAGMA_C_MAKE((a).x-(b).x, (a).y-(b).y)
    #define MAGMA_C_CNJG(a)       MAGMA_C_MAKE((a).x, -(a).y)
    
#else
    // generic types if no GPU code available
    typedef int magma_queue_t;
    typedef int magma_event_t;
    typedef int magma_device_t;
    
    #ifdef __cplusplus
        #include <complex>
        typedef std::complex<float>   magmaFloatComplex;
        typedef std::complex<double>  magmaDoubleComplex;
        #define MAGMA_Z_MAKE(r, i)  std::complex<double>(r,i)
        #define MAGMA_Z_REAL(x)     x.real()
        #define MAGMA_Z_IMAG(x)     x.imag()
        
        #define MAGMA_C_MAKE(r, i)  std::complex<float> (r,i)
        #define MAGMA_C_REAL(x)     x.real()
        #define MAGMA_C_IMAG(x)     x.imag()
    #else
        #include <complex.h>
        typedef float  complex /*_Complex*/ magmaFloatComplex;
        typedef double complex /*_Complex*/ magmaDoubleComplex;
        #define MAGMA_Z_MAKE(r, i)  ((r) + (i)*_Complex_I)
        #define MAGMA_Z_REAL(x)     creal(x)
        #define MAGMA_Z_IMAG(x)     cimag(x)
        
        #define MAGMA_C_MAKE(r, i)  ((r) + (i)*_Complex_I)
        #define MAGMA_C_REAL(x)     crealf(x)
        #define MAGMA_C_IMAG(x)     cimagf(x)
    #endif
#endif

#define MAGMA_D_MAKE(r,i)     (r)
#define MAGMA_D_REAL(x)       (x)
#define MAGMA_D_IMAG(x)       (0.0f)
#define MAGMA_D_SET2REAL(a,r) (a) = (r)
#define MAGMA_D_ADD(a, b)     ((a) + (b))
#define MAGMA_D_SUB(a, b)     ((a) - (b))
#define MAGMA_D_CNJG(a)       (a)

#define MAGMA_S_MAKE(r,i)     (r)
#define MAGMA_S_REAL(x)       (x)
#define MAGMA_S_IMAG(x)       (0.0)
#define MAGMA_S_SET2REAL(a,r) (a) = (r)
#define MAGMA_S_ADD(a, b)     ((a) + (b))
#define MAGMA_S_SUB(a, b)     ((a) - (b))
#define MAGMA_S_CNJG(a)       (a)

#define MAGMA_Z_ZERO              MAGMA_Z_MAKE( 0.0, 0.0)
#define MAGMA_Z_ONE               MAGMA_Z_MAKE( 1.0, 0.0)
#define MAGMA_Z_HALF              MAGMA_Z_MAKE( 0.5, 0.0)
#define MAGMA_Z_NEG_ONE           MAGMA_Z_MAKE(-1.0, 0.0)
#define MAGMA_Z_NEG_HALF          MAGMA_Z_MAKE(-0.5, 0.0)

#define MAGMA_C_ZERO              MAGMA_C_MAKE( 0.0, 0.0)
#define MAGMA_C_ONE               MAGMA_C_MAKE( 1.0, 0.0)
#define MAGMA_C_HALF              MAGMA_C_MAKE( 0.5, 0.0)
#define MAGMA_C_NEG_ONE           MAGMA_C_MAKE(-1.0, 0.0)
#define MAGMA_C_NEG_HALF          MAGMA_C_MAKE(-0.5, 0.0)

#define MAGMA_D_ZERO              ( 0.0)
#define MAGMA_D_ONE               ( 1.0)
#define MAGMA_D_HALF              ( 0.5)
#define MAGMA_D_NEG_ONE           (-1.0)
#define MAGMA_D_NEG_HALF          (-0.5)

#define MAGMA_S_ZERO              ( 0.0)
#define MAGMA_S_ONE               ( 1.0)
#define MAGMA_S_HALF              ( 0.5)
#define MAGMA_S_NEG_ONE           (-1.0)
#define MAGMA_S_NEG_HALF          (-0.5)


#if HAVE_AMDBLAS
    // OpenCL uses opaque memory references on GPU
    typedef cl_mem magma_ptr;
    typedef cl_mem magmaFloat_ptr;
    typedef cl_mem magmaDouble_ptr;
    typedef cl_mem magmaFloatComplex_ptr;
    typedef cl_mem magmaDoubleComplex_ptr;
    
    typedef cl_mem magma_const_ptr;
    typedef cl_mem magmaFloat_const_ptr;
    typedef cl_mem magmaDouble_const_ptr;
    typedef cl_mem magmaFloatComplex_const_ptr;
    typedef cl_mem magmaDoubleComplex_const_ptr;
#else
    // CUDA uses regular pointers on GPU
    typedef void               *magma_ptr;
    typedef float              *magmaFloat_ptr;
    typedef double             *magmaDouble_ptr;
    typedef magmaFloatComplex  *magmaFloatComplex_ptr;
    typedef magmaDoubleComplex *magmaDoubleComplex_ptr;
    
    typedef void               const *magma_const_ptr;
    typedef float              const *magmaFloat_const_ptr;
    typedef double             const *magmaDouble_const_ptr;
    typedef magmaFloatComplex  const *magmaFloatComplex_const_ptr;
    typedef magmaDoubleComplex const *magmaDoubleComplex_const_ptr;
#endif


// ----------------------------------------
#define MagmaMaxGPUs 4


// ----------------------------------------
// error codes
#define MAGMA_SUCCESS                    0
#define MAGMA_ERR_UNKNOWN            -1000
#define MAGMA_ERR_NOT_IMPLEMENTED    -1001
#define MAGMA_ERR_HOST_ALLOC         -1002
#define MAGMA_ERR_DEVICE_ALLOC       -1003


// ----------------------------------------
// parameter constants
// numbering is consistent with CBLAS and PLASMA; see plasma/include/plasma.h
#define MagmaRowMajor      101
#define MagmaColMajor      102
                           
#define MagmaNoTrans       111
#define MagmaTrans         112
#define MagmaConjTrans     113
                            
#define MagmaUpper         121
#define MagmaLower         122
#define MagmaFull          123
                            
#define MagmaNonUnit       131
#define MagmaUnit          132
                            
#define MagmaLeft          141
#define MagmaRight         142

#define MagmaOneNorm       171
#define MagmaRealOneNorm   172
#define MagmaTwoNorm       173
#define MagmaFrobeniusNorm 174
#define MagmaInfNorm       175
#define MagmaRealInfNorm   176
#define MagmaMaxNorm       177
#define MagmaRealMaxNorm   178

#define MagmaDistUniform   201
#define MagmaDistSymmetric 202
#define MagmaDistNormal    203

#define MagmaHermGeev      241
#define MagmaHermPoev      242
#define MagmaNonsymPosv    243
#define MagmaSymPosv       244

#define MagmaNoPacking     291
#define MagmaPackSubdiag   292
#define MagmaPackSupdiag   293
#define MagmaPackColumn    294
#define MagmaPackRow       295
#define MagmaPackLowerBand 296
#define MagmaPackUpeprBand 297
#define MagmaPackAll       298

#define MagmaNoVec         301
#define MagmaVec           302

#define MagmaForward       391
#define MagmaBackward      392

#define MagmaColumnwise    401
#define MagmaRowwise       402

// remember to update min/max when adding constants!
#define MagmaMinConst      101
#define MagmaMaxConst      402

// constants for lapack
// todo: use translators instead? lapack_trans_const( MagmaUpper )
#define MagmaRowMajorStr   "Row"
#define MagmaColMajorStr   "Col"

#define MagmaNoTransStr    "NoTrans"
#define MagmaTransStr      "Trans"
#define MagmaConjTransStr  "ConjTrans"
                         
#define MagmaUpperStr      "Upper"
#define MagmaLowerStr      "Lower"
#define MagmaFullStr       "Full"
#define MagmaUpperLowerStr "All"
                         
#define MagmaNonUnitStr    "NonUnit"
#define MagmaUnitStr       "Unit"
                         
#define MagmaLeftStr       "Left"
#define MagmaRightStr      "Right"

#define MagmaForwardStr    "Forward"
#define MagmaBackwardStr   "Backward"

#define MagmaColumnwiseStr "Columnwise"
#define MagmaRowwiseStr    "Rowwise"

// these could be enums, but that isn't portable in C++, e.g., if -fshort-enums is used
typedef int magma_trans_t;
typedef int magma_uplo_t;
typedef int magma_diag_t;
typedef int magma_side_t;


#ifdef __cplusplus
extern "C" {
#endif

// --------------------
// translators
int magma_trans_const( char lapack_char );
int magma_uplo_const ( char lapack_char );
int magma_diag_const ( char lapack_char );
int magma_side_const ( char lapack_char );

char        lapacke_const( int magma_const );
const char* lapack_const ( int magma_const );

#define lapacke_order_const(c) lapack_const(c)
#define lapacke_trans_const(c) lapack_const(c)
#define lapacke_side_const( c) lapack_const(c)
#define lapacke_diag_const( c) lapack_const(c)
#define lapacke_uplo_const( c) lapack_const(c)

#define lapack_order_const(c)  lapack_const(c)
#define lapack_trans_const(c)  lapack_const(c)
#define lapack_side_const( c)  lapack_const(c)
#define lapack_diag_const( c)  lapack_const(c)
#define lapack_uplo_const( c)  lapack_const(c)

#ifdef HAVE_AMDBLAS
int                  amdblas_const      ( int magma_const );
clAmdBlasOrder       amdblas_order_const( int magma_const );
clAmdBlasTranspose   amdblas_trans_const( int magma_const );
clAmdBlasSide        amdblas_side_const ( int magma_const );
clAmdBlasDiag        amdblas_diag_const ( int magma_const );
clAmdBlasUplo        amdblas_uplo_const ( int magma_const );
#endif

#ifdef HAVE_CUBLAS
int                  cublas_const       ( int magma_const );
cublasOperation_t    cublas_trans_const ( int magma_const );
cublasSideMode_t     cublas_side_const  ( int magma_const );
cublasDiagType_t     cublas_diag_const  ( int magma_const );
cublasFillMode_t     cublas_uplo_const  ( int magma_const );
#endif

#ifdef HAVE_CBLAS
#include "cblas.h"
enum CBLAS_ORDER     cblas_order_const  ( int magma_const );
enum CBLAS_TRANSPOSE cblas_trans_const  ( int magma_const );
enum CBLAS_SIDE      cblas_side_const   ( int magma_const );
enum CBLAS_DIAG      cblas_diag_const   ( int magma_const );
enum CBLAS_UPLO      cblas_uplo_const   ( int magma_const );
#endif

// todo: above functions should all be inlined or macros for
// efficiency. Here's an example.
// In C99, static inline potentially wastes some space by
// emitting multiple definitions, but is portable.
static inline int cblas_const( int magma_const ) {
    assert( magma_const >= MagmaMinConst );
    assert( magma_const <= MagmaMaxConst );
    return magma_const;
}

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_TYPES_H
