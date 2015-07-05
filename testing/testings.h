#ifndef TESTINGS_H
#define TESTINGS_H

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
    if ( MAGMA_SUCCESS != magma_malloc_host( (void**)&ptr, (size)*sizeof(type) ) ) { \
        fprintf( stderr, "!!!! magma_malloc_host failed for: %s\n", #ptr );               \
        exit(-1);                                                                    \
    }

#define TESTING_MALLOC_DEV( ptr, type, size )                                   \
    if ( MAGMA_SUCCESS != magma_malloc( &ptr, (size)*sizeof(type) ) ) { \
        fprintf( stderr, "!!!! magma_malloc_dev failed for: %s\n", #ptr );     \
		exit(-1);                                                               \
    }

#define TESTING_FREE( ptr ) free( ptr )

#define TESTING_FREE_HOST( ptr ) magma_free_host( ptr )

#define TESTING_FREE_DEV( ptr ) magma_free( ptr )

// need double precision type that does NOT get converted to float
// by the codegen.py script, to store times without overflow.
typedef double real_Double_t;

#endif /* TESTINGS_H */
