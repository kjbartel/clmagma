/*
 *   -- clMAGMA (version 0.2.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 */

#include <stdlib.h>
#include <stdio.h>

#include "magma.h"

#ifdef HAVE_CUBLAS

// ========================================
// initialization
magma_err_t
magma_init()
{
    return MAGMA_SUCCESS;
}

// --------------------
magma_err_t
magma_finalize()
{
    return MAGMA_SUCCESS;
}


// ========================================
// memory allocation
magma_err_t
magma_malloc( magma_ptr* ptrPtr, size_t size )
{
    cudaError_t err;
    err = cudaMalloc( ptrPtr, size );
    return err;
}

// --------------------
magma_err_t
magma_free( magma_ptr ptr )
{
    cudaError_t err;
    err = cudaFree( ptr );
    return err;
}

// --------------------
magma_err_t
magma_malloc_host( void** ptrPtr, size_t size )
{
    cudaError_t err;
    err = cudaMallocHost( ptrPtr, size );
    return err;
}

// --------------------
magma_err_t
magma_free_host( void* ptr )
{
    cudaError_t err;
    err = cudaFreeHost( ptr );
    return err;
}


// ========================================
// device & queue support
// --------------------
magma_err_t
magma_get_devices(
	magma_device_t* devices,
	magma_int_t     size,
	magma_int_t*    numPtr )
{
    // todo, this should query for the number of devices
    printf( "%s not implemented\n", __func__ );
    for( int i = 0; i < size; ++i ) {
        devices[i] = i;
    }
    *numPtr = size;
    return 0;
}

// --------------------
magma_err_t
magma_queue_create( magma_device_t device, magma_queue_t* queuePtr )
{
    cudaStream_t   stream;
    cublasStatus_t stat;
    cudaError_t    err;
    err  = cudaSetDevice( device );
    stat = cublasCreate( queuePtr );
    err  = cudaStreamCreate( &stream );
    stat = cublasSetStream( *queuePtr, stream );
    return err;
}

// --------------------
magma_err_t
magma_queue_destroy( magma_queue_t  queue )
{
    cudaStream_t   stream;
    cublasStatus_t stat;
    cudaError_t    err;
    stat = cublasGetStream( queue, &stream );
    stat = cublasDestroy( queue );
    err  = cudaStreamDestroy( stream );
    return err;
}

// --------------------
magma_err_t
magma_queue_sync( magma_queue_t queue )
{
    cudaStream_t   stream;
    cublasStatus_t stat;
    cudaError_t    err;
    stat = cublasGetStream( queue, &stream );
    err  = cudaStreamSynchronize( stream );
    return err;
}


// ========================================
// event support
magma_err_t
magma_event_create( magma_event_t* event )
{
    printf( "%s not implemented\n", __func__ );
    return 0;
}

magma_err_t
magma_event_destroy( magma_event_t event )
{
    printf( "%s not implemented\n", __func__ );
    return 0;
}

magma_err_t
magma_event_record( magma_event_t event, magma_queue_t queue )
{
    printf( "%s not implemented\n", __func__ );
    return 0;
}

magma_err_t
magma_event_query( magma_event_t event )
{
    printf( "%s not implemented\n", __func__ );
    return 0;
}

magma_err_t
magma_event_sync( magma_event_t event )
{
    cudaError_t err;
    err  = cudaEventSynchronize( event );
    return err;
}

#endif // HAVE_CUBLAS
