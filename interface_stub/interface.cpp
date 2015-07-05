/*
 *   -- clMAGMA (version 0.2.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      April 2012
 *
 * @author Mark Gates
 */

// defines stub functions, if no GPU implementation is available.
#include <stdio.h>
#include <stdlib.h>

#include "magma.h"

#if !defined(HAVE_CUBLAS) && !defined(HAVE_clAmdBlas) && !defined(HAVE_CBLAS)

// ========================================
// initialization
magma_err_t
magma_init( void )
{
    printf( "%s()\n", __func__ );
    return 0;
}

// --------------------
magma_err_t
magma_finalize( void )
{
    printf( "%s()\n", __func__ );
    return 0;
}


// ========================================
// memory allocation
magma_err_t
magma_malloc( magma_ptr* ptrPtr, size_t size )
{
    printf( "%s( ptr, %ld )\n", __func__, size );
    return 0;
}

// --------------------
magma_err_t
magma_free( magma_ptr ptr )
{
    printf( "%s( ptr )\n", __func__ );
    return 0;
}

// --------------------
// this actually allocates memory, else lapack may segfault
magma_err_t
magma_malloc_host( void** ptrPtr, size_t size )
{
    *ptrPtr = malloc( size );
    printf( "%s( ptr, %ld )\n", __func__, size );
    return 0;
}

// --------------------
// this actually frees memory.
magma_err_t
magma_free_host( void* ptr )
{
    free( ptr );
    printf( "%s( ptr )\n", __func__ );
    return 0;
}


// ========================================
// device & queue support
magma_err_t
magma_get_devices(
	magma_device_t* devices,
	magma_int_t     size,
	magma_int_t*    numPtr )
{
    printf( "%s( ptr, %d, numPtr )\n", __func__, size );
    return 0;
}

magma_err_t
magma_queue_create( magma_device_t device, magma_queue_t* queuePtr )
{
    *queuePtr = 1;
    printf( "%s( %d, %d )\n", __func__, device, *queuePtr );
    return 0;
}

magma_err_t
magma_queue_destroy( magma_queue_t  queue )
{
    printf( "%s( %d )\n", __func__, queue );
    return 0;
}

magma_err_t
magma_queue_sync( magma_queue_t queue )
{
    printf( "%s( %d )\n", __func__, queue );
    return 0;
}


// ========================================
// event support
magma_err_t
magma_event_create( magma_event_t* eventPtr )
{
    *eventPtr = 2;
    printf( "%s( %d )\n", __func__, *eventPtr );
    return 0;
}

magma_err_t
magma_event_destroy( magma_event_t event )
{
    printf( "%s( %d )\n", __func__, event );
    return 0;
}

magma_err_t
magma_event_record( magma_event_t event, magma_queue_t queue )
{
    printf( "%s( %d )\n", __func__, event );
    return 0;
}

magma_err_t
magma_event_query( magma_event_t event )
{
    printf( "%s( %d )\n", __func__, event );
    return 0;
}

magma_err_t
magma_event_sync( magma_event_t event )
{
    printf( "%s( %d )\n", __func__, event );
    return 0;
}

#endif // not HAVE_CUBLAS and not HAVE_clAmdBlas and not HAVE_CBLAS
