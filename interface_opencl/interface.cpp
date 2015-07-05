/*
 *   -- clMAGMA (version 1.0.0) --
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
#include "CL_MAGMA_RT.h"

#ifdef HAVE_clAmdBlas

// ========================================
// globals
cl_platform_id gPlatform;
cl_context     gContext;

// Run time global variable used for LU
CL_MAGMA_RT *rt;

// ========================================
// initialization
magma_err_t
magma_init()
{
    cl_int err;
    err = clGetPlatformIDs( 1, &gPlatform, NULL );
    assert( err == 0 );
    
    cl_device_id devices[ MagmaMaxGPUs ];
    cl_uint num;
    err = clGetDeviceIDs( gPlatform, CL_DEVICE_TYPE_GPU, MagmaMaxGPUs, devices, &num );
    assert( err == 0 );
    
    cl_context_properties properties[3] =
        { CL_CONTEXT_PLATFORM, (cl_context_properties) gPlatform, 0 };
    gContext = clCreateContext( properties, num, devices, NULL, NULL, &err );
    assert( err == 0 );
    
    err = clAmdBlasSetup();
    assert( err == 0 );

    // Initialize kernels related to LU
	rt = CL_MAGMA_RT::Instance();
    rt->Init(gPlatform, gContext);
    
    return err;
}

// --------------------
magma_err_t
magma_finalize()
{
    cl_int err;
    clAmdBlasTeardown();
    err = clReleaseContext( gContext );

    // quit the RT
    rt->Quit();

    return err;
}


// ========================================
// memory allocation
// #include "CL/cl_ext.h"
magma_err_t
magma_malloc( magma_ptr* ptrPtr, size_t size )
{
    cl_int err;
    *ptrPtr = clCreateBuffer( gContext, CL_MEM_READ_WRITE, size, NULL, &err );
    // *ptrPtr = clCreateBuffer( gContext, CL_MEM_READ_WRITE | CL_MEM_USE_PERSISTENT_MEM_AMD, size, NULL, &err );
    return err;
}

// --------------------
magma_err_t
magma_free( magma_ptr ptr )
{
    cl_int err = clReleaseMemObject( ptr );
    return err;
}

// --------------------
magma_err_t
magma_malloc_host( void** ptrPtr, size_t size )
{
    *ptrPtr = malloc( size );
    if ( *ptrPtr == NULL ) {
        return MAGMA_ERR_HOST_ALLOC;
    }
    else {
        return MAGMA_SUCCESS;
    }
}

// --------------------
magma_err_t
magma_free_host( void* ptr )
{
    free( ptr );
    return MAGMA_SUCCESS;
}


// ========================================
// device & queue support
magma_err_t
magma_get_devices(
	magma_device_t* devices,
	magma_int_t     size,
	magma_int_t*    numPtr )
{
    cl_int err;
    //err = clGetDeviceIDs( gPlatform, CL_DEVICE_TYPE_GPU, 1, size, devices, num );
    size_t n;
    err = clGetContextInfo(
        gContext, CL_CONTEXT_DEVICES,
        size*sizeof(magma_device_t), devices, &n );
    *numPtr = n / sizeof(magma_device_t);
    return err;
}

// --------------------
magma_err_t
magma_queue_create( magma_device_t device, magma_queue_t* queuePtr )
{
    assert( queuePtr != NULL );
    cl_int err;
    *queuePtr = clCreateCommandQueue( gContext, device, 0, &err );
    return err;
}

// --------------------
magma_err_t
magma_queue_destroy( magma_queue_t  queue )
{
    cl_int err = clReleaseCommandQueue( queue );
    return err;
}

// --------------------
magma_err_t
magma_queue_sync( magma_queue_t queue )
{
    cl_int err = clFinish( queue );
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
    cl_int err = clWaitForEvents(1, &event);
    return err;
}

#endif // HAVE_clAmdBlas
