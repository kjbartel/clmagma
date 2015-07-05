/*
 *   -- clMAGMA (version 1.3.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      @date November 2014
 *
 * @author Mark Gates
 */

#include <stdlib.h>
#include <stdio.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(MAGMA_WITH_MKL)
#include <mkl_service.h>
#endif

#if defined(MAGMA_WITH_ACML)
// header conflicts with magma's lapack prototypes, so declare function directly
// #include <acml.h>
extern "C"
void acmlversion(int *major, int *minor, int *patch);
#endif

#include "clmagma_runtime.h"
#include "common_magma.h"
#include "error.h"

#ifdef HAVE_clBLAS

// ========================================
// globals
cl_context     gContext;
magma_event_t  *g_event;

const char* clmagma_kernels = "clmagma_kernels.co";


// ========================================
// initialization
// --------------------
extern "C" magma_int_t
magma_init()
{
    g_runtime.init();
    g_runtime.load_kernels( 1, &clmagma_kernels );
    gContext = g_runtime.get_context();
    
    cl_int err = clblasSetup();
    check_error( err );
    
    g_event = NULL;

    return err;
}

// --------------------
//extern "C" magma_int_t
//magma_init_opencl( cl_platform_id platform, cl_context context, bool setup_clBlas )
//{
//    g_runtime.init( platform, context );
//    g_runtime.load_kernels( 1, &clmagma_kernels );
//    gPlatform = platform;
//    gContext  = Context;
//
//    cl_int err = 0;
//    if ( setup_clBlas ) {
//        err = clblasSetup();
//    }
//
//    g_event = NULL;
//
//    return err;
//}

// --------------------
extern "C" magma_int_t
magma_finalize()
{
    clblasTeardown();
    g_runtime.quit();
    return MAGMA_SUCCESS;
}

// --------------------
//extern "C" magma_int_t
//magma_finalize_opencl( bool finalize_clBlas )
//{
//    if ( finalize_clBlas ) {
//        clblasTeardown();
//    }
//    g_runtime->quit();
//
//    return MAGMA_SUCCESS;
//}

// --------------------
// Print the available GPU devices. Used in testing.
extern "C" void
magma_print_environment()
{
    magma_int_t major, minor, micro;
    magma_version( &major, &minor, &micro );
    printf( "%% clMAGMA %d.%d.%d %s\n",
            (int) major, (int) minor, (int) micro, MAGMA_VERSION_STAGE );

    // CUDA, OpenCL, OpenMP, MKL, ACML versions all printed on same line
    char device_name[1024], driver[1024];
    clGetPlatformInfo( g_runtime.get_platform(), CL_PLATFORM_VERSION, sizeof(device_name), device_name, NULL );
    printf( "%% OpenCL platform %s.", device_name );
    
#if defined(_OPENMP)
    int omp_threads = 0;
    #pragma omp parallel
    {
        omp_threads = omp_get_num_threads();
    }
    printf( " OpenMP threads %d.", omp_threads );
#else
    printf( " MAGMA not compiled with OpenMP." );
#endif

#if defined(MAGMA_WITH_MKL)
    MKLVersion mkl_version;
    mkl_get_version( &mkl_version );
    printf( " MKL %d.%d.%d, MKL threads %d.",
            mkl_version.MajorVersion,
            mkl_version.MinorVersion,
            mkl_version.UpdateVersion,
            mkl_get_max_threads() );
#endif
    
#if defined(MAGMA_WITH_ACML)
    int acml_major, acml_minor, acml_patch;
    acmlversion( &acml_major, &acml_minor, &acml_patch );
    printf( " ACML %d.%d.%d.", acml_major, acml_minor, acml_patch );
#endif

    printf( "\n" );
    
    // print devices
    int ndevices = g_runtime.get_num_devices();
    cl_device_id* devices = g_runtime.get_devices();
    cl_ulong mem_size, alloc_size;
    for( int dev=0; dev < ndevices; ++dev ) {
        clGetDeviceInfo( devices[dev], CL_DEVICE_NAME,               sizeof(device_name), device_name, NULL );
        clGetDeviceInfo( devices[dev], CL_DEVICE_GLOBAL_MEM_SIZE,    sizeof(mem_size),    &mem_size,   NULL );
        clGetDeviceInfo( devices[dev], CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(alloc_size),  &alloc_size, NULL );
        clGetDeviceInfo( devices[dev], CL_DRIVER_VERSION,            sizeof(driver),      driver,      NULL );
        printf( "%% Device: %s, %.1f MiB memory, max allocation %.1f MiB, driver  %s\n",
                device_name, mem_size/(1024.*1024.), alloc_size/(1024.*1024.), driver );
    }
}


// ========================================
// device support
extern "C" magma_int_t
magma_getdevices(
    magma_device_t* devices,
    magma_int_t     size,
    magma_int_t*    numPtr )
{
    // TODO just copy from g_runtime.get_devices()
    cl_int err;
    //err = clGetDeviceIDs( gPlatform, CL_DEVICE_TYPE_GPU, 1, size, devices, num );
    size_t n;
    err = clGetContextInfo(
        g_runtime.get_context(), CL_CONTEXT_DEVICES,
        size*sizeof(magma_device_t), devices, &n );
    *numPtr = n / sizeof(magma_device_t);
    check_error( err );
    return err;
}

// --------------------
extern "C" magma_int_t
magma_num_gpus( void )
{
    const char *ngpu_str = getenv("MAGMA_NUM_GPUS");
    cl_uint ngpu = 1;
    if ( ngpu_str != NULL ) {
        char* endptr;
        ngpu = strtol( ngpu_str, &endptr, 10 );

        cl_uint ndevices = g_runtime.get_num_devices();

        if ( ngpu < 1 || *endptr != '\0' ) {
            ngpu = 1;
            fprintf( stderr, "$MAGMA_NUM_GPUS='%s' is an invalid number; using %d GPU.\n",
                     ngpu_str, (int) ngpu );
        }
        else if ( ngpu > MagmaMaxGPUs || ngpu > ndevices ) {
            ngpu = min( ndevices, MagmaMaxGPUs );
            fprintf( stderr, "$MAGMA_NUM_GPUS='%s' exceeds MagmaMaxGPUs=%d or available GPUs=%d; using %d GPUs.\n",
                     ngpu_str, MagmaMaxGPUs, ndevices, (int) ngpu );
        }
        assert( 1 <= ngpu && ngpu <= ndevices );
    }
    return (magma_int_t)ngpu;
}

// --------------------
extern "C" magma_int_t
magma_queue_meminfo( magma_queue_t queue )
{
    cl_device_id dev;
    clGetCommandQueueInfo(queue, CL_QUEUE_DEVICE, sizeof(cl_device_id), &dev, NULL);
  
    cl_ulong mem_size;
    clGetDeviceInfo(dev, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &mem_size, NULL);
    mem_size /= sizeof(magmaDoubleComplex);

    return mem_size;
}


// ========================================
// queue support

// --------------------
extern "C" magma_int_t
magma_queue_create( magma_device_t device, magma_queue_t* queuePtr )
{
    assert( queuePtr != NULL );
    cl_int err;
    #ifdef TRACING
    *queuePtr = clCreateCommandQueue( g_runtime.get_context(), device, CL_QUEUE_PROFILING_ENABLE, &err );
    #else
    *queuePtr = clCreateCommandQueue( g_runtime.get_context(), device, 0, &err );
    #endif
    check_error( err );
    return err;
}

// --------------------
extern "C" magma_int_t
magma_queue_destroy( magma_queue_t  queue )
{
    cl_int err = clReleaseCommandQueue( queue );
    check_error( err );
    return err;
}

// --------------------
extern "C" magma_int_t
magma_queue_sync( magma_queue_t queue )
{
    cl_int err = clFinish( queue );
    clFlush( queue );
    check_error( err );
    return err;
}


// ========================================
// event support
// --------------------
extern "C" magma_int_t
magma_event_create( magma_event_t* event )
{
    printf( "%s not implemented\n", __func__ );
    return 0;
}

// --------------------
extern "C" magma_int_t
magma_event_destroy( magma_event_t event )
{
    printf( "%s not implemented\n", __func__ );
    return 0;
}

// --------------------
extern "C" magma_int_t
magma_event_record( magma_event_t event, magma_queue_t queue )
{
    printf( "%s not implemented\n", __func__ );
    return 0;
}

// --------------------
extern "C" magma_int_t
magma_event_query( magma_event_t event )
{
    printf( "%s not implemented\n", __func__ );
    return 0;
}

// --------------------
// blocks CPU until event occurs
extern "C" magma_int_t
magma_event_sync( magma_event_t event )
{
    cl_int err = clWaitForEvents(1, &event);
    check_error( err );
    return err;
}

// --------------------
// blocks queue (but not CPU) until event occurs
// In OpenCL, the event should be passed directly to the kernel, instead of a separate function call.
// Make an empty kernel? Seems a waste.
extern "C" void
magma_queue_wait_event( magma_queue_t queue, magma_event_t event )
{
    printf( "%s not implemented\n", __func__ );
}

// --------------------
// Sets global event used for tracing.
// TODO put this into tracing code, not in interface.cpp?
extern "C" magma_int_t
magma_setevent( magma_event_t* event )
{
    #ifdef TRACING
    g_event = event;
    #else
    printf("%s not implemented\n", __func__ );
    #endif

    return 0;
}


// ========================================
// complex support

// --------------------
extern "C" double
magma_cabs(magmaDoubleComplex z)
{
    double __x = z.x;
    double __y = z.y;

    double __s = max(abs(__x), abs(__y));
    if(__s == 0.0)
        return __s;
    __x /= __s;
    __y /= __s;
    return __s * sqrt(__x * __x + __y * __y);
}

// --------------------
extern "C" float
magma_cabsf(magmaFloatComplex z)
{
    float __x = z.x;
    float __y = z.y;

    float __s = max(abs(__x), abs(__y));
    if(__s == 0.0)
        return __s;
    __x /= __s;
    __y /= __s;
    return __s * sqrt(__x * __x + __y * __y);
}

#endif // HAVE_clBLAS
