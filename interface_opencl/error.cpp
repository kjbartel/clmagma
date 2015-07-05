#include <stdio.h>

#include "error.h"


// ----------------------------------------
// C++ function is overloaded for different error types,
// which depends on error types being enums to be differentiable.
void magma_xerror( cl_int err, const char* func, const char* file, int line )
{
    if ( err != CL_SUCCESS ) {
        fprintf( stderr, "OpenCL runtime error: %s (%d) in %s at %s:%d\n",
                 magma_clGetErrorString( err ), err, func, file, line );
    }
}


// --------------------
void magma_xerror( clblasStatus err, const char* func, const char* file, int line )
{
    if ( err != clblasSuccess ) {
        fprintf( stderr, "clBLAS error: %s (%d) in %s at %s:%d\n",
                 magma_clblasGetErrorString( err ), err, func, file, line );
    }
}


// --------------------
// Can't differentiate OpenCL errors (cl_int) from MAGMA errors (magma_int_t).
//void magma_xerror( magma_int_t err, const char* func, const char* file, int line )
//{
//    if ( err != MAGMA_SUCCESS ) {
//        fprintf( stderr, "MAGMA error: %s (%d) in %s at %s:%d\n",
//                 magma_strerror( err ), (int) err, func, file, line );
//    }
//}


// ----------------------------------------
extern "C"
const char* magma_clGetErrorString( cl_int error )
{
    switch( error ) {
        case CL_SUCCESS                                  : return "CL_SUCCESS";
        case CL_DEVICE_NOT_FOUND                         : return "CL_DEVICE_NOT_FOUND";
        case CL_DEVICE_NOT_AVAILABLE                     : return "CL_DEVICE_NOT_AVAILABLE";
        case CL_COMPILER_NOT_AVAILABLE                   : return "CL_COMPILER_NOT_AVAILABLE";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE            : return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case CL_OUT_OF_RESOURCES                         : return "CL_OUT_OF_RESOURCES";
        case CL_OUT_OF_HOST_MEMORY                       : return "CL_OUT_OF_HOST_MEMORY";
        case CL_PROFILING_INFO_NOT_AVAILABLE             : return "CL_PROFILING_INFO_NOT_AVAILABLE";
        case CL_MEM_COPY_OVERLAP                         : return "CL_MEM_COPY_OVERLAP";
        case CL_IMAGE_FORMAT_MISMATCH                    : return "CL_IMAGE_FORMAT_MISMATCH";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED               : return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
        case CL_BUILD_PROGRAM_FAILURE                    : return "CL_BUILD_PROGRAM_FAILURE";
        case CL_MAP_FAILURE                              : return "CL_MAP_FAILURE";
        case CL_MISALIGNED_SUB_BUFFER_OFFSET             : return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
        case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
        
        case CL_INVALID_VALUE                            : return "CL_INVALID_VALUE";
        case CL_INVALID_DEVICE_TYPE                      : return "CL_INVALID_DEVICE_TYPE";
        case CL_INVALID_PLATFORM                         : return "CL_INVALID_PLATFORM";
        case CL_INVALID_DEVICE                           : return "CL_INVALID_DEVICE";
        case CL_INVALID_CONTEXT                          : return "CL_INVALID_CONTEXT";
        case CL_INVALID_QUEUE_PROPERTIES                 : return "CL_INVALID_QUEUE_PROPERTIES";
        case CL_INVALID_COMMAND_QUEUE                    : return "CL_INVALID_COMMAND_QUEUE";
        case CL_INVALID_HOST_PTR                         : return "CL_INVALID_HOST_PTR";
        case CL_INVALID_MEM_OBJECT                       : return "CL_INVALID_MEM_OBJECT";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR          : return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
        case CL_INVALID_IMAGE_SIZE                       : return "CL_INVALID_IMAGE_SIZE";
        case CL_INVALID_SAMPLER                          : return "CL_INVALID_SAMPLER";
        case CL_INVALID_BINARY                           : return "CL_INVALID_BINARY";
        case CL_INVALID_BUILD_OPTIONS                    : return "CL_INVALID_BUILD_OPTIONS";
        case CL_INVALID_PROGRAM                          : return "CL_INVALID_PROGRAM";
        case CL_INVALID_PROGRAM_EXECUTABLE               : return "CL_INVALID_PROGRAM_EXECUTABLE";
        case CL_INVALID_KERNEL_NAME                      : return "CL_INVALID_KERNEL_NAME";
        case CL_INVALID_KERNEL_DEFINITION                : return "CL_INVALID_KERNEL_DEFINITION";
        case CL_INVALID_KERNEL                           : return "CL_INVALID_KERNEL";
        case CL_INVALID_ARG_INDEX                        : return "CL_INVALID_ARG_INDEX";
        case CL_INVALID_ARG_VALUE                        : return "CL_INVALID_ARG_VALUE";
        case CL_INVALID_ARG_SIZE                         : return "CL_INVALID_ARG_SIZE";
        case CL_INVALID_KERNEL_ARGS                      : return "CL_INVALID_KERNEL_ARGS";
        case CL_INVALID_WORK_DIMENSION                   : return "CL_INVALID_WORK_DIMENSION";
        case CL_INVALID_WORK_GROUP_SIZE                  : return "CL_INVALID_WORK_GROUP_SIZE";
        case CL_INVALID_WORK_ITEM_SIZE                   : return "CL_INVALID_WORK_ITEM_SIZE";
        case CL_INVALID_GLOBAL_OFFSET                    : return "CL_INVALID_GLOBAL_OFFSET";
        case CL_INVALID_EVENT_WAIT_LIST                  : return "CL_INVALID_EVENT_WAIT_LIST";
        case CL_INVALID_EVENT                            : return "CL_INVALID_EVENT";
        case CL_INVALID_OPERATION                        : return "CL_INVALID_OPERATION";
        case CL_INVALID_GL_OBJECT                        : return "CL_INVALID_GL_OBJECT";
        case CL_INVALID_BUFFER_SIZE                      : return "CL_INVALID_BUFFER_SIZE";
        case CL_INVALID_MIP_LEVEL                        : return "CL_INVALID_MIP_LEVEL";
        case CL_INVALID_GLOBAL_WORK_SIZE                 : return "CL_INVALID_GLOBAL_WORK_SIZE";
        case CL_INVALID_PROPERTY                         : return "CL_INVALID_PROPERTY";
        
        default:
            return "unknown OpenCL error code";
    }
}


// ----------------------------------------
extern "C"
const char* magma_clblasGetErrorString( clblasStatus error )
{
    switch( error ) {        
        case clblasNotImplemented     : return "clblasNotImplemented";
        case clblasNotInitialized     : return "clblasNotInitialized";
        case clblasInvalidMatA        : return "clblasInvalidMatA";
        case clblasInvalidMatB        : return "clblasInvalidMatB";
        case clblasInvalidMatC        : return "clblasInvalidMatC";
        case clblasInvalidVecX        : return "clblasInvalidVecX";
        case clblasInvalidVecY        : return "clblasInvalidVecY";
        case clblasInvalidDim         : return "clblasInvalidDim";
        case clblasInvalidLeadDimA    : return "clblasInvalidLeadDimA";
        case clblasInvalidLeadDimB    : return "clblasInvalidLeadDimB";
        case clblasInvalidLeadDimC    : return "clblasInvalidLeadDimC";
        case clblasInvalidIncX        : return "clblasInvalidIncX";
        case clblasInvalidIncY        : return "clblasInvalidIncY";
        case clblasInsufficientMemMatA: return "clblasInsufficientMemMatA";
        case clblasInsufficientMemMatB: return "clblasInsufficientMemMatB";
        case clblasInsufficientMemMatC: return "clblasInsufficientMemMatC";
        case clblasInsufficientMemVecX: return "clblasInsufficientMemVecX";
        case clblasInsufficientMemVecY: return "clblasInsufficientMemVecY";
        
        default: return "unknown clblas error code";
    }
}


// ----------------------------------------
extern "C"
const char* magma_strerror( magma_int_t error )
{
    // LAPACK-compliant errors
    if ( error > 0 ) {
        return "function-specific error, see documentation";
    }
    else if ( error < 0 && error > MAGMA_ERR ) {
        return "invalid argument";
    }
    // MAGMA-specific errors
    switch( error ) {
        case MAGMA_SUCCESS:
            return "success";
        
        case MAGMA_ERR:
            return "unknown error";
        
        case MAGMA_ERR_NOT_INITIALIZED:
            return "not initialized";
        
        case MAGMA_ERR_REINITIALIZED:
            return "reinitialized";
        
        case MAGMA_ERR_NOT_SUPPORTED:
            return "not supported";
        
        case MAGMA_ERR_ILLEGAL_VALUE:
            return "illegal value";
        
        case MAGMA_ERR_NOT_FOUND:
            return "not found";
        
        case MAGMA_ERR_ALLOCATION:
            return "allocation";
        
        case MAGMA_ERR_INTERNAL_LIMIT:
            return "internal limit";
        
        case MAGMA_ERR_UNALLOCATED:
            return "unallocated error";
        
        case MAGMA_ERR_FILESYSTEM:
            return "filesystem error";
        
        case MAGMA_ERR_UNEXPECTED:
            return "unexpected error";
        
        case MAGMA_ERR_SEQUENCE_FLUSHED:
            return "sequence flushed";
        
        case MAGMA_ERR_HOST_ALLOC:
            return "cannot allocate memory on CPU host";
        
        case MAGMA_ERR_DEVICE_ALLOC:
            return "cannot allocate memory on GPU device";
        
        case MAGMA_ERR_CUDASTREAM:
            return "CUDA stream error";
        
        case MAGMA_ERR_INVALID_PTR:
            return "invalid pointer";
        
        default:
            return "unknown MAGMA error code";
    }
}
