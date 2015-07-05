#ifndef CLMAGMA_RUNTIME_H
#define CLMAGMA_RUNTIME_H

#include <string>
#include <map>
#include <vector>

#include "common_magma.h"  // includes OpenCL, etc.
#include "error.h"


// ------------------------------------------------------------
class clmagma_runtime
{
public:
    // ------------------------------
    const static int MAX_DEVICES = 8;
    
    // ------------------------------
    clmagma_runtime():
        m_context      ( NULL ),
        m_num_devices  ( 0 )
    {
        for( int dev=0; dev < MAX_DEVICES; ++dev ) {
            m_devices[dev] = NULL;
        }
    }

    // ------------------------------
    ~clmagma_runtime()
    {
        quit();
    }
    
    // ------------------------------
    void init( bool require_double=true );
    void quit();
    int compile_file( const char* infile, const char* outfile );
    void save_programs( std::vector< cl_program >& programs, const char* filename );
    void load_programs( int nfiles, const char* const* infiles, std::vector< cl_program >& programs );
    void archive_files( int nfiles, const char* const* infiles, const char* outfile );
    void load_kernels( int nfiles, const char* const* infiles );
    
    // ------------------------------
    cl_kernel get_kernel( const char* name )
    {
        //printf( "kernel: %s\n", name );
        cl_kernel k = m_kernels[ name ];
        if ( k == NULL ) {
            fprintf( stderr, "Error: kernel '%s' not found\n", name );
            return NULL;
        }
        return k;
    }
    
    // ------------------------------
    cl_platform_id get_platform()     const { return m_platform;    }
    cl_context     get_context()      const { return m_context;     }
    int            get_num_devices()  const { return m_num_devices; }
    cl_device_id*  get_devices()            { return m_devices;     }
    
    // ==============================
private:
    cl_platform_id   m_platform;
    cl_context       m_context;
    cl_uint          m_num_devices;
    cl_device_id     m_devices[ MAX_DEVICES ];
    std::map< std::string, cl_kernel > m_kernels;
};


// ------------------------------------------------------------
// global runtime
extern clmagma_runtime g_runtime;

#endif        //  #ifndef CLMAGMA_RUNTIME_H
