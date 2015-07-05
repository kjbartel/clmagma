#include <stdio.h>
#include <sys/time.h>    // gettimeofday
#include <assert.h>
#include <errno.h>
#include <string.h>      // strerror_r

#include "trace.h"

// define TRACING to compile these functions, e.g.,
// gcc -DTRACING -c trace.cpp
#ifdef TRACING

// set TRACE_METHOD = 2 to record start time as
// later of CPU time and previous event's end time.
// set TRACE_METHOD = 1 to record start time using CUDA event.

//#define TRACE_METHOD 2

extern cl_context gContext;


// ----------------------------------------
struct event_log
{
    int    ncore;
    int    cpu_id   [ MAX_CORES ];
    double cpu_first;
    double cpu_start[ MAX_CORES ][ MAX_EVENTS ];
    double cpu_end  [ MAX_CORES ][ MAX_EVENTS ];
    char   cpu_tag  [ MAX_CORES ][ MAX_EVENTS ][ MAX_LABEL_LEN ];
    char   cpu_label[ MAX_CORES ][ MAX_EVENTS ][ MAX_LABEL_LEN ];
    
    int          ngpu;
    int          nqueue;
    magma_queue_t queues  [ MAX_GPU_QUEUES ];
    int          gpu_id   [ MAX_GPU_QUEUES ];
    magma_event_t  gpu_first[ MAX_GPU_QUEUES ];
#if TRACE_METHOD == 2
    double       gpu_start[ MAX_GPU_QUEUES ][ MAX_EVENTS ];
#else
//  magma_event_t  gpu_start[ MAX_GPU_QUEUES ][ MAX_EVENTS ];
#endif
    magma_event_t  gpu_end  [ MAX_GPU_QUEUES ][ MAX_EVENTS ];
    char         gpu_tag  [ MAX_GPU_QUEUES ][ MAX_EVENTS ][ MAX_LABEL_LEN ];
    char         gpu_label[ MAX_GPU_QUEUES ][ MAX_EVENTS ][ MAX_LABEL_LEN ];
};

// global log object
struct event_log glog;


// ----------------------------------------
void trace_init( int ncore, int ngpu, int nqueue, magma_queue_t* queues )
{
    if ( ncore > MAX_CORES ) {
        fprintf( stderr, "Error in trace_init: ncore %d > MAX_CORES %d\n",
                 ncore, MAX_CORES );
        exit(1);
    }
    if ( ngpu*nqueue > MAX_GPU_QUEUES ) {
        fprintf( stderr, "Error in trace_init: (ngpu=%d)*(nqueue=%d) > MAX_GPU_QUEUES=%d\n",
                 ngpu, nqueue, MAX_GPU_QUEUES );
        exit(1);
    }
    
    glog.ncore   = ncore;
    glog.ngpu    = ngpu;
    glog.nqueue = nqueue;
    
    // initialize ID = 0
    for( int core = 0; core < ncore; ++core ) {
        glog.cpu_id[core] = 0;
    }
    for( int dev = 0; dev < ngpu; ++dev ) {
        for( int s = 0; s < nqueue; ++s ) {
            int t = dev*glog.nqueue + s;
            glog.gpu_id[t] = 0;
            glog.queues[t] = queues[t];
        }
		/* In OpenCL, queues are assocaited on different devices
        cudaSetDevice( dev );
        cudaDeviceSynchronize();
		*/
    }
    // now that all GPUs are sync'd, record start time
	// using clEnqueueCopyBuffer as a GPU dummy point
	cl_mem dA = clCreateBuffer(gContext, CL_MEM_READ_WRITE, 10*sizeof(int), NULL, NULL);
	cl_mem dB = clCreateBuffer(gContext, CL_MEM_READ_WRITE, 10*sizeof(int), NULL, NULL);

    for( int dev = 0; dev < ngpu; ++dev ) {
        for( int s = 0; s < nqueue; ++s ) {
            int t = dev*glog.nqueue + s;
            clEnqueueCopyBuffer(glog.queues[t], dA, dB, 0, 0, sizeof(int), 0, NULL, &(glog.gpu_first[t]));
        }
    }
    // sync again
    for( int dev = 0; dev < ngpu; ++dev ) {
		for(int s = 0; s < nqueue; ++s){
            int t = dev*glog.nqueue + s;
			clFinish(glog.queues[t]);
		}
    }
    glog.cpu_first = get_time();
}


// ----------------------------------------
void trace_cpu_start( int core, const char* tag, const char* lbl )
{
    int id = glog.cpu_id[core];
    glog.cpu_start[core][id] = get_time();
    magma_strlcpy( glog.cpu_tag  [core][id], tag, MAX_LABEL_LEN );
    magma_strlcpy( glog.cpu_label[core][id], lbl, MAX_LABEL_LEN );
}


// ----------------------------------------
void trace_cpu_end( int core )
{
    int id = glog.cpu_id[core];
    glog.cpu_end[core][id] = get_time();
    if ( id+1 < MAX_EVENTS ) {
        glog.cpu_id[core] = id+1;
    }
}


// ----------------------------------------
void trace_gpu_start( int dev, int s, const char* tag, const char* lbl )
{
    int t = dev*glog.nqueue + s;
    int id = glog.gpu_id[t];
#if TRACE_METHOD == 2
    glog.gpu_start[t][id] = get_time();
#else
	//glog.gpu_start[t][id] = NULL;	
	/*
    cudaEventCreate( &glog.gpu_start[t][id] );
    cudaEventRecord(  glog.gpu_start[t][id], glog.streams[t] );
	*/
#endif
    magma_strlcpy( glog.gpu_tag  [t][id], tag, MAX_LABEL_LEN );
    magma_strlcpy( glog.gpu_label[t][id], lbl, MAX_LABEL_LEN );
}


// ----------------------------------------
// pass returned event into kernels to profile
magma_event_t* trace_gpu_end( int dev, int s )
{
    int t = dev*glog.nqueue + s;
    int id = glog.gpu_id[t];
	magma_event_t* re_ptr = &(glog.gpu_end[t][id]);
    if ( id+1 < MAX_EVENTS ) {
        glog.gpu_id[t] = id+1;
    }
	return re_ptr;
}


// ----------------------------------------
void trace_finalize( const char* filename, const char* cssfile )
{
    double xscale = 1e3;
    double height = 20.;  // of each row
    double margin =  5.;  // page margin
    double space  =  1.;  // between rows
    char buf[ 1024 ];
    magma_event_t prof_event;

    // sync devices
	/* sync in src
    for( int dev = 0; dev < glog.ngpu; ++dev ) {
        for( int s = 0; s < glog.nqueue; ++s ) {
            int t = dev*glog.nqueue + s;
			clFinish(glog.queues[t]);
        }
    }
    */
    FILE* trace_file = fopen( filename, "w" );
    if ( trace_file == NULL ) {
        strerror_r( errno, buf, sizeof(buf) );
        fprintf( stderr, "Can't open file '%s': %s (%d)\n", filename, buf, errno );
        return;
    }
    fprintf( stderr, "writing trace to '%s'\n", filename );
    
    int h = (int)( (glog.ncore + glog.ngpu*glog.nqueue)*(height + space) - space + 2*margin );
    int w = (int)( (get_time() - glog.cpu_first) * xscale + 2*margin );
    fprintf( trace_file,
             "<?xml version=\"1.0\" standalone=\"no\"?>\n"
             "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n"
             "    \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\n"
             "<svg version=\"1.1\" baseProfile=\"full\"\n"
             "    xmlns=\"http://www.w3.org/2000/svg\"\n"
             "    xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\"\n"
             "    viewBox=\"0 0 %d %d\" width=\"%d\" height=\"%d\" preserveAspectRatio=\"none\">\n\n",
             w, h, w, h );
    
    // Inkscape does not currently (Jan 2012) support external CSS;
    // see http://wiki.inkscape.org/wiki/index.php/CSS_Support
    // So embed CSS file here
    FILE* css_file = fopen( cssfile, "r" );
    if ( css_file == NULL ) {
        strerror_r( errno, buf, sizeof(buf) );
        fprintf( stderr, "Can't open file '%s': %s; skipping CSS\n", cssfile, buf );
    }
    else {
        fprintf( trace_file, "<style type=\"text/css\">\n" );
        while( fgets( buf, sizeof(buf), css_file ) != NULL ) {
            fputs( buf, trace_file );
        }
        fclose( css_file );
        fprintf( trace_file, "</style>\n\n" );
    }
    
    // format takes: x, y, width, height, class (tag), id (label)
    const char* format = 
        "<rect x=\"%8.3f\" y=\"%4.0f\" width=\"%8.3f\" height=\"%2.0f\" class=\"%-8s\" inkscape:label=\"%s\"/>\n";
    
    // output CPU events
    for( int core = 0; core < glog.ncore; ++core ) {
        if ( glog.cpu_id[core] > MAX_EVENTS ) {
            glog.cpu_id[core] += 1;  // count last event
            fprintf( stderr, "WARNING: trace on core %d, reached limit of %d events; output will be truncated.\n",
                     core, glog.cpu_id[core] );
        }
        fprintf( trace_file, "<!-- core %d, nevents %d -->\n", core, glog.cpu_id[core] );
        fprintf( trace_file, "<g inkscape:groupmode=\"layer\" inkscape:label=\"core %d\">\n", core );
        for( int i = 0; i < glog.cpu_id[core]; ++i ) {
            double start  = glog.cpu_start[core][i] - glog.cpu_first;
            double end    = glog.cpu_end  [core][i] - glog.cpu_first;
            fprintf( trace_file, format,
                     margin + start * xscale,
                     margin + core * (height + space),
                     (end - start) * xscale,
                     height,
                     glog.cpu_tag[core][i],
                     glog.cpu_label[core][i] );
        }
        fprintf( trace_file, "</g>\n\n" );
    }
    
    // output GPU events
    cl_ulong start_time, end_time, first_time;
	double exec_time;
	size_t return_bytes;
	cl_int err;
    //try {
        for( int dev = 0; dev < glog.ngpu; ++dev ) {
		for( int s = 0; s < glog.nqueue; ++s ) {
            int t = dev*glog.nqueue + s;
            if ( glog.gpu_id[t] >= MAX_EVENTS-1 ) {
                glog.gpu_id[t] += 1;  // count last event
                fprintf( stderr, "WARNING: trace on gpu %d/queue %d reached limit of %d events; output will be truncated.\n",
                         dev, s, glog.gpu_id[t] );
            }
            fprintf( trace_file, "<!-- gpu %d, queue %d (t %d), nevents %d -->\n", dev, s, t, glog.gpu_id[t] );
            fprintf( trace_file, "<g inkscape:groupmode=\"layer\" inkscape:label=\"gpu %d queue %d\">\n", dev, s );
            //cudaSetDevice( dev );
			prof_event = glog.gpu_first[t];    
			err = clGetEventProfilingInfo(prof_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &first_time, &return_bytes);
			
			double start, end;

            for( int i = 0; i < glog.gpu_id[t]; ++i ) {
				prof_event = glog.gpu_end[t][i];
				err = clGetEventProfilingInfo(prof_event, CL_PROFILING_COMMAND_START, 
											sizeof(cl_ulong), &start_time, &return_bytes);
				err = clGetEventProfilingInfo(prof_event, CL_PROFILING_COMMAND_END, 
											sizeof(cl_ulong), &end_time, &return_bytes);
				// nanoseconds
				exec_time = (double)(start_time - first_time);
				start = exec_time * 1e-9;

				exec_time = (double)(end_time - first_time);
				end = exec_time * 1e-9;

				//fprtinf start to end
                fprintf( trace_file, format,
                         margin + start * xscale,
                         margin + (dev*glog.nqueue + s + glog.ncore) * (height + space),
                         (end - start) * xscale,
                         height,
                         glog.gpu_tag[t][i],
                         glog.gpu_label[t][i] );
            }
            fprintf( trace_file, "</g>\n\n" );
        }}
    //}
    //catch( cudaError_t err ) {
    //    fprintf( stderr, "CUDA error: %s (%d)\n", cudaGetErrorString( err ), err );
    //}
    
    fprintf( trace_file, "</svg>\n" );
    
    fclose( trace_file );
}

#endif // TRACING
