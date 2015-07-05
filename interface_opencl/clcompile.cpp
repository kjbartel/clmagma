#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "clmagma_runtime.h"


const char* usage =
    "Usage: %s {-a|-c} [-k] [-o output_file] input_files\n"
    "  -a  archive multiple .co files into one .co file.\n"
    "  -c  compile .cl files into .co files.\n"
    "  -k  keep going; continue after errors.\n"
    "  -o  specify output file.\n";

typedef enum {
    mode_none,
    mode_archive,
    mode_compile
} clcompile_mode_t;


// ------------------------------------------------------------
int main( int argc, char** argv )
{
    // parse command line arguments
    clcompile_mode_t mode        = mode_none;
    bool             keep_going  = false;
    const char*      output_file = NULL;
    
    const char* cmd = argv[0];
    
    int opt;
    while( (opt = getopt( argc, argv, "acko:" )) != -1 ) {
        switch( opt ) {
            case 'a': mode = mode_archive;   break;
            case 'c': mode = mode_compile;   break;
            case 'k': keep_going = true;     break;
            case 'o': output_file = optarg;  break;
            
            case '?':
            default:
                fprintf( stderr, usage, cmd );
                exit(1);
                break;
        }
    }
    argc -= optind;
    argv += optind;
    
    int err = 0;
    if ( mode == mode_compile ) {
        // compile one or more .cl files into .co files
        if ( argc == 0 ) {
            fprintf( stderr, "Error: no input files\n" );
            exit(1);
        }
        if ( output_file != NULL && argc > 1 ) {
            fprintf( stderr, "Error: cannot specify -o when generating multiple output files\n" );
            exit(1);
        }
        g_runtime.init( true );
        for( int i=0; i < argc; ++i ) {
            if ( g_runtime.compile_file( argv[i], output_file ) != 0 ) {
                err = 1;
                if ( ! keep_going ) {
                    break;
                }
            }
        }
        g_runtime.quit();
    }
    else if ( mode == mode_archive ) {
        // archive multiple .co files into one combined .co file
        if ( argc == 0 ) {
            fprintf( stderr, "Error: no archive members specified\n" );
            exit(1);
        }
        if ( output_file == NULL ) {
            fprintf( stderr, "Error: no output file given\n" );
            exit(1);
        }
        g_runtime.init( true );
        g_runtime.archive_files( argc, argv, output_file );
        g_runtime.quit();
    }
    else {
        fprintf( stderr, usage, cmd );
        exit(1);
    }
    
    return err;
}
