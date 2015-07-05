#include <stdio.h>
#include <stdlib.h>

#include "CL_MAGMA_RT.h"

int main( int argc, char** argv )
{
    if ( argc != 2 ) {
        printf( "Usage: %s file.cl\n", argv[0] );
        exit(1);
    }
    
    CL_MAGMA_RT *runtime = CL_MAGMA_RT::Instance();
    runtime->Init();
    runtime->CompileFile( argv[1] );
    runtime->Quit();

    return 0;
}
