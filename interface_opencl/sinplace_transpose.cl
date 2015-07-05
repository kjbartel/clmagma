/*
    -- clMAGMA (version 0.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
        Univ. of Colorado, Denver
       April 2012

       @generated s Fri Nov 11 17:02:30 2011

*/
// #include "common_magma.h"
// #define PRECISION_s
// #include "commonblas.h"


#define NB 16

#define NBhalf (NB/2)
#define __mul24( x, y )  ((x)*(y))

//#if NB == 32
#define COPY_1D_TO_2D( a, lda, b, inx, iny )        \
    b[iny][inx]			= a[0];                        \
    b[iny+NBhalf][inx]	= a[NBhalf*lda];

#define COPY_2D_TO_1D( b, inx, iny, a, lda )        \
    a[0]			= b[inx][iny];                        \
    a[NBhalf*lda]	= b[inx][iny+NBhalf];

//#else

/*
#define COPY_1D_TO_2D( a, lda, b, inx, iny )        \
    b[iny][inx] = a[0];

#define COPY_2D_TO_1D( b, inx, iny, a, lda )        \
    a[0]        = b[inx][iny];

#endif
*/



__kernel void sinplace_T_even_kernel(__global float *matrix, int offset, int lda, int half )
{     
	__local float a[NB][NB+1];
	__local float b[NB][NB+1];

	int inx = get_local_id(0);
	int iny = get_local_id(1);

	bool bottom = ( get_group_id(0) > get_group_id(1) );
	int ibx = bottom ? (get_group_id(0) - 1) : (get_group_id(1) + half);
	int iby = bottom ? get_group_id(1)       : (get_group_id(0) + half);

	ibx *= NB;
	iby *= NB;

	matrix += offset;

	__global float *A = matrix + ibx + inx + (iby + iny)*lda;
	COPY_1D_TO_2D( A, lda, a, inx, iny);

	if( ibx == iby )
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		COPY_2D_TO_1D( a, inx, iny, A, lda);
	}
	else
	{
		__global float *B = matrix + iby + inx + __mul24( ibx + iny, lda );

		COPY_1D_TO_2D( B, lda, b, inx, iny);
		barrier(CLK_LOCAL_MEM_FENCE);

		COPY_2D_TO_1D( b, inx, iny, A, lda);
		COPY_2D_TO_1D( a, inx, iny, B, lda);
	}
} 

__kernel void sinplace_T_odd_kernel (__global float *matrix, int offset, int lda, int half )
{        
        __local float a[NB][NB+1];
        __local float b[NB][NB+1];
        
        int inx = get_local_id(0);
        int iny = get_local_id(1);

        bool bottom = ( get_group_id(0) >= get_group_id(1) );
        int ibx = bottom ? get_group_id(0)  : (get_group_id(1) + half - 1);
        int iby = bottom ? get_group_id(1)  : (get_group_id(0) + half);

        ibx *= NB;
        iby *= NB;

		matrix += offset;

        __global float *A = matrix + ibx + inx + __mul24( iby + iny, lda );

        COPY_1D_TO_2D( A, lda, a, inx, iny);

        if( ibx == iby )
        {
                barrier(CLK_LOCAL_MEM_FENCE);
                COPY_2D_TO_1D( a, inx, iny, A, lda);
        }
        else
        {
                __global float *B = matrix + iby + inx + __mul24( ibx + iny, lda );

                COPY_1D_TO_2D( B, lda, b, inx, iny);
                barrier(CLK_LOCAL_MEM_FENCE);

                COPY_2D_TO_1D( b, inx, iny, A, lda);
                COPY_2D_TO_1D( a, inx, iny, B, lda);
        }
} 
