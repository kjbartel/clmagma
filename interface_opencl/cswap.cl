/*
 *   -- clMAGMA (version 1.0.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      August 2012
 *
 * @generated c Wed Oct 24 00:32:59 2012
 */

#define PRECISION_c
#define BLOCK_SIZE 64
#define __mul24( x, y )  ((x)*(y))

#if defined(PRECISION_c) || defined(PRECISION_z)
typedef float2 magmaFloatComplex;
#endif

typedef struct {
    int n, offset_dA1, lda1, offset_dA2, lda2;
} magmagpu_cswap_params_t;

__kernel void magmagpu_cswap(__global magmaFloatComplex *dA1, __global magmaFloatComplex *dA2, magmagpu_cswap_params_t params )
{
    unsigned int x = get_local_id(0) + __mul24(get_local_size(0), get_group_id(0));
	unsigned int offset1 = __mul24( x, params.lda1);
	unsigned int offset2 = __mul24( x, params.lda2);
	if( x < params.n ){
		__global magmaFloatComplex *A1  = dA1 + params.offset_dA1 + offset1;
		__global magmaFloatComplex *A2  = dA2 + params.offset_dA2 + offset2;
		magmaFloatComplex temp = *A1;
		*A1 = *A2;
		*A2 = temp;
	}
}
