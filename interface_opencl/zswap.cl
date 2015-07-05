/*
 *   -- clMAGMA (version 1.0.0) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      August 2012
 *
 * @precisions normal z -> s d c 
 */

#define PRECISION_z
#define BLOCK_SIZE 64
#define __mul24( x, y )  ((x)*(y))

#if defined(PRECISION_c) || defined(PRECISION_z)
typedef double2 magmaDoubleComplex;
#endif

typedef struct {
    int n, offset_dA1, lda1, offset_dA2, lda2;
} magmagpu_zswap_params_t;

__kernel void magmagpu_zswap(__global magmaDoubleComplex *dA1, __global magmaDoubleComplex *dA2, magmagpu_zswap_params_t params )
{
    unsigned int x = get_local_id(0) + __mul24(get_local_size(0), get_group_id(0));
	unsigned int offset1 = __mul24( x, params.lda1);
	unsigned int offset2 = __mul24( x, params.lda2);
	if( x < params.n ){
		__global magmaDoubleComplex *A1  = dA1 + params.offset_dA1 + offset1;
		__global magmaDoubleComplex *A2  = dA2 + params.offset_dA2 + offset2;
		magmaDoubleComplex temp = *A1;
		*A1 = *A2;
		*A2 = temp;
	}
}
