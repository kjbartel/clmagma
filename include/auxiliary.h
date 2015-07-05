/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014
*/

#ifndef MAGMA_AUXILIARY_H
#define MAGMA_AUXILIARY_H

#include "magma_types.h"

/* ------------------------------------------------------------
 *   -- MAGMA Auxiliary structures and functions
 * --------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

real_Double_t magma_wtime( void );
real_Double_t magma_sync_wtime( magma_queue_t queue );

size_t magma_strlcpy(char *dst, const char *src, size_t siz);

magma_int_t magma_num_gpus( void );

double magma_cabs(magmaDoubleComplex x);
float  magma_cabsf(magmaFloatComplex x);

void magma_print_environment();

void swp2pswp(magma_trans_t trans, magma_int_t n, magma_int_t *ipiv, magma_int_t *newipiv);

#ifdef __cplusplus
}
#endif

#endif  // MAGMA_AUXILIARY_H
