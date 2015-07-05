/*
    -- clMAGMA (version 1.1.0-beta2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2013
*/

#ifndef MAGMA_AUXILIARY_H
#define MAGMA_AUXILIARY_H

#include "magma_types.h"

typedef struct magma_timestr_s
{
  unsigned int sec;
  unsigned int usec;
} magma_timestr_t;

#ifdef __cplusplus
extern "C" {
#endif

magma_timestr_t get_current_time(void);
double GetTimerValue(magma_timestr_t time_1, magma_timestr_t time_2);
double get_time(void);

void magma_print_devices();

void spanel_to_q(int uplo, int ib, float *a, int lda, float *work);
void sq_to_panel(int uplo, int ib, float *a, int lda, float *work);

void swp2pswp   (int trans, int n, int *ipiv, int *newipiv);

void cpanel_to_q(int uplo, int ib, magmaFloatComplex *a, int lda, magmaFloatComplex *work);
void cq_to_panel(int uplo, int ib, magmaFloatComplex *a, int lda, magmaFloatComplex *work);

void dpanel_to_q(int uplo, int ib, double *a, int lda, double *work);
void dq_to_panel(int uplo, int ib, double *a, int lda, double *work);

void zpanel_to_q(int uplo, int ib, magmaDoubleComplex *a, int lda, magmaDoubleComplex *work);
void zq_to_panel(int uplo, int ib, magmaDoubleComplex *a, int lda, magmaDoubleComplex *work);

double magma_cabs(magmaDoubleComplex x);
float magma_cabsf(magmaFloatComplex x);

#ifdef __cplusplus
}
#endif

#endif  // MAGMA_AUXILIARY_H
