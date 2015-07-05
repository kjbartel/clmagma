/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014
       
       @author Stan Tomov
       @author Mark Gates
       @author Azzam Haidar
*/

#include "magma.h"
#include "common_magma.h"

#ifdef __cplusplus
extern "C" {
#endif

// ==== Definition of blocking sizes for AMD Tahiti cards
#ifdef HAVE_clBLAS

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for potrf based on m
*/
magma_int_t magma_get_spotrf_nb( magma_int_t m )
{
    if      (m <= 1024) return 128;
    else                return 320;
}

magma_int_t magma_get_dpotrf_nb( magma_int_t m )
{
    if      (m <= 4256) return 128;
    else                return 256;
}

magma_int_t magma_get_cpotrf_nb( magma_int_t m )
{
    return 128;
}

magma_int_t magma_get_zpotrf_nb( magma_int_t m )
{
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for geqp3 based on m
*/
magma_int_t magma_get_sgeqp3_nb( magma_int_t m )
{
    return 32;
}

magma_int_t magma_get_dgeqp3_nb( magma_int_t m )
{
    return 32;
}

magma_int_t magma_get_cgeqp3_nb( magma_int_t m )
{
    return 32;
}

magma_int_t magma_get_zgeqp3_nb( magma_int_t m )
{
    return 32;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for geqrf based on m
*/
magma_int_t magma_get_sgeqrf_nb( magma_int_t m )
{
    if      (m <  2000) return 128;
    else                return 128;
}

magma_int_t magma_get_dgeqrf_nb( magma_int_t m )
{
    if      (m <= 2048) return 64;
    else                return 128;
}

magma_int_t magma_get_cgeqrf_nb( magma_int_t m )
{
    if      (m <= 2048) return 32;
    else if (m <= 4032) return 64;
    else                return 128;
}

magma_int_t magma_get_zgeqrf_nb( magma_int_t m )
{
    if      (m <= 2048) return 32;
    else if (m <= 4032) return 64;
    else                return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for geqlf based on m
*/
magma_int_t magma_get_sgeqlf_nb( magma_int_t m )
{
    return magma_get_sgeqrf_nb(m);
}

magma_int_t magma_get_dgeqlf_nb( magma_int_t m )
{
    return magma_get_dgeqrf_nb(m);
}

magma_int_t magma_get_cgeqlf_nb( magma_int_t m )
{
    if      (m <= 2048) return 32;
    else if (m <= 4032) return 64;
    else                return 128;
}

magma_int_t magma_get_zgeqlf_nb( magma_int_t m )
{
    if      (m <= 1024) return 64;
    else                return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for gelqf based on m
*/
magma_int_t magma_get_sgelqf_nb( magma_int_t m )
{
    return magma_get_sgeqrf_nb(m);
}

magma_int_t magma_get_dgelqf_nb( magma_int_t m )
{
    return magma_get_dgeqrf_nb(m);
}

magma_int_t magma_get_cgelqf_nb( magma_int_t m )
{
    if      (m <= 2048) return 32;
    else if (m <= 4032) return 64;
    else                return 128;
}

magma_int_t magma_get_zgelqf_nb( magma_int_t m )
{
    if      (m <= 1024) return 64;
    else                return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for getrf based on m
*/
magma_int_t magma_get_sgetrf_nb( magma_int_t m )
{
    if      (m <= 3200) return 128;
    else if (m <  9000) return 256;
    else                return 320;
}

magma_int_t magma_get_dgetrf_nb( magma_int_t m )
{
    if      (m <= 2048) return 64;
    else if (m <  7200) return 192;
    else                return 256;
}

magma_int_t magma_get_cgetrf_nb( magma_int_t m )
{
    if      (m <= 2048) return 64;
    else                return 128;
}

magma_int_t magma_get_zgetrf_nb( magma_int_t m )
{
    if      (m <= 3072) return 32;
    else if (m <= 9024) return 64;
    else                return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for gehrd based on m
*/
magma_int_t magma_get_sgehrd_nb( magma_int_t m )
{
    if      (m <= 1024) return 32;
    else                return 96;
}

magma_int_t magma_get_dgehrd_nb( magma_int_t m )
{
    if      (m <= 2048) return 32;
    else                return 64;
}

magma_int_t magma_get_cgehrd_nb( magma_int_t m )
{
    if      (m <= 1024) return 32;
    else                return 64;
}

magma_int_t magma_get_zgehrd_nb( magma_int_t m )
{
    if      (m <= 2048) return 32;
    else                return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sytrd based on m
*/
magma_int_t magma_get_ssytrd_nb( magma_int_t m )
{
    return 32;
}

magma_int_t magma_get_dsytrd_nb( magma_int_t m )
{
    return 32;
}

magma_int_t magma_get_chetrd_nb( magma_int_t m )
{
    return 32;
}

magma_int_t magma_get_zhetrd_nb( magma_int_t m )
{
    return 32;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sytrf based on m
  */
magma_int_t magma_get_zhetrf_nb( magma_int_t m )
{
    return 256;
}

magma_int_t magma_get_chetrf_nb( magma_int_t m )
{
    return 256;
}

magma_int_t magma_get_dsytrf_nb( magma_int_t m )
{
    return 96;
}

magma_int_t magma_get_ssytrf_nb( magma_int_t m )
{
    return 256;
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for gebrd based on m
*/
magma_int_t magma_get_sgebrd_nb( magma_int_t m )
{
    return 32;
}

magma_int_t magma_get_dgebrd_nb( magma_int_t m )
{
    return 32;
}

magma_int_t magma_get_cgebrd_nb( magma_int_t m )
{
    return 32;
}

magma_int_t magma_get_zgebrd_nb( magma_int_t m )
{
    return 32;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sygst based on m
*/
magma_int_t magma_get_ssygst_nb( magma_int_t m )
{
    return 64;
}

magma_int_t magma_get_dsygst_nb( magma_int_t m )
{
    return 64;
}

magma_int_t magma_get_chegst_nb( magma_int_t m )
{
    return 64;
}

magma_int_t magma_get_zhegst_nb( magma_int_t m )
{
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for getri based on m
*/
magma_int_t magma_get_sgetri_nb( magma_int_t m )
{
    return 64;
}

magma_int_t magma_get_dgetri_nb( magma_int_t m )
{
    return 64;
}

magma_int_t magma_get_cgetri_nb( magma_int_t m )
{
    return 64;
}

magma_int_t magma_get_zgetri_nb( magma_int_t m )
{
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for gesvd based on m
*/
magma_int_t magma_get_sgesvd_nb( magma_int_t m )
{
    return magma_get_sgebrd_nb(m);
}

magma_int_t magma_get_dgesvd_nb( magma_int_t m )
{
    return magma_get_dgebrd_nb(m);
}

magma_int_t magma_get_cgesvd_nb( magma_int_t m )
{
    return magma_get_cgebrd_nb(m);
}

magma_int_t magma_get_zgesvd_nb( magma_int_t m )
{
    return magma_get_zgebrd_nb(m);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return smlsiz for the divide and conquewr routine dlaex0 dstedx zstedx
*/
magma_int_t magma_get_smlsize_divideconquer()
{
    return 128;
}

#endif  // HAVE_clBLAS

#ifdef __cplusplus
} // extern "C"
#endif
