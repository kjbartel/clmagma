/*
    -- clMAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014
*/

#include "common_magma.h"

// -------------------------
// Returns version of MAGMA, as defined by
// MAGMA_VERSION_MAJOR, MAGMA_VERSION_MINOR, MAGMA_VERSION_MICRO constants.
extern "C"
void magma_version( magma_int_t* major, magma_int_t* minor, magma_int_t* micro )
{
    if ( major != NULL && minor != NULL && micro != NULL ) {
        *major = MAGMA_VERSION_MAJOR;
        *minor = MAGMA_VERSION_MINOR;
        *micro = MAGMA_VERSION_MICRO;
    }
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Auxiliary function: ipiv(i) indicates that row i has been swapped with
      ipiv(i) from top to bottom. This function rearranges ipiv into newipiv
      where row i has to be moved to newipiv(i). The new pivoting allows for
      parallel processing vs the original one assumes a specific ordering and
      has to be done sequentially.
*/
extern "C"
void swp2pswp( magma_trans_t trans, magma_int_t n, magma_int_t *ipiv, magma_int_t *newipiv)
{
  magma_int_t i, newind, ind;

  for(i=0; i<n; i++)
    newipiv[i] = -1;
  
  if ( trans == MagmaNoTrans ){
    for(i=0; i<n; i++){
      newind = ipiv[i] - 1;
      if (newipiv[newind] == -1) {
        if (newipiv[i]==-1){
          newipiv[i] = newind;
          if (newind>i)
            newipiv[newind]= i;
        }
        else
          {
            ind = newipiv[i];
            newipiv[i] = newind;
            if (newind>i)
              newipiv[newind]= ind;
          }
      }
      else {
        if (newipiv[i]==-1){
          if (newind>i){
            ind = newipiv[newind];
            newipiv[newind] = i;
            newipiv[i] = ind;
          }
          else
            newipiv[i] = newipiv[newind];
        }
        else{
          ind = newipiv[i];
          newipiv[i] = newipiv[newind];
          if (newind > i)
            newipiv[newind] = ind;
        }
      }
    }
  } else {
    for(i=n-1; i>=0; i--){
      newind = ipiv[i] - 1;
      if (newipiv[newind] == -1) {
        if (newipiv[i]==-1){
          newipiv[i] = newind;
          if (newind>i)
            newipiv[newind]= i;
        }
        else
          {
            ind = newipiv[i];
            newipiv[i] = newind;
            if (newind>i)
              newipiv[newind]= ind;
          }
      }
      else {
        if (newipiv[i]==-1){
          if (newind>i){
            ind = newipiv[newind];
            newipiv[newind] = i;
            newipiv[i] = ind;
          }
          else
            newipiv[i] = newipiv[newind];
        }
        else{
          ind = newipiv[i];
          newipiv[i] = newipiv[newind];
          if (newind > i)
            newipiv[newind] = ind;
        }
      }
    }
  }
}
