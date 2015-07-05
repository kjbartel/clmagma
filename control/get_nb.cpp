/*
    -- clMAGMA (version 0.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       April 2012
*/

// ==== Definition of blocking sizes ===============================
#ifdef HAVE_CUBLAS
  #ifdef (GPUSHMEM < 200)
    // ==== Definitions for Tesla 
    #include "get_nb_tesla.cpp"
  #else
    // ==== Definitions for Fermi
    #include "get_nb_fermi.cpp"
  #endif
#endif

#ifdef HAVE_AMDBLAS
  // ==== Definitions for AMD Tahiti
  #include "get_nb_tahiti.cpp"
#endif

