#//////////////////////////////////////////////////////////////////////////////
#   -- clMAGMA (version 0.3.0) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      April 2012
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

INC += -I$(MAGMA_DIR)/control

ZSRC = \
	zpotrf_gpu.cpp	\
	zlarfb_gpu.cpp  \
	zgeqrf2_gpu.cpp	\
	zgetrf_gpu.cpp	\
	zunmqr_gpu.cpp	\
	zgeqrf_gpu.cpp  \
        \
	zpotrs_gpu.cpp	\
	zposv_gpu.cpp	\
	zgetrs_gpu.cpp	\
	zgesv_gpu.cpp	\
	zgeqrs_gpu.cpp	\
	zgels_gpu.cpp	\
        \
	zgehrd.cpp	\
	zlahru.cpp	\
	zlahr2.cpp	\
	zunghr.cpp	\
	zungqr.cpp	\
        \
	zlatrd.cpp	\
	zhetrd.cpp	\
        \
        zgebrd.cpp      \
        zlabrd_gpu.cpp  \

-include Makefile.src

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.cpp=.o) $(ALLSRC:.f=.o)

all: $(LIBMAGMA)

$(LIBMAGMA): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

.cpp.o:
	$(CC) $(OPTS) $(INC) -c $< -o $@

.f.o:
	$(FORT) $(OPTS) $(INC) -c $< -o $@

clean:
	rm -f *.o *~ *.linkinfo

cleanall: clean
	rm -f *.a
	rm -f $(LIBMAGMA)
