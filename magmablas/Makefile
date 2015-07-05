#//////////////////////////////////////////////////////////////////////////////
#   -- clMAGMA (version 1.0.0) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      April 2012
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

INC += -I$(MAGMA_DIR)/control

ZHDR  =

ZSRC  =

-include Makefile.src

ALLHDR = $(ZHDR) $(CHDR) $(DHDR) $(SHDR)

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC) \
	dummy.cpp

ifeq (${GPU_TARGET}, Tesla)
GPUOPT = ${TESLAOPT}
else ifeq (${GPU_TARGET}, Fermi)
GPUOPT = ${FERMIOPT}
endif

ALLOBJ = $(ALLSRC:.cpp=.o)

all: $(LIBMAGMABLAS)

$(LIBMAGMABLAS): $(ALLHDR) $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

%.o : %.cpp
	$(CC) $(OPTS) $(INC) -c $< -o $@

%.cu_o: %.cu
	$(NVCC) $(NVOPTS) $(GPUOPT) $(INC) -c $< -o $@

clean:
	rm -f *.cu_o *.o *~ *.linkinfo

cleanall: clean
	rm -f *.a
	rm -f $(LIBMAGMABLAS)
