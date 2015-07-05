#//////////////////////////////////////////////////////////////////////////////
#   -- clMAGMA (version 0.1) --
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

ALLSRC = dummy.cpp $(ZSRC) $(CSRC) $(DSRC) $(SSRC)

ifeq (${GPU_TARGET}, 0)
GPUOPT = ${TESLAOPT}
endif

ifeq (${GPU_TARGET}, 1)
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
