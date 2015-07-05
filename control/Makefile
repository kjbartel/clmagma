#//////////////////////////////////////////////////////////////////////////////
#   -- clMAGMA (version 1.1.0-beta2) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date November 2013
#//////////////////////////////////////////////////////////////////////////////

DIR = control
MAGMA_DIR = ..
include ../Makefile.internal

ZSRC := \
	zpanel_to_q.cpp	\
	zprint.cpp	\

ZSRC += \
	magma_zf77.cpp		\
	magma_param.f90		\
	magma_zfortran.f90	\
	magma.f90

SRC := \
	auxiliary.cpp	\
	constants.cpp	\
	get_nb_tahiti.cpp	\
	xerbla.cpp	\
	trace.cpp	\
	strlcpy.cpp	\

-include Makefile.local
-include Makefile.src

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC) $(SRC)

ALLOBJ := $(ALLSRC:.cpp=.o)
ALLOBJ := $(ALLOBJ:.f90=.o)
ALLOBJ := $(ALLOBJ:.F90=.o)

.PHONY: all clean cleanall

.DEFAULT_GOAL := all
all: $(LIBMAGMA)

$(LIBMAGMA): stamp
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@
	cp *.mod ../include

# this forces LIBMAGMA to be archived
stamp: $(ALLOBJ)
	touch $@

clean:
	rm -f *.cu_o *.o *~ *.linkinfo *.mod ../include/*.mod sizeptr stamp

cleanall: clean
	rm -f $(LIBMAGMA)

# Extra dependencies for .mod files
magma.o: magma_zfortran.o magma_cfortran.o magma_dfortran.o magma_sfortran.o 
magma_zfortran.o: magma_param.o
magma_cfortran.o: magma_param.o
magma_dfortran.o: magma_param.o
magma_sfortran.o: magma_param.o
