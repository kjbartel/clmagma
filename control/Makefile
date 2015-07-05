#//////////////////////////////////////////////////////////////////////////////
#   -- clMAGMA (version 0.2.0) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      April 2012
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

ZSRC := \
	zpanel_to_q.cpp	\
	zprint.cpp	\

-include Makefile.src

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC) \
	auxiliary.cpp	\
	constants.cpp	\
	get_nb_tesla.cpp	\
	get_nb_fermi.cpp	\
	get_nb_tahiti.cpp	\
	xerbla.cpp	\

ALLOBJ  = $(ALLSRC:.cpp=.o)
ALLOBJ := $(ALLOBJ:.f90=.o)

all: $(LIBMAGMA)

$(LIBMAGMA): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

#-cp *.mod ../include  # use if has Fortran files

.cpp.o:
	$(CC) $(OPTS) $(INC) -c $< -o $@

%.o : %.f90 
	$(FORT) $(FOPTS) $(INC) -c $< -o $@

clean:
	rm -f *.cu_o *.o *~ *.linkinfo *.mod

cleanall: clean
	rm -f *.a
	rm -f $(LIBMAGMA)
