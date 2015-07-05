#//////////////////////////////////////////////////////////////////////////////
#   -- clMAGMA (version 0.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      April 2012
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

ZSRC = \
	testing_zpotrf_gpu.cpp	\
	testing_zgeqrf_gpu.cpp  \
	testing_zgetrf_gpu.cpp	\
        \
        testing_zgemm.cpp

#ifeq ($(FORT), pgfortran)
#	FOBJ := fortran_thunking.o
#	ZSRCF = testing_zgetrf_gpu_f.cuf
#else
#	FOBJ := fortran.o
#	ZSRCF = testing_zgetrf_gpu_f.f90	\
#		testing_zgetrf_f.f90
#endif

ZSRC += $(ZSRCF)
-include Makefile.src

INC += -I../interface_opencl/cl_magma_runtime/

ALLSRC  = $(ZSRC) $(CSRC) $(DSRC) $(SSRC) \
	testing_constants.cpp
ALLOBJF  = $(ALLSRC:.f90=.o)
ALLOBJF := $(ALLOBJF:.cuf=.o)
ALLOBJ   = $(ALLOBJF:.cpp=.o)
EXE      = $(ALLOBJ:.o=)

all: $(EXE)

#LIBLAPACKTEST = lin/liblapacktest.a
#
#$(LIBLAPACKTEST):
#	cd lin && make liblapacktest.a

$(EXE): $(LIBMAGMA) $(LIBMAGMABLAS) #$(LIBLAPACKTEST)

clean:
	rm -f *.o *~

cleanall: clean
	rm -f $(EXE)

#ifeq ($(FORT), pgfortran)
#testing_%.o : testing_%.cuf
#	$(FORT) $(FOPTS) $(INC) -c $< -o $@ 
#else
#testing_%.o : testing_%.f90
#	$(FORT) $(FOPTS) $(INC) -c $< -o $@ 
#endif

%.o: %.cpp
	$(CC) $(OPTS) $(INC) -c $< -o $@

#
# MAGMA Fortran module. Once finalized, the generated magma.mod
# has to be moved in the include directory
#
# lin/liblapacktest.a
testing_%_f : testing_%_f.o $(FOBJ)
	$(FORT) $(OPTS) $(LDOPTS) $< $(FOBJ) -o $@ \
	-L../lib -lclmagma -lclmagmablas $(LIBDIR) $(LIB) 

# lin/liblapacktest.a
testing_% : testing_%.o
	$(CC) $(OPTS) $(LDOPTS) $(NOMAIN) $< $(OBJ) -o $@ \
	-L../lib -lclmagma -lclmagmablas $(LIBDIR) $(LIB) 

#
# These two files are provided by cuda
#
#fortran.o : $(CUDADIR)/src/fortran.c
#	$(CC) $(OPTS) $(INC) -c $< -o $@
#
#fortran_thunking.o : $(CUDADIR)/src/fortran_thunking.c
#	$(CC) $(OPTS) $(INC) -c $< -o $@

# keep intermediate object files, to speed up re-compiling when library changes
.SECONDARY: $(ALLOBJ)