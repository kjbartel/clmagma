#//////////////////////////////////////////////////////////////////////////////
#   -- clMAGMA (version 1.3.0) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date November 2014
#//////////////////////////////////////////////////////////////////////////////

DIR = interface_opencl
MAGMA_DIR = ..
include $(MAGMA_DIR)/Makefile.internal

# alphabetic order by base name (ignoring precision)
ZSRC := \
	blas_z.cpp		\
	zset_get.cpp		\

SRC := \
	alloc.cpp		\
	error.cpp		\
	interface.cpp		\
	set_get.cpp		\
	clmagma_runtime.cpp	\

-include Makefile.local
-include Makefile.src

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC) $(SRC)

ALLOBJ := $(ALLSRC:.cpp=.o)

echo:
	@echo ALLSRC $(ALLSRC)
	@echo ALLOBJ $(ALLOBJ)
	@echo CL_SRC $(CL_SRC)
	@echo CL_OBJ $(CL_OBJ)
	@echo CL_COMPILER $(CL_COMPILER)

.PHONY: all clean cleanall

.DEFAULT_GOAL := all
all: $(LIBMAGMA) $(CL_COMPILER)

$(CL_COMPILER): clcompile.o clmagma_runtime.o error.o
	$(CC) $(LDOPTS) -o $@ $^ $(LIBDIR) $(LIB)

$(LIBMAGMA): $(ALLOBJ) stamp
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

# this forces LIBMAGMA to be archived
stamp: $(ALLOBJ)
	touch $@

clean:
	rm -f *.o *~ *.linkinfo stamp

cleanall: clean
	rm -f $(LIBMAGMA)
	rm -f $(CL_COMPILER)
