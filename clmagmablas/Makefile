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

ZHDR += \
	clag2z.h		\
	zcaxpycp.h		\
	zlacpy.h		\
	zlag2c.h		\
	zlanhe.h		\
	zlaset.h		\
	zlaswp.h		\
	zlat2c.h		\
	zswap.h			\
	ztranspose.h		\
	ztranspose_inplace.h	\

# alphabetic order by base name (ignoring precision)
# OpenCL sources need to be here for precision generation,
# but later we separate them out into CL_SRC
ZSRC += \
	clag2z.cl		\
	clag2z.cpp		\
	empty.cl		\
	empty.cpp		\
	zcaxpycp.cl		\
	zcaxpycp.cpp		\
	zlacpy.cl		\
	zlacpy.cpp		\
	zlag2c.cl		\
	zlag2c.cpp		\
	zlanhe.cl		\
	zlanhe.cpp		\
	zlaset.cl		\
	zlaset.cpp		\
	zlaswp.cl		\
	zlaswp.cpp		\
	zlat2c.cl		\
	zlat2c.cpp		\
	zswap.cl		\
	zswap.cpp		\
	ztranspose.cl		\
	ztranspose.cpp		\
	ztranspose_inplace.cl	\
	ztranspose_inplace.cpp	\

# to re-convert from CUDA
ZSRC += \
	dznrm2.cl		\
	dznrm2.cpp		\
	zbcyclic.cpp		\
	zgemm_reduce.cl		\
	zgemm_reduce.cpp	\
	zlacpy_cnjg.cl		\
	zlacpy_cnjg.cpp		\
	zlarfbx.cl		\
	zlarfbx.cpp		\
	zlarfgx-v2.cl		\
	zlarfgx-v2.cpp		\
	zlarfx.cl		\
	zlascl_2x2.cl		\
	zlascl_2x2.cpp		\

-include Makefile.local
-include Makefile.src

ALLHDR := $(ZHDR) $(CHDR) $(DHDR) $(SHDR) $(HDR)

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC) $(SRC)

ALLOBJ := $(ALLSRC:.cpp=.o)
ALLOBJ := $(filter-out %.cl, $(ALLOBJ)) # don't put .cl files into lib

CL_SRC := $(filter %.cl, $(ALLSRC))
CL_OBJ := $(CL_SRC:.cl=.co)

OPTS += -I ../interface_opencl

$(CL_COMPILER):
	cd ../interface_opencl && $(MAKE) $(CL_COMPILER)

.PHONY: all clean cleanall

.DEFAULT_GOAL := all
all: $(LIBMAGMA) $(MAGMA_KERNELS)

$(ALLOBJ): $(ALLHDR)
$(CL_OBJ): $(ALLHDR)

%.co: %.cl $(CL_COMPILER)
	$(CL_COMPILER) -c $<

$(MAGMA_KERNELS): $(CL_OBJ)
	$(CL_COMPILER) -a -o $@ $^

$(LIBMAGMA): $(ALLOBJ) stamp
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

# this forces LIBMAGMA to be archived
stamp: $(ALLOBJ)
	touch $@

clean:
	rm -f *.o *~ *.linkinfo stamp
	rm -f *.co

cleanall: clean
	rm -f $(LIBMAGMA)
