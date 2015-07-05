#//////////////////////////////////////////////////////////////////////////////
#   -- clMAGMA (version 1.1.0-beta2) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      @date November 2013
#//////////////////////////////////////////////////////////////////////////////

DIR = interface_opencl
MAGMA_DIR = ..
include ../Makefile.internal

ZSRC := \
	magmablas_z.cpp		\
	zpermute.cpp		\
	zinplace_transpose.cpp	\
	ztranspose.cpp		\
	ztranspose2.cpp		\
	zauxiliary.cpp		\
	zlacpy.cpp		\
	zswap.cpp		\
	dznrm2.cpp		\
	zbcyclic.cpp		\
	zgemm_reduce.cpp	\
	zlarfbx.cpp		\
	zlarfgx-v2.cpp		\
	\
	zinplace_transpose.cl	\
	zpermute-v2.cl		\
	ztranspose.cl		\
	ztranspose-v2.cl	\
	zauxiliary.cl		\
	zlacpy.cl		\
	zswap.cl		\
	dznrm2.cl		\
	zgemm_reduce.cl		\
	zlarfbx.cl		\
	zlarfx.cl		\
	zlarfgx-v2.cl		\

SRC := \
	interface.cpp	\
	strerror.cpp	\
	CL_MAGMA_RT.cpp	\

-include Makefile.local
-include Makefile.src

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC) $(SRC)

ALLOBJ := $(ALLSRC:.cpp=.o)

CL_SRCS = \
	sinplace_transpose.cl	\
	spermute-v2.cl		\
	stranspose.cl		\
	stranspose-v2.cl	\
	sauxiliary.cl		\
	slacpy.cl		\
	sswap.cl		\
	snrm2.cl		\
	\
	dinplace_transpose.cl	\
	dpermute-v2.cl		\
	dtranspose.cl		\
	dtranspose-v2.cl	\
	dauxiliary.cl		\
	dlacpy.cl		\
	dswap.cl		\
	dnrm2.cl		\
	\
	cinplace_transpose.cl	\
	cpermute-v2.cl		\
	ctranspose.cl		\
	ctranspose-v2.cl	\
	cauxiliary.cl		\
	clacpy.cl		\
	cswap.cl		\
	scnrm2.cl		\
	\
	zinplace_transpose.cl	\
	zpermute-v2.cl		\
	ztranspose.cl		\
	ztranspose-v2.cl	\
	zauxiliary.cl		\
	zlacpy.cl		\
	zswap.cl		\
	\
	dznrm2.cl		\
	sgemm_reduce.cl		\
	dgemm_reduce.cl		\
	cgemm_reduce.cl		\
	zgemm_reduce.cl		\
	slarfbx.cl		\
	dlarfbx.cl		\
	clarfbx.cl		\
	zlarfbx.cl		\
	slarfx.cl		\
	dlarfx.cl		\
	clarfx.cl		\
	zlarfx.cl		\
	slarfgx-v2.cl		\
	dlarfgx-v2.cl		\
	clarfgx-v2.cl		\
	zlarfgx-v2.cl		\

CL_OBJS = $(CL_SRCS:.cl=.co)

CL_COMPILER = clcompile

.PHONY: all clean cleanall

.DEFAULT_GOAL := all
all: $(LIBMAGMA) $(CL_OBJS)

$(CL_COMPILER): clcompile.o CL_MAGMA_RT.o
	$(CC) $(LDOPTS) -o $@ $^ \
	$(LIBDIR) $(LIB)

%.co: %.cl $(CL_COMPILER)
	./$(CL_COMPILER) $<

$(LIBMAGMA): $(ALLOBJ) stamp
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

# this forces LIBMAGMA to be archived
stamp: $(ALLHDR) $(ALLOBJ)
	touch $@

clean:
	rm -f *.o *~ *.linkinfo stamp
	rm -f *.co $(CL_COMPILER)

cleanall: clean
	rm -f $(LIBMAGMA)
