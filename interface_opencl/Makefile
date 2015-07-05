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

ZSRC  = magmablas_z.cpp

-include Makefile.src
-include Makefile.local

ALLSRC = \
	interface.cpp \
	strerror.cpp \
	sinplace_transpose.cpp \
	spermute.cpp     \
	stranspose.cpp   \
	stranspose2.cpp  \
	CL_MAGMA_RT.cpp  \
	$(ZSRC) $(CSRC) $(DSRC) $(SSRC)

ALLOBJ = $(ALLSRC:.cpp=.o)

CL_SRCS = \
	sinplace_transpose.cl	\
	spermute-v2.cl		\
	stranspose.cl		\
	stranspose-v2.cl	\

CL_OBJS = $(CL_SRCS:.cl=.co)

CL_COMPILER = clcompile

all: $(LIBMAGMA) $(CL_OBJS) $(CL_COMPILER)

$(CL_COMPILER): clcompile.o $(LIBMAGMA)
	$(CC) $(LDOPTS) $< -o $@ \
	-L../lib -lclmagma $(LIBDIR) $(LIB)

%.co: %.cl $(CL_COMPILER)
	./$(CL_COMPILER) $<

$(LIBMAGMA): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

.cpp.o:
	$(CC) $(OPTS) $(INC) -c $< -o $@

.f.o:
	$(FORT) $(OPTS) $(INC) -c $< -o $@

clean:
	rm -f *.o *~ *.linkinfo $(CL_COMPILER)

cleanall: clean
	rm -f *.a
	rm -f $(LIBMAGMA)
