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
	zungqr.cpp	\
        \
	zlatrd.cpp	\
	zhetrd.cpp	\
        \
        zgebrd.cpp      \
        zlabrd_gpu.cpp  \

# Matrix Inverse
ZSRC += \
	ztrtri_gpu.cpp  \
	zgetri_gpu.cpp  \
	zpotri_gpu.cpp  \
	zlauum_gpu.cpp  \

# SVD
ZSRC += \
    dgesvd.cpp  \
	zgesvd.cpp  \

# symmetric eigenvalues, CPU interface
ZSRC += \
	dsyevd.cpp      \
	zheevd.cpp      \
		\
	zstedx.cpp      \
	dstedx.cpp      \
	dlaex0.cpp      \
	dlaex1.cpp      \
	dlaex3.cpp      \
	zunmtr.cpp      \
	zunmqr.cpp      \
	zunmql.cpp      \

# non-symmetric eigenvalues
ZSRC += \
	dgeev.cpp   \
	zgeev.cpp   \
	zgehrd.cpp  \
	zlahru.cpp  \
	zlahr2.cpp  \
	zunghr.cpp  \

-include Makefile.src

ALLSRC_F = dlaed2.f	\
			slaed2.f

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.cpp=.o) $(ALLSRC_F:.f=.o)

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
