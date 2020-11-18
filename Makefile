#CC=/hpc/packages/minerva-common/gcc/7.3.0/bin/g++

CC=g++

DIC=$(PWD)
CFLAGS= -c -Wall -g -O3  -I $(DIC)








LDFLAGS= -DWITH_LAPACK -m64 -static  -lgsl  -lzlib  -latas -lboost  -lopenblas    -I  $(DIC)/armadillo-9.400.3/include    -lgslcblas   -llapack  -lgfortran   -lblas  -lz   -fopenmp   -lquadmath




SOURCES=conditional_function.cpp  MeQTLPolyG.cpp  caviar_PostCal.cpp PostCal.cpp  Util.cpp TopKSNP.cpp   gemma.cpp  gemma_gzstream.cpp gemma_io.cpp gemma_lapack.cpp gemma_lmm.cpp gemma_lm.cpp  gemma_mathfunc.cpp gemma_param.cpp  gemma_debug.cpp  gemma_eigenlib.cpp  gemma_fastblas.cpp prepare_input.cpp

EXECUT=MMQTL23

$(EXECUT): $(SOURCES)
        $(CC) $(SOURCES)   $(LDFLAGS) -o $@
clean:
        rm MMQTL23
