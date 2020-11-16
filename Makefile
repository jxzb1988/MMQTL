CC=/hpc/packages/minerva-common/gcc/7.3.0/bin/g++

#CC=g++

DIC=$(PWD)
CFLAGS= -c -Wall -g -O3  -I $(DIC)  





LDFLAGS= -DWITH_LAPACK -m64 -static  -I gsl-2.5/myObj/include/ -L lib-1.2./zlib-1.2.11/lib/  -I zlib-1.2.11/include/  -I ATLAS/myobj/include/  -I boost_1_58_0/include    -L openblas-0.2./lib/  -L boost_1_58_0/lib  -I  $(DIC)/armadillo-9.400.3/include    -L gsl-2.5/myObj/lib -lgsl -lgslcblas     lapack/liblapack.a -lgfortran   ATLAS/myobj/lib/libatlas.a blas/libblas.a  -lz   -fopenmp   -lquadmath  



 
SOURCES=conditional_function.cpp  MeQTLPolyG.cpp  caviar_PostCal.cpp PostCal.cpp  Util.cpp TopKSNP.cpp   gemma.cpp  gemma_gzstream.cpp gemma_io.cpp gemma_lapack.cpp gemma_lmm.cpp gemma_lm.cpp  gemma_mathfunc.cpp gemma_param.cpp  gemma_debug.cpp  gemma_eigenlib.cpp  gemma_fastblas.cpp prepare_input.cpp 

EXECUT=MMQTL23
	
$(EXECUT): $(SOURCES) 
	$(CC) $(SOURCES)   $(LDFLAGS) -o $@
clean:
	rm MMQTL23
