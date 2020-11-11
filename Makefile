CC=/hpc/packages/minerva-common/gcc/7.3.0/bin/g++

#CC=g++

DIC=$(PWD)
CFLAGS= -c -Wall -g -O3  -I $(DIC)  








LDFLAGS= -DWITH_LAPACK -m64 -static  -I/hpc/users/zengb02/work/software/DPR/gsl-2.5/myObj/include/ -L/hpc/users/zengb02/work/software/lib/zlib-1.2./zlib-1.2.11/lib/  -I/hpc/users/zengb02/work/software/lib/zlib-1.2./zlib-1.2.11/include/  -I/hpc/users/zengb02/work/software/DPR/ATLAS/ATLAS/myobj/include/  -I/hpc/users/zengb02/work/software/Sinai_Package/MMeQTL/MMeQTL_v2/boost_1_58_0/include    -L/hpc/users/zengb02/work/software/lib/openblas-0.2./lib/  -L/hpc/users/zengb02/work/software/lib/  -L/hpc/users/zengb02/work/software/Sinai_Package/MMeQTL/MMeQTL_v2/boost_1_58_0/lib  -I  $(DIC)/armadillo-9.400.3/include    -L/hpc/users/zengb02/work/software/DPR/gsl-2.5/myObj/lib -lgsl -lgslcblas     /hpc/users/zengb02/work/software/DPR/lapack/liblapack.a -lgfortran   /hpc/users/zengb02/work/software/DPR/ATLAS/ATLAS/myobj/lib/libatlas.a /hpc/users/zengb02/work/software/DPR/blas/libblas.a  -lz   -fopenmp   -lquadmath  



 
SOURCES=conditional_function.cpp  MeQTLPolyG.cpp  caviar_PostCal.cpp PostCal.cpp  Util.cpp TopKSNP.cpp   gemma.cpp  gemma_gzstream.cpp gemma_io.cpp gemma_lapack.cpp gemma_lmm.cpp gemma_lm.cpp  gemma_mathfunc.cpp gemma_param.cpp  gemma_debug.cpp  gemma_eigenlib.cpp  gemma_fastblas.cpp prepare_input.cpp 

EXECUT=MMQTL23
	
$(EXECUT): $(SOURCES) 
	$(CC) $(SOURCES)   $(LDFLAGS) -o $@
clean:
	rm MMQTL23
