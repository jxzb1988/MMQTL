#CC=/hpc/packages/minerva-common/gcc/7.3.0/bin/g++

CC=g++

DIC=$(PWD)
CFLAGS= -c -Wall -g -O3  -I $(DIC)






gsl_dir = /hpc/users/zengb02/work/software/DPR/gsl-2.5/myObj
Zlib_dir = /hpc/users/zengb02/work/software/lib/zlib-1.2./zlib-1.2.11
atlas_dir = /hpc/users/zengb02/work/software/DPR/ATLAS/ATLAS/myobj
boost_dir = /hpc/users/zengb02/work/software/Sinai_Package/MMeQTL/MMeQTL_v2/boost_1_58_0
openblas_dir = /hpc/users/zengb02/work/software/lib/openblas-0.2.
lapack_dir =  /hpc/users/zengb02/work/software/DPR/lapack
blas_dir = /hpc/users/zengb02/work/software/DPR/blas




LDFLAGS= -DWITH_LAPACK -m64 -static  -I${gsl_dir}/include/ -L${Zlib_dir}/lib/  -I${Zlib_dir}/include/  -I${atlas_dir}/include/  -I${boost_dir}/include    -L${openblas_dir}/lib/    -L${boost_dir}/lib  -I  $(DIC)/armadillo-9.400.3/include    -L${gsl_dir}/lib -lgsl -lgslcblas     ${lapack_dir}/liblapack.a -lgfortran   ${atlas_dir}/lib/libatlas.a ${blas_dir}/libblas.a  -lz   -fopenmp   -lquadmath



SOURCES=conditional_function.cpp  MeQTLPolyG.cpp  caviar_PostCal.cpp PostCal.cpp  Util.cpp TopKSNP.cpp   gemma.cpp  gemma_gzstream.cpp gemma_io.cpp gemma_lapack.cpp gemma_lmm.cpp gemma_lm.cpp  gemma_mathfunc.cpp gemma_param.cpp  gemma_debug.cpp  gemma_eigenlib.cpp  gemma_fastblas.cpp prepare_input.cpp

EXECUT=MMQTL23

$(EXECUT): $(SOURCES)
        $(CC) $(SOURCES)   $(LDFLAGS) -o $@
clean:
        rm MMQTL23
