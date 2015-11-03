
PROJECT(peps_itensor)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_BUILD_TYPE Release)

#set(CMAKE_C_COMPILER "/usr/public/intel/composer_xe_2013_sp1.3.174/bin/intel64/icc")
#set(CMAKE_CXX_COMPILER "/usr/public/intel/composer_xe_2013_sp1.3.174/bin/intel64/icpc")
#set(CMAKE_INCLUDE_PATH "/usr/public/intel/composer_xe_2013_sp1.3.174/compiler/include")
#set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/usr/public/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64")
set(CMAKE_CXX_COMPILER "icpc")

#set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
#set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#add_definitions(-std=c++11)

SET(BLAS_LAPACK_INCLUDE_DIR "/usr/public/intel/composer_xe_2013_sp1.3.174/mkl/include") #only use for mkl
SET(BLAS_LAPACK_LIB_DIR "/usr/public/intel/composer_xe_2013_sp1.3.174/mkl/lib/intel64" )
#SET(BLAS_LIB_DIR "/usr/public/blas/1.0")
#SET(LAPACK_LIB_DIR "/usr/public/lapack/lapack-3.4.2/lib")
SET(ITENSOR_INCLUDE_DIR "/home/jiangsb/lib/itensor/include")
SET(ITENSOR_LIB_DIR "/home/jiangsb/lib/itensor/lib")
SET(ITPP_INCLUDE_DIR "/home/jiangsb/opt/itpp_intel/include")
SET(ITPP_LIB_DIR "/home/jiangsb/opt/itpp_intel/lib")

#include directories storing head files
INCLUDE_DIRECTORIES(peps lattice peps_utilities peps_ansatz simple_update ${ITPP_INCLUDE_DIR} ${ITENSOR_INCLUDE_DIR} ${WIGNERSYMBOL_INCLUDE_DIR} ${BLAS_LAPACK_INCLUDE_DIR})

#add directories for source file generate libraries
ADD_SUBDIRECTORY(lattice)
ADD_SUBDIRECTORY(peps)
ADD_SUBDIRECTORY(peps_utilities)
ADD_SUBDIRECTORY(peps_ansatz)
ADD_SUBDIRECTORY(simple_update)
#ADD_SUBDIRECTORY(boundary_theory)

#include directories for library
LINK_DIRECTORIES(peps lattice peps_utilities peps_ansatz simple_update ${ITPP_LIB_DIR} ${ITENSOR_LIB_DIR} ${WIGNERSYMBOL_LIB_DIR} ${BLAS_LAPACK_LIB_DIR})

SET(ITENSOR_LIB itensor matrix utilities)
SET(NUMERIC_LIB itpp_static gsl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -liomp5 -lpthread)

#add directories for source file generating executable file
ADD_SUBDIRECTORY(test)
#ADD_SUBDIRECTORY(execute)
