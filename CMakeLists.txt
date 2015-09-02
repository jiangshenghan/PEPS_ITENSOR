
PROJECT(peps_itensor)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

SET(BLAS_LAPACK_INCLUDE_DIR "") #only use for mkl
SET(BLAS_LIB_DIR "/usr/public/blas/1.0")
SET(LAPACK_LIB_DIR "/usr/public/lapack/lapack-3.4.2/lib")
SET(ITENSOR_INCLUDE_DIR "/home/jiangsb/lib/itensor/include")
SET(ITENSOR_LIB_DIR "/home/jiangsb/lib/itensor/lib")

#include directories storing head files
INCLUDE_DIRECTORIES(peps lattice utilities ${ITENSOR_INCLUDE_DIR} ${BLAS_LAPACK_INCLUDE_DIR})

#add directories for source file
ADD_SUBDIRECTORY(lattice)
ADD_SUBDIRECTORY(peps)

#include directories for library
LINK_DIRECTORIES(peps lattice ${ITENSOR_LIB_DIR} ${LAPACK_LIB_DIR} ${BLAS_LIB_DIR})

#add directories for source file generating executable file
ADD_SUBDIRECTORY(test)
