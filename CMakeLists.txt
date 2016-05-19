
#PROJECT(peps_itensor)
#CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#SET(CMAKE_VERBOSE_MAKEFILE ON)
#SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_BUILD_TYPE Release)

#set(CMAKE_C_COMPILER "icc")
#set(CMAKE_CXX_COMPILER "icpc")

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -mkl")

FIND_PACKAGE(MPI REQUIRED)

#SET(ITENSOR_INCLUDE_DIR "/home/jiangsb/tn_ying/tensor_network/ITensor/include")
#SET(ITENSOR_LIB_DIR "/home/jiangsb/tn_ying/tensor_network/ITensor/utilities")
#SET(ITPP_INCLUDE_DIR "/home/jiangsb/opt/itpp_intel/include")
#SET(ITPP_LIB_DIR "/home/jiangsb/opt/itpp_intel/lib")

#include directories storing head files
INCLUDE_DIRECTORIES(peps lattice peps_utilities peps_ansatz simple_update full_update tensor_update tensor_rg tensor_vmc ${ITPP_INCLUDE_DIR} ${ITENSOR_INCLUDE_DIR} ${BLAS_LAPACK_INCLUDE_DIR} ${MPI_INCLUDE_PATH})

#add directories for source file generate libraries
ADD_SUBDIRECTORY(lattice)
ADD_SUBDIRECTORY(peps)
ADD_SUBDIRECTORY(peps_utilities)
ADD_SUBDIRECTORY(peps_ansatz)
ADD_SUBDIRECTORY(simple_update)
ADD_SUBDIRECTORY(full_update)
ADD_SUBDIRECTORY(tensor_update)
ADD_SUBDIRECTORY(tensor_rg)
ADD_SUBDIRECTORY(tensor_vmc)
#ADD_SUBDIRECTORY(corner_transfer_matrix)
#ADD_SUBDIRECTORY(boundary_theory)

#include directories for library
LINK_DIRECTORIES(peps lattice peps_utilities peps_ansatz full_update simple_update tensor_update tensor_rg tensor_vmc ${ITPP_LIB_DIR} ${BLAS_LAPACK_LIB_DIR} )

#SET(ITENSOR_LIB utility itensor matrix)
#SET(NUMERIC_LIB itpp_static gsl)
#SET(NUMERIC_LIB itpp_static gsl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -liomp5 -lpthread)
#SET(NUMERIC_LIB itpp_static openblas lapack gsl)

#add directories for source file generating executable file
ADD_SUBDIRECTORY(test)
ADD_SUBDIRECTORY(execute)

#ADD_EXECUTABLE(dmrg dmrg.cc)
#TARGET_LINK_LIBRARIES(dmrg ${ITENSOR_LIB} ${NUMERIC_LIB})

