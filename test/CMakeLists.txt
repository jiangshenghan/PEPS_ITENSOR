
#SET(ITENSOR_LIBFLAGS "itensor matrix utilities lapack blas")

ADD_EXECUTABLE(lattice_test.exe lattice_test.cc)
ADD_EXECUTABLE(peps_test.exe peps_test.cc)

TARGET_LINK_LIBRARIES(lattice_test.exe lattice)
TARGET_LINK_LIBRARIES(peps_test.exe lattice peps itensor matrix utilities lapack blas)
