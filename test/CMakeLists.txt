
#SET(ITENSOR_LIBFLAGS "itensor matrix utilities lapack blas")

ADD_EXECUTABLE(lattice_test.exe lattice_test.cc)
#ADD_EXECUTABLE(peps_test.exe peps_test.cc)
#ADD_EXECUTABLE(cgtensor_test.exe cgtensor_test.cc)
#ADD_EXECUTABLE(utilities_test.exe utilities_test.cc)
#ADD_EXECUTABLE(simple_update_test.exe simple_update_test.cc)

TARGET_LINK_LIBRARIES(lattice_test.exe lattice)
#TARGET_LINK_LIBRARIES(peps_test.exe lattice peps itensor matrix utilities lapack blas)
#TARGET_LINK_LIBRARIES(cgtensor_test.exe peps_utilities wignerSymbols itensor matrix utilities lapack blas)
#TARGET_LINK_LIBRARIES(utilities_test.exe peps_utilities wignerSymbols itensor matrix utilities lapack blas)
#TARGET_LINK_LIBRARIES(simple_update_test.exe simple_update peps lattice peps_utilities wignerSymbols itensor matrix utilities lapack blas)
