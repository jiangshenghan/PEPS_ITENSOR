
#ADD_EXECUTABLE(peps_test.exe peps_test.cc)
ADD_EXECUTABLE(utilities_test.exe utilities_test.cc)
#ADD_EXECUTABLE(simple_update_test.exe simple_update_test.cc)
#ADD_EXECUTABLE(full_update_test.exe full_update_test.cc)

#TARGET_LINK_LIBRARIES(peps_test.exe peps_ansatz peps lattice peps_utilities ${ITENSOR_LIB} ${NUMERIC_LIB})
TARGET_LINK_LIBRARIES(utilities_test.exe peps_itebd peps_ansatz peps lattice peps_utilities ${ITENSOR_LIB} ${NUMERIC_LIB} ${MPI_LIBRARIES})
#TARGET_LINK_LIBRARIES(simple_update_test.exe simple_update peps_ansatz peps lattice peps_utilities ${ITENSOR_LIB} ${NUMERIC_LIB})
#TARGET_LINK_LIBRARIES(full_update_test.exe full_update simple_update peps_ansatz peps lattice peps_utilities tnetwork ${ITENSOR_LIB} ${NUMERIC_LIB} ${MPI_LIBRARIES})

