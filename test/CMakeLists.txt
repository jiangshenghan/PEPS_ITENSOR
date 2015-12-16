

#ADD_EXECUTABLE(lattice_test.exe lattice_test.cc)
#ADD_EXECUTABLE(peps_test.exe peps_test.cc)
#ADD_EXECUTABLE(utilities_test.exe utilities_test.cc)
#ADD_EXECUTABLE(simple_update_test.exe simple_update_test.cc)
#ADD_EXECUTABLE(small_sample_test.exe small_sample_test.cc)
ADD_EXECUTABLE(ctm_test.exe ctm_test.cc)

#TARGET_LINK_LIBRARIES(lattice_test.exe lattice)
#TARGET_LINK_LIBRARIES(peps_test.exe lattice peps peps_utilities ${ITENSOR_LIB} ${NUMERIC_LIB})
#TARGET_LINK_LIBRARIES(utilities_test.exe peps_ansatz peps lattice peps_utilities ${ITENSOR_LIB} ${NUMERIC_LIB})
#TARGET_LINK_LIBRARIES(simple_update_test.exe peps_ansatz peps lattice peps_utilities simple_update ${ITENSOR_LIB} ${NUMERIC_LIB})
#TARGET_LINK_LIBRARIES(small_sample_test.exe lattice peps peps_utilities wignerSymbols itensor matrix utilities -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -liomp5 -lpthread)
TARGET_LINK_LIBRARIES(ctm_test.exe peps_ansatz peps lattice peps_utilities simple_update corner_transfer_matrix ${ITENSOR_LIB} ${NUMERIC_LIB})
