
#ADD_EXECUTABLE(lattice_test.exe lattice_test.cc)
ADD_EXECUTABLE(peps_test.exe peps_test.cc)
ADD_EXECUTABLE(utilities_test.exe utilities_test.cc)
ADD_EXECUTABLE(simple_update_test.exe simple_update_test.cc)
#ADD_EXECUTABLE(patch_test.exe patch_test.cc)
#ADD_EXECUTABLE(ctm_test.exe ctm_test.cc)

#TARGET_LINK_LIBRARIES(lattice_test.exe peps lattice ${ITENSOT_LIB} ${NUMERIC_LIB})
TARGET_LINK_LIBRARIES(peps_test.exe peps_ansatz peps peps_utilities lattice ${ITENSOR_LIB} ${NUMERIC_LIB})
TARGET_LINK_LIBRARIES(utilities_test.exe simple_update peps_ansatz peps lattice peps_utilities ${ITENSOR_LIB} ${NUMERIC_LIB})
TARGET_LINK_LIBRARIES(simple_update_test.exe peps_ansatz peps lattice peps_utilities simple_update ${ITENSOR_LIB} ${NUMERIC_LIB})
#TARGET_LINK_LIBRARIES(patch_test.exe simple_update peps_ansatz peps lattice peps_utilities ${ITENSOR_LIB} ${NUMERIC_LIB})
#TARGET_LINK_LIBRARIES(ctm_test.exe peps_ansatz peps lattice peps_utilities simple_update corner_transfer_matrix ${ITENSOR_LIB} ${NUMERIC_LIB})
