
#ADD_EXECUTABLE(featureless_honeycomb_peps.exe featureless_honeycomb_peps.cc)
#ADD_EXECUTABLE(honeycomb_correlators.exe honeycomb_correlators.cc)
#ADD_EXECUTABLE(honeycomb_TPO_correlators.exe honeycomb_TPO_correlators.cc)
#ADD_EXECUTABLE(peps_change_size.exe peps_change_size.cc)
ADD_EXECUTABLE(peps_to_tnetwork_storage.exe peps_to_tnetwork_storage.cc)
ADD_EXECUTABLE(trg_vmc_measure.exe trg_vmc_measure.cc)

#TARGET_LINK_LIBRARIES(featureless_honeycomb_peps.exe peps_ansatz peps lattice peps_utilities ${ITENSOR_LIB} ${NUMERIC_LIB})
#TARGET_LINK_LIBRARIES(honeycomb_correlators.exe peps_ansatz peps lattice peps_utilities ${ITENSOR_LIB} ${NUMERIC_LIB})
#TARGET_LINK_LIBRARIES(honeycomb_TPO_correlators.exe peps_ansatz peps lattice peps_utilities itensor matrix utilities -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -liomp5 -lpthread)
#TARGET_LINK_LIBRARIES(peps_change_size.exe peps_ansatz peps lattice peps_utilities ${ITENSOR_LIB} ${NUMERIC_LIB})
TARGET_LINK_LIBRARIES(peps_to_tnetwork_storage.exe simple_update peps_ansatz peps lattice peps_utilities ${ITENSOR_LIB} ${NUMERIC_LIB})
TARGET_LINK_LIBRARIES(trg_vmc_measure.exe tensor_vmc tensor_rg peps_ansatz peps lattice peps_utilities ${ITENSOR_LIB} ${NUMERIC_LIB} ${MPI_LIBRARIES})
