
ADD_EXECUTABLE(featureless_honeycomb_peps.exe featureless_honeycomb_peps.cc)
ADD_EXECUTABLE(honeycomb_correlators.exe honeycomb_correlators.cc)

TARGET_LINK_LIBRARIES(featureless_honeycomb_peps.exe peps_ansatz peps lattice peps_utilities wignerSymbols itensor matrix utilities -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -liomp5 -lpthread)
TARGET_LINK_LIBRARIES(honeycomb_correlators.exe peps_ansatz peps lattice peps_utilities wignerSymbols itensor matrix utilities -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -liomp5 -lpthread)
