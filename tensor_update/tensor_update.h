
//This file provide general class for tensor update

#ifndef _TENSOR_UPDATE_H_
#define _TENSOR_UPDATE_H_

#include "full_update.h"


void kagome_normal_rvb_tensor_update(IQPEPS &kagome_rvb);

std::vector<double> heisenberg_energy_derivative(PEPSt_RDM<IQTensor> &kagome_normal_rdm, const std::vector<std::vector<IQTensor>> &symmetric_site_basis, const std::vector<double> &site_params);

#endif
