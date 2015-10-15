
//This file defines functions transform peps of other lattices to square peps

#ifndef _TRANSFER_TO_SQUARE_
#define _TRANSFER_TO_SQUARE_

#include "featureless_honeycomb.h"
#include "peps.h"

//Combine two site tensors and bond tensors connect them to make a single site tensor for square lattice
//the remaining two bond tensors are square bond tensors
IQPEPS spin_sym_square_peps_from_honeycomb_tensor_uc(const std::vector<IQTensor> &honeycomb_site_tensors_uc, const std::vector<IQTensor> &honeycomb_bond_tensors_uc, const Lattice_Base &square_lattice);

#endif
