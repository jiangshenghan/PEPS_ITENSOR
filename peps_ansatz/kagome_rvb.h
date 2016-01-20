
#ifndef _KAGOME_RVB_H_
#define _KAGOME_RVB_H_

#include "peps.h"

//The following parameters labels different class of symmetric wf on kagome lattice
namespace kagome_psg
{
    extern double mu_12;
    extern double mu_c6;
    extern double mu_sigma;
    extern double chi_sigmac6;
}

//Depends on values of mu_12 and mu_c6, there are four kinds of srvb on kagome lattice
IQPEPS kagome_cirac_srvb_peps(int Lx, int Ly);

void random_init_kagome_rvb_cirac_peps(IQPEPS &kagome_rvb, std::array<double,2> bond_param_norms={1.,1.});

void random_init_kagome_rvb_cirac_site_tensors(IQPEPS &kagome_rvb);
void random_init_kagome_rvb_cirac_bond_tensors(IQPEPS &kagome_rvb, std::array<double,2> bond_param_norms);

//TODO: generalize to site v,w?
//make T^u rotation symmetric
void rotation_symmetrize_kagome_rvb_site_tensor(IQTensor &site_tensor);
//TODO: generalize to plaquette q?
//make P_p rotation symmetric
void rotation_symmetrize_kagome_rvb_bond_tensor(IQTensor &bond_tensor);

//there are two kinds of plaquette tensors. One kind only contains integer spins, and the other kind contains two half-int spins and one int spin. However, in any configuration, the number of each kind is fixed. In other words, changing the ratio between two kinds tensors does nothing to physical wavefunction. But, it changes measurement use simple update env. We should always fix the ratio?
void fix_ratio_kagome_rvb_bond_tensor(IQTensor &bond_tensor, const Singlet_Tensor_Basis &bond_tensor_basis, std::array<double,2> bond_param_norms={1.,1.});

#endif
