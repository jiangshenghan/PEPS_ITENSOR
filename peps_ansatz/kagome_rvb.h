
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
IQPEPS kagome_srvb_cirac_peps(int Lx, int Ly);

//void random_init_kagome_rvb_cirac_peps(IQPEPS &kagome_rvb);
//
//void random_init_kagome_rvb_cirac_site_tensors(IQPEPS &kagome_rvb);
//void random_init_kagome_rvb_cirac_bond_tensors(IQPEPS &kagome_rvb);
//
//void rotation_symmetrize_kagome_rvb_site_tensors_uc(std::vector<IQTensor> &site_tensors_uc);
//void rotation_symmetrize_kagome_rvb_bond_tensors_uc(std::vector<IQTensor> &bond_tensors_uc);

#endif
