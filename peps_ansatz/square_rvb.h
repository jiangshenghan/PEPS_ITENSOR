#ifndef _SQUARE_RVB_H_
#define _SQUARE_RVB_H_

#include "peps.h"

//The following parameters labels different class of symmetric wf on square lattice
namespace square_psg
{
    extern double mu_12;
    extern double mu_t2c4;
    extern double mu_t1T;
    extern double chi_c4;
    extern double chi_T;
    extern Complex Theta_c4;
}

//namespace square_psg
//{
//    static const double mu_12=-1;
//    static const double mu_t2c4=-1;
//    static const double mu_t1T=1;
//    static const double chi_c4=1;
//    static const double chi_T=1;
//    static const Complex Theta_c4=1;
//}

IQPEPS square_srvb_peps(int Lx, int Ly);

void random_init_square_rvb_peps(IQPEPS &square_rvb);

void random_init_square_rvb_site_tensors(IQPEPS &square_rvb);

void init_square_rvb_bond_tensors(IQPEPS &square_rvb);

void rotation_symmetrize_square_rvb_site_tensor(IQTensor &site_tensor);
void spin_symmetrize_tensor(IQTensor &site_tensor, const Singlet_Tensor_Basis &site_tensor_basis);

#endif
