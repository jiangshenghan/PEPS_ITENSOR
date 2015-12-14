
#include "kagome_rvb.h"

namespace kagome_psg
{
    double mu_12=1;
    double mu_c6=-1;
    double mu_sigma=-1;
    double chi_sigmac6=1;
}

using namespace kagome_psg;

IQPEPS kagome_srvb_cirac_peps(int Lx, int Ly)
{
    Kagome_Cirac_Lattice_Torus kagome_torus({Lx,Ly});
    IQPEPS_IndexSet_SpinHalf indexset(3,kagome_torus);
    IQPEPS kagome_srvb(kagome_torus,indexset);

    //input kagome srvb ansatz
    Singlet_Tensor_Basis site_tensor_basis(kagome_srvb.site_tensors(0).indices()),
                         bond_tensor_basis(kagome_srvb.bond_tensors(0).indices());

    std::vector<Complex> site_tensor_params={mu_12*mu_c6>0?sqrt(2):Complex_i*sqrt(2),sqrt(2)},
                        bond_tensor_params={1,sqrt(2),-mu_12*mu_c6*sqrt(2),mu_12*sqrt(2)};

    Print(site_tensor_params[0]);

    IQTensor site_tensor=singlet_tensor_from_basis_params(site_tensor_basis,site_tensor_params),
             bond_tensor=singlet_tensor_from_basis_params(bond_tensor_basis,bond_tensor_params);
    PrintDat(site_tensor);
    PrintDat(bond_tensor);

    //TODO: finish this function

    return kagome_srvb;
}
