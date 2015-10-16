
#ifndef _FEATURELESS_HONEYCOMB_H_
#define _FEATURELESS_HONEYCOMB_H_

#include "singlet_tensor_basis.h"

//ansatz for featureless insulator on honeycomb lattice on one uc
//A1 and A2 are tunable parameters
inline void generate_featureless_honeycomb_ansatz_uc(double A1, double A2, std::vector<IQTensor> &honeycomb_site_tensors_uc, std::vector<IQTensor> &honeycomb_bond_tensors_uc)
{
    //construct phys legs and virt legs for site tensors and bond tensors of honeycomb lattice in one uc
    std::vector<IQIndex> phys_legs={Spin_leg(std::vector<int>{0,1},"phys_leg 1",Out,Site), Spin_leg(std::vector<int>{0,1},"phys_leg 2",Out,Site)};

    std::vector<IQIndex> virt_legs;
    for (int i=0; i<8; i++)
    {
        virt_legs.push_back(Spin_leg(std::vector<int>{0,2},nameint("virt_leg ",i)));
    }

    //obtain singlet basis for two site tensors, then construct symmetric site tensors using these basis
    //in our gauge, these two site tensors are the same
    std::vector<std::vector<IQIndex>> site_tensors_indices={{phys_legs[0],virt_legs[0],virt_legs[1],virt_legs[2]},{phys_legs[1],virt_legs[3],virt_legs[4],virt_legs[5]}};
    std::vector<Singlet_Tensor_Basis> site_tensors_basis={Singlet_Tensor_Basis(site_tensors_indices[0]),Singlet_Tensor_Basis(site_tensors_indices[1])};
    //PrintDat(site_tensors_basis[0]);
    //PrintDat(site_tensors_basis[1]);

    //In current convention \widetilde{T}^i_{\alpha\beta\gamma} transform to [8*(i-1)+BinarytoDecimal(\alpha\beta\gamma)]'s elem of T_tilde
    //Here, i=1,2 labels the fusion channel, \alpha=0,1 labels deg in flavor space
    //cout << "params: A1=" << A1 << ", A2=" << A2 << endl;

    std::array<double,2> norms={2.,std::sqrt(12.)};
    //peps suppose to have Z2 IGG and independ on A1/A2
    //std::vector<double> T_tilde={0,-0.5*A2*norms[0],-0.5*A2*norms[0],A1*norms[0],A2*norms[0],-0.5*A1*norms[0],-0.5*A1*norms[0],0,0,0.5*A2*norms[1],-0.5*A2*norms[1],0,0,-0.5*A1*norms[1],0.5*A1*norms[1],0};
    //trivial IGG peps
    std::vector<double> T_tilde={0,-0.5*A2*norms[0],-0.5*A2*norms[0],0,A2*norms[0],1.5*A1*norms[0],-1.5*A1*norms[0],0,0,0.5*A2*norms[1],-0.5*A2*norms[1],A1*norms[1],0,-0.5*A1*norms[1],-0.5*A1*norms[1],0};
    
    honeycomb_site_tensors_uc.push_back(singlet_tensor_from_basis_params(site_tensors_basis[0],T_tilde));
    honeycomb_site_tensors_uc.push_back(singlet_tensor_from_basis_params(site_tensors_basis[1],T_tilde));


    //obtain singlet basis for three bond tensors, and then construct bond tensors
    std::vector<std::vector<IQIndex>> bond_tensors_indices={{dag(virt_legs[0]),dag(virt_legs[3])},{dag(virt_legs[5]),dag(virt_legs[6])},{dag(virt_legs[4]),dag(virt_legs[7])}};

    std::vector<Singlet_Tensor_Basis> bond_tensors_basis={Singlet_Tensor_Basis(bond_tensors_indices[0]),Singlet_Tensor_Basis(bond_tensors_indices[1]),Singlet_Tensor_Basis(bond_tensors_indices[2])};

    std::vector<double> B_tilde={std::sqrt(2),0,0,std::sqrt(2)};

    for (const auto &basis : bond_tensors_basis)
    {
        honeycomb_bond_tensors_uc.push_back(singlet_tensor_from_basis_params(basis,B_tilde));
    }


    //cout << "\n========================================\n" << endl;
    //cout << "parameters= " << A1 << " " << A2 << endl << endl;
    //cout << "honeycomb site tensor:" << endl;
    //for (const auto &tensor : honeycomb_site_tensors_uc)
    //{
    //    PrintDat(tensor);
    //}
    //cout << "honeycomb bond tensor:" << endl;
    //for (const auto &tensor : honeycomb_bond_tensors_uc)
    //{
    //    PrintDat(tensor);
    //}
    //cout << "\n========================================\n" << endl;

}

#endif
