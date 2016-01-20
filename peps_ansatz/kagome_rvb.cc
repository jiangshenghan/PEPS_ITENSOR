
#include "kagome_rvb.h"

namespace kagome_psg
{
    double mu_12=1;
    double mu_c6=-1;
    double mu_sigma=-1;
    double chi_sigmac6=1;
}

using namespace kagome_psg;

IQPEPS kagome_cirac_srvb_peps(int Lx, int Ly)
{
    Kagome_Cirac_Lattice_Torus kagome_torus({Lx,Ly});
    IQPEPS_IndexSet_SpinHalf indexset(3,kagome_torus);
    IQPEPS kagome_srvb(kagome_torus,indexset);

    //input kagome srvb ansatz
    Singlet_Tensor_Basis site_tensor_basis(kagome_srvb.site_tensors(0).indices()),
                         bond_tensor_basis(kagome_srvb.bond_tensors(0).indices());

    std::vector<Complex> site_tensor_params={sqrt(2)*std::pow(Complex(mu_12*mu_c6),-0.5),sqrt(2)},
                         bond_tensor_params={sqrt(2),sqrt(2),-mu_12*mu_c6*sqrt(2),mu_12*sqrt(2)};

    //Print(mu_12);
    //Print(mu_c6);

    IQTensor site_tensor=singlet_tensor_from_basis_params(site_tensor_basis,site_tensor_params),
             bond_tensor=singlet_tensor_from_basis_params(bond_tensor_basis,bond_tensor_params);

    //PrintDat(site_tensor);
    //PrintDat(tensor_permutation({0,2,1},site_tensor));
    //PrintDat(bond_tensor);
    //PrintDat(tensor_permutation({2,0,1},bond_tensor));

    kagome_srvb.generate_site_tensors({site_tensor,site_tensor,tensor_permutation({0,2,1},site_tensor)});
    kagome_srvb.generate_bond_tensors({bond_tensor,tensor_permutation({2,0,1},bond_tensor)},mu_12);

    return kagome_srvb;
}



void random_init_kagome_rvb_cirac_peps(IQPEPS &kagome_rvb, std::array<double,2> bond_param_norms)
{
    random_init_kagome_rvb_cirac_site_tensors(kagome_rvb);
    random_init_kagome_rvb_cirac_bond_tensors(kagome_rvb,bond_param_norms);
}

void random_init_kagome_rvb_cirac_site_tensors(IQPEPS &kagome_rvb)
{
    //generate site tensors by random parameters
    //random number generator
    std::default_random_engine generator(std::time(0));
    std::uniform_real_distribution<double> distribution(-1.0,1.0);
    auto rand_param=std::bind(distribution,generator);

    //Init symmetric site tensor in one site, then generating all other sites
    Singlet_Tensor_Basis site_tensor_basis(kagome_rvb.site_tensors(0).indices());
    //consider case where S_a half-int and S_b int, and set the parameters in this case to be real
    std::vector<Complex> site_tensor_params(site_tensor_basis.dim());

    for (int basei=0; basei<site_tensor_basis.dim(); basei++)
    {
        const auto &spin_list=site_tensor_basis.spin_configs(basei);
        if (spin_list[1]%2==1) site_tensor_params[basei]=rand_param();

        //Print(basei);
        //Print(spin_list);
        //Print(site_tensor_basis[basei]);
        //Print(site_tensor_params[basei]);
    }

    auto site_tensor=singlet_tensor_from_basis_params(site_tensor_basis,site_tensor_params);

    //Print(site_tensor_params);
    //PrintDat(site_tensor);

    rotation_symmetrize_kagome_rvb_site_tensor(site_tensor);

    kagome_rvb.generate_site_tensors({site_tensor,site_tensor,tensor_permutation({0,2,1},site_tensor)});
}

void random_init_kagome_rvb_cirac_bond_tensors(IQPEPS &kagome_rvb, std::array<double,2> bond_param_norms)
{
    //generate site tensors by random parameters
    //random number generator
    std::default_random_engine generator(std::time(0)+1);
    std::uniform_real_distribution<double> distribution(-1.0,1.0);
    auto rand_param=std::bind(distribution,generator);

    //Init one symmetric plaquette tensor, then generate all other plaquette tensors
    Singlet_Tensor_Basis bond_tensor_basis(kagome_rvb.bond_tensors(0).indices());
    //There are two cases: 
    //1. three spins are all integers;
    //2. one spin is integer and the other two are half-int. We consider cases where S_a to be the integer spin.
    std::vector<double> bond_tensor_params(bond_tensor_basis.dim());
    for (int basei=0; basei<bond_tensor_basis.dim(); basei++)
    {
        const auto &spin_list=bond_tensor_basis.spin_configs(basei);
        if (spin_list[0]%2==0) 
        {
            bond_tensor_params[basei]=rand_param();
        }
    }

    auto bond_tensor=singlet_tensor_from_basis_params(bond_tensor_basis,bond_tensor_params);
    rotation_symmetrize_kagome_rvb_bond_tensor(bond_tensor);
    fix_ratio_kagome_rvb_bond_tensor(bond_tensor,bond_tensor_basis,bond_param_norms);

    kagome_rvb.generate_bond_tensors({bond_tensor,tensor_permutation({2,0,1},bond_tensor)},mu_12);
}

void rotation_symmetrize_kagome_rvb_site_tensor(IQTensor &site_tensor)
{
    //we set (T_{sym}^u)^i_{\alpha\beta}=1/2*[(T^u)^i_{\alpha\beta}+(\mu_12\mu_c6)^{1/2}(\eta_12\eta_c6)_{\beta\beta'}(T^u)^i_{\beta'\alpha}]
    auto site_tensor_permute=tensor_permutation({0,2,1},site_tensor);
    site_tensor=0.5*(site_tensor+std::pow((Complex)(mu_12*mu_c6),0.5)*tensor_after_eta_action(mu_12*mu_c6,site_tensor_permute,site_tensor.indices()[2]));
    site_tensor.clean();
}

void rotation_symmetrize_kagome_rvb_bond_tensor(IQTensor &bond_tensor)
{
    //we set (P^p_{sym})_{\alpha\beta\gamma}=1/3*[(P^p)_{\alpha\beta\gamma}+\eta_12(\betta)\eta_c6(\gamma)(P^p)_{\gamma\alpha\beta}+\eta_12(\gamma)\eta_c6(\alpha)(P^p)_{\beta\gamma\alpha}]
    std::vector<IQTensor> bond_tensors_permute={bond_tensor,tensor_permutation({2,0,1},bond_tensor),tensor_permutation({1,2,0},bond_tensor)};
    obtain_tensor_after_eta_action(mu_12,bond_tensors_permute[1],bond_tensor.indices()[1]);
    obtain_tensor_after_eta_action(mu_c6,bond_tensors_permute[1],bond_tensor.indices()[2]);
    obtain_tensor_after_eta_action(mu_12,bond_tensors_permute[2],bond_tensor.indices()[2]);
    obtain_tensor_after_eta_action(mu_c6,bond_tensors_permute[2],bond_tensor.indices()[0]);

    bond_tensor=1./3.*(bond_tensors_permute[0]+bond_tensors_permute[1]+bond_tensors_permute[2]);
}


void fix_ratio_kagome_rvb_bond_tensor(IQTensor &bond_tensor, const Singlet_Tensor_Basis &bond_tensor_basis, std::array<double,2> bond_param_norms)
{
    std::vector<Complex> bond_params;
    obtain_singlet_tensor_params(bond_tensor,bond_tensor_basis,bond_params);
    //classify singlet basis into two types
    std::array<std::vector<int>,2> classified_base_no;
    for (int basei=0; basei<bond_tensor_basis.dim(); basei++)
    {
        int type=0;
        for (auto spin : bond_tensor_basis.spin_configs(basei))
        {
            if (spin%2!=0)
            {
                type=1;
                break;
            }
        }
        classified_base_no[type].push_back(basei);
    }

    //obtain the original norm
    std::array<double,2> origin_bond_param_norms={0,0};
    for (int typei=0; typei<2; typei++)
    {
        for (int base_no : classified_base_no[typei]) origin_bond_param_norms[typei]+=std::abs(bond_params[base_no]*bond_params[base_no]);
        origin_bond_param_norms[typei]=sqrt(origin_bond_param_norms[typei]);
    }

    //rescale params
    for (int typei=0; typei<2; typei++)
    {
        for (int base_no: classified_base_no[typei]) bond_params[base_no]*=bond_param_norms[typei]/origin_bond_param_norms[typei];
    }
    
    bond_tensor=singlet_tensor_from_basis_params(bond_tensor_basis,bond_params);
}
