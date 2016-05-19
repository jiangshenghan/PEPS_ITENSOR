
#include "kagome_rvb.h"

namespace kagome_psg
{
    double mu_12=1;
    double mu_c6=1;
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

    //auto site_tensor_rotation=site_tensor;
    //rotation_symmetrize_kagome_rvb_cirac_site_tensor(site_tensor_rotation);
    //Print((site_tensor-site_tensor_rotation).norm()/site_tensor.norm());
    //auto bond_tensor_rotation=bond_tensor;
    //rotation_symmetrize_kagome_rvb_cirac_bond_tensor(bond_tensor_rotation);
    //Print((bond_tensor-bond_tensor_rotation).norm()/bond_tensor.norm());

    kagome_srvb.generate_site_tensors({site_tensor,site_tensor,tensor_permutation({0,2,1},site_tensor)});
    kagome_srvb.generate_bond_tensors({bond_tensor,tensor_permutation({2,0,1},bond_tensor)},mu_12);

    return kagome_srvb;
}

IQPEPS kagome_normal_srvb_peps(int Lx, int Ly)
{
    Kagome_Normal_Lattice_Torus kagome_torus({Lx,Ly});
    IQPEPS_IndexSet_SpinHalf indexset(3,kagome_torus);
    IQPEPS kagome_srvb(kagome_torus,indexset);

    //input kagome_srvb ansatz
    Singlet_Tensor_Basis site_tensor_basis(kagome_srvb.site_tensors(0).indices()),
                         bond_tensor_basis(kagome_srvb.bond_tensors(0).indices());
    //Print(site_tensor_basis);
    std::vector<Complex> site_tensor_params(site_tensor_basis.dim()),
                         bond_tensor_params={1,Cplx_i*sqrt(2)};
    site_tensor_params[0]=sqrt(2)*std::pow(Complex(mu_sigma),-0.5);
    site_tensor_params[1]=sqrt(2)*std::pow(Complex(mu_12*mu_c6),-0.5)*std::pow(Complex(mu_sigma),-0.5)*mu_12*mu_c6;
    site_tensor_params[2]=sqrt(2)*std::pow(Complex(mu_12*mu_c6),-0.5);
    site_tensor_params[5]=sqrt(2);
    IQTensor site_tensor=singlet_tensor_from_basis_params(site_tensor_basis,site_tensor_params),
             bond_tensor=singlet_tensor_from_basis_params(bond_tensor_basis,bond_tensor_params);

    kagome_srvb.generate_site_tensors({site_tensor,site_tensor,site_tensor});
    kagome_srvb.generate_bond_tensors({bond_tensor,tensor_after_eta_action(mu_c6,tensor_permutation({0,1},bond_tensor),bond_tensor.indices()[0]),bond_tensor,tensor_permutation({0,1},bond_tensor),bond_tensor,bond_tensor},mu_12);
    //kagome_srvb.generate_bond_tensors({bond_tensor,chi_sigmac6*tensor_after_eta_action(mu_c6*mu_sigma,bond_tensor,bond_tensor.indices()[0]),bond_tensor,chi_sigmac6*tensor_after_eta_action(mu_sigma,bond_tensor,bond_tensor.indices()[0]),bond_tensor,bond_tensor},mu_12);

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
        //if (spin_list[1]%2==1) site_tensor_params[basei]=0.5;

        //Print(basei);
        //Print(spin_list);
        //Print(site_tensor_basis[basei]);
        //Print(site_tensor_params[basei]);
    }

    auto site_tensor=singlet_tensor_from_basis_params(site_tensor_basis,site_tensor_params);

    //Print(site_tensor_params);
    //PrintDat(site_tensor);

    rotation_symmetrize_kagome_rvb_cirac_site_tensor(site_tensor);

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
    //1. three spins are all integers. 
    //2. one spin is integer and the other two are half-int. 
    //We consider cases where S_a to be the integer spin, and we can generate the whole tensor by rotation symmetry
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
    rotation_symmetrize_kagome_rvb_cirac_bond_tensor(bond_tensor);
    fix_ratio_kagome_rvb_bond_tensor(bond_tensor,bond_tensor_basis,bond_param_norms);

    kagome_rvb.generate_bond_tensors({bond_tensor,tensor_permutation({2,0,1},bond_tensor)},mu_12);
}

void rotation_symmetrize_kagome_rvb_cirac_site_tensor(IQTensor &site_tensor)
{
    //we set (T_{sym}^u)^i_{\alpha\beta}=1/2*[(T^u)^i_{\alpha\beta}+(\mu_12\mu_c6)^{1/2}(\eta_12\eta_c6)_{\beta\beta'}(T^u)^i_{\beta'\alpha}]
    auto site_tensor_permute=tensor_permutation({0,2,1},site_tensor);
    site_tensor=0.5*(site_tensor+std::pow((Complex)(mu_12*mu_c6),0.5)*tensor_after_eta_action(mu_12*mu_c6,site_tensor_permute,site_tensor.indices()[2]));
    site_tensor.clean();
}

void rotation_symmetrize_kagome_rvb_cirac_bond_tensor(IQTensor &bond_tensor)
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
    std::vector<double> bond_params;
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

    //Print(bond_tensor_basis);
    //Print(classified_base_no[0]);
    //Print(classified_base_no[1]);
    Print(bond_params);

}

void random_init_kagome_rvb_normal_peps(IQPEPS &kagome_rvb)
{
    random_init_kagome_rvb_normal_site_tensors(kagome_rvb);
    init_kagome_rvb_normal_bond_tensors(kagome_rvb);
}

void random_init_kagome_rvb_normal_site_tensors(IQPEPS &kagome_rvb)
{
    //Init symmetric site tensor in one site, then generating all other sites
    Singlet_Tensor_Basis site_tensor_basis(kagome_rvb.site_tensors(0).indices());
    //PrintDat(site_tensor_basis);
    assert(site_tensor_basis.indices()[0].type()==Site);
    //the free params for basis satisfying
    //1. S_a half int, others int
    //2. S_a int, others half int
    //knowing above, the symmetric tensor can be generated by rotation
    std::vector<Complex> site_tensor_params(site_tensor_basis.dim());

    for (int basei=0; basei<site_tensor_basis.dim(); basei++)
    {
        auto spin_list=site_tensor_basis.spin_configs(basei);
        //auto flavor_list=site_tensor_basis.flavor_configs(basei);
        //int fusion_channel=site_tensor_basis.fusion_channel(basei);
        //we only generate params with spin_oddness of leg_a differ from other three virt legs.
        std::vector<int> spin_oddness;
        for (auto S : spin_list) spin_oddness.push_back(S%2);
        if (spin_oddness[1]==spin_oddness[2] || spin_oddness[1]==spin_oddness[3] || spin_oddness[1]==spin_oddness[4]) continue;
        site_tensor_params[basei]=5*rand_gen();

        //cout << "Spins: " << spin_list << endl
        //     << "Flavors: " << flavor_list << endl
        //     << "Fusion Channel: " << fusion_channel << endl
        //     << "Params: " << site_tensor_params[basei] << endl;
    }
    //Print(site_tensor_params);

    auto site_tensor=singlet_tensor_from_basis_params(site_tensor_basis,site_tensor_params);
    rotation_reflection_symmetrize_kagome_rvb_normal_site_tensor(site_tensor);
    kagome_rvb.generate_site_tensors({site_tensor,site_tensor,site_tensor});
}

void init_kagome_rvb_normal_bond_tensors(IQPEPS &kagome_rvb)
{
    Singlet_Tensor_Basis bond_tensor_basis(kagome_rvb.bond_tensors(0).indices());
    //PrintDat(bond_tensor_basis);
    std::vector<Complex> bond_tensor_params(bond_tensor_basis.dim(),0);

    for (int basei=0; basei<bond_tensor_basis.dim(); basei++)
    {
        auto spin_list=bond_tensor_basis.spin_configs(basei);
        auto flavor_list=bond_tensor_basis.flavor_configs(basei);
        //we choose gauge such that bond tensor in flavor space is either diag or ~\ii\sigma^y\otimes\mathrm{I}_{n/2}
        int flavor_dim=bond_tensor_basis.flavor_deg(0)[spin_list[0]];
        if (flavor_list[0]!=flavor_list[1] && (flavor_dim%2==1 || std::abs(flavor_list[0]-flavor_list[1])!=flavor_dim/2)) continue;

        //singlet form by two integer spins
        //chi_sigmac6=1: real sym
        //chi_sigmac6=-1: real antisym
        if (spin_list[0]%2==0)
        {
            if (std::abs(chi_sigmac6-1)<EPSILON)
            {
                if (flavor_list[0]!=flavor_list[1]) continue;
                //TODO: consider the case with -1 on diag
                bond_tensor_params[basei]=1.*sqrt(spin_list[0]+1.);

                //set spin 1 part to be minus
                //if (spin_list[0]==2) bond_tensor_params[basei]*=-1.;
            }
            if (std::abs(chi_sigmac6+1)<EPSILON)
            {
                assert(flavor_dim%2==0);
                if (flavor_list[0]-flavor_list[1]==flavor_dim/2) bond_tensor_params[basei]=-1.*sqrt(spin_list[0]+1.);
                if (flavor_list[1]-flavor_list[0]==flavor_dim/2) bond_tensor_params[basei]=1.*sqrt(spin_list[0]+1);
            }
            continue;
        }

        //singlet form by two half integer spins
        //mu_sigma=1: real
        //mu_sigma=-1: imag
        //mu_sigma.chi_sigmac6=-1: sym
        //mu_sigma.chi_sigmac6=+1: antisym
        if (spin_list[0]%2==1)
        {
            if (std::abs(mu_sigma*chi_sigmac6+1)<EPSILON)
            {
                if (flavor_list[0]!=flavor_list[1]) continue;
                //TODO: consider the case with -1 on diag
                bond_tensor_params[basei]=1.*sqrt(spin_list[0]+1.);
            }
            if (std::abs(mu_sigma*chi_sigmac6-1)<EPSILON)
            {
                assert(flavor_dim%2==0);
                if (flavor_list[0]-flavor_list[1]==flavor_dim/2) bond_tensor_params[basei]=-1.*sqrt(spin_list[0]+1.);
                if (flavor_list[1]-flavor_list[0]==flavor_dim/2)
                    bond_tensor_params[basei]=1.*sqrt(spin_list[0]+1.);
            }
            if (std::abs(mu_sigma+1)<EPSILON)
            {
                bond_tensor_params[basei]*=Complex_i;
            }
        }
    }

    auto bond_tensor=singlet_tensor_from_basis_params(bond_tensor_basis,bond_tensor_params);
    bond_tensor.clean();
    kagome_rvb.generate_bond_tensors({bond_tensor,tensor_after_eta_action(mu_c6,tensor_permutation({0,1},bond_tensor),bond_tensor.indices()[0]),bond_tensor,tensor_permutation({0,1},bond_tensor),bond_tensor,bond_tensor},mu_12);
}

void rotation_reflection_symmetrize_kagome_rvb_normal_site_tensor(IQTensor &site_tensor)
{
    //we set (T_sym)^i_{\alpha\beta\gamma\delta}=1/4*[(T^u)^i_{\alpha\beta\gamma\delta}+(mu_12*mu_c6)^{1/2}*(eta_12*eta_c6)(beta)*(\eta_12*eta_c6)(\delta)*(T^u)^i_{\beta\alpha\delta\gamma}+(mu_12*mu_c6*mu_sigma)^{1/2}*(eta_12*eta_c6)(\alpha)*(eta_12*eta_c6)(\beta)*eta_sigma(\gamma)*eta_sigma(\delta)*(T^u)^i_{\gamma\delta\alpha\beta}+(mu_sigma)^{1/2}*(eta_12*eta_c6)(\beta)*(eta_12*eta_c6*eta_sigma)*(\gamma)*(eta_sigma)(\delta)*(T^u)^i_{\delta\gamma\beta\alpha}]
    std::vector<IQTensor> site_tensor_permute={site_tensor,std::pow((Complex)(mu_12*mu_c6),0.5)*tensor_permutation({0,2,1,4,3},site_tensor),std::pow((Complex)(mu_12*mu_c6),0.5)*std::pow((Complex)mu_sigma,0.5)*tensor_permutation({0,3,4,1,2},site_tensor),std::pow((Complex)mu_sigma,0.5)*tensor_permutation({0,4,3,2,1},site_tensor)};

    obtain_tensor_after_eta_action(mu_12*mu_c6,site_tensor_permute[1],site_tensor.indices()[2]);
    obtain_tensor_after_eta_action(mu_12*mu_c6,site_tensor_permute[1],site_tensor.indices()[4]);

    obtain_tensor_after_eta_action(mu_12*mu_c6,site_tensor_permute[2],site_tensor.indices()[1]);
    obtain_tensor_after_eta_action(mu_12*mu_c6,site_tensor_permute[2],site_tensor.indices()[2]);
    obtain_tensor_after_eta_action(mu_sigma,site_tensor_permute[2],site_tensor.indices()[3]);
    obtain_tensor_after_eta_action(mu_sigma,site_tensor_permute[2],site_tensor.indices()[4]);

    obtain_tensor_after_eta_action(mu_12*mu_c6,site_tensor_permute[3],site_tensor.indices()[2]);
    obtain_tensor_after_eta_action(mu_12*mu_c6*mu_sigma,site_tensor_permute[3],site_tensor.indices()[3]);
    obtain_tensor_after_eta_action(mu_sigma,site_tensor_permute[3],site_tensor.indices()[4]);

    site_tensor=0.25*(site_tensor_permute[0]+site_tensor_permute[1]+site_tensor_permute[2]+site_tensor_permute[3]);
}


void obtain_kagome_rvb_normal_site_tensor_symmetric_basis(const Singlet_Tensor_Basis &site_singlet_basis, std::vector<IQTensor> &symmetric_singlet_basis)
{
    symmetric_singlet_basis.clear();
    for (int basei=0; basei<site_singlet_basis.dim(); basei++)
    {
        IQTensor base_tensor=site_singlet_basis[basei];
        rotation_reflection_symmetrize_kagome_rvb_normal_site_tensor(base_tensor);
        for (const auto &pre_base : symmetric_singlet_basis)
        {
            base_tensor-=(dag(pre_base)*base_tensor).toComplex()*pre_base;
        }
        if (base_tensor.norm()<1e-10) continue;

        base_tensor/=base_tensor.norm();
        symmetric_singlet_basis.push_back(base_tensor);

        //Print(basei);
        //Print(site_singlet_basis.spin_configs(basei));
    }
    Print(symmetric_singlet_basis.size());
}
