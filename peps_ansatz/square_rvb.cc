
#include "square_rvb.h"

using namespace square_psg;

void random_init_square_rvb_peps(IQPEPS &square_rvb)
{
    random_init_square_rvb_site_tensors(square_rvb);
    init_square_rvb_bond_tensors(square_rvb);
}

void random_init_square_rvb_site_tensors(IQPEPS &square_rvb)
{
    //generate site tensors by random parameters
    //random number generator
    std::default_random_engine generator(std::time(0));
    std::uniform_real_distribution<double> distribution(-1.0,1.0);
    auto rand_param=std::bind(distribution,generator);

    //Init a symmetric site, then generating all other sites
    Singlet_Tensor_Basis site_tensor_basis(square_rvb.site_tensors(0).indices());
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
        site_tensor_params[basei]=rand_param();

        //cout << "Spins: " << spin_list << endl
        //     << "Flavors: " << flavor_list << endl
        //     << "Fusion Channel: " << fusion_channel << endl
        //     << "Params: " << site_tensor_params[basei] << endl;
    }

    auto site_tensor=singlet_tensor_from_basis_params(site_tensor_basis,site_tensor_params);
    rotation_symmetrize_square_rvb_site_tensor(site_tensor);

    PrintDat(site_tensor);
    //auto sym_site_tensor=site_tensor;
    //spin_symmetrize_square_rvb_site_tensor(sym_site_tensor,site_tensor_basis);
    //rotation_symmetrize_square_rvb_site_tensor(sym_site_tensor);
    //sym_site_tensor*=0.25;
    //Print((sym_site_tensor-site_tensor).norm());

    square_rvb.generate_site_tensors({site_tensor});
}

void init_square_rvb_bond_tensors(IQPEPS &square_rvb)
{
    Singlet_Tensor_Basis bond_tensor_basis(square_rvb.bond_tensors(0).indices());
    std::vector<Complex> bond_tensor_params(bond_tensor_basis.dim(),0);

    for (int basei=0; basei<bond_tensor_basis.dim(); basei++)
    {
        auto spin_list=bond_tensor_basis.spin_configs(basei);
        auto flavor_list=bond_tensor_basis.flavor_configs(basei);
        //we choose gauge such that bond tensor in flavor space is either diag or ~\ii\sigma^y\otimes\mathrm{I}_{n/2}
        int flavor_dim=bond_tensor_basis.flavor_deg(0)[spin_list[0]];
        if (flavor_list[0]!=flavor_list[1] && (flavor_dim%2==1 || std::abs(flavor_list[0]-flavor_list[1])!=flavor_dim/2)) continue;

        //singlet form by two integer spins
        //chi_c4=1: real sym
        //chi_c4=-1: real antisym
        if (spin_list[0]%2==0)
        {
            if (std::abs(chi_c4-1)<EPSILON)
            {
                if (flavor_list[0]!=flavor_list[1]) continue;
                //TODO: consider the case with -1 on diag
                bond_tensor_params[basei]=1.;
            }
            if (std::abs(chi_c4+1)<EPSILON)
            {
                assert(flavor_dim%2==0);
                if (flavor_list[0]-flavor_list[1]==flavor_dim/2) bond_tensor_params[basei]=-1.;
                if (flavor_list[1]-flavor_list[0]==flavor_dim/2) bond_tensor_params[basei]=1.;
            }
            continue;
        }

        //singlet form by two half integer spins
        //mu_t1T=1: real
        //mu_t1T=-1: imag
        //mu_t2c4.chi_c4=-1: sym
        //mu_t2c4.chi_c4=+1: antisym
        if (spin_list[0]%2==1)
        {
            if (std::abs(mu_t2c4*chi_c4+1)<EPSILON)
            {
                if (flavor_list[0]!=flavor_list[1]) continue;
                //TODO: consider the case with -1 on diag
                bond_tensor_params[basei]=1.;
            }
            if (std::abs(mu_t2c4*chi_c4-1)<EPSILON)
            {
                assert(flavor_dim%2==0);
                if (flavor_list[0]-flavor_list[1]==flavor_dim/2) bond_tensor_params[basei]=-1.;
                if (flavor_list[1]-flavor_list[0]==flavor_dim/2)
                    bond_tensor_params[basei]=1.;
            }
            if (std::abs(mu_t1T+1)<EPSILON)
            {
                bond_tensor_params[basei]*=Complex_i;
            }
        }
    }//init params bond tensor

    auto bond_tensor=singlet_tensor_from_basis_params(bond_tensor_basis,bond_tensor_params);
    bond_tensor.clean();
    square_rvb.generate_bond_tensors({bond_tensor,chi_c4*bond_tensor},mu_12);
}


void rotation_symmetrize_square_rvb_site_tensor(IQTensor &site_tensor)
{
    IQTensor site_tensor_sym(site_tensor);
    //set all elems of site_tensor_sym to be zeor
    site_tensor_sym-=site_tensor;

    std::vector<int> max_val_list;
    for (const auto &leg : site_tensor.indices()) max_val_list.push_back(leg.m());

    for (int val_num=0; val_num<site_tensor.indices().dim(); val_num++)
    {
        auto val_list=list_from_num(val_num,max_val_list);
        auto val_list_rotated=val_list;
        Complex tensor_elem=0;
        for (int i=0; i<site_tensor.r()-1; i++)
        {
            tensor_elem+=std::pow(chi_c4*Theta_c4,i*1.)*site_tensor(site_tensor.indices()[0](val_list_rotated[0]+1),site_tensor.indices()[1](val_list_rotated[1]+1),site_tensor.indices()[2](val_list_rotated[2]+1),site_tensor.indices()[3](val_list_rotated[3]+1),site_tensor.indices()[4](val_list_rotated[4]+1));
            std::rotate(val_list_rotated.begin()+1,val_list_rotated.end()-1,val_list_rotated.end());
        }
        if (std::abs(tensor_elem)<EPSILON) continue;
        site_tensor_sym.set(site_tensor.indices()[0](val_list[0]+1),site_tensor.indices()[1](val_list[1]+1),site_tensor.indices()[2](val_list[2]+1),site_tensor.indices()[3](val_list[3]+1),site_tensor.indices()[4](val_list_rotated[4]+1),tensor_elem);
    }
    
    site_tensor=site_tensor_sym;
    site_tensor.clean();
}

void spin_symmetrize_square_rvb_site_tensor(IQTensor &site_tensor, const Singlet_Tensor_Basis &site_tensor_basis)
{
    std::vector<Complex> site_tensor_params;
    
    for (const auto &tensor_base : site_tensor_basis)
    {
        site_tensor_params.push_back((dag(tensor_base)*site_tensor).toComplex());
    }

    site_tensor=singlet_tensor_from_basis_params(site_tensor_basis,site_tensor_params);
    site_tensor.clean();
}