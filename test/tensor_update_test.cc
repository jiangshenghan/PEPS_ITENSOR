
#include "simple_update.h"
#include "tensor_update.h"

const dP imps_normalization_factor=1.E2;

int main()
{
    //control PSG parameters
    kagome_psg::mu_12=1; 
    kagome_psg::mu_c6=1;

    Print(kagome_psg::mu_12);
    Print(kagome_psg::mu_c6);


    //init params
    int Lx=8, Ly=8;
    //Kagome_Cirac_Lattice_Torus kagome_cirac_lattice({Lx,Ly});
    Kagome_Normal_Lattice_Torus kagome_normal_lattice({Lx,Ly});
    int D=6;

    //construct random peps
    // /*
    IQPEPS_IndexSet_SpinHalf index_set(D,kagome_normal_lattice);
    IQPEPS kagome_rvb(kagome_normal_lattice,index_set);
    //random_init_kagome_rvb_normal_peps(kagome_rvb);

    //construct good init state
    double init_sl_energy=0;
    do
    {
        random_init_kagome_rvb_normal_peps(kagome_rvb);
        std::array<IQIndex,2> siteuw_legs{kagome_rvb.phys_legs(0),kagome_rvb.phys_legs(2)};
        NN_Heisenberg_Hamiltonian hamiltonian_gate(siteuw_legs);
        std::array<std::vector<IQTensor>,2> env_tens;

        int comm_bond=kagome_rvb.lattice().comm_bond(0,2);
        IQTensor comm_bond_tensor=kagome_rvb.bond_tensors(comm_bond);

        get_env_tensor_minimization(kagome_rvb.site_tensors(0),kagome_rvb.site_tensors(2)*comm_bond_tensor,env_tens);

        std::array<IQTensor,2> site_env_tens{kagome_rvb.site_tensors(0),kagome_rvb.site_tensors(2)};
        for (int sitei=0; sitei<2; sitei++)
        {
            for (const auto &env_leg_tensor : env_tens[sitei])
            {
                site_env_tens[sitei]*=env_leg_tensor;
            }
            site_env_tens[sitei].noprime();
        }

        init_sl_energy=heisenberg_energy_from_site_env_tensors(site_env_tens,comm_bond_tensor,hamiltonian_gate);
        Print(init_sl_energy);
    }
    while (init_sl_energy>-0.15);
    // */

    //construct PEPS from file
    /*
    std::stringstream ss;
    //ss << "/home/jiangsb/tn_ying/tensor_network/result/peps_storage/peps_optimized/kagome_rvb_normal_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_step=1e-3_optimized";
    ss << "/home/jiangsb/tn_ying/tensor_network/result/peps_storage/kagome_reduced_update_optimized_peps/kagome_rvb_normal_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_fast_full_chi=100_step=0.1";

    std::string file_name=ss.str();
    IQPEPS kagome_rvb(kagome_normal_lattice);
    readFromFile(file_name,kagome_rvb);
    // */


    //tensor update
    kagome_normal_rvb_tensor_update(kagome_rvb);

    return 0;
}

