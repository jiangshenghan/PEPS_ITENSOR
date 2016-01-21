
//#include "peps.h"
//#include "square_rvb.h"
//#include "kagome_rvb.h"
#include "simple_update.h"

//using namespace square_psg;

int main()
{
    int D=6, Lx=8, Ly=8;

    //Sqaure Lattice:
    //Square_Lattice_Torus square_torus({Lx,Ly});
    //IQPEPS square_peps(square_torus);

    //std::stringstream ss;
    ////zero-flux state
    ////ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/optimized_peps/square_rvb_D=10_Lx=8_Ly=8_optimized_step_1e-4"; 
    ////pi-flux state
    //ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/optimized_peps/square_pi_rvb_D=6_Lx=8_Ly=8_patch=2x2_optimized_step_1e-4"; 
    //std::string file_name=ss.str();
    //readFromFile(file_name,square_peps);
    ////Print(square_peps.name());
    //PrintDat(square_peps.bond_tensors(0));
    //PrintDat(square_peps.bond_tensors(2));

    ////check simple update energy
    //std::array<std::vector<IQTensor>,2> env_tens;
    //auto combined_site_tens0=square_peps.site_tensors(0);
    //auto combined_site_tens1=square_peps.site_tensors(1);
    //for (int neighi=0; neighi<square_peps.n_bonds_to_one_site(); neighi++)
    //{
    //    int bondi=square_peps.lattice().site_neighbour_bonds(0,neighi);
    //    if (bondi==1) continue;
    //    combined_site_tens0*=square_peps.bond_tensors(bondi);
    //    bondi=square_peps.lattice().site_neighbour_bonds(1,neighi);
    //    if (bondi==1) continue;
    //    combined_site_tens1*=square_peps.bond_tensors(bondi);
    //}
    //get_env_tensor_minimization(combined_site_tens0*square_peps.bond_tensors(1),combined_site_tens1,env_tens);
    //std::array<IQTensor,2> site_env_tens{{combined_site_tens0,combined_site_tens1}};

    //for (int sitei=0; sitei<2; sitei++)
    //{
    //    for (const auto &env_leg_tensor : env_tens[sitei])
    //    {
    //        site_env_tens[sitei]*=env_leg_tensor;
    //    }
    //    site_env_tens[sitei].noprime();
    //}
    //Print(heisenberg_energy_from_site_env_tensors(site_env_tens,square_peps.bond_tensors(1),NN_Heisenberg_Hamiltonian({square_peps.phys_legs(0),square_peps.phys_legs(1)})));

    ////store as tnetwork_storage
    //Tnetwork_Storage<IQTensor> square_rvb_storage=peps_to_tnetwork_storage(square_peps);
    //std::stringstream ss_tnetwork;
    ////zero-flux state
    ////ss_tnetwork << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_rvb_D=10_Lx=8_Ly=8_optimized_1e-4";
    ////pi-flux state
    //ss_tnetwork << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_pi_rvb_D=6_Lx=8_Ly=8_patch=2x2_optimized_step_1e-4";
    //file_name=ss_tnetwork.str();
    //writeToFile(file_name,square_rvb_storage);

    //Kagome Cirac Lattice
    Kagome_Cirac_Lattice_Torus kagome_cirac_lattice({Lx,Ly});
    IQPEPS kagome_peps(kagome_cirac_lattice);
    kagome_psg::mu_12=1;
    kagome_psg::mu_c6=1;

    std::stringstream ss;
    //ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/optimized_peps/kagome_rvb_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_tree_shape_I_optimized_step_1e-1";
    ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/optimized_peps/kagome_rvb_D=3_Lx=" << Lx << "_Ly=" << Ly << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6;
    std::string file_name=ss.str();
    readFromFile(file_name,kagome_peps);

    Tnetwork_Storage<IQTensor> kagome_rvb_storage=peps_to_tnetwork_storage(kagome_peps);
    ss.str("");
    //ss << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/kagome_rvb_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_tree_shape_I_optimized_step_1e-1";
    ss << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/kagome_rvb_D=3_Lx=" << Lx << "_Ly=" << Ly << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6;
    file_name=ss.str();
    writeToFile(file_name,kagome_rvb_storage);
    
    return 0;
}
