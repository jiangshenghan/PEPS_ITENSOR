
#include "trotter_gate.h"
#include "simple_update_patch_general.h"

int main()
{
    square_psg::mu_12=1;

    int Lx=8, Ly=8;
    Square_Lattice_Torus square_torus({8,8});

    //short range rvb state
    //IQPEPS square_peps=square_srvb_peps(Lx,Ly);

    //read state from file
    IQPEPS square_peps(square_torus);
    std::stringstream ss;
    //zero-flux state
    ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_rvb_D=6_Lx=8_Ly=8_iter=5_step=9999"; 
    //pi-flux state
    //ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/optimized_peps/square_pi_rvb_D=6_Lx=8_Ly=8_optimized_step_1e-0"; 
    std::string file_name=ss.str();
    readFromFile(file_name,square_peps);
    Print(square_peps.name());
    //PrintDat(square_peps.bond_tensors(0));
    //PrintDat(square_peps.bond_tensors(2));


    std::array<std::vector<IQTensor>,2> env_tens;
    get_env_tensor_minimization(square_peps.site_tensors(0)*square_peps.bond_tensors(1),square_peps.site_tensors(1),env_tens);

    //Square_Patch_RDM square_RDM(square_peps,env_tens[0][0],{{0,1},{Lx,Lx+1}},{0,1});
    Square_Patch_RDM square_RDM(square_peps,env_tens[0][0],{{0,1},{Lx,Lx+1},{2*Lx,2*Lx+1}},{Lx,Lx+1});
    //Square_Patch_RDM square_RDM(square_peps,env_tens[0][0],{{0,1,2,3},{Lx,Lx+1,Lx+2,Lx+3},{2*Lx,2*Lx+1,2*Lx+2,2*Lx+3}},{Lx+1,Lx+2});
    Print(heisenberg_energy_from_RDM(square_RDM.two_sites_RDM()));

    //General_Patch_RDM<IQTensor> square_RDM_II("regular shape",square_peps, env_tens[0][0],{0,1},{0,1});
    //General_Patch_RDM<IQTensor> square_RDM_II("regular shape",square_peps, env_tens[0][0],{0,1,Lx,Lx+1},{0,1});
    General_Patch_RDM<IQTensor> square_RDM_II("regular shape",square_peps,env_tens[0][0],{0,1,Lx,Lx+1,2*Lx,2*Lx+1},{Lx,Lx+1});
    //General_Patch_RDM<IQTensor> square_RDM_II("regular shape", square_peps,env_tens[0][0],{0,1,2,3,Lx,Lx+1,Lx+2,Lx+3,2*Lx,2*Lx+1,2*Lx+2,2*Lx+3},{Lx+1,Lx+2});
    //General_Patch_RDM<IQTensor> square_RDM_II("special shape I", square_peps,env_tens[0][0],{1,2,Lx+1,Lx+2,2*Lx+1,2*Lx+2,Lx,Lx+3},{Lx+1,Lx+2});
    Print(heisenberg_energy_from_RDM(square_RDM_II.RDM()));

}
