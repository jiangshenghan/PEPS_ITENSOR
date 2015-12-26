
//#include "peps.h"
//#include "square_rvb.h"
//#include "kagome_rvb.h"
#include "simple_update.h"

using namespace square_psg;

int main()
{

    //Square_Lattice_Torus square_torus({8,8});

    //IQPEPS square_peps(square_torus);

    //std::stringstream ss;
    ////zero-flux state
    //ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_rvb_D=6_Lx=8_Ly=8_optimized"; 
    ////pi-flux state
    ////ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_pi_rvb_D=6_Lx=8_Ly=8_optimized"; 
    //std::string file_name=ss.str();
    //readFromFile(file_name,square_peps);
    //Print(square_peps.name());

    ////store as tnetwork_storage
    //Tnetwork_Storage<IQTensor> square_rvb_storage=peps_to_tnetwork_storage(square_peps);
    //std::stringstream ss_tnetwork;
    ////zero-flux state
    //ss_tnetwork << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_rvb_D=6_Lx=8_Ly=8_optimized";
    ////pi-flux state
    ////ss_tnetwork << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_pi_rvb_D=6_Lx=8_Ly=8_optimized";
    //file_name=ss_tnetwork.str();
    //writeToFile(file_name,square_rvb_storage);

    //check C4 symmetry
    //T^i_{abcd}=\chi_c4*theta_c4 T^i_{dabc}
    //IQTensor site_tensor=(square_peps.site_tensors(0)*square_peps.bond_tensors(0)*square_peps.bond_tensors(1)*square_peps.bond_tensors(15)*square_peps.bond_tensors(112));
    //IQIndex phys_leg=square_peps.phys_legs(0);
    ////PrintDat(site_tensor);
    //PrintDat(square_peps.bond_tensors(0));
    //PrintDat(square_peps.bond_tensors(1));
    //PrintDat(square_peps.bond_tensors(0)*(dag(square_peps.bond_tensors(0)).prime(dag(square_peps.virt_legs(35)))));

    //for (int indi=0; indi<1296; indi++)
    //{
    //    auto ind_list=list_from_num(indi,{6,6,6,6});
    //    double origin_elem=site_tensor(phys_leg(1),square_peps.virt_legs(30)(ind_list[0]+1),square_peps.virt_legs(35)(ind_list[1]+1),square_peps.virt_legs(4)(ind_list[2]+1),square_peps.virt_legs(225)(ind_list[3]+1)),
    //           rotated_elem=site_tensor(phys_leg(1),square_peps.virt_legs(30)(ind_list[3]+1),square_peps.virt_legs(35)(ind_list[0]+1),square_peps.virt_legs(4)(ind_list[1]+1),square_peps.virt_legs(225)(ind_list[2]+1));
    //    
    //    if (std::abs(origin_elem)>EPSILON && std::abs(rotated_elem)>EPSILON)
    //    {
    //        int eta_sign=std::abs(square_peps.virt_legs(35)(ind_list[1]+1).qn().sz()+square_peps.virt_legs(225)(ind_list[3]+1).qn().sz())%2;
    //        eta_sign=1-2*eta_sign;
    //        cout << eta_sign << endl;
    //        cout << origin_elem << " " << rotated_elem << endl;
    //        cout << origin_elem+eta_sign*rotated_elem << endl;
    //        cout << endl;
    //    }
    //}

    mu_12=1;
    IQPEPS square_peps=square_srvb_peps(8,8);

    std::array<std::vector<IQTensor>,2> env_tens;
    get_env_tensor_minimization(square_peps.site_tensors(0)*square_peps.bond_tensors(1),square_peps.site_tensors(1),env_tens);

    PrintDat(square_peps_two_sites_RDM_simple_update(square_peps,env_tens[0][0],{{0,1,2,3},{8,9,10,11},{16,17,18,19}},{9,10}));

    return 0;
}
