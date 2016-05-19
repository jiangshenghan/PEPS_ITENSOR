
#include "reduced_update.h"

const dP imps_normalization_factor=1.E2;

int main()
{
    //control PSG parameters
    kagome_psg::mu_12=-1; 
    kagome_psg::mu_c6=-1;

    Print(kagome_psg::mu_12);
    Print(kagome_psg::mu_c6);


    //init params
    int Lx=8, Ly=8;
    Kagome_Normal_Lattice_Torus kagome_normal_lattice({Lx,Ly});
    int D=6;

    //construct random peps
    // /*
    IQPEPS_IndexSet_SpinHalf index_set(D,kagome_normal_lattice);
    IQPEPS kagome_rvb(kagome_normal_lattice,index_set);
    random_init_kagome_rvb_normal_peps(kagome_rvb);
    // */

    //construct PEPS from file
    /*
    std::stringstream ss;
    ss << "/home/jiangsb/tn_ying/tensor_network/result/peps_storage/peps_optimized/kagome_rvb_normal_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_step=1e-3_optimized";
    //ss << "/home/jiangsb/tn_ying/tensor_network/result/peps_storage/kagome_rvb_D=6_Lx=8_Ly=8_mu12=1_muc6=1_iter=0_step=78";
    std::string file_name=ss.str();
    IQPEPS kagome_rvb(kagome_normal_lattice);
    readFromFile(file_name,kagome_rvb);
    // */


    //optimazation
    Evolution_Params su_params(4,{21,91,281,911},{1,1e-1,1e-2,1e-3});
    //Evolution_Params su_params(6,{37,97,297,697,1297,3997},{1,3e-1,1e-1,3e-2,1e-2,3e-3});
    //Evolution_Params su_params(4,{1297,3997},{1e-2,3e-3});
    //Evolution_Params su_params(1,{23},{1e-0});
    Print(su_params);

    //full update
    //kagome_normal_rvb_fast_full_update(kagome_rvb,su_params);

    //reduced update
    kagome_normal_rvb_reduced_update(kagome_rvb,su_params);

    return 0;
}
