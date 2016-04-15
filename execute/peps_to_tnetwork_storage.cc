
//#include "peps.h"
//#include "square_rvb.h"
//#include "kagome_rvb.h"
#include "simple_update.h"

//using namespace square_psg;

int main()
{
    int D=6, Lx=8, Ly=8;

    //Kagome Lattice
    // /*
    //Kagome_Cirac_Lattice_Torus kagome_cirac_lattice({Lx,Ly});
    Kagome_Normal_Lattice_Torus kagome_normal_lattice({Lx,Ly});
    IQPEPS kagome_peps(kagome_normal_lattice);
    kagome_psg::mu_12=1;
    kagome_psg::mu_c6=1;
    int su_env_option=1;

    std::stringstream ss;
    ss << "/home/jiangsb/tn_ying/tensor_network/result/peps_storage/peps_optimized/kagome_rvb_normal_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_cluster_triangle_step=1e-3_optimized";
    std::string file_name=ss.str();
    readFromFile(file_name,kagome_peps);

    Tnetwork_Storage<IQTensor> kagome_rvb_storage=peps_to_tnetwork_storage(kagome_peps);
    ss.str("");
    ss << "/home/jiangsb/tn_ying/tensor_network/result/tnetwork_storage/kagome_rvb_normal_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_cluster_triangle_step=1e-3_optimized";
    file_name=ss.str();
    writeToFile(file_name,kagome_rvb_storage);
    // */


    //Square Lattice
    /*
    Square_Lattice_Torus square_lattice({8,8});
    IQPEPS square_peps(square_lattice);
    std::stringstream ss;
    ss << "/home/jiangsb/tn_ying/tensor_network/result/peps_storage/peps_optimized/square_rvb_D=6_Lx=8_Ly=8_mu12=1_step=1e-3_optimized";
    std::string file_name=ss.str();
    readFromFile(file_name,square_peps);

    Tnetwork_Storage<IQTensor> square_rvb_storage=peps_to_tnetwork_storage(square_peps);
    ss.str("");
    ss << "/home/jiangsb/tn_ying/tensor_network/result/tnetwork_storage/square_rvb_D=6_Lx=8_Ly=8_mu12=1_step=1e-3_optimized";
    file_name=ss.str();
    writeToFile(file_name,square_rvb_storage);
    // */
    
    return 0;
}
