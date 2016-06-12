
#include "kagome_rvb.h"
#include "simple_update.h"

int main()
{
    kagome_psg::mu_12=-1;
    kagome_psg::mu_c6=-1;

    int D=6, Lx=8, Ly=8;
    Kagome_Normal_Lattice_Torus kagome_lattice({Lx,Ly});
    IQPEPS kagome_rvb(kagome_lattice);

    //reading PEPS from file
    std::stringstream ss;
    ss << "/home/jiangsb/tn_ying/tensor_network/result/peps_storage/kagome_simple_update/kagome_rvb_normal_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_step=1e-3";
    std::string file_name=ss.str();
    readFromFile(file_name,kagome_rvb);


    //change to another size of torus
    Lx=8;
    Ly=8;
    Kagome_Normal_Lattice_Torus new_kagome_lattice({Lx,Ly});
    IQPEPS_IndexSet_SpinHalf new_index_set(D,new_kagome_lattice);

    std::vector<IQTensor> site_tensors_uc, bond_tensors_uc;
    for (int sitei=0; sitei<kagome_lattice.n_sites_uc(); sitei++) site_tensors_uc.push_back(kagome_rvb.site_tensors(sitei));
    for (int bondi=0; bondi<kagome_rvb.n_bonds_uc(); bondi++) bond_tensors_uc.push_back(kagome_rvb.bond_tensors(bondi));

    IQPEPS new_kagome_rvb(new_kagome_lattice,new_index_set,site_tensors_uc,bond_tensors_uc,kagome_psg::mu_12);

    //stores new_kagome_rvb as tnetwork_storage
    Tnetwork_Storage<IQTensor> kagome_rvb_storage=peps_to_tnetwork_storage(new_kagome_rvb);
    ss.str(std::string());
    ss.clear();
    ss << "/home/jiangsb/tn_ying/tensor_network/result/tnetwork_storage/kagome_simple_update/kagome_rvb_normal_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_step=1e-3";
    file_name=ss.str();
    writeToFile(file_name,kagome_rvb_storage);

    return 0;
}
