
//#include "peps.h"
//#include "square_rvb.h"
//#include "kagome_rvb.h"
//#include "simple_update.h"
#include "simple_update_patch_general.h"

int main()
{
    kagome_psg::mu_12=-1;
    kagome_psg::mu_c6=1;

    auto kagome_rvb=kagome_cirac_srvb_peps(4,4);

    //Kagome_Cirac_Lattice_Torus kagome_cirac_lattice({4,4});
    //IQPEPS_IndexSet_SpinHalf index_set(6,kagome_cirac_lattice);
    //IQPEPS kagome_rvb(kagome_cirac_lattice,index_set);
    //random_init_kagome_rvb_cirac_peps(kagome_rvb);

    //IQTPO SdotS_kagome_cirac=SpinSpin_kagome_cirac({kagome_rvb.phys_legs(0),kagome_rvb.phys_legs(1),kagome_rvb.phys_legs(2)});
    //PrintDat(SdotS_kagome_cirac);
    //PrintDat(SdotS_kagome_cirac.site_tensors(0)*SdotS_kagome_cirac.bond_tensors(0)*SdotS_kagome_cirac.site_tensors(1)*SdotS_kagome_cirac.site_tensors(2));

    //IQTPO trotter_gate=trotter_gate_kagome_cirac({kagome_rvb.phys_legs(0),kagome_rvb.phys_legs(1),kagome_rvb.phys_legs(2)},0.1);
    //PrintDat(trotter_gate);
    //PrintDat(trotter_gate.site_tensors(0)*trotter_gate.bond_tensors(0)*trotter_gate.site_tensors(1)*trotter_gate.site_tensors(2));

    std::array<std::vector<IQTensor>,2> kagome_env_tens;
    IQTensor site_tensA=kagome_rvb.site_tensors(0)*kagome_rvb.bond_tensors(0)*kagome_rvb.site_tensors(1)*kagome_rvb.site_tensors(2);
    IQTensor site_tensB=kagome_rvb.site_tensors({1,0,0})*kagome_rvb.bond_tensors({0,-1,1})*kagome_rvb.site_tensors({1,-1,1});
    get_env_tensor_minimization(site_tensA,site_tensB,kagome_env_tens);

    General_Patch_RDM<IQTensor> kagome_tree_patch("tree shape I",kagome_rvb,kagome_env_tens[0][0],{0,1,2},{0,1,2});

    Print(kagome_tree_patch.wf_norm());
    Print(heisenberg_energy_from_RDM(kagome_tree_patch));

    return 0;
}

