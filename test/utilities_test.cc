
#include "kagome_rvb.h"
#include "full_update.h"
#include "cluster_env.h"


int main()
{
    kagome_psg::mu_12=1;
    kagome_psg::mu_c6=1;

    int Lx=8, Ly=8, D=3;
    Kagome_Normal_Lattice_Torus kagome_lattice({Lx,Ly});

    //IQPEPS_IndexSet_SpinHalf index_set(D,kagome_lattice);
    //IQPEPS kagome_rvb(kagome_lattice,index_set);
    //random_init_kagome_rvb_normal_peps(kagome_rvb);
    
    IQPEPS kagome_rvb=kagome_normal_srvb_peps(Lx,Ly);

    Cluster_Env kagome_cluster_env(kagome_rvb.name()+"triangle shape",{0,2},{kagome_rvb.site_tensors(0),kagome_rvb.site_tensors(1)*kagome_rvb.bond_tensors(0)*kagome_rvb.bond_tensors(2)*kagome_rvb.bond_tensors(3),kagome_rvb.site_tensors(2)*kagome_rvb.bond_tensors(1)*kagome_rvb.bond_tensors(4)});
    kagome_cluster_env.obtain_sl_env_iterative_nodeg();

    //obtain fake energy
    std::vector<int> cutting_sites={3*(Lx+1),3*(Lx+1)+2}, 
                     cutting_bonds={6*(Lx+1),6*(Lx+1)+1,6*(Lx+1)+4};

    std::array<std::vector<IQTensor>,2> su_env_mats;
    init_env_tensor(kagome_rvb.site_tensors(cutting_sites[0])*kagome_rvb.bond_tensors(cutting_bonds[0])*kagome_rvb.bond_tensors(cutting_bonds[1]),kagome_rvb.site_tensors(cutting_sites[1])*kagome_rvb.bond_tensors(cutting_bonds[2]),kagome_cluster_env.sl_env_tensor(),su_env_mats);
    PrintDat(su_env_mats[0]);
    PrintDat(su_env_mats[1]);
    //obtain env_tensors
    std::vector<IQTensor> env_tensors;
    obtain_kagome_normal_env_MPO(1,cutting_sites,cutting_bonds,su_env_mats,kagome_rvb,env_tensors);

    //init kagome RDM
    PEPSt_RDM<IQTensor> kagome_normal_rdm("two sites shape",cutting_sites,cutting_bonds,env_tensors,kagome_rvb,{0,0,0,0,1});

    Print(heisenberg_energy_from_RDM(kagome_normal_rdm));

    return 0;
}

