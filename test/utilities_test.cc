
#include "kagome_rvb.h"
#include "trotter_gate.h"
#include "peps_itebd.h"

int main()
{
    kagome_psg::mu_12=1;
    kagome_psg::mu_c6=1;
    
    Print(kagome_psg::mu_12);
    Print(kagome_psg::mu_c6);

    int Lx=4, Ly=4, D=3;

    IQPEPS kagome_rvb=kagome_normal_srvb_peps(Lx,Ly);
    //Kagome_Normal_Lattice_Torus kagome_lattice({Lx,Ly});
    //IQPEPS_IndexSet_SpinHalf index_set(3,kagome_lattice);
    //IQPEPS kagome_rvb(kagome_lattice,index_set);
    //random_init_kagome_rvb_normal_peps(kagome_rvb);

    Tnetwork_Storage<IQTensor> kagome_rvb_storage=peps_to_tnetwork_storage(kagome_rvb);

    PEPSt_iTEBD<IQTensor> kagome_itebd(kagome_rvb_storage,{"Maxm",30,"MaxIter",10,"MaxRestart",1,"Cutoff",1e-16,"ErrGoal",1e-16});
    kagome_itebd.env_tensors_from_itebd(1);

    NN_Heisenberg_Hamiltonian heisenberg_gate({findtype(kagome_rvb_storage._tensor_list(1),Site),findtype(kagome_rvb_storage._tensor_list(0),Site)});

    for (int measurei=0; measurei<30; measurei++)
    {
        //measure Heisenberg energy
        Print(measurei);
        Complex energy=kagome_itebd.expect_val_from_env_tensors({heisenberg_gate.site_tensors(0)*heisenberg_gate.bond_tensor(),heisenberg_gate.site_tensors(1)});
        Print(energy);
        kagome_itebd.env_tensors_from_itebd(1);
    }

    return 0;
}

