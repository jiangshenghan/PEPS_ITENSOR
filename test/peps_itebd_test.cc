
#include "kagome_rvb.h"
#include "peps_itebd.h"

int main()
{
    kagome_psg::mu_12=1;
    kagome_psg::mu_c6=1;
    
    Print(kagome_psg::mu_12);
    Print(kagome_psg::mu_c6);

    int Lx=8, Ly=8, D=3;

    IQPEPS kagome_rvb=kagome_normal_srvb_peps(Lx,Ly);
    Tnetwork_Storage<IQTensor> kagome_rvb_storage=peps_to_tnetwork_storage(kagome_rvb);

    PEPS_iTEBD<IQTensor> kagome_itebd(kagome_rvb_storage,{"Maxm",40,"MaxIter",10,"MaxRestart",1,"Cutoff",1e-16,"ErrGoal",1e-16});
    std::vector<IQTensor> env_tensors=kagome_itebd.env_tensors_from_itebd(4);

    NN_Heisenberg_Hamiltonian heisenberg_gate{(findtype(kagome_rvb_storage._tensor_list(1),Site),findtype(kagome_rvb_storage._tensor_list(0),Site)});

    for (int measurei=0; measurei<30; measurei++)
    {
        //measure Heisenberg energy
        Complex energy=kagome_itebd.expect_val_from_env_tensors({heisenberg_gate.site_tensors(0)*heisenberg_gate.bond_tensor(),heisenberg_gate.site_tensors(1)});
        Print(energy);
        env_tensors=kagome_itebd.env_tensors_from_itebd(2);
    }

    return 0;
}
