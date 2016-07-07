
#include "kagome_rvb.h"
#include "trotter_gate.h"
#include "peps_itebd.h"
#include "mkl.h"

int main(int argc, char *argv[])
{
    //to run this file, type e.g.
    // ./itebd_measure.exe file_path/tnetwork_storage_file_name maxm normal (threads_num)

    //int threads_num=std::atoi(argv[4]);
    //mkl_set_num_threads(threads_num);

    Tnetwork_Storage<IQTensor> tnetwork_storage;
    readFromFile(argv[1],tnetwork_storage);
    cout << "Reading files successfullt!" << endl;
    cout << argv[1] << endl;

    int maxm=std::atoi(argv[2]);
    Args itebd_opts("Maxm",maxm,"MaxIter",10,"MaxRestart",0,"Cutoff",1e-16,"ErrGoal",1e-16,"ContractMethod",argv[3]);
    Print(itebd_opts);

    NN_Heisenberg_Hamiltonian heisenberg_gate({findtype(tnetwork_storage._tensor_list(1),Site),findtype(tnetwork_storage._tensor_list(0),Site)});

    //IQTensor
    /*
    PEPSt_iTEBD<IQTensor> peps_itebd(tnetwork_storage,itebd_opts);
    peps_itebd.env_tensors_from_itebd(1);
    for (int measurei=0; measurei<50; measurei++)
    {
        //measure Heisenberg energy
        Print(measurei);
        Complex energy=peps_itebd.expect_val_from_env_tensors({heisenberg_gate.site_tensors(0)*heisenberg_gate.bond_tensor(),heisenberg_gate.site_tensors(1)});
        Print(energy);
        peps_itebd.env_tensors_from_itebd(1);
    }
    // */

    //ITensor
    // /*
    Tnetwork_Storage<ITensor> itensor_tnetwork_storage=tnetwork_storage_ITensor_from_IQTensor(tnetwork_storage);
    PEPSt_iTEBD<ITensor> peps_itebd(itensor_tnetwork_storage,itebd_opts);
    peps_itebd.env_tensors_from_itebd(1);
    
    for (int measurei=0; measurei<50; measurei++)
    {
        Print(measurei);
        Complex energy=peps_itebd.expect_val_from_env_tensors({(heisenberg_gate.site_tensors(0)*heisenberg_gate.bond_tensor()).toITensor(),heisenberg_gate.site_tensors(1).toITensor()});
        Print(energy);
        peps_itebd.env_tensors_from_itebd(1);
    }
    // */

    return 0;
}
