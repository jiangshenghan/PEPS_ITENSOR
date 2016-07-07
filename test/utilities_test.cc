
#include "kagome_rvb.h"
#include "trotter_gate.h"
#include "peps_itebd.h"

int main()
{
    //kagome_psg::mu_12=1;
    //kagome_psg::mu_c6=1;
    //
    //Print(kagome_psg::mu_12);
    //Print(kagome_psg::mu_c6);

    int Lx=4, Ly=4, D=3;

    //IQPEPS kagome_rvb=kagome_normal_srvb_peps(Lx,Ly);
    //Kagome_Normal_Lattice_Torus kagome_lattice({Lx,Ly});
    //IQPEPS_IndexSet_SpinHalf index_set(3,kagome_lattice);
    //IQPEPS kagome_rvb(kagome_lattice,index_set);
    //random_init_kagome_rvb_normal_peps(kagome_rvb);

    //Tnetwork_Storage<IQTensor> kagome_rvb_storage=peps_to_tnetwork_storage(kagome_rvb);
    Tnetwork_Storage<IQTensor> kagome_rvb_storage;
    std::string file_name="/home/jiangsb/tn_ying/kagome_srvb/kagome_original_srvb_eta12_0_etaC6_0_4_by_4";
    readFromFile(file_name,kagome_rvb_storage);
    Print(file_name);
    cout << "reading successfully" << endl;

    Args itebd_opts("Maxm",30,"MaxIter",10,"MaxRestart",1,"Cutoff",1e-16,"ErrGoal",1e-16,"ContractMethod","normal","AbsorbBond","right_bond");
    Print(itebd_opts);
    NN_Heisenberg_Hamiltonian heisenberg_gate({findtype(kagome_rvb_storage._tensor_list(1),Site),findtype(kagome_rvb_storage._tensor_list(0),Site)});

    //IQTensor
    /*
    PEPSt_iTEBD<IQTensor> kagome_itebd(kagome_rvb_storage,itebd_opts);
    kagome_itebd.env_tensors_from_itebd(1);
    for (int measurei=0; measurei<30; measurei++)
    {
        //measure Heisenberg energy
        Print(measurei);
        Complex energy=kagome_itebd.expect_val_from_env_tensors({heisenberg_gate.site_tensors(0)*heisenberg_gate.bond_tensor(),heisenberg_gate.site_tensors(1)});
        Print(energy);
        kagome_itebd.env_tensors_from_itebd(1);
    }
    // */

    //ITensor
    // /*
    Tnetwork_Storage<ITensor> kagome_peps_storage=tnetwork_storage_ITensor_from_IQTensor(kagome_rvb_storage);
    PEPSt_iTEBD<ITensor> kagome_itebd(kagome_peps_storage,itebd_opts);
    kagome_itebd.env_tensors_from_itebd(1);
    for (int measurei=0; measurei<30; measurei++)
    {
        //measure Heisenberg energy
        Print(measurei);
        Complex energy=kagome_itebd.expect_val_from_env_tensors({(heisenberg_gate.site_tensors(0)*heisenberg_gate.bond_tensor()).toITensor(),heisenberg_gate.site_tensors(1).toITensor()});
        Print(energy);
        kagome_itebd.env_tensors_from_itebd(1);
    }
    // */


    //test arnoldi of imps
    /*
    PEPSt_iTEBD<IQTensor> kagome_itebd(kagome_rvb_storage,itebd_opts);
    std::vector<DL_iMPOt<IQTensor>> dl_impos_test=kagome_itebd.lcols_dl_impos();
    DL_iMPSt<IQTensor> lb_dl_imps_test=DL_iMPSt<IQTensor>(dl_impos_test.back().ket_tensors(),dl_impos_test.back().ket_outgoing_inds(),dl_impos_test.back().ket_virt_inds()),
                       rb_dl_imps_test=lb_dl_imps_test;
    
    itebd_opts.add("AbsorbBond","left_bond");
    contract_dl_impo_imps(lb_dl_imps_test,dl_impos_test[0],itebd_opts);
    contract_dl_impo_imps(lb_dl_imps_test,dl_impos_test[1],itebd_opts);
    Print(lb_dl_imps_test);
    itebd_opts.add("AbsorbBond","right_bond");
    contract_dl_impo_imps(rb_dl_imps_test,dl_impos_test[0],itebd_opts);
    contract_dl_impo_imps(rb_dl_imps_test,dl_impos_test[1],itebd_opts);
    Print(rb_dl_imps_test);
    */
    

    return 0;
}

