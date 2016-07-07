
#include "peps_itebd.h"

template <class TensorT>
PEPSt_iTEBD<TensorT>::PEPSt_iTEBD(const Tnetwork_Storage<TensorT> &peps_storage, const Args &itebd_opts):
    peps_storage_(peps_storage),
    itebd_opts_(itebd_opts)
{
    //init cutting_tensors
    if (peps_storage_._tnetwork_type==8)
    {
        cutting_tensors_.push_back(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(1)));
        cutting_tensors_.push_back(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(0)));
    }

    init_impo();
}
template
PEPSt_iTEBD<ITensor>::PEPSt_iTEBD(const Tnetwork_Storage<ITensor> &peps_storage, const Args &itebd_opts);
template
PEPSt_iTEBD<IQTensor>::PEPSt_iTEBD(const Tnetwork_Storage<IQTensor> &peps_storage, const Args &itebd_opts);


template <class TensorT>
void PEPSt_iTEBD<TensorT>::init_impo()
{
    std::string contract_method=itebd_opts_.getString("ContractMethod","normal");
    //we consider case where Lx,Ly>=4
    int Lx=peps_storage_._Lx,
        Ly=peps_storage_._Ly;
    const auto &tensor_list=peps_storage_._tensor_list;
    const auto &coor_to_siteind=peps_storage_._coor_to_siteind;

    std::vector<std::vector<TensorT>> multicols_ket_tensors;
    std::vector<std::vector<IndexT>> multicols_linds, multicols_rinds, multicols_udinds;
    //implement kagome lattice
    //we set one step to be four cols
    if (peps_storage_._tnetwork_type==8 && contract_method.find("normal")!=std::string::npos)
    {
        //init ket tensors 
        multicols_ket_tensors.push_back(std::vector<TensorT>{tensor_list(coor_to_siteind(0,0)(1)),tensor_list(coor_to_siteind(0,0)(0))});
        multicols_ket_tensors.push_back(std::vector<TensorT>{tensor_list(coor_to_siteind(0,0)(2))});
        multicols_ket_tensors.push_back(std::vector<TensorT>{tensor_list(coor_to_siteind(1,0)(1)),tensor_list(coor_to_siteind(1,0)(0))});
        multicols_ket_tensors.push_back(std::vector<TensorT>{tensor_list(coor_to_siteind(1,0)(2))});

        //Print(multicols_ket_tensors);

        //init indices
        multicols_linds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(0,0)(1)),tensor_list(coor_to_siteind(Lx-1,1)(2))),commonIndex(tensor_list(coor_to_siteind(0,0)(0)),tensor_list(coor_to_siteind(Lx-1,0)(2)))});
            multicols_linds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(0,0)(2)),tensor_list(coor_to_siteind(0,0)(1))),commonIndex(tensor_list(coor_to_siteind(0,0)(2)),tensor_list(coor_to_siteind(0,0)(0)))});
            multicols_linds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(1,0)(1)),tensor_list(coor_to_siteind(0,1)(2))),commonIndex(tensor_list(coor_to_siteind(1,0)(0)),tensor_list(coor_to_siteind(0,0)(2)))});
            multicols_linds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(1,0)(2)),tensor_list(coor_to_siteind(1,0)(1))),commonIndex(tensor_list(coor_to_siteind(1,0)(2)),tensor_list(coor_to_siteind(1,0)(0)))});

            multicols_rinds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(0,0)(1)),tensor_list(coor_to_siteind(0,0)(2))),commonIndex(tensor_list(coor_to_siteind(0,0)(0)),tensor_list(coor_to_siteind(0,0)(2)))});
        multicols_rinds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(0,0)(2)),tensor_list(coor_to_siteind(1,0)(0))),commonIndex(tensor_list(coor_to_siteind(0,0)(2)),tensor_list(coor_to_siteind(1,Ly-1)(1)))});
        multicols_rinds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(1,0)(1)),tensor_list(coor_to_siteind(1,0)(2))),commonIndex(tensor_list(coor_to_siteind(1,0)(0)),tensor_list(coor_to_siteind(1,0)(2)))});
        multicols_rinds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(1,0)(2)),tensor_list(coor_to_siteind(2,0)(0))),commonIndex(tensor_list(coor_to_siteind(1,0)(2)),tensor_list(coor_to_siteind(2,Ly-1)(1)))});

        multicols_udinds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(0,0)(1)),tensor_list(coor_to_siteind(0,1)(0))),commonIndex(tensor_list(coor_to_siteind(0,0)(0)),tensor_list(coor_to_siteind(0,0)(1))),commonIndex(tensor_list(coor_to_siteind(0,Ly-1)(1)),tensor_list(coor_to_siteind(0,0)(0)))});
        multicols_udinds.push_back(std::vector<IndexT>{});
        multicols_udinds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(1,0)(1)),tensor_list(coor_to_siteind(1,1)(0))),commonIndex(tensor_list(coor_to_siteind(1,0)(0)),tensor_list(coor_to_siteind(1,0)(1))),commonIndex(tensor_list(coor_to_siteind(1,Ly-1)(1)),tensor_list(coor_to_siteind(1,0)(0)))});
        multicols_udinds.push_back(std::vector<IndexT>{});

        //Print(multicols_linds);
        //Print(multicols_rinds);
        //Print(multicols_udinds);

        //init impos
        lcols_dl_impos_.clear();
        lcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_one",multicols_ket_tensors[0],multicols_linds[0],multicols_rinds[0],multicols_udinds[0]));
        lcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_two",multicols_ket_tensors[1],multicols_linds[1],multicols_rinds[1]));
        lcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_one",multicols_ket_tensors[2],multicols_linds[2],multicols_rinds[2],multicols_udinds[2]));
        lcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_two",multicols_ket_tensors[3],multicols_linds[3],multicols_rinds[3]));

        rcols_dl_impos_.clear();
        rcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_one",multicols_ket_tensors[0],multicols_rinds[0],multicols_linds[0],multicols_udinds[0]));
        rcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_two",multicols_ket_tensors[3],multicols_rinds[3],multicols_linds[3]));
        rcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_one",multicols_ket_tensors[2],multicols_rinds[2],multicols_linds[2],multicols_udinds[2]));
        rcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_two",multicols_ket_tensors[1],multicols_rinds[1],multicols_linds[1]));

        //Print(lcols_dl_impos_);
        //Print(rcols_dl_impos_);
    }
    
    //implememt kagome lattice with one more extra delta tensor
    if (peps_storage_._tnetwork_type==8 && contract_method.find("extra_delta_tensor")!=std::string::npos)
    {
        //init indices
        multicols_linds.push_back(std::vector<IndexT>{isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(0,0)(1)),tensor_list(coor_to_siteind(Lx-1,1)(2))),"virt_leg_(0,0,v)-(-1,0,r)"),commonIndex(tensor_list(coor_to_siteind(0,0)(0)),tensor_list(coor_to_siteind(Lx-1,0)(2)))});
        multicols_linds.push_back(std::vector<IndexT>{isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(0,0)(2)),tensor_list(coor_to_siteind(0,0)(1))),"virt_leg_(0,0,r)-(0,0,v)"),commonIndex(tensor_list(coor_to_siteind(0,0)(2)),tensor_list(coor_to_siteind(0,0)(0)))});
        multicols_linds.push_back(std::vector<IndexT>{isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(1,0)(1)),tensor_list(coor_to_siteind(0,1)(2))),"virt_leg_(1,0,v)-(0,0,r)"),commonIndex(tensor_list(coor_to_siteind(1,0)(0)),tensor_list(coor_to_siteind(0,0)(2)))});
        multicols_linds.push_back(std::vector<IndexT>{isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(1,0)(2)),tensor_list(coor_to_siteind(1,0)(1))),"virt_leg_(1,0,r)-(1,0,v)"),commonIndex(tensor_list(coor_to_siteind(1,0)(2)),tensor_list(coor_to_siteind(1,0)(0)))});

        multicols_rinds.push_back(dag<IndexT>(multicols_linds[1]));
        multicols_rinds.push_back(dag<IndexT>(multicols_linds[2]));
        multicols_rinds.push_back(dag<IndexT>(multicols_linds[3]));
        multicols_rinds.push_back(std::vector<IndexT>{isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(1,1)(2)),tensor_list(coor_to_siteind(2,0)(1))),"virt_leg_(1,0,r)-(2,0,v)"),commonIndex(tensor_list(coor_to_siteind(1,0)(2)),tensor_list(coor_to_siteind(2,0)(0)))});

        multicols_udinds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(0,0)(1)),tensor_list(coor_to_siteind(0,1)(0))),commonIndex(tensor_list(coor_to_siteind(0,0)(0)),tensor_list(coor_to_siteind(0,0)(1))),commonIndex(tensor_list(coor_to_siteind(0,Ly-1)(1)),tensor_list(coor_to_siteind(0,0)(0)))});
        multicols_udinds.push_back(std::vector<IndexT>{isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(1,0)(1)),tensor_list(coor_to_siteind(0,1)(2))),"virt_leg_(0,0,r)-(0,1,w)"),isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(0,0)(2)),tensor_list(coor_to_siteind(0,0)(1))),"virt_leg_(0,0,w)-(0,0,r)"),isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(1,Ly-1)(1)),tensor_list(coor_to_siteind(0,0)(2))),"virt_leg_(0,-1,r)-(0,0,w)")});
        multicols_udinds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(1,0)(1)),tensor_list(coor_to_siteind(1,1)(0))),commonIndex(tensor_list(coor_to_siteind(1,0)(0)),tensor_list(coor_to_siteind(1,0)(1))),commonIndex(tensor_list(coor_to_siteind(1,Ly-1)(1)),tensor_list(coor_to_siteind(1,0)(0)))});
        multicols_udinds.push_back(std::vector<IndexT>{isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(2,0)(1)),tensor_list(coor_to_siteind(1,1)(2))),"virt_leg_(1,0,r)-(1,1,w)"),isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(1,0)(2)),tensor_list(coor_to_siteind(1,0)(1))),"virt_leg_(1,0,w)-(1,0,r)"),isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(2,Ly-1)(1)),tensor_list(coor_to_siteind(1,0)(2))),"virt_leg_(1,-1,r)-(1,0,w)")});

        //init ket tensors
        TensorT v_tensor=tensor_list(coor_to_siteind(0,0)(1));
        v_tensor.replaceIndex(commonIndex(v_tensor,tensor_list(coor_to_siteind(Lx-1,1)(2))),multicols_linds[0][0]);
        v_tensor.replaceIndex(commonIndex(v_tensor,tensor_list(coor_to_siteind(0,0)(2))),multicols_rinds[0][0]);
        multicols_ket_tensors.push_back(std::vector<TensorT>{v_tensor,tensor_list(coor_to_siteind(0,0)(0))});

        TensorT r_tensor=delta_tensor<TensorT>(multicols_linds[1][0],dag(multicols_udinds[1][1]))*delta_tensor<TensorT>(multicols_udinds[1][0],multicols_rinds[1][0]),
                w_tensor=tensor_list(coor_to_siteind(0,0)(2));
        w_tensor.replaceIndex(commonIndex(w_tensor,tensor_list(coor_to_siteind(0,0)(1))),multicols_udinds[1][1]);
        w_tensor.replaceIndex(commonIndex(w_tensor,tensor_list(coor_to_siteind(1,Ly-1)(1))),dag(multicols_udinds[1][2]));
        multicols_ket_tensors.push_back(std::vector<TensorT>{r_tensor,w_tensor});

        v_tensor=tensor_list(coor_to_siteind(1,0)(1));
        v_tensor.replaceIndex(commonIndex(v_tensor,tensor_list(coor_to_siteind(0,1)(2))),multicols_linds[2][0]);
        v_tensor.replaceIndex(commonIndex(v_tensor,tensor_list(coor_to_siteind(1,0)(2))),multicols_rinds[2][0]);
        multicols_ket_tensors.push_back(std::vector<TensorT>{v_tensor,tensor_list(coor_to_siteind(1,0)(0))});

        r_tensor=delta_tensor<TensorT>(multicols_linds[3][0],dag(multicols_udinds[3][1]))*delta_tensor<TensorT>(multicols_udinds[3][0],multicols_rinds[3][0]);
        w_tensor=tensor_list(coor_to_siteind(1,0)(2));
        w_tensor.replaceIndex(commonIndex(w_tensor,tensor_list(coor_to_siteind(1,0)(1))),multicols_udinds[3][1]);
        w_tensor.replaceIndex(commonIndex(w_tensor,tensor_list(coor_to_siteind(2,Ly-1)(1))),dag(multicols_udinds[3][2]));
        multicols_ket_tensors.push_back(std::vector<TensorT>{r_tensor,w_tensor});

        //init impos
        lcols_dl_impos_.clear();
        rcols_dl_impos_.clear();
        int total_cols=multicols_ket_tensors.size();
        for (int coli=0; coli<total_cols; coli++)
        {
            lcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_one",multicols_ket_tensors[coli],multicols_linds[coli],multicols_rinds[coli],multicols_udinds[coli]));
            rcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_one",multicols_ket_tensors[(total_cols-coli)%total_cols],multicols_rinds[(total_cols-coli)%total_cols],multicols_linds[(total_cols-coli)%total_cols],multicols_udinds[(total_cols-coli)%total_cols]));
            //rcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_one",multicols_ket_tensors[coli],multicols_rinds[coli],multicols_linds[coli],multicols_udinds[coli]));
        }
    }

    //implement kagome lattice with one more extra delta tensor contract for X dir for zero flux state
    if (peps_storage_._tnetwork_type==8 && contract_method.find("zero_slX")!=std::string::npos)
    {
        //init indices
        //
        //  \ / \ /
        //   u   r
        //  / \ / \ /
        //     v   w
        //    / \ / \
        // 
        multicols_linds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(0,0)(0)),tensor_list(coor_to_siteind(0,Ly-1)(1))),commonIndex(tensor_list(coor_to_siteind(0,0)(0)),tensor_list(coor_to_siteind(Lx-1,0)(2))),isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(Lx-1,0)(1)),tensor_list(coor_to_siteind(Lx-1,0)(2))),"virt_leg_(-1,0,r)-(-1,0,w)"),isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(Lx-1,0)(2)),tensor_list(coor_to_siteind(Lx-1,0)(1))),"virt_leg_(-1,0,r)-(-1,0,v)")});
        multicols_linds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(0,0)(1)),tensor_list(coor_to_siteind(0,0)(0))),isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(0,0)(1)),tensor_list(coor_to_siteind(Lx-1,1)(2))),"virt_leg_(0,0,v)-(-1,0,r)"),isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(Lx-1,1)(2)),tensor_list(coor_to_siteind(0,0)(1))),"virt_leg(-1,1,w)-(-1,0,r)")});

        multicols_rinds.push_back(std::vector<IndexT>{commonIndex(tensor_list(coor_to_siteind(0,0)(0)),tensor_list(coor_to_siteind(0,0)(2))),dag(multicols_linds[1][0]),dag(multicols_linds[1][1]),dag(multicols_linds[1][2])});
        multicols_rinds.push_back(std::vector<IndexT>{isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(0,0)(1)),tensor_list(coor_to_siteind(0,0)(2))),"virt_leg(0,0,v)-(0,0,r)"),commonIndex(tensor_list(coor_to_siteind(0,0)(1)),tensor_list(coor_to_siteind(0,1)(0))),commonIndex(tensor_list(coor_to_siteind(Lx-1,1)(2)),tensor_list(coor_to_siteind(0,1)(0))),isomorphic_legs(commonIndex(tensor_list(coor_to_siteind(Lx-1,1)(2)),tensor_list(coor_to_siteind(Lx-1,1)(1))),"virt_leg(-1,1,w)-(-1,1,r)")});

        //init ket tensors
        TensorT r_tensor=delta_tensor<TensorT>(multicols_linds[0][2],multicols_linds[0][3])*delta_tensor<TensorT>(multicols_rinds[0][2],multicols_rinds[0][3]);
        multicols_ket_tensors.push_back(std::vector<TensorT>{tensor_list(coor_to_siteind(0,0)(0)),r_tensor});

        TensorT v_tensor=tensor_list(coor_to_siteind(0,0)(1)),
                w_tensor=tensor_list(coor_to_siteind(Lx-1,1)(2));
        v_tensor.replaceIndex(commonIndex(v_tensor,w_tensor),multicols_linds[1][1]);
        v_tensor.replaceIndex(commonIndex(v_tensor,tensor_list(coor_to_siteind(0,0)(2))),multicols_rinds[1][0]);
        w_tensor.replaceIndex(commonIndex(w_tensor,v_tensor),multicols_linds[1][2]);
        w_tensor.replaceIndex(commonIndex(w_tensor,tensor_list(Lx-1,1)(1)),multicols_rinds[1][3]);
        multicols_ket_tensors.push_back(std::vector<TensorT>{v_tensor,w_tensor});

        //init impos
        lcols_dl_impos_.clear();
        rcols_dl_impos_.clear();
        int total_cols=multicols_ket_tensors.size();
        for (int coli=0; coli<total_cols; coli++)
        {
            lcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_two",multicols_ket_tensors[coli],multicols_linds[coli],multicols_rinds[coli]));
            rcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_two",multicols_ket_tensors[total_cols-coli-1],multicols_rinds[total_cols-coli-1],multicols_linds[total_cols-coli-1]));
        }
    }

    //TODO:implement kagome lattice with one more extra delta tensor contract for X dir for pi flux state

}
template 
void PEPSt_iTEBD<ITensor>::init_impo();
template
void PEPSt_iTEBD<IQTensor>::init_impo();


    template <class TensorT>
void PEPSt_iTEBD<TensorT>::env_tensors_from_itebd(int n_steps)
{
    //get left and right imps
    cout << endl << "------------------------------------------------------------------------" << endl;
    cout << "contraction left and right imps!" << endl;
    for (int stepi=0; stepi<n_steps; stepi++) 
    {
        Print(stepi);
        update_imps_one_step();
    }

    int Lx=peps_storage_._Lx,
        Ly=peps_storage_._Ly;

    const auto &tensor_list=peps_storage_._tensor_list;
    const auto &coor_to_siteind=peps_storage_._coor_to_siteind;

    cout << endl << "------------------------------------------------------------------------" << endl;
    cout << "contraction for up-down direction!" << endl;

    //obtain env_tensors from left and right imps
    //
    //kagome case
    //for kagome lattice, env tensors are
    //
    //      4
    //      |
    //   0--v--1
    //      |
    //   2--u--3
    //      |
    //      5
    //
    if (peps_storage_._tnetwork_type==8)
    {
        //get virt indices of cutting sites
        std::vector<IndexT> cutting_ket_inds;
        //left up
        cutting_ket_inds.push_back(commonIndex(tensor_list(coor_to_siteind(0,0)(1)),tensor_list(coor_to_siteind(Lx-1,1)(2))));
        //right up
        cutting_ket_inds.push_back(commonIndex(tensor_list(coor_to_siteind(0,0)(1)),tensor_list(coor_to_siteind(0,0)(2))));
        //left down
        cutting_ket_inds.push_back(commonIndex(tensor_list(coor_to_siteind(0,0)(0)),tensor_list(coor_to_siteind(Lx-1,0)(2))));
        //right down
        cutting_ket_inds.push_back(commonIndex(tensor_list(coor_to_siteind(0,0)(0)),tensor_list(coor_to_siteind(0,0)(2))));
        //up
        cutting_ket_inds.push_back(commonIndex(tensor_list(coor_to_siteind(0,0)(1)),tensor_list(coor_to_siteind(0,1)(0))));
        //down
        cutting_ket_inds.push_back(commonIndex(tensor_list(coor_to_siteind(0,0)(0)),tensor_list(coor_to_siteind(0,Ly-1)(1))));

        //obtain result for up-down contraction
        std::string contract_method=itebd_opts_.getString("ContractMethod","normal");
        if (contract_method.find("normal")!=std::string::npos || contract_method.find("extra_delta_tensor")!=std::string::npos)
        {
            //contract_tensors
            std::vector<TensorT> contract_tensors;
            int n_sites_uc=ldl_imps_.n_sites_uc();
            for (int sitei=0; sitei<n_sites_uc; sitei++)
            {
                TensorT limps_tensor=ldl_imps_.dl_site_tensors(sitei),
                        rimps_tensor=rdl_imps_.dl_site_tensors(sitei);
                limps_tensor.replaceIndex(ldl_imps_.ket_siteinds(sitei),dag(cutting_ket_inds[2*sitei]));
                limps_tensor.replaceIndex(ldl_imps_.bra_siteinds(sitei),prime(cutting_ket_inds[2*sitei]));
                rimps_tensor.replaceIndex(rdl_imps_.ket_siteinds(sitei),dag(cutting_ket_inds[2*sitei+1]));
                rimps_tensor.replaceIndex(rdl_imps_.bra_siteinds(sitei),prime(cutting_ket_inds[2*sitei+1]));

                contract_tensors.push_back(limps_tensor);
                contract_tensors.push_back(cutting_tensors_[sitei]);
                contract_tensors.push_back(dag(cutting_tensors_[sitei]).prime(Link));
                contract_tensors.push_back(rimps_tensor);
            }

            //contract_seq
            std::vector<int> VU_contract_seq, VD_contract_seq;
            VU_contract_seq.push_back(-1);
            VD_contract_seq.push_back(-1);
            for (int tensori=0; tensori<contract_tensors.size(); tensori++)
            {
                VU_contract_seq.push_back(tensori);
                VD_contract_seq.push_back(contract_tensors.size()-1-tensori);
            }

            //obtain up and down dominant eigenvector
            //
            //init using random tensor
            // /*
            //up_inds are input(output) inds for VU(VD), and down_inds are output(input) inds for VU(VD)
            std::vector<IndexT> up_inds, down_inds;
            up_inds.push_back(ldl_imps_.virt_inds().front());
            up_inds.push_back(cutting_ket_inds[4]);
            up_inds.push_back(dag(cutting_ket_inds[4]).prime());
            up_inds.push_back(rdl_imps_.virt_inds().front());

            down_inds.push_back(dag(ldl_imps_.virt_inds().back()));
            down_inds.push_back(cutting_ket_inds[5]);
            down_inds.push_back(dag(cutting_ket_inds[5]).prime());
            down_inds.push_back(dag(rdl_imps_.virt_inds().back()));

            std::vector<IndexT> up_inds_dag=dag<IndexT>(up_inds),
                down_inds_dag=dag<IndexT>(down_inds);
            TensorT VU(up_inds_dag),
                    VD(down_inds_dag);
            //we assume the case where dominant eigenvector has sz==0
            randTensor(VU);
            randTensor(VD);
            // */

            TensorT_Matrix_Arnoldi<TensorT> UMat(up_inds,down_inds,contract_tensors,VU_contract_seq),
                DMat(down_inds,up_inds,contract_tensors,VD_contract_seq);
            Args arnoldi_opts=itebd_opts_;
            Complex eta_U, eta_D;
            do
            {
                arnoldi_opts.add("MaxIter",arnoldi_opts.getInt("MaxIter")+10);
                eta_U=arnoldi(UMat,VU,itebd_opts_);
                eta_D=arnoldi(DMat,VD,itebd_opts_);
                cout << "MaxIter: " << arnoldi_opts.getInt("MaxIter") << endl;
                Print(eta_U);
                Print(eta_D);
            }
            while ((std::abs((eta_U-eta_D)/eta_U)>1e-8) && (arnoldi_opts.getInt("MaxIter")<=100));

            //get env_tensors_
            env_tensors_.clear();
            env_tensors_.push_back(contract_tensors[0]);
            env_tensors_.push_back(contract_tensors[3]);
            env_tensors_.push_back(contract_tensors[4]);
            env_tensors_.push_back(contract_tensors[7]);
            env_tensors_.push_back(VU);
            env_tensors_.push_back(VD);
        }

        /*
        if (contract_method.find("zero_slX")!=std::string::npos)
        {
            std::vector<TensorT> contract_tensors;
            int n_sites_uc=2;
            for (int sitei=0; sitei<n_sites_uc; sitei++)
            {
                TenosrT limps_tensor=ldl_imps_.dl_site_tensors(2*(n_sites_uc-sitei))*ldl_imps_.dl_site_tensors(2*(n_sites_uc-sitei)+1),
                        rimps_tensor=rdl_imps_.dl_site_tensors(2*(n_sites_uc-sitei))*rdl_imps_.dl_site_tensors(2*(n_sites_uc-sitei)+1);

                limps_tensor.replaceIndex(ldl_imps_.ket_siteinds(2*(n_sites_uc-sitei)),dag());
                    
                contract_tensors.push_back();
            }
        }
        */

        //TODO: implement pi_slX

        //get wf_norm_
        wf_norm_=(env_tensors_[4]*env_tensors_[0]*tensor_list(1)*dag(prime(tensor_list(1)).noprime(Site))*env_tensors_[1]*env_tensors_[2]*tensor_list(0)*dag(prime(tensor_list(0)).noprime(Site))*env_tensors_[3]*env_tensors_[5]).toComplex();
        Print(wf_norm_);
    }

}
template
void PEPSt_iTEBD<ITensor>::env_tensors_from_itebd(int n_steps);
template
void PEPSt_iTEBD<IQTensor>::env_tensors_from_itebd(int n_steps);


template <class TensorT>
void PEPSt_iTEBD<TensorT>::update_imps_one_step()
{
    std::string contract_method=itebd_opts_.getString("ContractMethod","normal");

    //kagome case
    if (peps_storage_._tnetwork_type==8 && contract_method.find("normal")!=std::string::npos)
    {
        //if default construct, then we init the imps using first col
        cout << endl << "------------------------------------------------------------------------" << endl;
        if (!ldl_imps_.valid())
        {
            cout << "init left imps!" << endl;
            ldl_imps_=DL_iMPSt<TensorT>(lcols_dl_impos_[0].ket_tensors(),lcols_dl_impos_[0].ket_outgoing_inds(),lcols_dl_impos_[0].ket_virt_inds());
        }
        else
        {
            cout << "contract for left first col!" << endl;
            contract_dl_impo_imps(ldl_imps_,lcols_dl_impos_[0],itebd_opts_);
        }
        if (!rdl_imps_.valid())
        {
            cout << "init right imps!" << endl;
            rdl_imps_=DL_iMPSt<TensorT>(rcols_dl_impos_[0].ket_tensors(),rcols_dl_impos_[0].ket_outgoing_inds(),rcols_dl_impos_[0].ket_virt_inds());
        }
        else
        {
            cout << "contract for right first col!" << endl;
            contract_dl_impo_imps(rdl_imps_,rcols_dl_impos_[0],itebd_opts_);
        }
        rdl_imps_.move_tensors();

        //contract for second col
        cout << endl << "------------------------------------------------------------------------" << endl;
        cout << "contract for left second col!" << endl;
        contract_dl_impo_imps(ldl_imps_,lcols_dl_impos_[1],itebd_opts_);
        ldl_imps_.move_tensors();
        cout << "contract for right second col!" << endl;
        contract_dl_impo_imps(rdl_imps_,rcols_dl_impos_[1],itebd_opts_);

        //contract for third col
        cout << endl << "------------------------------------------------------------------------" << endl;
        cout << "contract for left third col!" << endl;
        contract_dl_impo_imps(ldl_imps_,lcols_dl_impos_[2],itebd_opts_);
        cout << "contract for right third col!" << endl;
        contract_dl_impo_imps(rdl_imps_,rcols_dl_impos_[2],itebd_opts_);
        rdl_imps_.move_tensors();

        //contract for fourth col
        cout << endl << "------------------------------------------------------------------------" << endl;
        cout << "contract for left fourth col!" << endl;
        contract_dl_impo_imps(ldl_imps_,lcols_dl_impos_[3],itebd_opts_);
        ldl_imps_.move_tensors();
        cout << "contract for right fourth col!" << endl;
        contract_dl_impo_imps(rdl_imps_,rcols_dl_impos_[3],itebd_opts_);
    }

    if (peps_storage_._tnetwork_type==8 && contract_method.find("extra_delta_tensor")!=std::string::npos)
    {
        //if default construct, then we init the imps using first col
        cout << endl << "------------------------------------------------------------------------" << endl;
        int n_cols=lcols_dl_impos_.size();
        if (!ldl_imps_.valid())
        {
            cout << "init left imps!" << endl;
            ldl_imps_=DL_iMPSt<TensorT>(lcols_dl_impos_[n_cols-1].ket_tensors(),lcols_dl_impos_[n_cols-1].ket_outgoing_inds(),lcols_dl_impos_[n_cols-1].ket_virt_inds());
        }
        if (!rdl_imps_.valid())
        {
            cout << "init right imps!" << endl;
            rdl_imps_=DL_iMPSt<TensorT>(rcols_dl_impos_[n_cols-1].ket_tensors(),rcols_dl_impos_[n_cols-1].ket_outgoing_inds(),rcols_dl_impos_[n_cols-1].ket_virt_inds());
        }

        for (int coli=0; coli<n_cols; coli++)
        {
            cout << endl << "------------------------------------------------------------------------" << endl;
            Print(coli);
            cout << "left col: " << endl;
            contract_dl_impo_imps(ldl_imps_,lcols_dl_impos_[coli],itebd_opts_);
            cout << "right col: " << endl;
            contract_dl_impo_imps(rdl_imps_,rcols_dl_impos_[coli],itebd_opts_);
        }

    }
}
template
void PEPSt_iTEBD<ITensor>::update_imps_one_step();
template
void PEPSt_iTEBD<IQTensor>::update_imps_one_step();


template <class TensorT>
Complex PEPSt_iTEBD<TensorT>::expect_val_from_env_tensors(const std::vector<TensorT> &two_sites_mpo) const
{
    //kagome case
    if (peps_storage_._tnetwork_type==8)
    {
        TensorT result_tensor=env_tensors_[4]*env_tensors_[0]*peps_storage_._tensor_list(1)*two_sites_mpo[0]*dag(prime(peps_storage_._tensor_list(1)))*env_tensors_[1]*env_tensors_[2]*peps_storage_._tensor_list(0)*two_sites_mpo[1]*dag(prime(peps_storage_._tensor_list(0)))*env_tensors_[3]*env_tensors_[5];
        return result_tensor.toComplex()/wf_norm_;
    }
}
template
Complex PEPSt_iTEBD<ITensor>::expect_val_from_env_tensors(const std::vector<ITensor> &two_sites_mpo) const;
template
Complex PEPSt_iTEBD<IQTensor>::expect_val_from_env_tensors(const std::vector<IQTensor> &two_sites_mpo) const;
