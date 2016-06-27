
#include "peps_itebd.h"

template <class TensorT>
PEPSt_iTEBD<TensorT>::PEPSt_iTEBD(const Tnetwork_Storage<TensorT> &peps_storage, const Args &itebd_opts):
    peps_storage_(peps_storage),
    itebd_opts_(itebd_opts)
{
    init_impo();
}
template
PEPSt_iTEBD<ITensor>::PEPSt_iTEBD(const Tnetwork_Storage<ITensor> &peps_storage, const Args &itebd_opts);
template
PEPSt_iTEBD<IQTensor>::PEPSt_iTEBD(const Tnetwork_Storage<IQTensor> &peps_storage, const Args &itebd_opts);


template <class TensorT>
void PEPSt_iTEBD<TensorT>::init_impo()
{
    //we consider case where Lx,Ly>=4
    int Lx=peps_storage_._Lx,
        Ly=peps_storage_._Ly;

    std::vector<std::vector<TensorT>> multicols_ket_tensors;
    std::vector<std::vector<IndexT>> multicols_linds, multicols_rinds, multicols_udinds;
    //implement kagome lattice
    //we set one step to be four cols
    if (peps_storage_._tnetwork_type==8)
    {
        //init ket tensors 
        multicols_ket_tensors.push_back(std::vector<TensorT>{peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(1)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(0))});
        multicols_ket_tensors.push_back(std::vector<TensorT>{peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(2))});
        multicols_ket_tensors.push_back(std::vector<TensorT>{peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(1)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(0))});
        multicols_ket_tensors.push_back(std::vector<TensorT>{peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(2))});

        //init indices
        multicols_linds.push_back(std::vector<IndexT>{commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(1)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(Lx-1,1)(2))),commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(0)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(Lx-1,0)(2)))});
            multicols_linds.push_back(std::vector<IndexT>{commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(2)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(1))),commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(2)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(0)))});
            multicols_linds.push_back(std::vector<IndexT>{commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(1)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,1)(2))),commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(0)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(2)))});
            multicols_linds.push_back(std::vector<IndexT>{commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(2)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(1))),commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(2)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(0)))});

            multicols_rinds.push_back(std::vector<IndexT>{commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(1)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(2))),commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(0)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(2)))});
        multicols_rinds.push_back(std::vector<IndexT>{commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(2)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(0))),commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(2)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,Ly-1)(1)))});
        multicols_rinds.push_back(std::vector<IndexT>{commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(1)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(2))),commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(0)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(2)))});
        multicols_rinds.push_back(std::vector<IndexT>{commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(2)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(2,0)(0))),commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(2)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(2,Ly-1)(1)))});

        multicols_udinds.push_back(std::vector<IndexT>{commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(1)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,1)(0))),commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(0)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(1))),commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,Ly-1)(1)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(0)))});
        multicols_udinds.push_back(std::vector<IndexT>{});
        multicols_udinds.push_back(std::vector<IndexT>{commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(1)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,1)(0))),commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(0)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(1,0)(1))),commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,Ly-1)(1)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(0)))});
        multicols_udinds.push_back(std::vector<IndexT>{});

        //init impos
        lcols_dl_impos_.clear();
        lcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_one",multicols_ket_tensors[0],multicols_linds[0],multicols_rinds[0],multicols_udinds[0]));
        lcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_two",multicols_ket_tensors[1],multicols_linds[1],multicols_rinds[1]));
        lcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_one",multicols_ket_tensors[2],multicols_linds[2],multicols_rinds[2],multicols_udinds[2]));
        lcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_two",multicols_ket_tensors[3],multicols_linds[3],multicols_rinds[3]));

        rcols_dl_impos_.clear();
        rcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_one",multicols_ket_tensors[0],multicols_rinds[0],multicols_linds[0],multicols_udinds[0]));
        rcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_two",multicols_ket_tensors[1],multicols_rinds[1],multicols_linds[1]));
        rcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_one",multicols_ket_tensors[2],multicols_rinds[2],multicols_linds[2],multicols_udinds[2]));
        rcols_dl_impos_.push_back(DL_iMPOt<TensorT>("type_two",multicols_ket_tensors[3],multicols_rinds[3],multicols_linds[3]));
    }

}
template 
void PEPSt_iTEBD<ITensor>::init_impo();
template
void PEPSt_iTEBD<IQTensor>::init_impo();


template <class TensorT>
const std::vector<TensorT> &PEPSt_iTEBD<TensorT>::env_tensors_from_itebd(int n_steps)
{
    //get left and right imps
    for (int stepi=0; stepi<n_steps; stepi++) update_imps_one_step();

    int Lx=peps_storage_._Lx,
        Ly=peps_storage_._Ly;

    //obtain env_tensors from left and right imps
    //kagome case
    //for kagome lattice, env tensors are
    //
    //      4
    //      |
    //   0--v--2
    //      |
    //   1--u--3
    //      |
    //      5
    //
    if (peps_storage_._tnetwork_type==8)
    {
        //get virt indices of cutting sites
        std::vector<IndexT> cutting_ket_inds;
        //left up
        cutting_ket_inds.push_back(commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(1)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(Lx-1,1)(2))));
        //right up
        cutting_ket_inds.push_back(commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(1)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(2))));
        //left down
        cutting_ket_inds.push_back(commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(0)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(Lx-1,0)(2))));
        //right down
        cutting_ket_inds.push_back(commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(0)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(2))));
        //up
        cutting_ket_inds.push_back(commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(1)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,1)(0))));
        //down
        cutting_ket_inds.push_back(commonIndex(peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,0)(0)),peps_storage_._tensor_list(peps_storage_._coor_to_siteind(0,Ly-1)(1))));

        //construct member variables of TensorT_Matrix Class
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
            contract_tensors.push_back(peps_storage_._tensor_list(sitei));
            contract_tensors.push_back(dag(peps_storage_._tensor_list(sitei)).prime(Link));
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

        //up and down indices
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

        //obtain up and down dominant eigenvector
        std::vector<IndexT> up_inds_dag=dag<IndexT>(up_inds),
                            down_inds_dag=dag<IndexT>(down_inds);
        TensorT VU(up_inds_dag),
                VD(down_inds_dag);
        randTensor(VU);
        randTensor(VD);
        TensorT_Matrix_Arnoldi<TensorT> UMat(up_inds,down_inds,contract_tensors,VU_contract_seq),
                                        DMat(down_inds,up_inds,contract_tensors,VD_contract_seq);
        Complex eta_U=arnoldi(UMat,VU,itebd_opts_),
                eta_D=arnoldi(DMat,VD,itebd_opts_);
        Print(eta_U);
        Print(eta_D);

        //get env_tensors_
        env_tensors_.clear();
        env_tensors_.push_back(contract_tensors[0]);
        env_tensors_.push_back(contract_tensors[4]);
        env_tensors_.push_back(contract_tensors[3]);
        env_tensors_.push_back(contract_tensors[7]);
        env_tensors_.push_back(VU);
        env_tensors_.push_back(VD);
    }

    return env_tensors_;
}
template
const std::vector<ITensor> &PEPSt_iTEBD<ITensor>::env_tensors_from_itebd(int n_steps);
template
const std::vector<IQTensor> &PEPSt_iTEBD<IQTensor>::env_tensors_from_itebd(int n_steps);


template <class TensorT>
void PEPSt_iTEBD<TensorT>::update_imps_one_step()
{
    //kagome case
    if (peps_storage_._tnetwork_type==8)
    {
        //if default construct, then we init the imps using first col
        if (!ldl_imps_.valid())
        {
            ldl_imps_=DL_iMPSt<TensorT>(lcols_dl_impos_[0].ket_tensors(),lcols_dl_impos_[0].ket_outgoing_inds(),lcols_dl_impos_[0].ket_virt_inds());
        }
        else
        {
            contract_dl_impo_imps(ldl_imps_,lcols_dl_impos_[0],itebd_opts_);
        }
        if (!rdl_imps_.valid())
        {
            rdl_imps_=DL_iMPSt<TensorT>(rcols_dl_impos_[0].ket_tensors(),rcols_dl_impos_[0].ket_outgoing_inds(),rcols_dl_impos_[0].ket_virt_inds());
        }
        else
        {
            contract_dl_impo_imps(rdl_imps_,rcols_dl_impos_[0],itebd_opts_);
        }
        rdl_imps_.move_tensors();

        //contract for second col
        contract_dl_impo_imps(ldl_imps_,lcols_dl_impos_[1],itebd_opts_);
        ldl_imps_.move_tensors();
        contract_dl_impo_imps(rdl_imps_,rcols_dl_impos_[1],itebd_opts_);

        //contract for third col
        contract_dl_impo_imps(ldl_imps_,lcols_dl_impos_[2],itebd_opts_);
        contract_dl_impo_imps(rdl_imps_,rcols_dl_impos_[2],itebd_opts_);
        rdl_imps_.move_tensors();

        //contract for fourth col
        contract_dl_impo_imps(ldl_imps_,lcols_dl_impos_[3],itebd_opts_);
        ldl_imps_.move_tensors();
        contract_dl_impo_imps(rdl_imps_,rcols_dl_impos_[3],itebd_opts_);
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
        TensorT result_tensor=env_tensors_[4]*env_tensors_[0]*peps_storage_._tensor_list(1)*two_sites_mpo[0]*dag(prime(peps_storage_._tensor_list(1)))*env_tensors_[2]*env_tensors_[1]*peps_storage_._tensor_list(0)*two_sites_mpo[1]*dag(prime(peps_storage_._tensor_list(0)))*env_tensors_[3]*env_tensors_[5];
        return result_tensor.toComplex();
    }
}
template
Complex PEPSt_iTEBD<ITensor>::expect_val_from_env_tensors(const std::vector<ITensor> &two_sites_mpo) const;
template
Complex PEPSt_iTEBD<IQTensor>::expect_val_from_env_tensors(const std::vector<IQTensor> &two_sites_mpo) const;
