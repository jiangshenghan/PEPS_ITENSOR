
#include "imps.h"

//
//class iMPSt
//
template <class TensorT>
iMPSt<TensorT>::iMPSt(const std::vector<TensorT> &site_tensors, const IndexT &lind, const IndexT &rind):
    n_sites_uc_(site_tensors.size()),
    site_tensors_(site_tensors)
{
    //init site_inds_
    for (int sitei=0; sitei<n_sites_uc_; sitei++)
    {
        for (const auto &ind: site_tensors_[sitei].indices())
        {
            int ni=(sitei+1)%n_sites_uc_,
            pi=(sitei+n_sites_uc_-1)%n_sites_uc_;
            if (ind==commonIndex(site_tensors_[sitei],site_tensors_[ni]) || ind==commonIndex(site_tensors_[sitei],site_tensors_[pi])) continue;
            if (ind==lind || ind==rind) continue;
            site_inds_.push_back(ind);
            break;
        }
    }
    //init virt_inds_
    virt_inds_=std::vector<IndexT>(n_sites_uc_+1);
    virt_inds_[0]=lind;
    virt_inds_[n_sites_uc_]=dag(rind);
    for (int sitei=1; sitei<n_sites_uc_; sitei++)
    {
        virt_inds_[sitei]=commonIndex(site_tensors_[sitei],site_tensors_[sitei-1]);
    }
}
template 
iMPSt<ITensor>::iMPSt(const std::vector<ITensor> &site_tensors, const ITensor::IndexT &lind, const ITensor::IndexT &rind);
template
iMPSt<IQTensor>::iMPSt(const std::vector<IQTensor> &site_tensors, const IQTensor::IndexT &lind, const IQTensor::IndexT &rind);


template <class TensorT>
iMPSt<TensorT>::iMPSt(const std::vector<TensorT> &site_tensors, const std::vector<IndexT> &site_inds):
    n_sites_uc_(site_tensors.size()),
    site_inds_(site_inds),
    site_tensors_(site_tensors)
{
    //init virt_inds_
    virt_inds_=std::vector<IndexT>(n_sites_uc_+1);
    //init bulk_virt_inds_
    for (int sitei=1; sitei<n_sites_uc_; sitei++)
    {
        virt_inds_[sitei]=commonIndex(site_tensors_[sitei],site_tensors_[sitei-1]);
    }
    //init boundary virt_inds_
    for (const auto &ind: site_tensors_[0].indices())
    {
        if (ind==site_inds_[0]) continue;
        if (ind==dag(virt_inds_[1])) continue;
        virt_inds_[0]=ind;
    }
    for (const auto &ind: site_tensors_.back().indices())
    {
        if (ind==site_inds_.back()) continue;
        if (ind==virt_inds_[n_sites_uc_-1]) continue;
        virt_inds_.back()=dag(ind);
    }
}
template
iMPSt<ITensor>::iMPSt(const std::vector<ITensor> &site_tensors, const std::vector<ITensor::IndexT> &site_inds);
template
iMPSt<IQTensor>::iMPSt(const std::vector<IQTensor> &site_tensors, const std::vector<IQTensor::IndexT> &site_inds);


template <class TensorT>
void iMPSt<TensorT>::move_tensors(int movei)
{
    std::rotate(site_inds_.begin(),site_inds_.end()-movei,site_inds_.end());
    //pay attention to boundary virt legs
    std::vector<IndexT> new_virt_inds;
    for (int indi=0; indi<=n_sites_uc_; indi++)
    {
        new_virt_inds.push_back(isomorphic_legs(virt_inds_[(n_sites_uc_-movei+indi)%n_sites_uc_],nameint("virt_ind_",indi)));
    }

    for (int sitei=0; sitei<n_sites_uc_; sitei++)
    {
        site_tensors_[sitei].replaceIndex(virt_inds_[sitei],new_virt_inds[(sitei+movei)%n_sites_uc_]);
        if ((sitei+1+movei)==n_sites_uc_)
        {
            site_tensors_[sitei].replaceIndex(dag(virt_inds_[sitei+1]),dag(new_virt_inds[n_sites_uc_]));
        }
        else
        {
            site_tensors_[sitei].replaceIndex(dag(virt_inds_[sitei+1]),dag(new_virt_inds[(sitei+1+movei)%n_sites_uc_]));
        }
    }
    std::rotate(site_tensors_.begin(),site_tensors_.end()-movei,site_tensors_.end());
    virt_inds_=new_virt_inds;
}



//
//class DL_iMPSt
//
template <class TensorT>
DL_iMPSt<TensorT>::DL_iMPSt(const std::vector<IndexT> &ket_siteinds, const std::vector<IndexT> &bra_siteinds, std::vector<TensorT> site_tensors):
    ket_siteinds_(ket_siteinds), 
    bra_siteinds_(bra_siteinds)
{
    //cout << "construct dl_imps!" << endl;
    std::vector<IndexT> site_inds;
    for (int sitei=0; sitei<site_tensors.size(); sitei++)
    {
        //avoid the case where ketind and braind the same up to dag and prime
        if (ket_siteinds_[sitei]==(dag(bra_siteinds_[sitei]).noprime())) 
        {
            //cout << "replace ket and bra siteinds to avoid contraction" << endl;
            IndexT temp_ketind, temp_braind;
            temp_ketind=isomorphic_legs(ket_siteinds_[sitei],nameint("ket_siteind_",sitei));
            temp_braind=isomorphic_legs(bra_siteinds_[sitei],nameint("bra_siteind_",sitei));
            site_tensors[sitei].replaceIndex(ket_siteinds_[sitei],temp_ketind);
            site_tensors[sitei].replaceIndex(bra_siteinds_[sitei],temp_braind);
            ket_siteinds_[sitei]=temp_ketind;
            bra_siteinds_[sitei]=temp_braind;
        }
        CombinerT temp_combiner(ket_siteinds_[sitei],bra_siteinds_[sitei]);
        temp_combiner.init(nameint("siteind_",sitei));
        siteind_combiners_.push_back(temp_combiner);
        site_inds.push_back(siteind_combiners_[sitei].right());
        site_tensors[sitei]=site_tensors[sitei]*siteind_combiners_[sitei];
    }
    imps_=iMPSt<TensorT>(site_tensors,site_inds);
}
template
DL_iMPSt<ITensor>::DL_iMPSt(const std::vector<ITensor::IndexT> &ket_siteinds, const std::vector<ITensor::IndexT> &bra_siteinds, std::vector<ITensor> site_tensors);
template
DL_iMPSt<IQTensor>::DL_iMPSt(const std::vector<IQTensor::IndexT> &ket_siteinds, const std::vector<IQTensor::IndexT> &bra_siteinds, std::vector<IQTensor> site_tensors);


//template <class TensorT>
//DL_iMPSt<TensorT>::DL_iMPSt(std::vector<TensorT> site_tensors, const IndexT &lind, const IndexT &rind, const std::vector<IndexT> &ket_siteinds, std::vector<IndexT> &bra_siteinds): 
//    ket_siteinds_(ket_siteinds), 
//    bra_siteinds_(bra_siteinds)
//{
//    std::vector<IndexT> site_inds;
//    for (int sitei=0; sitei<site_tensors.size(); sitei++)
//    {
//        //avoid the case where ketind and braind the same up to dag and prime
//        if (ket_siteinds_[sitei]==(dag(bra_siteinds_[sitei]).noprime())) 
//        {
//            IndexT temp_ketind, temp_braind;
//            temp_ketind=isomorphic_legs(ket_siteinds_[sitei],nameint("ket_leg_",sitei));
//            temp_braind=isomorphic_legs(bra_siteinds_[sitei],nameint("bra_leg_",sitei));
//            site_tensors[sitei].replaceIndex(ket_siteinds_[sitei],temp_ketind);
//            site_tensors[sitei].replaceIndex(bra_siteinds_[sitei],temp_braind);
//            ket_siteinds_[sitei]=temp_ketind;
//            bra_siteinds_[sitei]=temp_braind;
//        }
//        siteind_combiners_.push_back(CombinerT(ket_siteinds_[sitei],bra_siteinds_[sitei]));
//        site_inds.push_back(siteind_combiners_[sitei].right());
//        site_tensors[sitei]=site_tensors[sitei]*siteind_combiners_[sitei];
//    }
//    imps_=iMPSt<TensorT>(site_tensors,lind,rind);
//}
//template
//DL_iMPSt<ITensor>::DL_iMPSt(std::vector<ITensor> site_tensors, const ITensor::IndexT &lind, const ITensor::IndexT &rind, const std::vector<ITensor::IndexT> &ket_siteinds, std::vector<ITensor::IndexT> &bra_siteinds);
//template
//DL_iMPSt<IQTensor>::DL_iMPSt(std::vector<IQTensor> site_tensors, const IQTensor::IndexT &lind, const IQTensor::IndexT &rind, const std::vector<IQTensor::IndexT> &ket_siteinds, std::vector<IQTensor::IndexT> &bra_siteinds);

template <class TensorT>
DL_iMPSt<TensorT>::DL_iMPSt(std::vector<TensorT> ket_site_tensors, const std::vector<IndexT> &ket_siteinds, const std::vector<IndexT> &sl_virt_inds)
{
    //Print(ket_site_tensors);
    //Print(ket_siteinds);
    //Print(sl_virt_inds);

    //init ket and bra site_tensors
    std::vector<TensorT> bra_site_tensors;
    for (int sitei=0; sitei<ket_site_tensors.size(); sitei++)
    {
        bra_siteinds_.push_back(dag(isomorphic_legs(ket_siteinds[sitei],nameint("bra_siteind_",sitei))));
        bra_site_tensors.push_back(dag(ket_site_tensors[sitei]));
        bra_site_tensors[sitei].prime(dag(sl_virt_inds[sitei]));
        bra_site_tensors[sitei].prime(sl_virt_inds[sitei+1]);
        bra_site_tensors[sitei].replaceIndex(dag(ket_siteinds[sitei]),bra_siteinds_[sitei]);
        //modidy ket_site_leg
        ket_siteinds_.push_back(isomorphic_legs(ket_siteinds[sitei],nameint("ket_siteind_",sitei)));
        ket_site_tensors[sitei].replaceIndex(ket_siteinds[sitei],ket_siteinds_[sitei]);
        //init siteind_combiners_
        CombinerT temp_combiner(ket_siteinds_[sitei],bra_siteinds_[sitei]);
        temp_combiner.init(nameint("siteind_",sitei));
        siteind_combiners_.push_back(temp_combiner);
    }
    //Print(ket_site_tensors);
    //Print(bra_site_tensors);
    //Print(ket_siteinds_);
    //Print(bra_siteinds_);
    //Print(siteind_combiners_);

    //obtain dl_site_tensors
    std::vector<CombinerT> virt_leg_combiners;
    for (int legi=0; legi<sl_virt_inds.size(); legi++)
    {
        CombinerT temp_combiner(sl_virt_inds[legi],dag(sl_virt_inds[legi]).prime());
        temp_combiner.init(nameint("virt_ind_",legi));
        virt_leg_combiners.push_back(temp_combiner);
    }
    //Print(virt_leg_combiners);
    std::vector<IndexT> dl_siteinds;
    std::vector<TensorT> dl_site_tensors;
    for (int sitei=0; sitei<ket_site_tensors.size(); sitei++)
    {
        dl_siteinds.push_back(siteind_combiners_[sitei].right());
        dl_site_tensors.push_back(ket_site_tensors[sitei]*bra_site_tensors[sitei]*virt_leg_combiners[sitei]*dag(virt_leg_combiners[sitei+1])*siteind_combiners_[sitei]);
    }
    //Print(dl_site_tensors);
    //Print(dl_siteinds);

    imps_=iMPSt<TensorT>(dl_site_tensors,dl_siteinds);
}
template
DL_iMPSt<ITensor>::DL_iMPSt(std::vector<ITensor> ket_site_tensors, const std::vector<ITensor::IndexT> &ket_siteinds, const std::vector<ITensor::IndexT> &sl_virt_inds);
template
DL_iMPSt<IQTensor>::DL_iMPSt(std::vector<IQTensor> ket_site_tensors, const std::vector<IQTensor::IndexT> &ket_siteinds, const std::vector<IQTensor::IndexT> &sl_virt_inds);


template <class TensorT>
void DL_iMPSt<TensorT>::move_tensors(int movei)
{
    std::rotate(ket_siteinds_.begin(),ket_siteinds_.end()-movei,ket_siteinds_.end());
    std::rotate(bra_siteinds_.begin(),bra_siteinds_.end()-movei,bra_siteinds_.end());
    std::rotate(siteind_combiners_.begin(),siteind_combiners_.end()-movei,siteind_combiners_.end());
    imps_.move_tensors(movei);
}
template
void DL_iMPSt<ITensor>::move_tensors(int movei);
template
void DL_iMPSt<IQTensor>::move_tensors(int movei);



//
//class DL_iMPOt
//
template <class TensorT>
DL_iMPOt<TensorT>::DL_iMPOt(std::string type_name, const std::vector<TensorT> &init_ket_tensors, const std::vector<IndexT> &init_ket_incoming_inds, const std::vector<IndexT> &init_ket_outgoing_inds, std::vector<IndexT> init_ket_virt_inds):
    type_name_(type_name),
    n_tensors_uc_(init_ket_tensors.size())
{
    //reconstruct tensors for different types separately
    if (type_name_.find("type_one")!=std::string::npos)
    {
        //get the incoming and outgoing inds
        for (int tensori=0; tensori<n_tensors_uc_; tensori++)
        {
            ket_incoming_inds_.push_back(isomorphic_legs(init_ket_incoming_inds[tensori],nameint("ket_incoming_ind_",tensori)));
            ket_outgoing_inds_.push_back(isomorphic_legs(init_ket_outgoing_inds[tensori],nameint("ket_outgoing_ind_",tensori)));
            bra_incoming_inds_.push_back(isomorphic_legs(dag(init_ket_incoming_inds[tensori]),nameint("bra_incoming_ind_",tensori)));
            bra_outgoing_inds_.push_back(isomorphic_legs(dag(init_ket_outgoing_inds[tensori]),nameint("bra_outgoing_ind_",tensori)));
        }
        //obtain virt indices
        for (int indi=0; indi<=n_tensors_uc_; indi++)
        {
            ket_virt_inds_.push_back(isomorphic_legs(init_ket_virt_inds[indi],nameint("ket_virt_ind_",indi)));
            bra_virt_inds_.push_back(isomorphic_legs(dag(init_ket_virt_inds[indi]),nameint("bra_virt_ind_",indi)));
        }

        //construct ket and bra tensor
        for (int tensori=0; tensori<n_tensors_uc_; tensori++)
        {
            TensorT temp_ket_tensor=init_ket_tensors[tensori], 
                    temp_bra_tensor=dag(temp_ket_tensor);

            temp_ket_tensor.replaceIndex(init_ket_incoming_inds[tensori],ket_incoming_inds_[tensori]);
            temp_ket_tensor.replaceIndex(init_ket_outgoing_inds[tensori],ket_outgoing_inds_[tensori]);
            temp_ket_tensor.replaceIndex(init_ket_virt_inds[tensori],ket_virt_inds_[tensori]);
            temp_ket_tensor.replaceIndex(dag(init_ket_virt_inds[tensori+1]),dag(ket_virt_inds_[tensori+1]));
            ket_tensors_.push_back(temp_ket_tensor);
            //Print(temp_ket_tensor);

            temp_bra_tensor.replaceIndex(dag(init_ket_incoming_inds[tensori]),bra_incoming_inds_[tensori]);
            temp_bra_tensor.replaceIndex(dag(init_ket_outgoing_inds[tensori]),bra_outgoing_inds_[tensori]);
            temp_bra_tensor.replaceIndex(dag(init_ket_virt_inds[tensori]),bra_virt_inds_[tensori]);
            temp_bra_tensor.replaceIndex(init_ket_virt_inds[tensori+1],dag(bra_virt_inds_[tensori+1]));
            bra_tensors_.push_back(temp_bra_tensor);
            //Print(temp_bra_tensor);
        }
    }

    if (type_name_.find("type_two")!=std::string::npos)
    {
        //get the incoming and outgoing inds
        for (int tensori=0; tensori<2*n_tensors_uc_; tensori++)
        {
            ket_incoming_inds_.push_back(isomorphic_legs(init_ket_incoming_inds[tensori],nameint("ket_incoming_ind_",tensori)));
            ket_outgoing_inds_.push_back(isomorphic_legs(init_ket_outgoing_inds[tensori],nameint("ket_outgoing_ind_",tensori)));
            bra_incoming_inds_.push_back(isomorphic_legs(dag(init_ket_incoming_inds[tensori]),nameint("bra_incoming_ind_",tensori)));
            bra_outgoing_inds_.push_back(isomorphic_legs(dag(init_ket_outgoing_inds[tensori]),nameint("bra_outgoing_ind_",tensori)));
        }
        for (int tensori=0; tensori<n_tensors_uc_; tensori++)
        {
            TensorT temp_ket_tensor=init_ket_tensors[tensori],
                    temp_bra_tensor=dag(temp_ket_tensor);

            temp_ket_tensor.replaceIndex(init_ket_incoming_inds[2*tensori],ket_incoming_inds_[2*tensori]);
            temp_ket_tensor.replaceIndex(init_ket_incoming_inds[2*tensori+1],ket_incoming_inds_[2*tensori+1]);
            temp_ket_tensor.replaceIndex(init_ket_outgoing_inds[2*tensori],ket_outgoing_inds_[2*tensori]);
            temp_ket_tensor.replaceIndex(init_ket_outgoing_inds[2*tensori+1],ket_outgoing_inds_[2*tensori+1]);
            ket_tensors_.push_back(temp_ket_tensor);
            //Print(temp_ket_tensor);

            temp_bra_tensor.replaceIndex(dag(init_ket_incoming_inds[2*tensori]),bra_incoming_inds_[2*tensori]);
            temp_bra_tensor.replaceIndex(dag(init_ket_incoming_inds[2*tensori+1]),bra_incoming_inds_[2*tensori+1]);
            temp_bra_tensor.replaceIndex(dag(init_ket_outgoing_inds[2*tensori]),bra_outgoing_inds_[2*tensori]);
            temp_bra_tensor.replaceIndex(dag(init_ket_outgoing_inds[2*tensori+1]),bra_outgoing_inds_[2*tensori+1]);
            bra_tensors_.push_back(temp_bra_tensor);
            //Print(temp_bra_tensor);
        }
    }

    //TODO: implement other two types
}
template
DL_iMPOt<ITensor>::DL_iMPOt(std::string type_name, const std::vector<ITensor> &init_ket_tensors, const std::vector<ITensor::IndexT> &init_ket_incoming_inds, const std::vector<ITensor::IndexT> &init_ket_outgoing_inds, std::vector<ITensor::IndexT> init_ket_virt_inds);
template
DL_iMPOt<IQTensor>::DL_iMPOt(std::string type_name, const std::vector<IQTensor> &init_ket_tensors, const std::vector<IQTensor::IndexT> &init_ket_incoming_inds, const std::vector<IQTensor::IndexT> &init_ket_outgoing_inds, std::vector<IQTensor::IndexT> init_ket_virt_inds);


//Methods to update imps
template <class TensorT>
void contract_dl_impo_imps(DL_iMPSt<TensorT> &dl_imps, const DL_iMPOt<TensorT> &dl_impo, Args contract_opts)
{
    using IndexT=typename TensorT::IndexT;
    using IndexValT=typename TensorT::IndexValT;
    using CombinerT=typename TensorT::CombinerT;

    //we should do the contraction/compression for different types separately
    //left and right dominant eigenvalue, should be real and equal if no degenerate
    Complex eta_L, eta_R;
    TensorT VR, VL;
    //tensors_uc are for imps truncation
    std::vector<TensorT> tensors_uc;
    TensorT_Matrix_Arnoldi<TensorT> LMat, RMat;

    Print(dl_impo.type_name());
    //Print(dl_impo);
    //Print(dl_imps);
    
    if (dl_impo.type_name().find("type_one")!=std::string::npos)
    {
        //init contract_tensors and tensors_uc
        std::vector<TensorT> contract_tensors;
        tensors_uc.clear();
        int n_tensors_uc=dl_impo.n_tensors_uc();
        for (int sitei=0; sitei<n_tensors_uc; sitei++)
        {
            //decombine the siteind of imps tensor
            TensorT imps_tensor=dl_imps.dl_site_tensors(sitei);
            //replace the incoming inds of impo tensor with siteind of imps tensor for contraction
            TensorT impo_ket_tensor=dl_impo.ket_tensors(sitei),
                    impo_bra_tensor=dl_impo.bra_tensors(sitei);
            impo_ket_tensor.replaceIndex(dl_impo.ket_incoming_inds(sitei),dag(dl_imps.ket_siteinds(sitei)));
            impo_bra_tensor.replaceIndex(dl_impo.bra_incoming_inds(sitei),dag(dl_imps.bra_siteinds(sitei)));

            contract_tensors.push_back(imps_tensor);
            contract_tensors.push_back(impo_ket_tensor);
            contract_tensors.push_back(impo_bra_tensor);
            contract_tensors.push_back(dag(impo_ket_tensor).prime().noprime(prime(dl_impo.ket_outgoing_inds(sitei))));
            contract_tensors.push_back(dag(impo_bra_tensor).prime().noprime(prime(dl_impo.bra_outgoing_inds(sitei))));
            contract_tensors.push_back(dag(imps_tensor).prime());

            tensors_uc.push_back(imps_tensor);
            tensors_uc.push_back(impo_ket_tensor);
            tensors_uc.push_back(impo_bra_tensor);

        }

        //contract_seq
        std::vector<int> VL_contract_seq, VR_contract_seq;
        VL_contract_seq.push_back(-1);
        VR_contract_seq.push_back(-1);
        for (int tensori=0; tensori<contract_tensors.size(); tensori++) 
        {
            VL_contract_seq.push_back(tensori);
            VR_contract_seq.push_back(contract_tensors.size()-1-tensori);
        }

        //leg_combiners, combining virt_inds of both imps and impo
        std::vector<CombinerT> leg_combiners;
        //we will not use 0th leg_combiner
        leg_combiners.push_back(CombinerT(dl_imps.virt_inds(0)));
        //combiners of virt_inds of both imps and impo
        Print(n_tensors_uc);
        for (int indi=0; indi<=n_tensors_uc; indi++)
        {
            CombinerT temp_combiner(dl_imps.virt_inds(indi),dl_impo.ket_virt_inds(indi),dl_impo.bra_virt_inds(indi));
            temp_combiner.init();
            leg_combiners.push_back(temp_combiner);
            leg_combiners.push_back(prime(dag(temp_combiner)));
        }
        
        //left and right indices 
        //left_inds are input(output) inds for VL(VR), and right_inds are output(input) inds for VL(VR)
        std::vector<IndexT> left_inds, right_inds;
        for (const auto &ind: leg_combiners[1].left()) left_inds.push_back(dag(ind));
        for (const auto &ind: leg_combiners[2].left()) left_inds.push_back(dag(ind));
        for (const auto &ind: leg_combiners[2*n_tensors_uc+1].left()) right_inds.push_back(ind);
        for (const auto &ind: leg_combiners.back().left()) right_inds.push_back(ind);
        //Print(left_inds);
        //Print(right_inds);

        //combiner_seq
        std::vector<std::vector<int>> VL_combiner_seq, VR_combiner_seq;
        //left combiner seq
        VL_combiner_seq.push_back(std::vector<int>{});
        for (int tensori=0; tensori<n_tensors_uc; tensori++)
        {
            VL_combiner_seq.back().push_back(-(2*tensori+2));
            VL_combiner_seq.push_back(std::vector<int>{});
            VL_combiner_seq.push_back(std::vector<int>{});
            VL_combiner_seq.push_back(std::vector<int>{-(2*(tensori+1)+1),2*tensori+2});
            VL_combiner_seq.push_back(std::vector<int>{});
            VL_combiner_seq.push_back(std::vector<int>{});
            VL_combiner_seq.push_back(std::vector<int>{2*(tensori+1)+1});
        }
        //right combiner seq
        VR_combiner_seq.push_back(std::vector<int>{});
        for (int tensori=n_tensors_uc-1; tensori>=0; tensori--)
        {
            VR_combiner_seq.back().push_back(2*(tensori+1)+1);
            VR_combiner_seq.push_back(std::vector<int>{});
            VR_combiner_seq.push_back(std::vector<int>{});
            VR_combiner_seq.push_back(std::vector<int>{2*tensori+2,-(2*(tensori+1)+1)});
            VR_combiner_seq.push_back(std::vector<int>{});
            VR_combiner_seq.push_back(std::vector<int>{});
            VR_combiner_seq.push_back(std::vector<int>{-(2*tensori+2)});
        }

        //Print(left_inds);
        //Print(right_inds);
        //Print(VL_contract_seq);
        //Print(VR_contract_seq);
        //Print(leg_combiners);
        //Print(VL_combiner_seq);
        //Print(VR_combiner_seq);

        //implement arnoldi method
        LMat=TensorT_Matrix_Arnoldi<TensorT>(left_inds,right_inds,contract_tensors,VL_contract_seq,leg_combiners,VL_combiner_seq);
        RMat=TensorT_Matrix_Arnoldi<TensorT>(right_inds,left_inds,contract_tensors,VR_contract_seq,leg_combiners,VR_combiner_seq);

        //TODO: missing the case where the dominant eigenvector has sz!=0
        //using identity matrix to init VL and VR
        std::vector<IndexT> VL_in_inds, VL_out_inds,
                            VR_in_inds, VR_out_inds;
        for (const auto &ind: leg_combiners[1].left()) VL_in_inds.push_back(ind);
        for (const auto &ind: leg_combiners[2].left()) VL_out_inds.push_back(ind);
        for (const auto &ind: leg_combiners[2*n_tensors_uc+1].left()) VR_in_inds.push_back(dag(ind));
        for (const auto &ind: leg_combiners.back().left()) VR_out_inds.push_back(dag(ind));
        VL=delta_tensor<TensorT>(VL_in_inds,VL_out_inds);
        VR=delta_tensor<TensorT>(VR_in_inds,VR_out_inds);

        //using random matrix to init VL and VR
        //std::vector<IndexT> left_inds_dag=dag<IndexT>(left_inds),
        //                    right_inds_dag=dag<IndexT>(right_inds);
        //VL=TensorT(left_inds_dag);
        //VR=TensorT(right_inds_dag);
        //randTensor(VL);
        //randTensor(VR);

    }

    //FIXME: debug this case
    if (dl_impo.type_name().find("type_two")!=std::string::npos)
    {
        //init contract_tensors and tensors_uc
        std::vector<TensorT> contract_tensors;
        tensors_uc.clear();
        int impo_n_tensors_uc=dl_impo.n_tensors_uc(),
            imps_n_tensors_uc=2*impo_n_tensors_uc;
        for (int tensori=0; tensori<impo_n_tensors_uc; tensori++)
        {
            //contract two imps tensor and decombine the siteind
            TensorT imps_tensor=dl_imps.dl_site_tensors(tensori*2)*dl_imps.dl_site_tensors(tensori*2+1);
            //replace the incoming inds of imps tensor with siteind of imps tensor for contraction
            TensorT impo_ket_tensor=dl_impo.ket_tensors(tensori),
                    impo_bra_tensor=dl_impo.bra_tensors(tensori);
            impo_ket_tensor.replaceIndex(dl_impo.ket_incoming_inds(tensori*2),dag(dl_imps.ket_siteinds(tensori*2)));
            impo_ket_tensor.replaceIndex(dl_impo.ket_incoming_inds(tensori*2+1),dag(dl_imps.ket_siteinds(tensori*2+1)));
            impo_bra_tensor.replaceIndex(dl_impo.bra_incoming_inds(tensori*2),dag(dl_imps.bra_siteinds(tensori*2)));
            impo_bra_tensor.replaceIndex(dl_impo.bra_incoming_inds(tensori*2+1),dag(dl_imps.bra_siteinds(tensori*2+1)));

            contract_tensors.push_back(imps_tensor);
            contract_tensors.push_back(impo_ket_tensor);
            contract_tensors.push_back(impo_bra_tensor);
            contract_tensors.push_back(dag(impo_ket_tensor).prime(dl_imps.ket_siteinds(tensori*2)).prime(dl_imps.ket_siteinds(tensori*2+1)).prime(Site));
            contract_tensors.push_back(dag(impo_bra_tensor).prime(dl_imps.bra_siteinds(tensori*2)).prime(dl_imps.bra_siteinds(tensori*2+1)).prime(Site));
            contract_tensors.push_back(prime(dag(imps_tensor)));

            tensors_uc.push_back(imps_tensor);
            tensors_uc.push_back(impo_ket_tensor);
            tensors_uc.push_back(impo_bra_tensor);
        }

        //contract_seq
        std::vector<int> VL_contract_seq, VR_contract_seq;
        VL_contract_seq.push_back(-1);
        VR_contract_seq.push_back(-1);
        for (int tensori=0; tensori<contract_tensors.size(); tensori++)
        {
            VL_contract_seq.push_back(tensori);
            VR_contract_seq.push_back(contract_tensors.size()-1-tensori);
        }

        //left and right indices 
        //left_inds are input(output) inds for VL(VR), and right_inds are output(input) inds for VL(VR)
        std::vector<IndexT> left_inds, right_inds;
        left_inds.push_back(dl_imps.virt_inds().front());
        left_inds.push_back(dag(dl_imps.virt_inds().front()).prime());
        right_inds.push_back(dag(dl_imps.virt_inds().back()));
        right_inds.push_back(prime(dl_imps.virt_inds().back()));

        //implement arnoldi method
        LMat=TensorT_Matrix_Arnoldi<TensorT>(left_inds,right_inds,contract_tensors,VL_contract_seq);
        RMat=TensorT_Matrix_Arnoldi<TensorT>(right_inds,left_inds,contract_tensors,VR_contract_seq);

        //TODO: missing the case where the dominant eigenvector has sz!=0
        //using identity matrix to init VL and VR
        VL=delta_tensor<TensorT>(dag(dl_imps.virt_inds().front()),prime(dl_imps.virt_inds().front()));
        VR=delta_tensor<TensorT>(dl_imps.virt_inds().back(),dag(prime(dl_imps.virt_inds().back())));

        //using random matrix to init VL and VR
        //std::vector<IndexT> left_inds_dag=dag<IndexT>(left_inds),
        //                    right_inds_dag=dag<IndexT>(right_inds);
        //VL=TensorT(left_inds_dag);
        //VR=TensorT(right_inds_dag);
        //randTensor(VL);
        //randTensor(VR);

    }

    //TODO: implement other types

    bool check_herm=false;
    do
    {
        cout << "MaxIter: " << contract_opts.getInt("MaxIter") << endl;
        eta_L=arnoldi(LMat,VL,contract_opts);
        eta_R=arnoldi(RMat,VR,contract_opts);
        contract_opts.add("MaxIter",contract_opts.getInt("MaxIter")+10);
        TensorT VL_dag=dag(VL).prime().mapprime(2,0),
                VR_dag=dag(VR).prime().mapprime(2,0);
        if ((VL-VL_dag).norm()/VL.norm()<1e-8 && (VR-VR_dag).norm()/VR.norm()<1e-8) check_herm=true;

        Print(eta_L);
        Print(eta_R);
        Print((eta_L-eta_R)/eta_L);
        Print((VL-VL_dag).norm()/VL.norm());
        Print((VR-VR_dag).norm()/VR.norm());
    }
    while ((std::abs((eta_L-eta_R)/eta_L)>1e-8 || (!check_herm)) && (contract_opts.getInt("MaxIter")<30));

    clean(VL);
    clean(VR);
    //Print(eta_L);
    //Print(eta_R);

    //update imps by truncation
    dl_imps=dl_imps_from_truncation(tensors_uc,dl_impo.ket_outgoing_inds(),dl_impo.bra_outgoing_inds(),VL,VR,contract_opts);
    cout << "Finish contraction and compression!" << endl;
}
template
void contract_dl_impo_imps(DL_iMPSt<ITensor> &dl_imps, const DL_iMPOt<ITensor> &dl_impo, Args contract_opts);
template
void contract_dl_impo_imps(DL_iMPSt<IQTensor> &dl_imps, const DL_iMPOt<IQTensor> &dl_impo, Args contract_opts);


template <class TensorT>
DL_iMPSt<TensorT> dl_imps_from_truncation(const std::vector<TensorT> &tensors_uc, const std::vector<typename TensorT::IndexT> &ket_siteinds, const std::vector<typename TensorT::IndexT> &bra_siteinds, const TensorT &VL, const TensorT &VR, Args trunc_opts)
{
    using IndexT=typename TensorT::IndexT;
    using IndexValT=typename TensorT::IndexValT;
    using CombinerT=typename TensorT::CombinerT;

    //Print(tensors_uc);
    //Print(VL);
    //Print(VR);

    //set svd options
    trunc_opts.add("DoRelCutoff",true);
    trunc_opts.add("AbsoluteCutoff",false);

    //eigendecompose VL and VR
    TensorT X,Y;
    Print(eigen_factor(VR,X));
    Print(eigen_factor(VL,Y));

    //combine the primed legs of X/Y
    CombinerT X_primed_combiner, Y_primed_combiner;
    for (const auto &ind: X.indices())
    {
        if (ind.primeLevel()!=0) X_primed_combiner.addleft(ind);
    }
    for (const auto &ind: Y.indices())
    {
        if (ind.primeLevel()!=0) Y_primed_combiner.addleft(ind);
    }
    X=X*X_primed_combiner; clean(X);
    Y=Y*Y_primed_combiner; clean(Y);

    IndexT Uind=uniqueIndex(Y,VL,Link), Vind=uniqueIndex(X,VR,Link);
    cout << "check eigen-factor:" << endl;
    Print((X*dag(prime(X).noprime(prime(Vind)))-VR).norm());
    Print((Y*dag(prime(Y).noprime(prime(Uind)))-VL).norm());

    //singular value decompose of Y*X
    std::vector<IndexT> Ycontract_inds, Xcontract_inds;
    for (const auto &ind: Y.indices())
    {
        if (ind==Uind) continue;
        Ycontract_inds.push_back(ind);
    }
    for (const auto &ind: X.indices())
    {
        if (ind==Vind) continue;
        Xcontract_inds.push_back(ind);
    }
    TensorT YX=tensor_contraction<TensorT,IndexT>(Y,X,Ycontract_inds,Xcontract_inds);
    TensorT U(Uind),D,V;
    //Print(tensor_svd(YX,U,D,V,trunc_opts));
    Print(svd(YX,U,D,V,trunc_opts));
    D/=D.norm();
    //Print(Y.norm());
    //Print(X.norm());
    //Print(YX.norm());
    //Print(D.diag());

    //update union_tensor
    int site_indi=0;
    std::vector<CombinerT> siteind_combiners;
    TensorT union_tensor=dag(U)*Y;
    for (const auto &tens: tensors_uc) 
    {
        union_tensor*=tens;
        while (hasindex(union_tensor,ket_siteinds[site_indi]) && hasindex(union_tensor,bra_siteinds[site_indi]))
        {
            siteind_combiners.push_back(CombinerT(ket_siteinds[site_indi],bra_siteinds[site_indi]));
            int combiner_last=siteind_combiners.size()-1;
            if (combiner_last>0) siteind_combiners.back().addleft(siteind_combiners[combiner_last-1].right());
            union_tensor=union_tensor*siteind_combiners.back();
            //Print(siteind_combiners.back());
            //Print(union_tensor.indices());
            site_indi++;
        }
    }
    union_tensor*=X*dag(V);
    for (int combineri=siteind_combiners.size()-1; combineri>=0; combineri--)
    {
        union_tensor=union_tensor*dag(siteind_combiners[combineri]);
    }
    //very important to do clean?
    clean(union_tensor);
    Print(union_tensor.indices());

    //update site_tensors_ by svd. 
    //for example for n=2 case
    //
    // --union_tensor-- = --P--nD--Q--
    //       | |            |      |
    //
    //Then, we have
    //
    //  --tens0-- = --P--
    //      |         |
    //
    //  --tens1-- = --nD--Q--Dinv--
    //      |             |
    //
    std::vector<TensorT> site_tensors;
    //choose the way to absorb bond tensor (entanglement spectrum)
    std::string absorb_dir=trunc_opts.getString("AbsorbBond","left_bond");
    //absorb LEFT bond tensor to site tensor
    if (absorb_dir.find("left_bond")!=std::string::npos)
    {
        IndexT lind=commonIndex(union_tensor,U);
        for (int sitei=0; sitei<ket_siteinds.size()-1; sitei++)
        {
            //decompose using svd method
            TensorT P(ket_siteinds[sitei],bra_siteinds[sitei],lind),nD,Q;
            Print(P.indices());
            //Print(tensor_svd(union_tensor,P,nD,Q,trunc_opts));
            Print(svd(union_tensor,P,nD,Q,trunc_opts));
            nD/=nD.norm();
            //Print(union_tensor.norm());
            //Print(nD.diag());
            site_tensors.push_back(P);
            union_tensor=nD*Q;
            lind=commonIndex(union_tensor,nD);

            //decompose using denmatDecomp method
            //TensorT P(ket_siteinds[sitei],bra_siteinds[sitei],lind),nDQ;
            //Print(denmatDecomp(union_tensor,P,nDQ,Fromleft,trunc_opts));
            //site_tensors.push_back(P);
            //lind=commonIndex(nDQ,P);
            //union_tensor=nDQ;
        }
        TensorT Dinv=dag(D);
        Dinv.pseudoInvert();
        IndexT rind=commonIndex(Dinv,U);
        Dinv.replaceIndex(rind,isomorphic_legs(rind,"rind"));
        site_tensors.push_back(union_tensor*Dinv);
    }

    //absorb RIGHT bond tensor to site tensor
    if (absorb_dir.find("right_bond")!=std::string::npos)
    {
        IndexT rind=commonIndex(union_tensor,V);
        for (int sitei=ket_siteinds.size()-1; sitei>0; sitei--)
        {
            //decompose using svd method
            TensorT P,nD,Q(ket_siteinds[sitei],bra_siteinds[sitei],rind);
            //Print(tensor_svd(union_tensor,P,nD,Q,trunc_opts));
            Print(svd(union_tensor,P,nD,Q,trunc_opts));
            nD/=nD.norm();
            //Print(union_tensor.norm());
            //Print(nD.diag());
            site_tensors.insert(site_tensors.begin(),Q);
            union_tensor=P*nD;
            rind=commonIndex(union_tensor,nD);
        }
        TensorT Dinv=dag(D);
        Dinv.pseudoInvert();
        IndexT lind=commonIndex(Dinv,V);
        Dinv.replaceIndex(lind,isomorphic_legs(lind,"lind"));
        site_tensors.insert(site_tensors.begin(),Dinv*union_tensor);
    }

    for (auto &site_tensor: site_tensors)
    {
        int tensor_dim=1;
        for (const auto &ind: site_tensor.indices()) tensor_dim*=ind.m();
        site_tensor=site_tensor/site_tensor.norm()*std::sqrt(1.*tensor_dim);
        //Print(site_tensor.indices());
        //Print(site_tensor.norm());
    }

    DL_iMPSt<TensorT> dl_imps(ket_siteinds,bra_siteinds,site_tensors);
    //Print(dl_imps);
    return dl_imps;
}
template
DL_iMPSt<ITensor> dl_imps_from_truncation(const std::vector<ITensor> &tensors_uc, const std::vector<typename ITensor::IndexT> &ket_siteinds, const std::vector<typename ITensor::IndexT> &bra_siteinds, const ITensor &VL, const ITensor &VR, Args trunc_opts);
template
DL_iMPSt<IQTensor> dl_imps_from_truncation(const std::vector<IQTensor> &tensors_uc, const std::vector<typename IQTensor::IndexT> &ket_siteinds, const std::vector<typename IQTensor::IndexT> &bra_siteinds, const IQTensor &VL, const IQTensor &VR, Args trunc_opts);
