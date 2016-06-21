
#include "imps.h"

//
//class iMPSt
//
template <class TensorT>
iMPSt<TenosrT>::iMPSt(const std::vector<TensorT> &site_tensors, const IndexT &lind, const IndexT &rind):
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
    virt_inds=std::vector<IndexT>(n_sites_uc_+1);
    virt_inds[0]=lind;
    virt_inds[n_sites_uc_]=dag(rind);
    for (int sitei=1; sitei<n_sites_uc_; sitei++)
    {
        virt_inds[sitei]=commonIndex(site_tensors_[sitei],site_tensors_[sitei-1]);
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
    virt_inds=std::vector<IndexT>(n_sites_uc_+1);
    //init bulk_virt_inds_
    for (int sitei=1; sitei<n_sites_uc_; sitei++)
    {
        virt_inds[sitei]=commonIndex(site_tensors_[sitei],site_tensors_[sitei-1]);
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


//
//class DL_iMPSt
//
template <class TensorT>
DL_iMPSt<TensorT>::DL_iMPSt(std::vector<TensorT> site_tensors, const std::vector<IndexT> &ket_siteinds, std::vector<IndexT> &bra_siteinds):
    ket_siteinds_(ket_siteinds), 
    bra_siteinds_(bra_siteinds)
{
    std::vector<IndexT> site_inds;
    for (int sitei=0; sitei<site_tensors.size(); sitei++)
    {
        if (ket_siteinds_[sitei]==prime(dag(bottom_siteinds_sitei]))) 
        {
            cout << "Should not set ket_siteinds=prime(dag(bra_siteinds))!" << endl;
            exit(1);
        }
        siteind_combiners_.push_back(CombinerT(ket_siteinds_[sitei],bottom_siteinds_[sitei]));
        site_inds.push_back(siteind_combiners_[sitei].right());
        site_tensors[sitei]=site_tensors[sitei]*siteind_combiners_[sitei];
    }
    imps_=iMPSt<TensorT>(site_tensors,site_inds);
}
template
DL_iMPSt<ITensor>::DL_iMPSt(std::vector<ITensor> site_tensors, const std::vector<ITensor::IndexT> &ket_siteinds, std::vector<ITensor::IndexT> &bra_siteinds);
template
DL_iMPSt<IQTensor>::DL_iMPSt(std::vector<IQTensor> site_tensors, const std::vector<IQTensor::IndexT> &ket_siteinds, std::vector<IQTensor::IndexT> &bra_siteinds);


template <class TensorT>
DL_iMPSt<TensorT>::DL_iMPSt(std::vector<TensorT> site_tensors, const IndexT &lind, const IndexT &rind, const std::vector<IndexT> &ket_siteinds, std::vector<IndexT> &bra_siteinds): 
    ket_siteinds_(ket_siteinds), 
    bra_siteinds_(bra_siteinds)
{
    std::vector<IndexT> site_inds;
    for (int sitei=0; sitei<site_tensors.size(); sitei++)
    {
        if (ket_siteinds_[sitei]==prime(dag(bottom_siteinds_sitei]))) 
        {
            cout << "Should not set ket_siteinds=prime(dag(bra_siteinds))!" << endl;
            exit(1);
        }
        siteind_combiners_.push_back(CombinerT(ket_siteinds_[sitei],bottom_siteinds_[sitei]));
        site_inds.push_back(siteind_combiners_[sitei].right());
        site_tensors[sitei]=site_tensors[sitei]*siteind_combiners_[sitei];
    }
    imps_=iMPSt<TensorT>(site_tensors,lind,rind);
}
template
DL_iMPSt<ITensor>::DL_iMPSt(std::vector<ITensor> site_tensors, const ITensor::IndexT &lind, const ITensor::IndexT &rind, const std::vector<ITensor::IndexT> &ket_siteinds, std::vector<ITensor::IndexT> &bra_siteinds);
template
DL_iMPSt<IQTensor>::DL_iMPSt(std::vector<IQTensor> site_tensors, const IQTensor::IndexT &lind, const IQTensor::IndexT &rind, const std::vector<IQTensor::IndexT> &ket_siteinds, std::vector<IQTensor::IndexT> &bra_siteinds);


//
//class DL_iMPOt
//
template <class TensorT>
DL_iMPOt<TensorT>::DL_iMPOt(std::string type_name, const std::vector<IndexT> &init_ket_tensors, const std::vector<IndexT> &init_incoming_inds, const std::vector<IndexT> &init_outgoing_inds, const std::vector<IndexT> init_boundary_inds=std::vector<IndexT>()):
    type_name_(type_name),
    n_tensors_uc_(init_ket_tensors.size())
{
    //get the incoming and outgoing inds
    for (int tensori=0; tensori<n_tensors_uc_; tensori++)
    {
        ket_incoming_inds_.push_back(isomorphic_legs(init_incoming_inds[tensori],nameint("ket_incoming_ind_",tensori)));
        ket_outgoing_inds_.push_back(isomorphic_legs(init_outgoing_inds[tensori],nameint("ket_outgoing_ind_",tensori)));
        bra_incoming_inds_.push_back(isomorphic_legs(dag(init_incoming_inds[tensori]),nameint("bra_incoming_ind_",tensori)));
        bra_outgoing_inds_.push_back(isomorphic_legs(dag(init_outgoing_inds[tensori]),nameint("bra_outgoing_ind_",tensori)));
    }

    //reconstruct tensors for different types separately
    if (type_name_.find("type_one")!=std::string::npos)
    {
        //obtain virt indices
        std::vector<IndexT> init_virt_inds(n_tensors_uc_+1);

        init_virt_inds[0]=init_boundary_inds[0];
        init_virt_inds[n_tensors_uc_]=dag(init_boundary_inds[1]);
        for (int tensori=1; tensori<n_tensors_uc_; tensori++)
        {
            init_virt_inds[tensori]=commonIndex(init_ket_tensors[tensori],init_ket_tensors[tensori-1]);
        }
        for (int indi=0; indi<=n_tensors_uc_; indi++)
        {
            ket_virt_inds_.push_back(isomorphic_legs(init_virt_inds[indi],nameinit("ket_virt_ind_",indi)));
            bra_virt_inds_.push_back(isomorphic_legs(dag(init_virt_inds[indi]),nameint("bra_virt_ind_",indi)));
        }

        //construct ket and bra tensor
        for (int tensori=0; tensori<n_tensors_uc_; tensori++)
        {
            TensorT temp_ket_tensor=init_ket_tensors[tensori], 
                    temp_bra_tensor=dag(temp_ket_tensor);

            temp_ket_tensor.replaceIndex(init_incoming_inds[tensori],ket_incoming_inds_[tensori]);
            temp_ket_tensor.replaceIndex(init_outgoing_inds[tensori],ket_outgoing_inds_[tensori]);
            temp_ket_tensor.replaceIndex(init_virt_inds[tensori],ket_virt_inds_[tensori]);
            temp_ket_tensor.replaceIndex(dag(init_virt_inds[tensori+1]),dag(ket_virt_inds_[tensori+1])));
            ket_tensors_.push_back(temp_ket_tensor);

            temp_bra_tensor.replaceIndex(dag(init_incoming_inds[tensori],bra_incoming_inds_[tensori]));
            temp_bra_tensor.replaceIndex(dag(init_outgoing_inds[tensori],bra_outgoing_inds_[tensori]));
            temp_bra_tensor.replaceIndex(dag(init_virt_inds[tensori]),bra_virt_inds_[tensori]);
            temp_bra_tensor.replaceIndex(init_virt_inds[tensori+1],bra_virt_inds_[tensori+1]);
            bra_tensors_.push_back(temp_bra_tensor);
        }
    }

    if (type_name_.find("type_two")!=std::string::npos)
    {
        for (int tensori=0; tensori<n_tensors_uc_; tensori++)
        {
            TensorT temp_ket_tensor=init_ket_tensors[tensori],
                    temp_bra_tensor=dag(temp_ket_tensor);

            temp_ket_tensor.replaceIndex(init_incoming_inds[2*tensori],ket_incoming_inds_[2*tensori]);
            temp_ket_tensor.replaceIndex(init_incoming_inds[2*tensori+1],ket_incoming_inds_[2*tensori+1]);
            temp_ket_tensor.replaceIndex(init_outgoing_inds[2*tensori],ket_outgoing_inds_[2*tensori]);
            temp_ket_tensor.replaceIndex(init_outgoing_inds[2*tensori+1],ket_outgoing_inds_[2*tensori+1]);
            ket_tensors_.push_back(temp_ket_tensor);

            temp_bra_tensor.replaceIndex(dag(init_incoming_inds[2*tensori]),bra_incoming_inds_[2*tensori]);
            temp_bra_tensor.replaceIndex(dag(init_incoming_inds[2*tensori+1]),bra_incoming_inds_[2*tensori+1]);
            temp_bra_tensor.replaceIndex(dag(init_outgoing_inds[2*tensori]),bra_outgoing_inds_[2*tensori]);
            temp_bra_tensor.replaceIndex(dag(init_outgoing_inds[2*tensori+1]),bra_outgoing_inds_[2*tensori+1]);
            bra_tensors_.push_back(temp_bra_tensor);
        }
    }

    //TODO: implement other two types
}
template
DL_iMPOt<ITensor>::DL_iMPOt(std::string type_name, const std::vector<ITensor> &init_ket_tensors, const std::vector<ITensor::IndexT> &init_incoming_inds, const std::vector<ITensor::IndexT> &init_outgoing_inds, const std::vector<ITensor::IndexT> init_boundary_inds=std::vector<ITensor::IndexT>());
template
DL_iMPOt<IQTensor>::DL_iMPOt(std::string type_name, const std::vector<IQTensor> &init_ket_tensors, const std::vector<IQTensor::IndexT> &init_incoming_inds, const std::vector<IQTensor::IndexT> &init_outgoing_inds, const std::vector<IQTensor::IndexT> init_boundary_inds=std::vector<ITensor::IndexT>());


//Methods to update imps
template <TensorT>
void contract_dl_impo_imps(DL_iMPSt<TensorT> &dl_imps, DL_iMPOt<TensorT> &dl_impo, const Args &contract_opts)
{
    using IndexT=typename TensorT::IndexT;
    using IndexValT=typename TensorT::IndexValT;
    using CombinerT=typename TensorT::CombinerT;

    //we should do the contraction/compression for different types separately
    //left and right dominant eigenvalue, should be real and equal if no degenerate
    Complex eta_L, eta_R;
    TensorT VR, VL;

    //init arnoldi_args_
    Args arnoldi_args;
    arnoldi_args.add("MaxIter",10);
    arnoldi_args.add("MaxRestart",1);
    arnoldi_args.add("ErrGoal",1e-16);

    if (type_name_.find("type_one")!=std::string::npos)
    {
        //construct member variables of TensorT_Matrix Class
        //contract_tensors
        std::vector<TensorT> contract_tensors;
        int n_tensors_uc=dl_impo.n_tensors_uc();
        for (int sitei=0; sitei<n_tensors_uc; sitei++)
        {
            //decombine the siteind of imps tensor
            TensorT imps_tensor=dl_imps.site_tensors(sitei)*dag(dl_imps.siteind_combiners(sitei));
            //replace the incoming inds of impo tensor with siteind of imps tensor for contraction
            TensorT impo_ket_tensor=dl_impo.ket_tensors(sitei),
                    impo_bra_tensor=dl_impo.bra_tensors(sitei);
            impo_ket_tensor.replaceIndex(dl_impo.ket_incoming_inds(sitei),dag(dl_imps.ket_siteinds(sitei)));
            impo_bra_tensor.replaceIndex(dl_impo.bra_incoming_inds(sitei),dag(dl_imps.bra_siteinds(sitei)));

            contract_tensors.push_back(imps_tensor);
            contract_tensors.push_back(impo_ket_tensor);
            contract_tensors.push_back(impo_bra_tensor);
            contract_tensors.push_back(dag(impo_ket_tensor).prime().noprime(dl_impo.ket_outgoing_inds(sitei)));
            contract_tensors.push_back(dag(impo_bra_tensor).prime().noprime(dl_impo.bra_outgoing_inds(sitei)));
            contract_tensors.push_back(dag(imps_tensor).prime());

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
        //set 0th combiner to be empty
        leg_combiners.push_back(CombinerT());
        //combiners of virt_inds of both imps and impo
        for (int indi=0; indi<=n_tensors_uc; indi++)
        {
            CombinerT temp_combiner(dl_imps.virt_inds(indi),dl_impo.ket_virt_inds(indi),dl_impo.bra_virt_inds(indi));
            leg_combiners.push_back(temp_combiner);
            leg_combiners.push_back(dag(temp_combiner).prime());
        }
        
        //FIXME: modify left and right indices as well as combine seq for VL and VR
        //left and right indices 
        //left_inds are input(output) inds for VL(VR), and right_inds are output(input) inds for VL(VR)
        std::vector<IndexT> left_inds, right_inds;
        for (const auto &ind: leg_combiners[1].left()) left_inds.push_back(ind);
        for (const auto &ind: leg_combiners[2].left()) left_inds.push_back(ind);
        for (const auto &ind: leg_combiners[2*n_tensors_uc+1].left()) right_inds.push_back(dag(ind));
        for (const auto &ind: leg_combiners.back().left()) right_inds.push_back(dag(ind));
        //combiner_seq
        std::vector<std::vector<int>> VL_combiner_seq, VR_combiner_seq;
        //left combiner seq
        for (int tensori=0; tensori<n_tensors_uc; tensori++)
        {
            VL_combiner_seq.push_back(std::vector<int>{-2*(tensori+2)};
            VL_combiner_seq.push_back(std::vector<int>{});
            VL_combiner_seq.push_back(std::vector<int>{});
            VL_combiner_seq.push_back(std::vector<int>{-(2*(tensori+1)+1),2*tensori+2});
            VL_combiner_seq.push_back(std::vector<int>{});
            VL_combiner_seq.push_back(std::vector<int>{});
            VL_combiner_seq.push_back(std::vector<int>{2*(tensori+1)+1});
        }
        //right combiner seq
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

        //implement arnoldi method
        TensorT_Matrix_Arnoldi<TensorT> LMat(left_inds,right_inds,contract_tensors,VL_contract_seq,leg_combiners,VL_combiner_seq),
                                        RMat(right_inds,left_inds,contract_tensors,VR_contract_seq,leg_combiners,VR_combiner_seq);
        VL=TensorT(dag<IndexT>(left_inds));
        VR=TensorT(dag<IndexT>(right_inds));
        randTensor(VL);
        randTensor(VR);
        eta_L=arnoldi(LMat,VL,arnoldi_args);
        eta_R=arnoldi(RMat,VR,arnoldi_args);
    }

    if (type_name_.find("type_two")!=std::string::npos)
    {
        //contract_tensors
        std::vector<TensorT> contract_tensors;
        int impo_n_tensors_uc=dl_impo.n_tensors_uc(),
            imps_n_tensors_uc=2*impo_n_tensors_uc;
        for (int tensori=0; tensori<impo_n_tensors_uc; tensori++)
        {
            //contract two imps tensor and decombine the siteind
            TensorT imps_tensor=(dl_imps.site_tensors(tensori*2)*dag(dl_imps.siteind_combiners(tensori*2)))*(dl_imps.site_tensors(tensori*2+1)*dag(dl_imps.siteind_combiners(tensori*2+1)));
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
            contract_tensors.push_back(dag(impo_ket_tensor).prime(dl_imps.ket_siteinds(tensori*2)).prime(dl_imps.ket_siteinds(tensori*2+1)));
            contract_tensors.push_back(dag(impo_bra_tensor).prime(dl_imps.bra_siteinds(tensori*2)).prime(dl_imps.bra_siteinds(tensori*2+1)));
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
        rigt_inds.push_back(dag(dl_imps.virt_inds().back()));
        right_inds.push_back(prime(dl_imps.virt_inds().back()));

        //implement arnoldi method
        TensorT_Matrix_Arnoldi<TensorT> LMat(left_inds,right_inds,contract_tensors,VL_contract_seq),
                                        RMat(right_inds,left_inds,contract_tensors,VR_contract_seq);
        VL=TensorT(dag<IndexT>(left_inds));
        VR=TensorT(dag<IndexT>(right_inds));
        randTensorT(VL);
        randTensorT(VR);
        eta_L=arnoldi(LMat,VL,contract_opts);
        eta_R=arnoldi(RMat,VR,contract_opts);
    }

    //TODO: implement other types

    Print(eta_L);
    Print(eta_R);

    //update the imps
    //we set the outgoing inds of dl_impo to be new siteinds
    dl_imps.update_siteinds(dl_impo.ket_outgoing_inds(),dl_impo.bra_outgoing_inds());
    dl_imps.truncate_tensors(VL,VR,contract_opts);
}
template
void contract_dl_impo_imps(DL_iMPSt<ITensor> &dl_imps, DL_iMPOt<ITensor> &dl_impo, const Args &contract_opts);
template
void contract_dl_impo_imps(DL_iMPSt<IQTensor> &dl_imps, DL_iMPOt<IQTensor> &dl_impo, const Args &contract_opts);


template <TensorT>
iMPSt<TensorT> dl_imps_from_truncation(const std::vector<TensorT> &tensors_uc, const std::vector<typename TensorT::IndexT> &ket_siteinds, const std::vector<typename TensorT::IndexT> &bra_siteinds, const TensorT &VL, const TensorT &VR, int maxm, double cutoff)
{
    using IndexT=typename TensorT::IndexT;
    using IndexValT=typename TensorT::IndexValT;
    using CombinerT=typename TensorT::CombinerT;

    //set svd options
    Args svd_args;
    svd_args.add("Cutoff",cutoff);
    svd_args.add("Maxm",maxm);
    svd_args.add("DoRelCutoff",true);
    svd_args.add("AbsoluteCutoff",false);

    //eigendecompose VL and VR
    TensorT X,Y;
    eigen_factor(VR,X);
    eigen_factor(VL,Y);

    //singular value decompose of Y*X
    IndexT Uind=uniqueIndex(Y,VL), Vind=uniqueIndex(X,VR);
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
    tensor_svd(YX,U,D,V,svd_args);
    D/=D.norm();
    Print(D.diag());

    //update union_tensor
    TensorT Dinv_sqrt=dag(D);
    Dinv_sqrt.pseudoInvert();
    Dinv_sqrt.mapElems([](double x){ return std::sqrt(abs(x)); });
    TensorT union_tensor=dag(U)*Y;
    for (const auto &tens: tensors_uc) union_tensor*=tens;
    union_tensor*=X*dag(V);

    //update site_tensors_ by svd. 
    //for example for n=2 case
    //
    // --union_tensor-- = --P--nlD--Q--
    //       | |            |       |
    //
    //Then, we have
    //
    //  --tens0-- = --lDinv_sqrt--P--nlD_sqrt--
    //      |                     |
    //
    //  --tens1-- = --nlD_sqrt--Q--rDinv_sqrt--
    //      |                   |
    //
    TensorT lDinv_sqrt=prime(Dinv_sqrt,commonIndex(Dinv_sqrt,V)),
            rDinv_sqrt=prime(Dinv_sqrt,commonIndex(Dinv_sqrt,U));
    std::vector<TensorT> site_tensors;
    for (int sitei=0; sitei<n_sites_uc_-1; sitei++)
    {
        TensorT P(ket_siteinds[sitei],bra_siteinds[sitei],dag(commonIndex(U,lDinv_sqrt))),nlD,Q;
        svd(union_tensor,P,nlD,Q,svd_args);
        nlD/=nlD.norm();
        Print(nlD.diag());

        TensorT nlD_sqrt=nlD;
        nlD_sqrt.mapElems([](double x){ return std::sqrt(abs(x)); });

        //update site_tensors_ and bulk virt_inds_
        site_tensors.push_back(lDinv_sqrt*P*nlD_sqrt);
        site_tensors.noprime();

        union_tensor=nlD*Q;
        lDinv_sqrt=dag(nlD);
        lDinv_sqrt.pseudoInvert();
        lDinv_sqrt.mapElems([](double x){ return std::sqrt(abs(x)); });
    }
    site_tensors.push_back(lDinv_sqrt*union_tensor*rDinv_sqrt);
    site_tensors.noprime();

    DL_iMPSt<TensorT> dl_imps(site_tensors,site_inds);
    return dl_imps;
}
