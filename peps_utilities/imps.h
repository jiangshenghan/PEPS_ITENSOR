
#ifndef _IMPS_H_
#define _IMPS_H_

#include "utilities.h"
#include "tensor_matrix.h"

//
//class iMPSt stores infinite MPS (with only site tensors)
//
template <class TensorT>
class iMPSt
{
    public:
        //
        //type alias
        //
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;

        //
        //Constructors
        //
        iMPSt() {}
        iMPSt(const std::vector<TensorT> &site_tensors, const IndexT &lind, const IndexT &rind):
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

            //init uc_tensor_
            uc_tensor_=site_tensors_[0];
            for (int sitei=1; sitei<n_sites_uc_; sitei++) uc_tensor_*=site_tensors_[sitei];
        }

        //TODO: remove following
        {
            //init arnoldi_args_
            arnoldi_args_.add("MaxIter",10);
            arnoldi_args_.add("MaxRestart",1);
            arnoldi_args_.add("ErrGoal",1e-16);

        }

        //truncate imps where the dominant eigenvectors are known
        //VL_indices={dag(lind),prime(lind)}
        //VR_indices={dag(rind),prime(rind)}
        void truncate_tensors(const TensorT &VL, const TensorT &VR, int maxm, double cutoff=1e-16)
        {
            //set svd options
            Args svd_args;
            svd_args.add("Cutoff",cutoff);
            svd_args.add("Maxm",maxm);
            svd_args.add("DoRelCutoff",true);
            svd_args.add("AbsoluteCutoff",false);

            //eigendecompose VL and VR
            TensorT X,Y;
            eigen_factor(VL,X);
            eigen_factor(VR,Y);

            //singular value decompose of Y*X
            TensorT YX=tensor_contraction<TensorT,IndexT>(Y,X,{dag(virt_inds_[0])},{virt_inds_.back()});
            TensorT U(commonIndex(YX,Y)),D,V;
            tensor_svd(YX,U,D,V,svd_args);
            D/=D.norm();
            Print(D.diag());

            //update imps
            TensorT D_sqrt=D, Dinv_sqrt=dag(D);
            D_sqrt.mapElems([](double x){ return std::sqrt(abs(x)); });
            Dinv_sqrt.pseudoInvert();
            Dinv_sqrt.mapElems([](double x){ return std::sqrt(abs(x)); });
            uc_tensor_=uc_tensor_*(Dinv_sqrt*dag(U)*Y)*(Dinv_sqrt*dag(V)*X);
            //update boundary virt_inds_
            virt_inds_[0]=commonIndex(Dinv_sqrt,V);
            virt_inds_[n_sites_uc_]=commonIndex(Dinv_sqrt,U);

            if (n_sites_uc_==1)
            {
                site_tensors_[0]=uc_tensor_;
                return;
            }
            //update site_tensors_ svd. Notice, we should also update virt_inds
            //for example for n=2 case
            //
            // --lD_sqrt--uc_tensor--rD_sqrt-- = --P--nlD--Q--
            //              | |                       | |
            //
            //Then, we have
            //
            //  --tens0-- = --lDinv_sqrt--P--nlD_sqrt--
            //      |                     |
            //
            //  --tens1-- = --nlD_sqrt--Q--rDinv_sqrt--
            //      |                   |
            //
            TensorT lD_sqrt=prime(D_sqrt,virt_inds_.back())
                    lDinv_sqrt=prime(Dinv_sqrt,virt_inds_.back()),
                    rD_sqrt=prime(D_sqrt,virt_inds_[0]),
                    rDinv_sqrt=prime(Dinv_sqrt,virt_inds_[0]);
            TensorT rtensor=lD_sqrt*uc_tensor_*rD_sqrt;
            IndexT llDind=prime(virt_inds_.back());
            for (int sitei=0; sitei<n_sites_uc_-1; sitei++)
            {
                TensorT P(site_inds_[sitei],llDind),nlD,Q;
                svd(rtensor,P,nlD,Q,svd_args);
                nlD/=nlD.norm();
                Print(nlD.diag());

                TensorT nlD_sqrt=nlD, nlDinv_sqrt=dag(nlD);
                nlD_sqrt.mapElems([](double x){ return std::sqrt(abs(x)); });
                nlDinv_sqrt.pseudoInvert();
                nlDinv_sqrt.mapElems([](double x){ return std::sqrt(abs(x)); });

                site_tensors_[sitei]=lDinv_sqrt*P*nlD_sqrt;
                virt_inds_[sitei+1]=commonIndex(Q,nlD_sqrt);
                rtensor=nlD*Q;
                lD_sqrt=nlD_sqrt;
                lDinv_sqrt=nlDinv_sqrt;
                llDind=commonIndex(nlD,P);
            }
            site_tensors_[n_sites_uc_-1]=lDinv_sqrt*rtensor*rDinv_sqrt;
        }

        //
        //Access Methods
        //
        int n_sites_uc() const { return n_sites_uc_; }
        const IndexT &virt_inds(int indi) const { return virt_inds_[indi]; }
        const IndexT &site_inds(int sitei) const { return site_inds_[sitei]; }
        const TensorT &site_tensors(int sitei) const { return site_tensors_[sitei]; }
        const TensorT &uc_tensor() const { return uc_tensor_; }

        //
        //Assign Methods
        //

    protected:
        int n_sites_uc_;
        //site_inds may not have type Site. We fix # of site indice of each tensor to be 1
        std::vector<IndexT> site_inds_;
        //virt legs for imp, we use the convention
        //
        //                   |
        //                site_inds
        //                   |
        // --virt_inds[i]--tensor[i]--dag(virt_inds[i+1])--
        //
        std::vector<IndexT> virt_inds_;
        std::vector<TensorT> site_tensors_;
        //uc_tensor_ is formed by multiplication of all site_tensors_ 
        TensorT uc_tensor_;
};


//
//class for imps from double layer tensors
//compare with conventional peps, the major difference is there are in fact two phys legs per site_tensor.
//
//       |      ket_ind
//       |     /       \
//  site_tensor         siteind_combiner--
//       |     \       /
//       |      bra_ind
//
//Notice, we should not set tind to be prime(dag(bind)) to avoid confusion
//
template <class TensorT>
class DL_iMPSt 
{
    public:
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;

        //
        //Constructors
        //
        DL_iMPSt() {}
        DL_iMPSt(std::vector<TensorT> site_tensors, const IndexT &lind, const IndexT &rind, const std::vector<IndexT> &ket_siteinds, std::vector<IndexT> &bra_siteinds): ket_siteinds_(ket_siteinds), bra_siteinds_(bra_siteinds)
        {
            std::vector<IndexT> site_inds;
            for (int sitei=0; sitei<site_tensors.size(); sitei++)
            {
                if (ket_siteinds_[sitei]==prime(dag(bottom_siteinds_sitei])) 
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

        //
        //Acess Methods
        //
        int n_sites_uc() const { return imps_.n_sites_uc(); }
        const IndexT &site_inds(int tensori) const { return imps_.site_inds(tensori); }
        const IndexT &virt_inds(int indi) const { return imps_.virt_inds(indi); }
        const TensorT &site_tensors(int tensori) const { return imps_.site_tensors(tensori); }
        const IndexT &ket_siteinds(int indi) const { return ket_siteinds_[indi]; }
        const IndexT &bra_siteinds(int indi) const { return bra_siteinds_[indi]; }
        const CombinerT &siteind_combiners(int combineri) const { return siteind_combiners_[combineri]; }
        const iMPSt<TensorT> &imps const { return imps_; }

        //
        //Methods
        //
        void truncate_tensors(const TensorT &VL, const TensorT &VR, const Args &args) 
        { 
            int maxm=args.getInt("Maxm"),
                cutoff=args.getReal("Cutoff");
            imps_.truncate_tensors(VL,VR,maxm,cutoff); 
        }

    private:
        iMPSt<TensorT> imps_;
        std::vector<IndexT> ket_siteinds_, bra_siteinds_;
        std::vector<CombinerT> siteind_combiners_;
};


//
//class for impo from double layer peps
//an impo site tensor is formed by contraction of ket and bra tensor, which share the same phys_legs
//
//          /
// i   --ttens--
// M       /|
// P        |/
// S   --btens--
//         /
//
//we will save the incoming and outgoing legs, where other legs can be obtained by commonIndex method
//We implement different types of impos, shown in Fig.5 arXiv:0711.3960
//type_one: one-one-one
//type_two: two-one-two
//type_three: two-one-one
//type_four: one-one-two
template <class TensorT>
class DL_iMPOt
{
    public:
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;

        //
        //Constructor
        //
        DL_iMPOt() {}
        //we will reconstruct the indices for all tensors
        DL_iMPOt(std::string type_name, const std::vector<IndexT> &init_ket_tensors, const std::vector<IndexT> &init_incoming_inds, const std::vector<IndexT> &init_outgoing_inds, const std::vector<IndexT> init_boundary_inds=std::vector<IndexT>()):
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

        //
        //Acess Method
        //
        std::string type_name() const { return type_name_; }
        int n_tensors_uc() const { return n_tensors_uc_; }
        const IndexT &bra_incoming_inds(int indi) const { return bra_incoming_inds_[indi]; }
        const IndexT &bra_outgoing_inds(int indi) const { return bra_outgoing_inds_[indi]; }
        const IndexT &ket_incoming_inds(int indi) const { return ket_incoming_inds_[indi]; }
        const IndexT &ket_outgoing_inds(int indi) const { return ket_outgoing_inds_[indi]; }
        const IndexT &bra_virt_inds(int indi) const { return bra_virt_inds_[indi]; }
        const IndexT &ket_virt_inds(int indi) const { return ket_virt_inds_[indi]; }
        const TensorT &ket_tensors(int tensori) const { return ket_tensors_[tensori]; }
        const TensorT &bra_tensors(int tensori) const { return bra_tensors_[tensori]; }


    private:
        //type_name_ follows fig.5 in arXiv:0711.3960
        std::string type_name_;
        //number of impo tensors in one unit cell
        //n_tensors_uc_ should be compatible with the acting imps
        int n_tensors_uc_;
        //four type of site legs (site legs are not physical legs)
        std::vector<IndexT> bra_incoming_inds_, ket_incoming_inds_, bra_outgoing_inds_, ket_outgoing_inds_;
        //virt legs for impo, we use the convention (for type_one)
        //
        //                   |
        //                outgoing
        //                   |
        // --virt_inds[i]--tensor[i]--dag(virt_inds[i+1])--
        //                   |
        //                 incoming
        //                   |
        //
        std::vector<IndexT> bra_virt_inds_, ket_virt_inds_;
        std::vector<TensorT> bra_tensors_, ket_tensors_;
};


//contract given dl_imps and dl_impo, and then do compression
//contract_opts:
//getInt: Maxm(for svd), MaxIter(for arnoldi), MaxRestart(for arnoldi)
//getReal: Cutoff(for svd), ErrorGoal(for arnoldi)
//getString: ContractMethod("single_layer")
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

    if (type_name_.find("type_one")!=std::string::npos)
    {
        //construct member variables of TensorT_Matrix Class
        //contract_tensors
        std::vector<TensorT> contract_tensors;
        int n_tensors_uc=dl_impo.n_tensors_uc();
        for (int sitei=0; sitei<n_tensors_uc; sitei++)
        {
            contract_tensors.push_back(dl_imps.site_tensors(sitei));
            contract_tensors.push_back(dl_impo.ket_tensors(sitei));
            contract_tensors.push_back(dl_impo.bra_tensors(sitei));
            contract_tensors.push_back(dag(dl_impo.ket_tensors(sitei)).prime().noprime(dl_impo.ket_outgoing_inds(sitei)));
            contract_tensors.push_back(dag(dl_impo.bra_tensors(sitei)).prime().noprime(dl_impo.bra_outgoing_inds(sitei)));
            contract_tensors.push_back(dag(dl_imps.site_tensors()));
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
        //leg_combiners, which combines virt_inds of both imps and impo
        std::vector<CombinerT> leg_combiners;
        //set 0th combiner to be empty
        leg_combiners.push_back(CombinerT());
        for (int indi=0; indi<n_tensors_uc; indi++)
        {
            CombinerT temp_combiner(dl_imps.virt_inds(indi),dl_impo.ket_virt_inds(indi),dl_impo.bra_virt_inds(indi));
            leg_combiners.push_back(temp_combiner);
            leg_combiners.push_back(dag(temp_combiner).prime());
        }
        //left and right indices 
        //left_inds are input(output) inds for VL(VR), and right_inds are output(input) inds for VL(VR)
        std::vector<IndexT> left_inds, right_inds;
        left_inds.push_back(leg_combiners[1].right());
        left_inds.push_back(leg_combiners[2].right());
        right_inds.push_back(dag(leg_combiners[2*n_tensors_uc+1].right()));
        right_inds.push_back(dag(leg_combiners.back()));
        //combiner_seq
        std::vector<std::vector<int>> VL_combiner_seq, VR_combiner_seq;
        VL_combiner_seq.push_back(std::vector<int>{});
        VR_combiner_seq.push_back(std::vector<int>{});
        for (int tensori=0; tensori<n_tensors_uc; tensori++)
        {
            VL_combiner_seq.back().push_back(2*tensori+1);
            VL_combiner_seq.push_back(std::vector<int>{});
            VL_combiner_seq.push_back(std::vector<int>{});
            VL_combiner_seq.push_back(std::vector<int>{-(2*(tensori+1)+1),2*tensori+2});
            VL_combiner_seq.push_back(std::vector<int>{});
            VL_combiner_seq.push_back(std::vector<int>{});
            VL_combiner_seq.push_back(std::vector<int>{-(2*(tensori+1)+2)});
        }
        for (int tensori=n_tensor_uc-1; tensori>=0; tensori--)
        {
            VR_combiner_seq.back().push_back(-(2*(tensori+1)+2));
            VR_combiner_seq.push_back(std::vector<int>{});
            VR_combiner_seq.push_back(std::vector<int>{});
            VR_combiner_seq.push_back(std::vector<int>{2*tensori+2,-(2*(tensori+1)+1)});
            VR_combiner_seq.push_back(std::vector<int>{});
            VR_combiner_seq.push_back(std::vector<int>{});
            VR_combiner_seq.push_back(std::vector<int>{2*tensori+1});
        }

        TensorT_Matrix_Arnoldi<TensorT> LMat(left_inds,right_inds,contract_tensors,VL_contract_seq,leg_combiners,VL_combiner_seq),
                                        RMat(right_inds,left_inds,contract_tensors,VR_contract_seq,leg_combiners,VR_combiner_seq);

        //implement arnoldi method
        VL=TensorT(dag<IndexT>(left_indices));
        VR=TensorT(dag<IndexT>(right_indices));
        randTensor(VL);
        randTensor(VR);
        eta_L=arnoldi(LMat,VL,contract_opts);
        eta_R=arnoldi(RMat,VR,contract_opts);
    }

    if (type_name_.find("type_two")!=std::string::npos)
    {
    }

    //TODO: implement other types

    Print(eta_L);
    Print(eta_R);
    //truncate the imps
    dl_imps.truncate_tensors(VL,VR,contract_opts);
}

#endif
