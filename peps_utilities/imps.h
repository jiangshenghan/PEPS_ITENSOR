
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
            lind_(lind), rind_(rind),
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
                    if (ind==lind_ || ind==rind_) continue;
                    site_inds_.push_back(ind);
                    break;
                }
            }

            //init uc_tensor_
            uc_tensor_=site_tensors_[0];
            for (int sitei=1; sitei<n_sites_uc_; sitei++) uc_tensor_*=site_tensors_[sitei];
        }

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
            TensorT YX=tensor_contraction<TensorT,IndexT>(Y,X,{dag(lind_)},{dag(rind_)});
            TensorT U(commonIndex(YX,Y)),D,V;
            tensor_svd(YX,U,D,V,svd_args);
            D/=D.norm();

            //update imps
            TensorT D_sqrt=D, Dinv_sqrt=dag(D);
            D_sqrt.mapElems([](double x){ return std::sqrt(abs(x)); });
            Dinv_sqrt.pseudoInvert();
            Dinv_sqrt.mapElems([](double x){ return std::sqrt(abs(x)); });

            uc_tensor_=uc_tensor_*(Dinv_sqrt*dag(U)*Y)*(Dinv_sqrt*dag(V)*X);
            lind_=commonIndex(Dinv_sqrt,V);
            rind_=commonIndex(Dinv_sqrt,U);

            if (n_sites_uc_==1)
            {
                site_tensors_[0]=uc_tensor_;
                return;
            }
            //update site_tensors_ by svd
            //TODO: keep check of the indices
            TensorT lD_sqrt=prime(D_sqrt,rind_)
                    lDinv_sqrt=prime(Dinv_sqrt,rind_),
                    rD_sqrt=prime(D_sqrt,lind_),
                    rDinv_sqrt=prime(Dinv_sqrt,lind_);
            TensorT rtensor=lD_sqrt*uc_tensor_*rD_sqrt;
            IndexT lDind=prime(rind_);
            for (int sitei=0; sitei<n_sites_uc_-1; sitei++)
            {
                TensorT P(site_inds_[sitei],lDind),nlD,Q;
                svd(rtensor,P,nlD,Q,svd_args);
                nlD/=nlD.norm();

                TensorT nlD_sqrt=nlD, nlDinv_sqrt=dag(nlD);
                nlD_sqrt.mapElems([](double x){ return std::sqrt(abs(x)); });
                nlDinv_sqrt.pseudoInvert();
                nlDinv_sqrt.mapElems([](double x){ return std::sqrt(abs(x)); });

                site_tensors_[sitei]=lDinv_sqrt*P*nlD_sqrt;
                rtensor=nlD*Q;
                lD_sqrt=nlD_sqrt;
                lDinv_sqrt=nlDinv_sqrt;
            }
            site_tensors_[n_sites_uc_-1]=lDinv_sqrt*rtensor*rDinv_sqrt;
        }

        //
        //Access Methods
        //
        const IndexT &lind() const { return lind_; }
        const IndexT &rind() const { return rind_; }
        const TensorT &site_tensors(int sitei) const { return site_tensors_[sitei]; }
        const TensorT &uc_tensor() const { return uc_tensor_; }

        //
        //Assign Methods
        //

    protected:
        int n_sites_uc_;
        //lind_/rind_ stores the leftmost/rightmost indices. They should has same dim and opposite direction
        //lr may also denote up/down depends on the orientation of imps
        IndexT lind_, rind_;
        //site_inds may not have type Site. We fix # of site indice of each tensor to be 1
        std::vector<IndexT> site_inds_;
        std::vector<TensorT> site_tensors_;
        //uc_tensor_ is formed by multiplication of all site_tensors_ 
        TensorT uc_tensor_;
};


//
//class for imps from double layer tensors
//compare with conventional peps, the major difference is there are in fact two phys legs per site_tensor.
//
//       |      tind
//       |     /    \
//  site_tensor      site_inds--
//       |     \    /
//       |      bind
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
        DL_iMPSt(std::vector<TensorT> site_tensors, const IndexT &lind, const IndexT &rind, const std::vector<IndexT> &top_siteinds, std::vector<IndexT> &bottom_siteinds): top_siteinds_(top_siteinds), bottom_site_inds_(bottom_siteinds)
        {
            std::vector<IndexT> site_inds;
            for (int sitei=0; sitei<site_tensors.size(); sitei++)
            {
                if (top_siteinds_[sitei]==prime(dag(bottom_siteinds_sitei])) 
                {
                    cout << "Should not set top_siteinds=prime(dag(bottom_siteinds))!" << endl;
                    exit(1);
                }
                siteind_combiners_.push_back(CombinerT(top_siteinds_[sitei],bottom_siteinds_[sitei]));
                site_inds.push_back(siteind_combiners_[sitei].right());
                site_tensors[sitei]=site_tensors[sitei]*siteind_combiners_[sitei];
            }
            imps_=iMPSt<TensorT>(site_tensors,lind,rind);
        }

    private:
        iMPSt<TensorT> imps_;
        std::vector<IndexT> top_siteinds_, bottom_site_inds_;
        std::vector<CombinerT> siteind_combiners_;
}


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
class DLPEPS_iMPOt
{
    public:
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;

        //
        //Constructor
        //
        DLPEPS_iMPOt() {}
        //we will reconstruct the indices for all tensors
        DLPEPS_iMPOt(std::string type_name, const std::vector<IndexT> &init_ket_tensors, const std::vector<IndexT> &init_incoming_inds, const std::vector<IndexT> &init_outgoing_inds, const std::vector<IndexT> init_boundary_inds=std::vector<IndexT>())
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

                init_virt_inds[0]=dag(init_boundary_inds[0]);
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


    private:
        //type_name_ follows fig.5 in arXiv:0711.3960
        std::string type_name_;
        //number of impo tensors in one unit cell
        int n_tensors_uc_;
        //four type of site legs (site legs are not physical legs)
        std::vector<IndexT> bra_incoming_inds_, ket_incoming_inds_, bra_outgoing_inds_, ket_outgoing_inds_;
        //virt legs for impo, we use the convention
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
}

#endif
