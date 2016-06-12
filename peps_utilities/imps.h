
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
//an impo site tensor is formed by contraction of top and bottom tensor, which share the same phys_legs
//
//          /
// i   --ttens--
// M       /|
// P        |/
// S   --btens--
//         /
//
//there are four site legs per tensor. Top tensor and bottom tensor share a phys leg (type==Site)
//we will save all site legs as well as four boundary legs 
//
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
        DLPEPS_iMPOt(const std::vector<TensorT> &btensors, const IndexT &blind, const IndexT &brind, const std::vector<IndexT> &inb_siteinds, const std::vector<IndexT> &outb_siteinds):
            n_sites_uc_(site_tensors.size())
        {
            bottom_lind_=isomorphic_legs(blind,"bottom_left_ind");
            bottom_rind_=isomorphic_legs(brind,"bottom_right_ind");
            top_lind_=isomorphic_legs(dag(blind),"top_left_ind");
            top_rind_=isomorphic_legs(dag(brind),"top_right_ind");

            for (int sitei=0; sitei<n_sites_uc_; sitei++)
            {
                inbottom_siteinds_.push_back(isomorphic_legs(inb_siteinds[sitei]),nameint("in_bottom_siteind_",sitei));
                outbottom_siteinds_.push_back(isomorphic_legs(outb_siteinds[sitei],nameint("out_bottom_siteind_",sitei)));
                intop_siteinds_.push_back(isomorphic_legs(dag(inb_siteinds[sitei])),nameint("in_top_siteind_",sitei));
                outtop_siteinds_.push_back(isomorphic_legs(dag(outb_siteinds[sitei]),nameint("out_top_siteind_",sitei)));
            }

            //TODO
            for ()
        }

    private:
        int n_sites_uc_;
        //four boundary legs
        IndexT top_lind_, bottom_lind_, top_rind_, bottom_rind_;
        //four type of site legs
        std::vector<IndexT> intop_siteinds_, inbottom_siteinds_, outtop_siteinds_, outbottom_siteinds_;
        std::vector<TensorT> top_tensors_, bottom_tensors_;
}

#endif
