
#ifndef _IMPS_H_
#define _IMPS_H_

#include "utilities.h"
#include "tensor_matrix.h"
#include "tensor_svd.h"

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
        iMPSt(const std::vector<TensorT> &site_tensors, const IndexT &lind, const IndexT &rind);
        iMPSt(const std::vector<TensorT> &site_tensors, const std::vector<IndexT> &site_inds);

        //
        //Methods
        //
        bool valid() { return !site_tensors_.empty(); }
        //change the order of tensors in one u.c. by moving tensors
        void move_tensors(int movei=1);
        

        //
        //Access Methods
        //
        int n_sites_uc() const { return n_sites_uc_; }
        const IndexT &virt_inds(int indi) const { return virt_inds_[indi]; }
        const std::vector<IndexT> &virt_inds() const { return virt_inds_; }
        const IndexT &site_inds(int sitei) const { return site_inds_[sitei]; }
        const std::vector<IndexT> &site_inds() const { return site_inds_; }
        const TensorT &site_tensors(int sitei) const { return site_tensors_[sitei]; }
        const std::vector<TensorT> &site_tensors() const { return site_tensors_; }


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
        DL_iMPSt(const std::vector<IndexT> &ket_siteinds, const std::vector<IndexT> &bra_siteinds, std::vector<TensorT> site_tensors);
        //DL_iMPSt(std::vector<TensorT> site_tensors, const IndexT &lind, const IndexT &rind, const std::vector<IndexT> &ket_siteinds, std::vector<IndexT> &bra_siteinds);
        //init using ket tensors, the site tensors are obtained from contraction ket tensors and bra tensors with ket_siteinds and sl_virt_inds uncontracted 
        //indices will be replaced to avoid incorrect contraction
        DL_iMPSt(std::vector<TensorT> ket_site_tensors, const std::vector<IndexT> &ket_siteinds, const std::vector<IndexT> &sl_virt_inds);

        //
        //Methods
        //
        //true if not default construct
        bool valid() { return imps_.valid(); }
        //change the order of tensors in one u.c. by moving tensors
        void move_tensors(int movei=1);

        //
        //Acess Methods
        //
        int n_sites_uc() const { return imps_.n_sites_uc(); }
        const IndexT &site_inds(int tensori) const { return imps_.site_inds(tensori); }
        const IndexT &virt_inds(int indi) const { return imps_.virt_inds(indi); }
        const std::vector<IndexT> &virt_inds() const { return imps_.virt_inds(); }
        const TensorT &site_tensors(int tensori) const { return imps_.site_tensors(tensori); }
        //get site_tensors with ket and bra indices separate
        TensorT dl_site_tensors(int tensori) const { return imps_.site_tensors(tensori)*dag(siteind_combiners_[tensori]); }
        const IndexT &ket_siteinds(int indi) const { return ket_siteinds_[indi]; }
        const std::vector<IndexT> &ket_siteinds() const { return ket_siteinds_; }
        const IndexT &bra_siteinds(int indi) const { return bra_siteinds_[indi]; }
        const std::vector<IndexT> &bra_siteinds() const { return bra_siteinds_; }
        const CombinerT &siteind_combiners(int combineri) const { return siteind_combiners_[combineri]; }
        const std::vector<CombinerT> &siteind_combiners() const { return siteind_combiners_; }
        const iMPSt<TensorT> &imps() const { return imps_; }


    private:
        std::vector<IndexT> ket_siteinds_, bra_siteinds_;
        std::vector<CombinerT> siteind_combiners_;
        iMPSt<TensorT> imps_;
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
        DL_iMPOt(std::string type_name, const std::vector<TensorT> &init_ket_tensors, const std::vector<IndexT> &init_ket_incoming_inds, const std::vector<IndexT> &init_ket_outgoing_inds, std::vector<IndexT> init_ket_virt_inds=std::vector<IndexT>());

        //
        //Acess Method
        //
        std::string type_name() const { return type_name_; }
        int n_tensors_uc() const { return n_tensors_uc_; }
        const IndexT &bra_incoming_inds(int indi) const { return bra_incoming_inds_[indi]; }
        const std::vector<IndexT> &bra_incoming_inds() const { return bra_incoming_inds_; }
        const IndexT &bra_outgoing_inds(int indi) const { return bra_outgoing_inds_[indi]; }
        const std::vector<IndexT> &bra_outgoing_inds() const { return bra_outgoing_inds_; }
        const IndexT &ket_incoming_inds(int indi) const { return ket_incoming_inds_[indi]; }
        const std::vector<IndexT> &ket_incoming_inds() const { return ket_incoming_inds_; }
        const IndexT &ket_outgoing_inds(int indi) const { return ket_outgoing_inds_[indi]; }
        const std::vector<IndexT> &ket_outgoing_inds() const { return ket_outgoing_inds_; }
        const IndexT &bra_virt_inds(int indi) const { return bra_virt_inds_[indi]; }
        const std::vector<IndexT> &bra_virt_inds() const { return bra_virt_inds_; }
        const IndexT &ket_virt_inds(int indi) const { return ket_virt_inds_[indi]; }
        const std::vector<IndexT> &ket_virt_inds() const { return ket_virt_inds_; }
        const TensorT &ket_tensors(int tensori) const { return ket_tensors_[tensori]; }
        const std::vector<TensorT> &ket_tensors() const { return ket_tensors_; }
        const TensorT &bra_tensors(int tensori) const { return bra_tensors_[tensori]; }
        const std::vector<TensorT> &bra_tensors() const { return bra_tensors_; }

        //
        //Methods
        //
        bool valid() { return !ket_tensors_.empty(); }


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


//
//Contraction Methods for imps
//
//contract_opts:
//getInt: Maxm(for svd), MaxIter(for arnoldi), MaxRestart(for arnoldi)
//getReal: Cutoff(for svd), ErrGoal(for arnoldi)
//getString: AbsorbBond("left_bond","right_bond")
//
//contract given imps and impo, and then do compression
template <class TensorT>
void contract_impo_imps(iMPSt<TenosrT> &imps, const iMPOt<TensorT> &imps, Args contract_opts);

//contract given dl_imps and dl_impo, and then do compression
template <class TensorT>
void contract_dl_impo_imps(DL_iMPSt<TensorT> &dl_imps, const DL_iMPOt<TensorT> &dl_impo, Args contract_opts);

//truncate imps where the dominant eigenvectors are known
//VL_indices={dag(lind),prime(lind)}
//VR_indices={dag(rind),prime(rind)}
//Here product of tensors_uc equals product of all site tensors in one uc (may not share the same virt legs)
//VL=X.X^dagger, VR=Y.Y^dagger
//union_tensor=U^dagger.Y.(prod of tensors_uc).X.V^dagger
//we then do svd on union_tensor to obtain site_tensors
//
//trunc_opts:
//getInt: Maxm
//getReal: Cutoff
//getString: AbsorbBond("left_bond","right_bond")
template <class TensorT>
DL_iMPSt<TensorT> dl_imps_from_truncation(const std::vector<TensorT> &tensors_uc, const std::vector<typename TensorT::IndexT> &ket_siteinds, const std::vector<typename TensorT::IndexT> &bra_siteinds, const TensorT &VL, const TensorT &VR, Args trunc_opts);

//print methods for above classes
template <class TensorT>
inline std::ostream &operator<<(std::ostream &s, const iMPSt<TensorT> imps)
{
    s << endl;
    s << "n_sites_uc=" << imps.n_sites_uc() << endl;
    s << "site_tensors:" << endl << imps.site_tensors() << endl;
    s << "site_inds:" << endl << imps.site_inds() << endl;
    s << "virt_inds:" << endl << imps.virt_inds() << endl;
    return s;
}

template <class TensorT>
inline std::ostream &operator<<(std::ostream &s, const DL_iMPSt<TensorT> dl_imps)
{
    s << endl;
    s << "imps:" << endl << dl_imps.imps() << endl;
    s << "ket_siteinds:" << endl << dl_imps.ket_siteinds() << endl;
    s << "bra_siteinds:" << endl << dl_imps.bra_siteinds() << endl;
    s << "siteind_combiners:" << endl << dl_imps.siteind_combiners() << endl;
}

template <class TensorT>
inline std::ostream &operator<<(std::ostream &s, const DL_iMPOt<TensorT> &dl_impo)
{
    s << endl;
    s << dl_impo.type_name() << endl;
    s << "n_tensors_uc=" << dl_impo.n_tensors_uc() << endl;
    s << "ket_tensors:" << endl << dl_impo.ket_tensors() << endl;
    s << "bra_tensors:" << endl << dl_impo.bra_tensors() << endl;
    s << "ket_incoming_inds:" << endl << dl_impo.ket_incoming_inds() << endl;
    s << "ket_outgoing_inds:" << endl << dl_impo.ket_outgoing_inds() << endl;
    s << "ket_virt_inds" << endl << dl_impo.ket_virt_inds() << endl;
    s << "bra_incoming_inds:" << endl << dl_impo.bra_incoming_inds() << endl;
    s << "bra_outgoing_inds:" << endl << dl_impo.bra_outgoing_inds() << endl;
    s << "bra_virt_inds:" << endl << dl_impo.bra_virt_inds() << endl;
    return s;
}


#endif
