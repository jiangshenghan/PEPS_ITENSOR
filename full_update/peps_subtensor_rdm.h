
#ifndef _PEPS_SUBTENSOR_RDM_H_
#define _PEPS_SUBTENSOR_RDM_H_

#include "square_rvb.h"
#include "kagome_rvb.h"
#include "trotter_gate.h"
#include "simple_update_env.h"

template <class TensorT>
class PEPSt_Subtensor_RDM
{
    public:
        //
        //type alias
        //
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;

        //
        //Constructor
        //
        PEPSt_Subtensor_RDM(const std::vector<TensorT> &cutting_site_tensors, const std::vector<TensorT> &cutting_bond_tensors, const std::vector<TensorT> &env_tensors, const std::vector<int> contract_seq=std::vector<int>());

        //
        //Acess Methods
        //
        int cutting_sites_no() const { return cutting_site_tensors_.size(); }
        int cutting_bonds_no() const { return cutting_bond_tensors_.size(); }

        const IndexT &cutting_phys_legs(int sitei) const
        {
            for (const auto &ind : cutting_site_tensors_[sitei].indices())
            {
                if (ind.type()==Site)
                    return ind;
            }
            return IndexT::Null();
        }

        const std::vector<TensorT> &cutting_site_tensors() const { return cutting_site_tensors_; }
        const TensorT &cutting_site_tensors(int sitei) const { return cutting_site_tensors_[sitei]; }
        const std::vector<TensorT> &cutting_bond_tensors() const { return cutting_bond_tensors_; }
        const TensorT &cutting_bond_tensors(int bondi) const { return cutting_bond_tensors_[bondi]; }
        const std::vector<TensorT> &sub_tensors() const { return sub_tensors_; }
        const TensorT &sub_tensors(int sitei) const { return sub_tensors_[sitei]; }
        const TensorT &sub_env_tensor() const { return sub_env_tensor_; }

        Complex wf_norm_sq() const { return wf_norm_sq_; }
        const TensorT &RDM() const { return RDM_; }


        //
        //Calculate env_tensor_ and RDM_
        //
        void obtain_RDM();
        void update_RDM(const std::vector<TensorT> &cutting_site_tensors, const std::vector<TensorT> &cutting_bond_tensors, const std::vector<TensorT> &env_tensors, const std::vector<int> &contract_seq=std::vector<int>());
        void obtain_sub_env_tensor();
        void obtain_sub_tensors();

        //Calculate expect val by replacing tensors
        Complex expect_val_from_replaced_tensors(std::array<std::vector<TensorT>,2> replaced_sub_tensors_ket_bra, std::array<std::vector<TensorT>,2> replaced_bond_tensors_ket_bra) const
        {
            //reverse dir and add prime to bra tensor
            IndexT ket_phys_leg=findtype(replaced_sub_tensors_ket_bra[0][0],Site), 
                   bra_phys_leg=findtype(replaced_sub_tensors_ket_bra[1][0],Site);
            if (ket_phys_leg.dir()==bra_phys_leg.dir()) 
            {
                for (auto &sub_tensor_bra : replaced_sub_tensors_ket_bra[1]) sub_tensor_bra.dag().prime(Link);
                for (auto &bond_tensor_bra : replaced_bond_tensors_ket_bra[1]) bond_tensor_bra.dag().prime();
            }

            TensorT expect_val_tensor=sub_env_tensor_;
            for (int sitei=0; sitei<cutting_site_tensors_.size(); sitei++) expect_val_tensor=expect_val_tensor*replaced_sub_tensors_ket_bra[0][sitei]*replaced_sub_tensors_ket_bra[1][sitei];
            for (int bondi=0; bondi<cutting_bond_tensors_.size(); bondi++) expect_val_tensor=expect_val_tensor*replaced_bond_tensors_ket_bra[0][bondi]*replaced_bond_tensors_ket_bra[1][bondi];
            return expect_val_tensor.toComplex();
        }

        Complex expect_val_from_replaced_tensors(const std::vector<TensorT> &replaced_sub_tensors_ket, const std::vector<TensorT> &replaced_bond_tensors_ket) const
        {
            std::vector<TensorT> replaced_sub_tensors_bra, replaced_bond_tensors_bra;
            for (const auto &site_tensor : replaced_sub_tensors_ket) replaced_sub_tensors_bra.push_back(dag(site_tensor).prime(Link));
            for (const auto &bond_tensor : replaced_bond_tensors_ket) replaced_bond_tensors_bra.push_back(dag(bond_tensor).prime());
            return expect_val_from_replaced_tensors({replaced_sub_tensors_ket,replaced_sub_tensors_bra},{replaced_bond_tensors_ket,replaced_bond_tensors_bra});
        }

    private:
        std::vector<int> contract_seq_;
        //Cutting tensors can be decomposed to two part. sub_tensors_ are tensors with physical legs, and left_tensors_ are tensors left without physical legs
        std::vector<TensorT> cutting_site_tensors_, sub_tensors_, left_tensors_;
        //cutting_bond_tensors_ are bulk bonds. We always absorb boundary bonds to cutting_site_tensors_
        std::vector<TensorT> cutting_bond_tensors_;
        //env_tensors_ stores MPO environment
        std::vector<TensorT> env_tensors_;
        //sub_env_tensor stores effective environment for sub_tensors (includes left_tensors_ to env_tensors)
        TensorT sub_env_tensor_, RDM_;
        Complex wf_norm_sq_;
};

double heisenberg_energy_from_RDM(const PEPSt_Subtensor_RDM<IQTensor> &peps_rdm);

#endif
