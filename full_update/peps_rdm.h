
#ifndef _PEPS_RDM_H_
#define _PEPS_RDM_H_

#include "square_rvb.h"
#include "kagome_rvb.h"
#include "trotter_gate.h"
#include "simple_update_env.h"

//this class stores reduced density matrix used for measurement and full update, where env is approximate by matrix product state
template <class TensorT>
class PEPSt_RDM
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
        PEPSt_RDM(const std::string &name, const std::vector<int> &cutting_sites, const std::vector<int> &cutting_bonds, const std::vector<TensorT> &env_tensors, const PEPSt<TensorT> &peps, std::vector<int> env_contract_seq=std::vector<int>());

        //
        //Access Method
        //
        const std::string &peps_name() const { return peps_.name(); }

        int cutting_sites_no() const { return cutting_sites_.size(); }
        const std::vector<int> &cutting_sites() const { return cutting_sites_; }
        const std::vector<int> &cutting_bonds() const { return cutting_bonds_; }

        const IndexT &cutting_phys_legs(int cuti) const { return peps_.phys_legs(cutting_sites_[cuti]); }

        std::vector<TensorT> cutting_site_tensors() const
        {
            std::vector<TensorT> site_tensors;
            for (int sitei : cutting_sites_) site_tensors.push_back(peps_.site_tensors(sitei));
            return site_tensors;
        }
        const TensorT &cutting_site_tensors(int cuti) const { return peps_.site_tensors(cutting_sites_[cuti]); }

        std::vector<TensorT> cutting_bond_tensors() const
        {
            std::vector<TensorT> bond_tensors;
            for (int bondi : cutting_bonds_) bond_tensors.push_back(peps_.bond_tensors(bondi));
            return bond_tensors;
        }
        const TensorT &cutting_bond_tensors(int cuti) const { return peps_.bond_tensors(cutting_bonds_[cuti]); }

        const TensorT &RDM() const { return RDM_; }

        Complex wf_norm_sq() const { return wf_norm_sq_; }

        //
        //update peps_rdm
        //
        void update_peps_rdm(const std::vector<IQTensor> &env_tensors, const PEPSt<TensorT> &peps)
        {
            peps_=peps;
            env_tensors_=env_tensors;
            obtain_RDM_and_wf_norm();
        }
        void update_peps_rdm(const std::vector<IQTensor> &env_tensors, const PEPSt<TensorT> &peps, std::vector<int> env_contract_seq)
        {
            env_contract_seq_=env_contract_seq;
            peps_=peps;
            env_tensors_=env_tensors;
            obtain_RDM_and_wf_norm();
        }

        //
        //Reduced Density Matrix and related
        //
        void obtain_RDM_and_wf_norm();

        Complex expect_val_from_replaced_tensors(std::array<std::vector<TensorT>,2> replaced_site_tensors_ket_bra, std::array<std::vector<TensorT>,2> replaced_bond_tensors_ket_bra);
        
        Complex expect_val_from_replaced_tensors(std::array<std::vector<TensorT>,2> replaced_site_tensors_ket_bra)
        {
            return expect_val_from_replaced_tensors(replaced_site_tensors_ket_bra,{cutting_bond_tensors(),cutting_bond_tensors()});
        }

        Complex expect_val_from_replaced_tensors(std::vector<TensorT> replaced_site_tensors_ket, std::vector<TensorT> replaced_bond_tensors_ket)
        {
            return expect_val_from_replaced_tensors({replaced_site_tensors_ket,replaced_site_tensors_ket},{replaced_bond_tensors_ket,replaced_bond_tensors_ket});
        }

        Complex expect_val_from_replaced_tensors(std::vector<TensorT> replaced_site_tensors_ket)
        {
            return expect_val_from_replaced_tensors({replaced_site_tensors_ket,replaced_site_tensors_ket});
        }

        //contract all tensors as well as env tensors
        TensorT tensor_from_contract_patch(const std::array<std::vector<TensorT>,2> &tensors_ket_bra);
        TensorT tensor_from_contract_patch(const std::vector<TensorT> &tensors_ket, bool contract_phys_leg=true)
        {
            std::vector<TensorT> tensors_bra;
            for (const auto &tensor_ket : tensors_ket)
            {
                if (contract_phys_leg)
                {
                    tensors_bra.push_back(dag(tensor_ket).prime(Link));
                }
                else
                {
                    tensors_bra.push_back(dag(prime(tensor_ket)));
                }
            }
            return tensor_from_contract_patch({tensors_ket,tensors_bra});
        }

        //combine site and bond tensors
        void combine_site_bond_tensors(const std::vector<TensorT> &site_tensors, const std::vector<TensorT> &bond_tensors, std::vector<TensorT> &combined_tensors);


    private:
        //name_ characterize the shape of cutting sites
        //1. Square lattice:
        //
        //2. Kagome cirac lattice:
        //  a. tree shape I: three cutting site tensors and one cutting plaquette tensor
        //     |
        //     1
        //     |
        //     3
        //    / \
        //   0   2
        //  /     \
        //
        //3. kagome normal lattice
        //  a. two sites shape: two cutting site tensors (u,v) and three cutting bond tensors
        //    |         
        //  B0|         
        //    |   B1    |   B4
        //  --+-- --- --+-- ---
        //    |         |
        std::string name_;
        std::vector<int> cutting_sites_, cutting_bonds_;
        //env_tensors_ are MPS whose site legs are the boundary legs of cutting sites
        std::vector<TensorT> env_tensors_;
        //env_contract_seq_[i]==0 means contract to single layer
        //env_contract_seq_[i]==1 means contract to double layer
        std::vector<int> env_contract_seq_;
        PEPSt<TensorT> peps_;
        //RDM_ are obtained by contracting env_tensors_ and double layer peps tensors at cutting sites
        TensorT RDM_;
        //wf_norm_sq_ is obtained by contracting phys_legs of RDM_
        //Notice, wf_norm_sq_ may be complex number
        Complex wf_norm_sq_;

};


double heisenberg_energy_from_RDM(PEPSt_RDM<IQTensor> peps_rdm);

#endif

