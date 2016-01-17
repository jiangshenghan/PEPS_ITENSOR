
#ifndef _SIMPLE_UPDATE_PATCH_GENERAL_H_
#define _SIMPLE_UPDATE_PATCH_GENERAL_H_

#include "simple_update_patch_square.h"

//
//this class is used for obtain reduced density matrix (RDM) from a small patch of tensors with (idenical) env_matrix acting on boundary legs
//The convention is the boundary tensors are always SITE tensors rather than bond tensors, and env_tens is formed by leg s and leg s' with different direction
//
//   --T---T--
//  /  |   |  \
// E           E
// |           |
// E           E
//  \  |   |  /
//   --T---T--
//
//
//The following class is Patch_RDM of general case (not limited to square patch of square lattice)
//The detail algorithm is as follows
//1. obtain single layer tensors with env-mat multiplied on boundary legs 
//2. obtain double layer tensors
//3. separate the patch to several parts. contract every part of the patch, and finally get RDM_
//4. one can replace tensors at cutting sites to get various parameters for optimization
//However, the contraction algorithm for different geometry and patches are different
//1. Square lattice: we use the convention that two cutting sites are in horizontal direction --0---1--. We seperate the patch to left and right part, and each one includes one cutting site
//  a. regular shape: the patch_sites_ is ordered as row by row, for example, in the following case, patch_dim_={4,4,4} 
//      |8 9 A B|       Here, we assume 5,6 are two cutting sites.
//      |4 5 6 7|       We contract col by col. 
//      |0 1 2 3|
//  b. special shape I: the patch is as follows
//        |4 5|       Here, 2,3 are two cutting sites.
//      |6 2 3 7|     
//        |0 1|
//
//2. kagome cirac geometry: We consider RDM for three sites
//  a. tree shape I: includes three site tensors and one plaquette tensor
//     |
//     1
//     |
//     3
//    / \
//   0   2
//  /     \
//
template <class TensorT>
class General_Patch_RDM
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
        General_Patch_RDM(const std::string &name, PEPSt<TensorT> &peps, const TensorT &env_tens, const std::vector<int> &patch_sites, const std::vector<int> &cutting_sites);

        //
        //Acess Method
        //
        int cutting_sites_no() const
        {
            return cutting_sites_.size();
        }
        const std::vector<int> &cutting_sites() const
        {
            return cutting_sites_;
        }
        const IndexT &cutting_phys_legs(int cuti) const
        {
            return peps_.phys_legs(cutting_sites_[cuti]);
        }
        const TensorT &cutting_site_tensors(int cuti) const
        {
            return peps_.site_tensors(cutting_sites_[cuti]);
        }
        const TensorT &cutting_bond_tensor() const
        {
            int comm_bond=peps_.lattice().comm_bond(cutting_sites_);
            return peps_.bond_tensors(comm_bond);
        }
        IndexT cutting_virt_legs(int cuti) const
        {
            return commonIndex(cutting_site_tensors(cuti),cutting_bond_tensor());
        }
        const TensorT &RDM() const
        {
            return RDM_;
        }
        double patch_norm() const
        {
            return patch_norm_;
        }
        double wf_norm() const
        {
            return patch_norm_;
        }
        std::string peps_name() const
        {
            return peps_.name();
        }
        std::string patch_name() const
        {
            if (name_.find("square")!=std::string::npos) 
            {
                if (name_.find("regular shape")!=std::string::npos)
                {
                    return nameint("patch=",patch_dim_[0])+nameint("x",patch_dim_.size());
                }
                if (name_.find("special shape I")!=std::string::npos)
                {
                    return "patch=speical_shape_I";
                }
            }
            if (name_.find("kagome cirac")!=std::string::npos)
            {
                if (name_.find("tree shape I")!=std::string::npos)
                {
                    return "tree_shape_I";
                }
            }
            return "unnamed_patch";
        }

        //
        //Reduced Density Matrix and related
        //
        void obtain_RDM_and_patch_norm();
        std::vector<TensorT> double_layer_tensors_from_replaced_cutting_sites_tensors(std::vector<std::array<TensorT,2>> replaced_tensors_ket_bra); 
        TensorT tensor_from_contract_patch(const std::vector<TensorT> &double_layer_tensors);
        //for square lattice, dir=0 means left half, dir=1 right half
        TensorT tensor_from_contract_part_patch(const std::vector<TensorT> &double_layer_tensors, int parti);
        //calculate expectation value of tensO
        Complex expect_val_from_RDM(const TensorT &tensO)
        {
            return (RDM_*tensO).toComplex();
        }
        //contract network obtained by replace cutting site and bond tensors
        Complex expect_val_from_replaced_tensors(const std::vector<std::array<TensorT,2>> &replaced_tensors_ket_bra)
        {
            auto double_layer_tensors=double_layer_tensors_from_replaced_cutting_sites_tensors(replaced_tensors_ket_bra);
            return tensor_from_contract_patch(double_layer_tensors).toComplex();
        }
        //for the case ket tensor equals bra tensor (up to dagger)
        Complex expect_val_from_replaced_tensors(const std::vector<TensorT> &replaced_tensors_ket)
        {
            std::vector<std::array<TensorT,2>> replaced_tensors_ket_bra;
            for (const auto &ket_tensor : replaced_tensors_ket)
            {
                replaced_tensors_ket_bra.push_back({ket_tensor,ket_tensor});
            }
            return expect_val_from_replaced_tensors(replaced_tensors_ket_bra);
        }

    private:
        //name_ labels the lattice geometry as well as the patch geometry, which determines the order of contraction
        std::string name_;
        std::vector<int> patch_dim_, patch_sites_, cutting_sites_;
        //we should change indices of env_tens_ such to avoid wrong contraction 
        TensorT env_tens_;
        const PEPSt<TensorT> &peps_;
        //combine top and bottom legs to avoid index no overflow
        std::vector<std::vector<CombinerT>> virt_legs_combiners_;
        std::vector<CombinerT> phys_legs_combiners_;
        //all physical legs are contracted except ones on the cutting_sites_
        std::vector<TensorT> patch_double_layer_tensors_;
        //RDM_ are obtained by contracting patch_double_layer_tensors_
        //the contraction order is depend on patch geometry, and is chosen to minimize time cost
        TensorT RDM_;
        //patch_norm_ is obtained by contraction phys_legs of RDM_
        double patch_norm_;

        //
        //Constructor Helpers
        //
        void init_patch_dim();
        //after applying this function, boundary legs of patch_site_tensors and patch_site_tensors_dag can be directly contracted
        void init_env_tens();
        void init_legs_combiners();
        //TODO: for bulk sites on square lattice, it is very time costing to construct a double layer tensors with combiners (~D^10), modify algorithm in those case (for example, special shape I case)
        void init_patch_double_layer_tensors();
        //absorb env tens to boundary legs for ket and bra tensors locate at sitei
        //if the tensors is inside patch bulk, then leave it invariant
        void absorb_env_tens(int sitei, std::array<TensorT,2> &ket_bra_tensors);

        //absorb neighbouring bonds of sitei into double layer tensor, where we disgard forbid bond
        void absorb_neigh_bonds(int sitei, TensorT &double_layer_tensor, int forbid_bond=-1);
};



//we use small patch to do simple update for square PEPS
void spin_square_peps_patch_simple_update(IQPEPS &square_peps, const Evolution_Params &square_su_params, std::vector<int> patch_sites, std::array<int,2> evolved_sites, std::string patch_name);

//patch simple update for kagome cirac PEPS
void spin_kagome_cirac_peps_patch_simple_update(IQPEPS &kagome_rvb, const Evolution_Params &su_params, std::vector<int> patch_sites, std:vector<int> evolved_sites, std::string patch_name);


//obtain leg gate params of SQUARE peps with a given general RDM
bool obtain_spin_sym_leg_gates_params_minimization_from_RDM(General_Patch_RDM<IQTensor> &square_RDM, const Trotter_Gate &trotter_gate, const std::array<Singlet_Tensor_Basis,2> &leg_gates_basis, std::vector<double> &leg_gate_params, double cutoff=1E-5);
//obtain leg gates params for kagome cirac peps (both site and plaquette)
//0/1 labels site leg gates and plaquette leg_gates
bool obtain_kagome_cirac_leg_gates_params_minimization(General_Patch_RDM<IQTensor> &kagome_patch_RDM, const IQTPO &evolve_gate, const std::array<std::vector<Singlet_Tensor_Basis>,2> &leg_gates_basis, std::array<std::vector<double>,2> &leg_gates_params, double cutoff=1E-5);

//the following functions provides distance_sq for kagome cirac 
double kagome_cirac_wf_distance_f(const gsl_vector *x, void *params);
//TODO:improve the numerical derivative?
void kagome_cirac_wf_distance_df(const gsl_vector *x, void *params, gsl_vector *df);
void kagome_cirac_wf_distance_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);


//measure heisenberg energy using two sites RDM
double heisenberg_energy_from_RDM(const General_Patch_RDM<IQTensor> &patch_rdm);

#endif
