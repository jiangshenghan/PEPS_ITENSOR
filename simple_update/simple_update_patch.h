
#ifndef _SIMPLE_UPDATE_PATCH_H_
#define _SIMPLE_UPDATE_PATCH_H_

#include "simple_update.h"

//
//this class is used for obtain reduced density matrix (RDM) from a small patch of tensors with (idenical) env_matrix acting on boundary legs
//Here, we store sites of patch in order by a matrix of int, so the first and last rows as well as first and last cols are boundary sites, which should multiply E on their down(1st row), up(last row), left(1st col), right(last col) legs
//The convention is the boundary tensors are always site tensors rather than bond tensors, and env_tens is formed by In leg s and Out leg s'
//the patch is at least 2 by 2 size
//we assume the two sites connected by horizontal bond
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
class Square_Patch_RDM
{
    public:
        //
        //Constructor
        //
        Square_Patch_RDM(const IQPEPS &square_peps, const IQTensor &env_tens, const std::vector<std::vector<int>> &patch_sites, const std::array<int,2> cutting_sites);

        //
        //Acess Methods
        //
        const IQTensor &two_sites_RDM() const 
        { 
            return two_sites_RDM_;
        }

        //
        //Reduce Density Matrix and related
        //
        void obtain_two_sites_RDM();
        //Obtain expectation value with tensor at cutting sites replaced by some other tensors Si. 
        //Notice there will be no bond tensors between the replaced tensors. they can be directly contracted
        //However, the other legs should be kept the same as legs in original peps 
        //unlike RDM, the physical legs will be contracted
        //
        //   --S2-S3--
        //  /  |   |  \
        // E   |   |   E
        // |   |   |   |
        // E   |   |   E
        //  \  |   |  /
        //   --S0-S1--
        //
        //
        Complex expect_val_from_replaced_tensors(std::array<IQTensor,2> replaced_tensors, std::array<IQTensor,2> replaced_tensors_dag);

        //
        //Constructor Helpers
        //
        //obtain coordinate of cutting_site in the patch
        void init_patch_cutting_coords();

        //init patch tensors (dag) using peps tensors
        //every tensor is combined of site tens and bond tens
        //we take care of boundary tens in function modify_boundary_patch_tensors
        void init_patch_tensors();
        //get boundary patch_tensors and patch_tensors_dag by absorbing env matrix on boundary leg
        //after applying this function, boundary legs of patch_site_tensors and patch_site_tensors_dag can be directly contracted
        void modify_boundary_patch_tensors();

        //init left,down and physical legs combiners
        void init_legs_combiners();

        //init double layer tensors 
        void init_patch_double_layer_tensors();



    private:
        //patch_dim_[0] for total row # and patch_dim_[1] for total col #
        std::array<int,2> patch_dim_;
        //patch_sites_[i][j] stores site # in row i and col j of patch
        std::vector<std::vector<int>> patch_sites_;
        std::array<int,2> cutting_sites_;
        //patch_cutting_coords_ are position of cutting_sites_ at patch
        std::array<std::array<int,2>,2> patch_cutting_coords_;
        IQTensor env_tens_;
        IQPEPS square_peps_;
        //tensors multiply with env_tens_ on boundary legs
        //right and up bond tensors are absorbed
        //for tensors at cutting_sites_, physical legs are different
        std::vector<std::vector<IQTensor>> patch_tensors_, patch_tensors_dag_, patch_double_layer_tensors_;
        std::vector<std::vector<IQCombiner>> left_legs_combiners_, down_legs_combiners_;
        std::vector<IQCombiner> phys_legs_combiners_;
        IQTensor two_sites_RDM_;
};

//get tensors multiply env tens but leaves the original indices invariant
//we always assume contraction of noprime leg with boundary leg
void obtain_env_dressed_tensor(IQTensor &dressed_tens, const IQTensor &env_tens, const IQIndex &boundary_leg);

/*
//here we use small patch to do simple update
void spin_square_peps_patch_simple_update(IQPEPS &square_peps, const Evolution_Params &square_su_params, std::vector<std::vector<int>> patch_sites, std::array<int,2> evolved_sites={0,1});


//obtain leg gate params with a given two_sites_RDM
bool obtain_spin_sym_leg_gates_params_minimization_from_RDM(const IQTensor &two_sites_RDM, const Trotter_Gate &trotter_gate, const std::array<Singlet_Tensor_Basis,2> &leg_gates_basis, std::vector<double> &leg_gate_params, double cutoff=1E-5);
*/


//measure heisenberg energy using two sites RDM
double heisenberg_energy_from_RDM(const IQTensor &two_sites_RDM);

#endif
