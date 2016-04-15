
//This file provide general class for various update algorithm

#ifndef _FULL_UPDATE_H_
#define _FULL_UPDATE_H_

#include "peps_rdm.h"
#include "full_update_itebd.h"

//
//fast full update for kagome cirac geometry
//
//provide params for wf_distance_params fo kagome_cirac 
struct Kagome_Cirac_Wf_Distance_Params_FU
{
    public:
        Kagome_Cirac_Wf_Distance_Params_FU(double wf_norm, PEPSt_RDM<IQTensor> &kagome_patch_RDM, const std::vector<IQTensor> &site_tensors_evolved, const IQTensor &bond_tensor_evolved, const std::vector<IQCombiner> &legs_combiners, const std::vector<Singlet_Tensor_Basis> &legs_basis, const Singlet_Tensor_Basis &bond_basis, const std::vector<int> &params_index): 
            evolved_wf_norm(wf_norm), 
            kagome_cirac_rdm(kagome_patch_RDM), 
            evolved_site_tensors(site_tensors_evolved), 
            evolved_bond_tensor(bond_tensor_evolved), 
            evolve_legs_combiners(legs_combiners), 
            leg_gates_basis(legs_basis), 
            bond_tensor_basis(bond_basis),
            bond_params_index(params_index) {}

        double evolved_wf_norm;
        PEPSt_RDM<IQTensor> &kagome_cirac_rdm;
        std::vector<IQTensor> evolved_site_tensors;
        IQTensor evolved_bond_tensor;
        std::vector<IQCombiner> evolve_legs_combiners;
        std::vector<Singlet_Tensor_Basis> leg_gates_basis; 
        Singlet_Tensor_Basis bond_tensor_basis;
        std::vector<int> bond_params_index;
};

//fast full update algorithm
void kagome_cirac_rvb_fast_full_update(IQPEPS &kagome_rvb, const Evolution_Params &su_params, const std::array<double,2> &bond_param_norms={1,1});

//obtain MPO form of cluster env tens, where the cluster is formed by three honeycomb
//std::vector<IQTensor> cluster_env_tensors_for_kagome_cirac_geometry(const std::vector<int> &cutting_sites, const std::vector<int> &cutting_bonds, const PEPSt<IQTensor> &kagome_rvb, std::array<std::vector<IQTensor>,2> &env_mats);

bool obtain_kagome_cirac_site_leg_gate_bond_params_minimization(PEPSt_RDM<IQTensor> &kagome_cirac_rdm, const IQTPO &evolve_gate, const std::vector<Singlet_Tensor_Basis> &leg_gates_basis, const Singlet_Tensor_Basis &cutting_bond_basis, std::vector<double> &leg_gate_params, std::vector<double> &bond_params, const std::vector<int> &bond_params_index, double cutoff=1e-5);

double kagome_cirac_wf_distance_sq_fu_f(const gsl_vector *x, void *params);
void kagome_cirac_wf_distance_sq_fu_df(const gsl_vector *x, void *params, gsl_vector *df);
void kagome_cirac_wf_distance_sq_fu_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);
void kagome_cirac_wf_distance_sq_fu_df_check(const gsl_vector *x, void *params, gsl_vector *df);


//
//fast full update for kagome normal geometry
//
//
struct Kagome_Normal_Wf_Distance_Params_FU
{
    public:
        Kagome_Normal_Wf_Distance_Params_FU(Complex wf_norm_sq_evolved, PEPSt_RDM<IQTensor> &kagome_patch_RDM, const std::vector<IQTensor> &site_tensors_evolved, const IQTensor &bond_tensor_evolved, const std::vector<IQCombiner> &legs_combiners, const std::vector<Singlet_Tensor_Basis> &legs_basis): 
            evolved_wf_norm_sq(wf_norm_sq_evolved), 
            kagome_normal_rdm(kagome_patch_RDM), 
            evolved_site_tensors(site_tensors_evolved), 
            evolved_bond_tensor(bond_tensor_evolved), 
            evolve_legs_combiners(legs_combiners), 
            leg_gates_basis(legs_basis) {}

        Complex evolved_wf_norm_sq;
        PEPSt_RDM<IQTensor> &kagome_normal_rdm;
        std::vector<IQTensor> evolved_site_tensors;
        IQTensor evolved_bond_tensor;
        std::vector<IQCombiner> evolve_legs_combiners;
        std::vector<Singlet_Tensor_Basis> leg_gates_basis; 
};

//fast full update algorithm (can also be used to simple update)
void kagome_normal_rvb_fast_full_update(IQPEPS &kagome_rvb, const Evolution_Params &su_params);

bool obtain_kagome_normal_leg_gate_params_minimization(PEPSt_RDM<IQTensor> &kagome_normal_rdm, const IQTPO &evolve_gate, const std::vector<Singlet_Tensor_Basis> &leg_gates_basis, std::vector<double> &leg_gate_params, double cutoff=1e-5);

double kagome_normal_wf_distance_sq_fu_f(const gsl_vector *x, void *params);
void kagome_normal_wf_distance_sq_fu_df(const gsl_vector *x, void *params, gsl_vector *df);
void kagome_normal_wf_distance_sq_fu_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);

//get MPO environment from simple update env for kagome normal lattice with two sites rdm
//env_option labels different configuration:
// env_option=0 : simple update env
// env_option=1 : triangle shape
// env_option=2 : large tree shape
void obtain_kagome_normal_env_MPO(int env_option, std::vector<int> cutting_sites, std::vector<int> cutting_bonds, const std::array<std::vector<IQTensor>,2> &su_env_mats, const IQPEPS &kagome_rvb, std::vector<IQTensor> &env_tensors);

//get tensors from contracting double layer site tensors and su_env_mat
IQTensor env_site_combined_tensor(const IQTensor &env_mat, const IQTensor &site_tensor, const std::vector<IQIndex> &indices_contract, bool contract_phys_ind=true);

//obtain double layer tensor as
//      \ /
//       v
//      / \
//     /   \
//    u-----w
//where other legs of u and w are contracted by env_mat
//order of u,v,w can be changed, the last site_tensors served as v
//Notice we should include bond tensors in prior
IQTensor env_kagome_normal_uc_combined_tensor(const IQTensor &env_mat, const std::vector<IQTensor> &site_tensors, int bulk_no);

#endif
