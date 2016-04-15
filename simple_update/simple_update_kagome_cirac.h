
#ifndef _SIMPLE_UPDATE_KAGOME_CIRAC_H_
#define _SIMPLE_UPDATE_KAGOME_CIRAC_H_

#include "simple_update_patch_general.h"

//provide params for wf_distance_params fo kagome_cirac 
struct Kagome_Cirac_Wf_Distance_Params
{
    public:
        Kagome_Cirac_Wf_Distance_Params(double wf_norm, const General_Patch_RDM<IQTensor> &patch_RDM, const std::vector<IQTensor> &site_tensors_evolved, const IQTensor &bond_tensor_evolved, const std::vector<IQCombiner> &legs_combiners, const std::array<std::vector<Singlet_Tensor_Basis>,2> &basis): 
            evolved_wf_norm(wf_norm), 
            kagome_patch_RDM(patch_RDM), 
            evolved_site_tensors(site_tensors_evolved), 
            evolved_bond_tensor(bond_tensor_evolved), 
            evolve_legs_combiners(legs_combiners), 
            leg_gates_basis(basis) {}

        double evolved_wf_norm;
        General_Patch_RDM<IQTensor> kagome_patch_RDM;
        std::vector<IQTensor> evolved_site_tensors;
        IQTensor evolved_bond_tensor;
        std::vector<IQCombiner> evolve_legs_combiners;
        std::array<std::vector<Singlet_Tensor_Basis>,2> leg_gates_basis; 
};

struct Kagome_Cirac_Wf_Distance_Params_NBLG
{
    public:
        Kagome_Cirac_Wf_Distance_Params_NBLG(double wf_norm, const General_Patch_RDM<IQTensor> &patch_RDM, const std::vector<IQTensor> &site_tensors_evolved, const IQTensor &bond_tensor_evolved, const std::vector<IQCombiner> &legs_combiners, const std::vector<Singlet_Tensor_Basis> &legs_basis, const Singlet_Tensor_Basis &bond_basis, const std::vector<int> &params_index): 
            evolved_wf_norm(wf_norm), 
            kagome_patch_RDM(patch_RDM), 
            evolved_site_tensors(site_tensors_evolved), 
            evolved_bond_tensor(bond_tensor_evolved), 
            evolve_legs_combiners(legs_combiners), 
            leg_gates_basis(legs_basis), 
            bond_tensor_basis(bond_basis),
            bond_params_index(params_index) {}

        double evolved_wf_norm;
        General_Patch_RDM<IQTensor> kagome_patch_RDM;
        std::vector<IQTensor> evolved_site_tensors;
        IQTensor evolved_bond_tensor;
        std::vector<IQCombiner> evolve_legs_combiners;
        std::vector<Singlet_Tensor_Basis> leg_gates_basis; 
        Singlet_Tensor_Basis bond_tensor_basis;
        std::vector<int> bond_params_index;
};



//patch simple update for kagome cirac PEPS
void spin_kagome_cirac_peps_patch_simple_update(IQPEPS &kagome_rvb, const Evolution_Params &su_params, std::vector<int> patch_sites, std::vector<int> evolved_sites, std::string patch_name, std::array<double,2> bond_param_norms={1.,1.});
//update kagome cirac PEPS without plaquette leg gates
void spin_kagome_cirac_peps_patch_simple_update_no_bond_leg_gates(IQPEPS &kagome_rvb, const Evolution_Params &su_params, std::vector<int> patch_sites, std::vector<int> evolved_sites, std::string patch_name, std::array<double,2> bond_param_norms={1.,1.});

//obtain leg gates params for kagome cirac peps (both site and plaquette)
//0/1 labels site leg gates and plaquette leg_gates
bool obtain_kagome_cirac_leg_gates_params_minimization(General_Patch_RDM<IQTensor> &kagome_patch_RDM, const IQTPO &evolve_gate, const std::array<std::vector<Singlet_Tensor_Basis>,2> &leg_gates_basis, std::array<std::vector<double>,2> &leg_gates_params, double cutoff=1E-5);
//obtain site leg gate and bond tensors params
bool obtain_kagome_cirac_site_leg_gate_bond_params_minimization(General_Patch_RDM<IQTensor> &kagome_patch_RDM, const IQTPO &evolve_gate, const std::vector<Singlet_Tensor_Basis> &leg_gates_basis, const Singlet_Tensor_Basis &comm_bond_basis, std::vector<double> &leg_gate_params, std::vector<double> &bond_params, const std::vector<int> &bond_params_index, double cutoff=1E-5);

//the following functions provides distance_sq for kagome cirac 
double kagome_cirac_wf_distance_sq_f(const gsl_vector *x, void *params);
//TODO:improve the numerical derivative?
void kagome_cirac_wf_distance_sq_df(const gsl_vector *x, void *params, gsl_vector *df);
void kagome_cirac_wf_distance_sq_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);

double kagome_cirac_wf_distance_sq_nblg_f(const gsl_vector *x, void *params);
void kagome_cirac_wf_distance_sq_nblg_df(const gsl_vector *x, void *params, gsl_vector *df);
void kagome_cirac_wf_distance_sq_nblg_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);
void kagome_cirac_wf_distance_sq_nblg_df_check(const gsl_vector *x, void *params, gsl_vector *df);

#endif
