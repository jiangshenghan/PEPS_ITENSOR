
//This file provide general class for reduced update algorithm

#ifndef _REDUCED_UPDATE_H_
#define _REDUCED_UPDATE_H_

#include "full_update.h"


//
//reduced update for kagome normal geometry
//
struct Wf_Distance_Params_RU
{
    public:
        Wf_Distance_Params_RU(Complex wf_norm_sq_evolved, PEPSt_Subtensor_RDM<IQTensor> &patch_RDM, const std::vector<IQTensor> &sub_tensors_evolved, const IQTensor &bond_tensor_evolved, const std::vector<Singlet_Tensor_Basis> &legs_basis): 
            evolved_wf_norm_sq(wf_norm_sq_evolved), 
            rdm(patch_RDM), 
            evolved_sub_tensors(sub_tensors_evolved), 
            evolved_bond_tensor(bond_tensor_evolved), 
            leg_gates_basis(legs_basis) {}

        Complex evolved_wf_norm_sq;
        PEPSt_Subtensor_RDM<IQTensor> &rdm;
        std::vector<IQTensor> evolved_sub_tensors;
        IQTensor evolved_bond_tensor;
        std::vector<Singlet_Tensor_Basis> leg_gates_basis; 
};

void kagome_normal_rvb_reduced_update(IQPEPS &kagome_rvb, const Evolution_Params &su_params);

bool obtain_kagome_normal_leg_gate_params_minimization(PEPSt_Subtensor_RDM<IQTensor> &kagome_normal_rdm, const IQTPO &evolve_gate, const std::vector<Singlet_Tensor_Basis> &leg_gates_basis, std::vector<double> &leg_gate_params, double cutoff=1e-5);

double kagome_normal_wf_distance_sq_ru_f(const gsl_vector *x, void *params);
void kagome_normal_wf_distance_sq_ru_df(const gsl_vector *x, void *params, gsl_vector *df);
void kagome_normal_wf_distance_sq_ru_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);

#endif
