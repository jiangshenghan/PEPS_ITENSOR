
#ifndef _SIMPLE_UPDATE_ENV_H_
#define _SIMPLE_UPDATE_ENV_H_

#include "peps.h"
//#include "square_rvb.h"

//Params to do obtain optimal env_tens
//1. spin and flavor inf of virtual leg
//2. two "site" tensors to do SVD
//3. env_tens (we only need the indices)
struct Env_Tens_Params
{
    public:
        std::vector<int> flavor_deg;
        std::vector<int> flavor_accumulate_deg;
        std::vector<Spin_Basis> spin_basis;
        std::array<IQTensor,2> site_tens;
        std::array<std::vector<IQTensor>,2> env_tens;
};

//solve the env tensor for simple update by iterative method
//1. input tensors and initial env_tensor, env_tensor are tensors with two indices (s,s')
//2. multiply them, get combined tensors
//3. using combined tensors to obtain new env tensor
//4. view preivous steps as a function with input tensors as params, diag elements of env_tensor as auguments and f=updated_env-env, solve f=0
void get_env_tensor_iterative(const IQTensor &site_tensA, const IQTensor &site_tensB, std::array<std::vector<IQTensor>,2> &env_tens);

//TODO: consider flavor deg case
//solve the env tensor by minimization method
void get_env_tensor_minimization(const IQTensor &site_tensA, const IQTensor &site_tensB, std::array<std::vector<IQTensor>,2> &env_tens);


//init the env_tensors according to their spin and flavor qn
void init_env_tensor(const IQTensor &site_tensA, const IQTensor &site_tensB, std::array<std::vector<IQTensor>,2> &env_tens);

//update environment for symmetric tensor where all legs have no extra deg
//In this case, we can simply transfer the two tensors to two matrix:
//          ---A--B---
//Notice that we already include the previous environment matrix in A & B
//Then QR(LQ) decomposition is trivial in the non-deg case 
//Namely, diagonal elements of R(L) are just norms of cols of A (rows of B)
//and the non-diagonal elements of R(L) vanishes
//the env matrix = R.L is a diagonal matrix. We stores the diagonal part on a vector
std::vector<double> nondeg_spin_sym_env_updated(const IQTensor &tens_A, const IQTensor &tens_B);


//transfer between env_elems and env_tens
//env_elems is stored according to flavor and spin qns 
std::vector<double> env_elems_from_env_tens(const std::vector<int> &flavor_accumulate_deg, const std::vector<Spin_Basis> &spin_basis, const std::array<std::vector<IQTensor>,2> &env_tens);
void obtain_env_tens_from_env_elems(const std::vector<int> &flavor_accumulate_deg, const std::vector<Spin_Basis> &spin_basis, const std::vector<double> env_elems, std::array<std::vector<IQTensor>,2> &env_tens);

//The following provides f,df,fdf of gsl_minimization
double updated_env_tens_diff_f(const gsl_vector *x, void *params);
void updated_env_tens_diff_df(const gsl_vector *x, void *params, gsl_vector *df);
void updated_env_tens_diff_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);

#endif
