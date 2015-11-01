
#ifndef _SIMPLE_UPDATE_ENV_H_
#define _SIMPLE_UPDATE_ENV_H_

#include "peps.h"
//#include "square_rvb.h"

//function to solve the env tensor for simple update
//1. input tensors and initial env_tensor, env_tensor are tensors with two indices (s,s')
//2. multiply them, get combined tensors
//3. using combined tensors to obtain new env tensor
//4. view preivous steps as a function with input tensors as params, diag elements of env_tensor as auguments and f=updated_env-env, solve f=0
void get_env_tensor(const IQTensor &site_tensA, const IQTensor &site_tensB, std::array<std::vector<IQTensor>,2> &env_tens);

//update environment for symmetric tensor where all legs have no extra deg
//In this case, we can simply transfer the two tensors to two matrix:
//          ---A--B---
//Notice that we already include the previous environment matrix in A & B
//Then QR(LQ) decomposition is trivial in the non-deg case 
//Namely, diagonal elements of R(L) are just norms of cols of A (rows of B)
//and the non-diagonal elements of R(L) vanishes
//the env matrix = R.L is a diagonal matrix. We stores the diagonal part on a vector
std::vector<double> nondeg_spin_sym_env_updated(const IQTensor &tens_A, const IQTensor &tens_B);


#endif
