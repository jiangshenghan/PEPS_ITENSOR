
#ifndef _SIMPLE_UPDATE_ENV_H_
#define _SIMPLE_UPDATE_ENV_H_

#include "peps.h"

//update environment for symmetric tensor where all legs have no extra deg
//In this case, we can simply transfer the two tensors to two matrix:
//          ---A--B---
//Then QR(LQ) decomposition is trivial in the non-deg case 
//Namely, diagonal elements of R(L) are just norms of cols of A (rows of B)
//and the non-diagonal elements of R(L) vanishes
//the env matrix = R.L, and we store it on tensors formed by common index shared by A and B (with some suitable prime)
const IQTensor &nondeg_spin_sym_env(const IQTensor &tens_A, const IQTensor &tens_B);

#endif
