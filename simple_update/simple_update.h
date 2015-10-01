
#ifndef _SIMPLE_UPDATE_H_
#define _SIMPLE_UPDATE_H_

#include "trotter_gate.h"
#include "simple_update_env.h"

//Using simple update algorithm to obtain spin symmetric peps with the optimal energy. We list the algorithm as follows
//1. input init PEPS as well as the imaginary time evolution params. We should ensure that init PEPS satisfy all lattice symmetry as well as spin rotation symmetry
//2. Obtain the effective environment. In simple update, we appox the env by direct product operators
//3. Acting trotter gate, and approx it by leg_gates u,v
//4. Multiple leg gates on all virtual legs, and using lattice symmetry to symmetrize peps
//5. repeat step 2,3,4
//
//this is for square peps without extra deg in virtual bonds
//we use site 0 and site 1 to update
void spin_square_peps_simple_update(IQPEPS &square_peps, const Evolution_Params &square_su_params);

//Obtain the leg gate u,v to approximate trotter gate using iterative method
//Notice, we already include the environment in the site_tensors
void obtain_spin_sym_leg_gates_params(const std::array<IQTensor,2> &site_tensors, const IQTensor &bond_tensor, Trotter_Gate &trotter_gate, const std::array<Singlet_Tensor_Basis,2> &leg_gates_basis, std::vector<double> &leg_gate_params);


//Obtain the C4 symmetrized singlet tensor for certain class
//we assume there are four virtual legs for the site_tensor in this function
void C4_symmetrized_tensor(IQTensor &site_tensor);


#endif
