
#ifndef _SIMPLE_UPDATE_H_
#define _SIMPLE_UPDATE_H_

#include "trotter_gate.h"
#include "simple_update_env.h"
#include "square_rvb.h"
#include "kagome_rvb.h"

//provide struct for wavefunction overlap params (M_ijkl=\langle\phi_ij|\phi_lk\rangle and w_ij=\langle\phi_ij|\psi\rangle+\langle\psi|\phi_ij\rangle, as well as N_leg_basis) to do multiminimization
struct Wf_Distance_Params
{
    public:
        Wf_Distance_Params(int N, double wf_norm, std::vector<double> M_init, std::vector<double> w_init): N_leg_basis(N), evolved_wf_norm(wf_norm), M(M_init), w(w_init) {}
        int N_leg_basis;
        double evolved_wf_norm;
        std::vector<double> M;
        std::vector<double> w;
};


//Using simple update algorithm to obtain spin symmetric peps with the optimal energy. The algorithm is as follows
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
void obtain_spin_sym_leg_gates_params_iterative(const std::array<IQTensor,2> &site_tensors, const IQTensor &bond_tensor, const Trotter_Gate &trotter_gate, const std::array<Singlet_Tensor_Basis,2> &leg_gates_basis, std::vector<double> &leg_gate_params);

//Obtain the leg gate by miinmizing function d=||\phi-\psi||, where \psi is the evolved two-site wavefunction and \phi is obtained from applying leg gates to site tensors
//Notice, we already include the environment in the site_tensors
bool obtain_spin_sym_leg_gates_params_minimization(const std::array<IQTensor,2> &site_tensors, const IQTensor &bond_tensor, const Trotter_Gate &trotter_gate, const std::array<Singlet_Tensor_Basis,2> &leg_gates_basis, std::vector<double> &leg_gate_params, double cutoff=1E-5);

//The following provides f,df,fdf of gsl_minimization
double wf_distance_f(const gsl_vector *x, void *params);
void wf_distance_df(const gsl_vector *x, void *params, gsl_vector *df);
void wf_distance_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);
//check wf_distance_func
void wf_distance_func_check(const std::array<IQTensor,2> &site_tensors, const IQTensor &bond_tensor,const Trotter_Gate &trotter_gate, const std::array<Singlet_Tensor_Basis,2> &leg_gates_basis, std::vector<double> &leg_gate_params);

//measure heisenberg energy using two sites, env_tens already included in tensA and tensB
double heisenberg_energy_from_site_env_tensors(const std::array<IQTensor,2> &site_env_tens, const IQTensor &comm_bond_tensor, const NN_Heisenberg_Hamiltonian &hamiltonian_gate);
double heisenberg_energy_from_site_env_tensors(const std::array<IQTensor,2> &site_env_tens, const NN_Heisenberg_Hamiltonian &hamiltonian_gate);

#endif
