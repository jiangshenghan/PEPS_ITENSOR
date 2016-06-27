
#ifndef _PEPS_ITEBD_H_
#define _PEPS_ITEBD_H_

#include "imps.h"
#include "peps.h"

//
//this class is to obtain two_sites effective environment of PEPS from itebd method
//
template <class TensorT>
class PEPSt_iTEBD
{
    public:
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;

        //
        //Constructor
        //
        PEPSt_iTEBD(const Tnetwork_Storage<TensorT> &peps_storage, const Args &itebd_opts);

        //
        //Methods
        //
        //obtain the two sites env using left and right imps
        //1. for kagome lattice, env tensors are
        //      4
        //      |
        //   0--v--2
        //      |
        //   1--u--3
        //      |
        //      5
        //
        const std::vector<TensorT> &env_tensors_from_itebd(int n_steps);
        //update imps by contracting cols of impos
        void update_imps_one_step();
        //obtain expectation value of two-site mpo operator
        Complex expect_val_from_env_tensors(const std::vector<TensorT> &two_sites_mpo) const;

        friend void contract_dl_impo_imps(DL_iMPSt<TensorT> &dl_imps, const DL_iMPOt<TensorT> &dl_impo, const Args &contract_opts);

        //
        //Acess Method
        //
        const std::vector<TensorT> &env_tensors() const { return env_tensors_; }
        
        //
        //Constructor Helpers
        //
        void init_impo();


    private:
        //peps_storage_ are imposed in periodic BC. We use one u.c. (or two u.c. for pi-flux case) to generate iPEPS
        const Tnetwork_Storage<TensorT> &peps_storage_;
        //stores contraction results of left/right cols
        DL_iMPSt<TensorT> ldl_imps_, rdl_imps_;
        std::vector<DL_iMPOt<TensorT>> lcols_dl_impos_, rcols_dl_impos_;
        std::vector<TensorT> env_tensors_;
        //itebd_options:
        //getInt: Maxm(for svd), MaxIter(for arnoldi), MaxRestart(for arnoldi)
        //getReal: Cutoff(for svd), ErrGoal(for arnoldi)
        Args itebd_opts_; 

};



#endif
