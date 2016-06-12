
#ifndef _PEPS_ITEBD_H_
#define _PEPS_ITEBD_H_

#include "imps.h"
#include "peps.h"

//
//this class is to obtain effective environment of PEPS from itebd method
//TODO: Now, we only implement two-sites env. Consider env for more sites
//
template <class TensorT>
class PEPSt_iTEBD
{
    public:
        //
        //Constructor
        //
        PEPSt_iTEBD(const std::vector<int> &cut_sites_, const Tnetwork_Storage<TensorT> &peps_storage, const Args &args);

        //
        //Methods
        //
        

    private:
        std::vector<int> cut_sites_;
        //peps_storage_ are imposed in periodic BC. We use one u.c. (or two u.c. for pi-flux case) to generate iPEPS
        const Tnetwork_Storage<TensorT> &peps_storage_;
        //stores contraction results of left/right cols
        std::vector<DL_iMPS<TensorT>> peps_imps_;
        std::vector<TensorT> env_tensors_;
        //itebd_options:
        //getInt: Maxm 
        //getRead: Cutoff
        //getString: ContractMethod(single_layer)
        Args itebd_opts_; 
}

void update_peps_imps(int stepno)

#endif
