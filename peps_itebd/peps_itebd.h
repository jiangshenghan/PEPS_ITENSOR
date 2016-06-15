
#ifndef _PEPS_ITEBD_H_
#define _PEPS_ITEBD_H_

#include "imps.h"
#include "peps.h"

//
//this class is to obtain effective environment of PEPS from itebd method
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
        DL_iMPS<TensorT> left_dl_imps_, right_dl_imps;
        std::vector<TensorT> env_tensors_;
        //itebd_options:
        //getInt: Maxm 
        //getReal: Cutoff
        //getString: ContractMethod(single_layer), ContractGeometry(normal)
        Args itebd_opts_; 
};


#endif
