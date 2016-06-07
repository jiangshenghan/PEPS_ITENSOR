
#ifndef _IMPS_H_
#define _IMPS_H_

#include "utilities.h"
#include "tensor_matrix.h"

//
//class iMPSt stores infinite MPS in canonical form
//for canonical form, we have
//
// --G--l          |       l--G--         |
//   |   \  = eta. | ,    /   |    = eta. |
//   |   /         |      \   |           |
// --G--l                  l--G--          
//
//
template <class TensorT>
class iMPSt
{
    public:
        //
        //type alias
        //
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;

        //
        //Constructors
        //
        iMPSt(const std::vector<TensorT> &init_site_tensors, const std::vector<TensorT> &init_bond_tensors, TensorT &init_VL=TensorT(), TensorT &init_VR=TensorT());

        //
        //Constructor Helpers
        //
        void init_options();
        void obtain_dominant_eigensystems();

    private:
        int n_sites_uc_;
        //Here, bond_tensors_ are "Diag" matrix, which stores single values
        //Notice, bond_tensors_.size()=site_tensors.size()+1. Namely, we have
        // --bond[N]--site[0]--bond[0]--...--site[N-1]--bond[N-1]--
        std::vector<TensorT> site_tensors_, bond_tensors_;
        //VL_/VR_ contain left/right dominant eigenvectors, which should be positive hermitian
        //VL_.indices()={bond_virt_indices_[0][0],dag(prime(bond_virt_indices_[0][0]))}
        //VR_.indices()={bond_virt_indices_[0][1],dag(prime(bond_virt_indices_[0][1]))}
        TensorT VL_, VR_;
        Complex eta_L_, eta_R_;
        Args arnoldi_args_, svd_args_;
};


#endif
