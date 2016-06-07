
#include "imps.h"

template <class TensorT>
iMPSt<TensorT>::iMPSt(const std::vector<TensorT> &init_site_tensors, const std::vector<TensorT> &init_bond_tensors, TensorT &init_VL, TensorT &init_VR):
    n_sites_uc_(site_tensors.size()),
    site_tensors_(init_site_tensors),
    bond_tensors_(init_bond_tensors),
    VL_(init_VL),
    VR_(init_VR)
{
    init_options();
    //TODO: consider the case where VL_/VR_ are empty tensors
    obtain_dominant_eigensystems();
}


template <class TensorT>
void iMPSt<TensorT>::init_options()
{
    //init arnoldi_args_
    arnoldi_args_.add("MaxIter",10);
    arnoldi_args_.add("MaxRestart",1);
    arnoldi_args_.add("ErrGoal",1e-16);

    //TODO: init svd_args_
}
template
void iMPSt<ITensor>::init_options();
template 
void iMPSt<IQTensor>::init_options();

template <class TensorT>
void iMPSt<TensorT>::obtain_dominant_eigensystems()
{
    TensorT Gamma_tensor=site_tensors_[0];
    for (int sitei=1; sitei<n_sites_uc_; sitei++)
    {
        Gamma_tensor*=bond_tensors_[sitei-1]*site_tensors[1];
    }
    TensorT Gamma_tensor_dag=Gamma_tensor;
    Gamma_tensor_dag.prime(Link).dag();

    TensorT_Matrix L_tensor(Gamma_tensor*Gamma_tensor_dagbond_tensors_[n_sites_uc_],)
    eta_L_=arnoldi();
}
template
void iMPSt<ITensor>::obtain_dominant_eigensystems();
template
void iMPSt<IQTensor>::obtain_dominant_eigensystems();


