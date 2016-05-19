
#include "peps_subtensor_rdm.h"

template <class TensorT>
PEPSt_Subtensor_RDM<TensorT>::PEPSt_Subtensor_RDM(const std::vector<TensorT> &cutting_site_tensors, const std::vector<TensorT> &cutting_bond_tensors, const std::vector<TensorT> &env_tensors, const std::vector<int> contract_seq):
    contract_seq_(contract_seq),
    cutting_site_tensors_(cutting_site_tensors),
    cutting_bond_tensors_(cutting_bond_tensors),
    env_tensors_(env_tensors)
{
    if (contract_seq_.empty())
        contract_seq_=std::vector<int>(env_tensors_.size(),0);
    obtain_RDM();
}
template
PEPSt_Subtensor_RDM<ITensor>::PEPSt_Subtensor_RDM(const std::vector<ITensor> &cutting_site_tensors, const std::vector<ITensor> &cutting_bond_tensors, const std::vector<ITensor> &env_tensors, const std::vector<int> contract_seq);
template
PEPSt_Subtensor_RDM<IQTensor>::PEPSt_Subtensor_RDM(const std::vector<IQTensor> &cutting_site_tensors, const std::vector<IQTensor> &cutting_bond_tensors, const std::vector<IQTensor> &env_tensors, const std::vector<int> contract_seq);


template <class TensorT>
void PEPSt_Subtensor_RDM<TensorT>::obtain_RDM()
{
    obtain_sub_env_tensor();

    //get RDM_
    RDM_=sub_env_tensor_;
    for (const auto sub_tensor : sub_tensors_)
    {
        RDM_=RDM_*sub_tensor*prime(dag(sub_tensor));
    }
    for (const auto bond_tensor : cutting_bond_tensors_)
    {
        RDM_=RDM_*bond_tensor*prime(dag(bond_tensor));
    }
    clean(RDM_);

    //get wf_norm_sq_
    auto RDM_trace=RDM_;
    for (int sitei=0; sitei<cutting_site_tensors_.size(); sitei++)
    {
        IndexT phys_leg=this->cutting_phys_legs(sitei);
        RDM_trace.trace(phys_leg,dag(prime(phys_leg)));
        //Print(RDM_trace.norm());
    }
    //PrintDat(RDM_trace);
    wf_norm_sq_=RDM_trace.toComplex();
}
template
void PEPSt_Subtensor_RDM<ITensor>::obtain_RDM();
template 
void PEPSt_Subtensor_RDM<IQTensor>::obtain_RDM();

template <class TensorT>
void PEPSt_Subtensor_RDM<TensorT>::update_RDM(const std::vector<TensorT> &cutting_site_tensors, const std::vector<TensorT> &cutting_bond_tensors, const std::vector<TensorT> &env_tensors, const std::vector<int> &contract_seq)
{
    if (!contract_seq.empty())
        contract_seq_=contract_seq;
    cutting_site_tensors_=cutting_site_tensors;
    env_tensors_=env_tensors;
    obtain_RDM();
}
template
void PEPSt_Subtensor_RDM<ITensor>::update_RDM(const std::vector<ITensor> &cutting_tensor, const std::vector<ITensor> &cutting_bond_tensors, const std::vector<ITensor> &env_tensors, const std::vector<int> &contract_seq);
template
void PEPSt_Subtensor_RDM<IQTensor>::update_RDM(const std::vector<IQTensor> &cutting_tensor, const std::vector<IQTensor> &env_tensors, const std::vector<IQTensor> &cutting_bond_tensors, const std::vector<int> &contract_seq);


template <class TensorT>
void PEPSt_Subtensor_RDM<TensorT>::obtain_sub_env_tensor()
{
    obtain_sub_tensors();

    //obtain sub_env_tensor_ from contracting env_tensors_ with left_tensors_ 
    std::vector<bool> env_contracted(env_tensors_.size(),false);
    for (int sitei=0; sitei<left_tensors_.size(); sitei++)
    {
        TensorT double_layer_tensor=left_tensors_[sitei];
        //contract some env tensors with single layer tensors
        for (int envi=0; envi<env_tensors_.size(); envi++)
        {
            const auto &env_tensor=env_tensors_[envi];
            if (env_contracted[envi]) continue;
            if (contract_seq_[envi]==1) continue;
            if (commonIndex(env_tensor,left_tensors_[sitei])!=IndexT::Null())
            {
                double_layer_tensor*=env_tensor;
                env_contracted[envi]=true;
            }
        }

        double_layer_tensor*=prime(dag(left_tensors_[sitei]));
        //Print(double_layer_tensor.indices());

        //contract other env tensors with double layer tensors
        for (int envi=0; envi<env_tensors_.size(); envi++)
        {
            const auto &env_tensor=env_tensors_[envi];
            if (env_contracted[envi]) continue;
            if (contract_seq_[envi]==0) continue;
            if (commonIndex(env_tensor,left_tensors_[sitei])!=IndexT::Null())
            {
                double_layer_tensor*=env_tensor;
                env_contracted[envi]=true;
            }
        }
        //Print(double_layer_tensor.indices());

        if (sitei==0)
        {
            sub_env_tensor_=double_layer_tensor;
        }
        else
        {
            sub_env_tensor_*=double_layer_tensor;
        }
    }
    //PrintDat(sub_env_tensor_);

    //Make sub_env_tensor_ positive hermitian
    TensorT sub_env_tensor_herm=(sub_env_tensor_+conj(swapPrime(sub_env_tensor_,0,1)))/2.;
    TensorT U,D;
    diagHermitian(sub_env_tensor_herm,U,D);
    D.mapElems([](double x){return (std::abs(x));});
    TensorT sub_env_tensor_pos_herm=dag(U)*D*prime(U);

    Print((sub_env_tensor_herm-sub_env_tensor_pos_herm).norm()/sub_env_tensor_herm.norm());
    //Print((sub_env_tensor_herm+sub_env_tensor_pos_herm).norm()/sub_env_tensor_herm.norm());

    sub_env_tensor_=sub_env_tensor_pos_herm;
}
template
void PEPSt_Subtensor_RDM<ITensor>::obtain_sub_env_tensor();
template
void PEPSt_Subtensor_RDM<IQTensor>::obtain_sub_env_tensor();

template <class TensorT>
void PEPSt_Subtensor_RDM<TensorT>::obtain_sub_tensors()
{
    left_tensors_.clear();
    sub_tensors_.clear();

    for (int sitei=0; sitei<cutting_site_tensors_.size(); sitei++)
    {
        //pick out the left indices
        std::vector<IndexT> left_indices;
        for (const auto &ind : cutting_site_tensors_[sitei].indices())
        {
            if (ind.type()==Site) continue;
            for (const auto &env_tensor : env_tensors_)
            {
                if (hasindex(env_tensor,dag(ind)))
                {
                    left_indices.push_back(ind);
                    break;
                }
            }
        }

        left_tensors_.push_back(TensorT(left_indices));
        sub_tensors_.push_back(TensorT());
        //decompose cutting_site_tensors_
        denmatDecomp(cutting_site_tensors_[sitei],left_tensors_[sitei],sub_tensors_[sitei],Fromleft);

        //PrintDat(left_tensors_[sitei]);
        //PrintDat(sub_tensors_[sitei]);
        //Print((cutting_site_tensors_[sitei]-left_tensors_[sitei]*sub_tensors_[sitei]).norm());
    }
}
template
void PEPSt_Subtensor_RDM<ITensor>::obtain_sub_tensors();
template
void PEPSt_Subtensor_RDM<IQTensor>::obtain_sub_tensors();


double heisenberg_energy_from_RDM(const PEPSt_Subtensor_RDM<IQTensor> &peps_rdm)
{
    std::vector<IQIndex> phys_legs;
    for (int sitei=0; sitei<peps_rdm.cutting_sites_no(); sitei++)
    {
        phys_legs.push_back(peps_rdm.cutting_phys_legs(sitei));
    }

    Complex energy=0;

    NN_Heisenberg_Hamiltonian heisenberg_gate({phys_legs[0],phys_legs[1]});
    energy=(peps_rdm.RDM()*(heisenberg_gate.site_tensors(0)*heisenberg_gate.bond_tensor()*heisenberg_gate.site_tensors(1))).toComplex();

    return (energy/peps_rdm.wf_norm_sq()).real();
}

