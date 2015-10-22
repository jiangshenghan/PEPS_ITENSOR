
#ifndef _TENSOR_PRODUCT_OPERATOR_H_
#define _TENSOR_PRODUCT_OPERATOR_H_

#include "predef.h"

//class for tensor product operator
template <class TensorT>
class TPOt
{
    public:

        //
        //Type Alias
        //
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;

        //
        //Constructors
        //
        TPOt() {}
        TPOt(int n_sites, int n_bonds): site_tensors_(n_sites), bond_tensors_(n_bonds) {}
        

        //
        //Access Method (const and non-const)
        //
        int n_sites() const { return site_tensors_.size(); }
        int n_bonds() const { return bond_tensors_.size(); }
        const TensorT &site_tensors(int i) const { return site_tensors_[i]; }
        TensorT &site_tensors(int i) { return site_tensors_[i]; }
        const std::vector<TensorT> &site_tensors() const { return site_tensors_; }
        std::vector<TensorT> &site_tensors() { return site_tensors_; }

        const TensorT &bond_tensors(int i) const { return bond_tensors_[i]; }
        TensorT &bond_tensors(int i) { return bond_tensors_[i]; }
        const std::vector<TensorT> &bond_tensors() const { return bond_tensors_; }

        //return phys_legs with no prime
        IndexT phys_legs(int i)
        {
            IndexT leg=findtype(site_tensors_[i],Site);
            if (leg.primeLevel()!=0) leg.noprime().dag();
            return leg;
        }

    private:
        //site tensors denote tensors with phys_legs while bond tensors without
        std::vector<TensorT> site_tensors_, bond_tensors_;

};
using TPO=TPOt<ITensor>;
using IQTPO=TPOt<IQTensor>;

#endif
