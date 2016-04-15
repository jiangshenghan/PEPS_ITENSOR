
#ifndef _TENSOR_PRODUCT_OPERATOR_H_
#define _TENSOR_PRODUCT_OPERATOR_H_

#include "singlet_tensor_basis.h"

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

        //return phys_legs with no prime and point Out if IQIndex
        IndexT phys_legs(int i)
        {
            IndexT leg=findtype(site_tensors_[i],Site);
            if (leg.primeLevel()!=0) leg.noprime();
            if (leg.dir()==In) leg.dag();
            return leg;
        }

    private:
        //site tensors denote tensors with phys_legs while bond tensors without
        std::vector<TensorT> site_tensors_, bond_tensors_;

};
using TPO=TPOt<ITensor>;
using IQTPO=TPOt<IQTensor>;

//
//various correlator for half spins
//
//SpinSpin correlator
IQTPO SpinSpin(); 
IQTPO SpinSpin(const std::array<IQIndex,2> &phys_legs);

//\epsilon_{ijk}S_iS_jS_k
//IQTPO scalarSpinChirality();

//\epsilon_{ijk}S_jS_k
//This operator breaks Sz rotation symmetry, so strictly speaking, we should use ITensor instead of IQTensor
std::array<IQTPO,3> vectorSpinChirality();
std::array<IQTPO,3> vectorSpinChirality(const std::array<IQIndex,2> &phys_legs);

//Heisenberg Hamilltonian for kagome cirac lattice, which involves three sites
IQTPO SpinSpin_kagome_cirac(const std::vector<IQIndex> &phys_legs);
//NN_Heisenberg trotter gate for kagome cirac lattice
IQTPO trotter_gate_kagome_cirac(const std::vector<IQIndex> &phys_legs, double t=1.);
//NN Heisenberg trotter gate for two site
IQTPO trotter_gate_NN_Heisenberg(const std::vector<IQIndex> &phys_legs, double t=1.);


template <class TensorT>
std::ostream &operator<<(std::ostream &s, const TPOt<TensorT> &T)
{
    s << "n_sites=" << T.n_sites() << ", n_bonds=" << T.n_bonds() << endl;
    s << "site tensors:" << endl;
    for (int i=0; i<T.n_sites(); i++) s << T.site_tensors(i) << endl;
    s << "bond tensors:" << endl;
    for (int i=0; i<T.n_bonds(); i++) s << T.bond_tensors(i) << endl;
    return s;
}


#endif
