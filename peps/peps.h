
#ifndef _PEPS_H_
#define _PEPS_H_

#include "peps_indexset.h"

//
//class PEPSt_Torus
//the lower t is for template
//
//using PEPS_Torus for ITensor while IQPEPS_Torus for IQTensor
//
template <class TensorT>
class PEPSt_Torus
{
    public:

        //
        //type alias
        //
        using IndexT=typename TensorT::IndexT;

        //
        //Constructors
        //
        //Assign site tensors with random values
        PEPSt_Torus(const Lattice_Torus_Base &lattice, const PEPSt_IndexSet_Base<IndexT> &index_set);

        //Initialize site tensors by tensor data in one unit cell.
        //Thus, the PEPS is translationally invariant.
        PEPSt_Torus(const Lattice_Torus_Base &lattice, const PEPSt_IndexSet_Base<IndexT> &index_set, std::vector<TensorT> &site_tensors_uc);

        //
        //Access Methods
        //
        inline const int &d() const
        {
            return d_;
        }

        inline int &D() const
        {
            return D_;
        }

        inline const int &n_sites_uc() const
        {
            return lattice_.n_sites_uc();
        }

        inline const int &n_bonds_uc() const
        {
            return lattice_.n_bonds_uc();
        }

        inline const int &n_bonds_to_one_site() const
        {
            return lattice_.n_bonds_to_one_site();
        }

        inline const int &n_sites_total() const
        {
            return lattice_.n_sites_total();
        }

        inline const int &n_bonds_total() const
        {
            return lattice_.n_bonds_total();
        }

        inline const Lattice_Torus_Base &lattice() const
        {
            return lattice_;
        }

        inline TensorT &site_tensors(int i) const
        {
            return site_tensors_[i];
        }

        inline std::vector<TensorT> &site_tensors()
        {
            return site_tensors_;
        }

        //
        //Methods
        //


        //
        //Constructor Helpers
        //

        //Assign value of tensor TB to tensor TA without changing the index.
        //Hilbert space of TB should be morphism to that of TA
        //For IQTensor, arrows should also match
        void tensor_assignment(TensorT &TA, TensorT &TB);
        //According to the library, ITensor is constructed by IndexSet<Index>(indexset), while IQTensor is constructed by indexset directly
        void construct_tensor(TensorT &tensor, std::vector<IndexT> &indexset);
        //construct new tensors
        void new_site_tensors();
        //using random initial site tensors
        //for IQTensor we need to make sure at least one block is created
        void random_site_tensors();


    private:
        //
        //Data Member
        //
        //d_ is physical index dimension while D_ is virtual leg dim.
        int d_, D_;

        const Lattice_Torus_Base &lattice_;

        const PEPSt_IndexSet_Base<IndexT> &index_set_;

        std::vector<TensorT> site_tensors_; 

        //divergence of each site tensors
        //used for IQPEPS
        //default setting to 0
        std::vector<QN> tensors_div_;
};
using PEPS_Torus=PEPSt_Torus<ITensor>;
using IQPEPS_Torus=PEPSt_Torus<IQTensor>;

#endif
