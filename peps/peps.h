
#ifndef _PEPS_H_
#define _PEPS_H_

#include "singlet_tensor_basis.h"
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
        PEPSt_Torus() {}
        //Assign site tensors with random values
        PEPSt_Torus(const Lattice_Base &lattice, const PEPSt_IndexSet_Base<IndexT> &index_set);

        //Initialize site tensors by tensor data in one unit cell.
        //Thus, the PEPS is translationally invariant.
        PEPSt_Torus(const Lattice_Base &lattice, const PEPSt_IndexSet_Base<IndexT> &index_set, std::vector<TensorT> &site_tensors_uc, std::vector<TensorT> &bond_tensors_uc);

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

        inline const Lattice_Base &lattice() const
        {
            return lattice_;
        }

        inline const IndexT &phys_legs(const int &site_i) const
        {
            return index_set_.phys_legs(site_i);
        }

        inline const IndexT &virt_legs(const int &leg_i) const
        {
            return index_set_.virt_legs_(leg_i);
        }

        inline const TensorT &site_tensors(int i) const
        {
            return site_tensors_[i];
        }

        TensorT &site_tensors(int i)
        {
            return site_tensors_[i];
        }

        inline const std::vector<TensorT> &site_tensors() const
        {
            return site_tensors_;
        }

        inline const TensorT &bond_tensors(int i) const
        {
            return bond_tensors_[i];
        }

        TensorT &bond_tensors(int i)
        {
            return bond_tensors_[i];
        }

        inline const std::vector<TensorT> &bond_tensors() const
        {
            return bond_tensors_;
        }

        //
        //Methods
        //
        //given tensors in one uc, generate all translation related tensors
        void generate_site_tensors(std::vector<TensorT> site_tensors_uc);
        void generate_bond_tensors(std::vector<TensorT> bond_tensors_uc);


        //
        //Constructor Helpers
        //

        //Assign value of tensor TB to tensor TA without changing the index.
        //Hilbert space of TB should be morphism to that of TA
        //For IQTensor, arrows should also match
        void tensor_assignment(TensorT &TA, const TensorT &TB);
        //According to the itensor library, ITensor is constructed by IndexSet<Index>(indexset), while IQTensor is constructed by indexset directly
        void construct_tensor(TensorT &tensor, std::vector<IndexT> &indexset);
        //construct new tensors
        void new_site_tensors();
        void new_bond_tensors();
        //using random initial site tensors
        //for IQTensor we need to make sure at least one block is created
        //void random_site_tensors();


    private:
        //
        //Data Member
        //
        //d_ is physical index dimension while D_ is virtual leg dim.
        int d_, D_;

        const Lattice_Base &lattice_;

        const PEPSt_IndexSet_Base<IndexT> &index_set_;

        std::vector<TensorT> site_tensors_, bond_tensors_; 

        //name stores information about lattice, spin symmetry, extra degenerate
        std::string name_;

        //divergence of each site tensors
        //used for IQPEPS
        //default setting to 0
        //std::vector<QN> tensors_div_;
};
using PEPS_Torus=PEPSt_Torus<ITensor>;
using IQPEPS_Torus=PEPSt_Torus<IQTensor>;


void randomize_spin_sym_square_peps(IQPEPS_Torus &square_peps);


#endif
