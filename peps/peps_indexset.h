
#ifndef _PEPS_INDEXSET_H_
#define _PEPS_INDEXSET_H_

#include "lattice.h"

//
//PEPSt_IndexSet_Base class
//the lower t is for "template"
//
//restore the indices of physical ones and virtual ones for certain lattice
//
//
template <class IndexT>
class PEPSt_IndexSet_Base
{
    public:
        PEPSt_IndexSet_Base() {}
        PEPSt_IndexSet_Base(const int &d, const int &D, const int &n_sites_total, const int &n_bonds_to_one_site);

        //
        //Access methods
        //
        inline const int &d() const
        {
            return d_;
        }

        inline const int &D() const
        {
            return D_;
        }
        
        inline const IndexT &phys_legs(const int &site_i) const
        {
            return phys_legs_[site_i];
        }

        inline const IndexT &virt_legs(const int &leg_i) const
        {
            return virt_legs_[leg_i];
        }

        const std::string &name() const { return name_; }

    protected:
        //dimension of single physical/virtual leg.
        int d_, D_;
        std::vector<IndexT> phys_legs_, virt_legs_;

        std::string name_;
};

//
//Indices for general PEPS
//using Index
//
class PEPS_IndexSet : public PEPSt_IndexSet_Base<Index>
{
    public:
        //
        //Constructor
        //
        PEPS_IndexSet() {}
        PEPS_IndexSet(const int &d, const int &D, const Lattice_Base &lattice);

        //
        //Constructor Helper
        //
        void init_phys_legs();
        void init_virt_legs();
};

//
//Indices for spin one half PEPS
//using IQIndex
//
class IQPEPS_IndexSet_SpinHalf : public PEPSt_IndexSet_Base<IQIndex>
{
    public:
        //
        //Constructor
        //
        IQPEPS_IndexSet_SpinHalf() {}
        //D=sum_{i=0}^{n}{2*i/2+1}. Namely, a virtual leg is formed by
        //0\circplus 1/2\circplus 1 ...
        IQPEPS_IndexSet_SpinHalf(const int &D, const Lattice_Base &lattice);
        //specify the quantum number in virt_leg_spin: 
        //virt_leg_spin[i] = # of spin i/2
        //should be translated to Sz quantum number
        IQPEPS_IndexSet_SpinHalf(const int &D, const std::vector<int> &virt_leg_spin, const Lattice_Base &lattice);

        //
        //Constructor Helpers
        //Initialize physical/virtual legs
        //
        void init_phys_legs();
        void init_virt_legs(const int &spin_dim, const std::vector<int> &virt_indqn_deg);

    private:

};

#endif
