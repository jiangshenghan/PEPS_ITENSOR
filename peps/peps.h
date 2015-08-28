
#ifndef _IPEPS_H_
#define _IPEPS_H_

#include "ilattice.h"

//
//class PEPSt_Torus
//(the lower "t" stands for "template")
//
//Use PEPS_Torus for ITensors,
//while IQPEPS_Torus for IQTensors
//
class iPEPS
{
    public:
        //
        //Constructors
        //
        iPEPS(int &d, int &D, const Lattice_uc_Base &lattice_uc);

        //
        //Destructors
        //
        ~iPEPS();

        //
        //Access Methods
        //
        inline int &d() const
        {
            return d_;
        }

        inline int &D() const
        {
            return D_;
        }

        inline const Lattice_uc_Base &lattice_uc() const
        {
            return lattice_uc_;
        }

        inline ITensor &site_tensors(int i) const
        {
            return site_tensors_[i];
        }

        inline ITensor &link_tensors(int i) const
        {
            return link_tensors_[i];
        }

    private:
        //
        //Data Member
        //
        //d_ is physical index dimension while D_ is virtual leg dim.
        int d_, D_;
        const Lattice_uc_Base& lattice_uc_;
        //site_legs_[i] stores indice for physical leg i, where i labels the sublattice. One should be able to generate legs (x,y,i) by site_legs_[i].
        vector<Index> site_legs_;
        //link_legs_[i] stores indice for two indices associated with link i, where i is the sublattice.
        vector< array<Index,2> > link_legs_;
        //site_tensors_[i] is an abstract tensor for site with sublattice i. Similar for link_tensors[i]
        //we set link tensor as diagonal matrix
        vector<ITensor> site_tensors_, link_tensors_;
};

#endif
