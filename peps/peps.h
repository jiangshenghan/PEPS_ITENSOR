
#ifndef _PEPS_H_
#define _PEPS_H_

#include "lattice.h"

//
//peps class
//
//
class PEPS
{
    public:
        //
        //PEPS Constructors
        //
        PEPS(int d, int D, const Lattice_Base& lattice);

        //
        //PEPS Destructors
        //
        ~PEPS();

        //
        //Access Methods
        //

    private:
        //
        //Data Member
        //
        int d_, D_;
        Lattice_Base& lattice_;
        vector<Index> site_legs_, bulk_link_legs_, boundary_link_legs_;
        vector< IndexSet<Index> > site_indexsets_;
        vector<ITensor> site_tensors_;
};


#endif
