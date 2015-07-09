
#ifndef _PEPS_LATTICE_H_
#define _PEPS_LATTICE_H_

#include "predef.h"
#include "boost/multi_array.hpp"


//
//Lattice_Base
//
//Represent the structure of a 2D lattice
//

class Lattice_Base
{
    public:
        //
        //Constructors
        //

        Lattice_Base();

        //
        //Accessor Methods
        //
        inline int n_site_to_links() const
        {
            return n_site_to_links_;
        }

        inline int n_sites_() const
        {
            return n_sites_;
        }

        inline int n_bulk_links() const
        {
            return n_bulk_links_;
        }

        inline int n_boundary_links() const
        {
            return n_boundary_links_;
        }

    protected:
        //
        //Data Member
        //
        std::array<Boundary,2> boundary_condition_;

        //n_site_to_links stores number of links connect to a single site 
        int n_site_to_links_, n_sites_, n_bulk_links_, n_boundary_links_;
        
        //neighbour_sites_[i] stores sites connected to site i
        std::vector<std::vector<int>> neighbour_sites_, site_to_bulk_links_, site_to_boundary_links_;
        std::vector<std::array<int,2>> bulk_link_to_sites_;
        std::vector<int> boundary_link_to_sites_;
};


//
//Square_Lattice 
//
class Square_Lattice : public Lattice_Base
{
    public:
        //
        //Constructor
        //
        Square_Lattice(int& lx, int& ly, std::array<Boundary,2>& boundary_condition={Open,Open});

        //
        //Access Method
        //

    private:
        //
        //Data Member
        //
        int lx_, ly_;
};



#endif
