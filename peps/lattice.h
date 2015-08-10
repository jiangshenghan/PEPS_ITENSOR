
#ifndef _LATTICE_H_
#define _LATTICE_H_

#include "predef.h"


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
        Lattice_Base(int& n_site_to_links, int &n_sites, int& n_bulk_links, int& n_boundary_links);

        //
        //pure virtual destructor to make a abstract class
        //
        virtual ~Lattice_Base()=0;

        //
        //Accessor Methods
        //
        inline int n_site_to_links() const
        {
            return n_site_to_links_;
        };

        inline int n_sites_() const
        {
            return n_sites_;
        };

        inline int n_bulk_links() const
        {
            return n_bulk_links_;
        };

        inline int n_boundary_links() const
        {
            return n_boundary_links_;
        };

        inline int neighbour_sites(int& i, int &j) const
        {
            return neighbour_sites_[i][j];
        };

        inline int site_to_bulk_links(int& i, int& j) const
        {
            return site_to_bulk_links_[i][j];
        };

        inline int site_to_boundary_links(int& i, int& j) const
        {
            return site_to_boundary_links_[i][j];
        };

        inline int bulk_link_to_sites(int& i, int& j) const
        {
            return bulk_link_to_sites_[i][j];
        };

        inline int boundary_link_to_sites(int& i) const
        {
            return boundary_link_to_sites_[i];
        };

    protected:
        //
        //Data Member
        //
        //n_site_to_links stores number of links connect to a single site 
        int n_site_to_links_, n_sites_, n_bulk_links_, n_boundary_links_;
        //neighbour_sites_[i] stores sites connected to site i, i=1,...,n_sites_
        vector< vector<int> > neighbour_sites_, site_to_bulk_links_, site_to_boundary_links_;
        //bulk_link_to_sites implicitly define direction for links
        //useful for quantum number implementation
        vector< array<int,2> > bulk_link_to_sites_;
        vector<int> boundary_link_to_sites_;
        Boundary boundary_conditions_;
        
};


//
//Square_Lattice 
//
class Square_Lattice_Open : public Lattice_Base
{
    public:
        //
        //Constructor
        //
        Square_Lattice(int& lx, int& ly)

        //
        //Access Method
        //

    private:
        //
        //Data Member
        //
        int lx_, ly_;
        vector< array<int,2> > site_to_coord_list_;
        vector< vector<int> > coord_to_site_list_;
};



#endif
