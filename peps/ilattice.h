
#ifndef _ILATTICE_H_
#define _ILATTICE_H_

//
//iLattice_Base
//
//Represent the structure of a 2D lattice with infinite size
//

class iLattice_Base
{
    public:
        //
        //Constructors
        //
        iLattice_Base(int& n_sites_uc, int& n_links_uc);

        //
        //pure virtual destructor to make an abstract class
        //
        virtual ~iLattice_Base()=0;

        //
        //Accessor Methods
        //
        inline int n_sites_uc() const
        {
            return n_sites_uc_;
        };

        inline int n_links_uc() const
        {
            return n_links_uc_;
        };

        inline int n_site_to_links() const
        {
            return n_site_to_links_;
        };


    protected:
        //
        //Data member
        //
        //n_sites_uc and n_links_uc store number of sites and bond in one unit cell
        //n_site_to_links refer to number of links connect to one site
        int n_sites_uc_, n_links_uc_, n_site_to_links;
        //neighbour_sites_[i] stores neighbour sites of site (0,0,i), where i=0,...,n_sites_uc_-1
        //site_to_links[i] stores links connected to site (0,0,i), where i=0,...,n_sites_uc_-1
        vector< vector<Coordinate> > neighbour_sites_, site_to_links_;
        //link_to_sites_[i] stores the two sites connected to link (0,0,i), where i=0,...,n_links_uc_=1
        vector< array<Coordinate,2> > link_to_sites_;

}

//
//infinite square lattice with two sites per unit cell as following
//
//      |    |
//    --b----a--
//      |    |
//     1|    |
//    0-a-2--b--
//     3|    |
//
//a,b and the four legs around a form a unit cell. The two vector is (1,-1) and (1,1).
//
class iSquare_Lattice : public iLattice_Base
{
    public:
        //
        //Constructor
        //
        iSquare_Lattice();

         Neighbour {Left, Up, Right, Down};

    private:
        //
        //Data Members
        //
};

#endif
