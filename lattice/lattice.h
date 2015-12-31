
#ifndef _LATTICE_H_
#define _LATTICE_H_

#include "utilities.h"

//
//Lattice_Base
//generate various lattices, which can be defined on torus, cylinder and open boundary
//For torus, there is no boundary
//For cylinder, there are left boundary (labeled as -1) and right boundary (labelede as -2)
//Here "bonds" have more general meaning. For example, sometimes, we also call plaquette as bond

class Lattice_Base
{
    public:
        //
        //Constructors
        //
        //Lattice_Base() {n_sites_total_=0;}

        Lattice_Base(int n_sites_uc, int n_bonds_uc, const std::array<int,2> &n_uc, int n_boundary_legs=0, int n_sites_to_one_bond=2);

        //Constructor helpers
        void init_list_coord();

        //lattice is empty if use default constructor, thus n_sites_total==0
        //bool is_empty() const { return (n_sites_total_==0); }

        //
        //Accessor Methods
        //
        int n_sites_uc() const
        {
            return n_sites_uc_;
        }

        int n_bonds_uc() const
        {
            return n_bonds_uc_;
        }

        int n_bonds_to_one_site() const
        {
            return n_bonds_to_one_site_;
        }

        int n_sites_to_one_bond() const
        {
            return n_sites_to_one_bond_;
        }

        int n_sites_total() const
        {
            return n_sites_total_;
        }

        int n_bonds_total() const
        {
            return n_bonds_total_;
        }

        inline const std::array<int,2> &n_uc() const
        {
            return n_uc_;
        }

        int n_boundary_legs() const
        {
            return n_boundary_legs_;
        }

        //get jth bond connected to site_i
        int site_neighbour_bonds(const int &site_i, const int &j) const
        {
            return site_neighbour_bonds_[site_i][j];
        }
        const std::vector<int> &site_neighbour_bonds(const int &site_i) const
        {
            return site_neighbour_bonds_[site_i];
        }

        //get jth site connect bond_i
        int bond_neighbour_sites(const int &bond_i, const int &j) const
        {
            return bond_neighbour_sites_[bond_i][j];
        }
        inline const std::vector<int> &bond_neighbour_sites(const int &bond_i) const
        {
            return bond_neighbour_sites_[bond_i];
        }

        int boundary_neighbour_site(int boundary_i) const
        {
            return boundary_neighbour_site_[boundary_i];
        }
        int boundary_neighbour_bond(int boundary_i) const
        {
            return boundary_neighbour_bond_[boundary_i];
        }

        int site_neighbour_boundaries(int site_i, int j) const
        {
            return site_neighbour_boundaries_[site_i][j];
        }
        const std::vector<int> &site_neighbour_boundaries(int site_i) const
        {
            return site_neighbour_boundaries_[site_i];
        }

        int bond_neighbour_boundaries(int bond_i, int j) const
        {
            return bond_neighbour_boundaries_[bond_i][j];
        }
        const std::vector<int> &bond_neighbour_boundaries(int bond_i) const
        {
            return bond_neighbour_boundaries_[bond_i];
        }

        Coordinate site_list_to_coord(int site_i) const
        {
            if (name_.find("torus")==std::string::npos && (site_i<0 || site_i>n_sites_total_))
            {
                return Coordinate{site_i,site_i,site_i};
            }
            
            //ensure 0<=site_i<n_sites_total_ for torus
            site_i=(site_i%n_sites_total_+n_sites_total_)%n_sites_total_;
            return site_list_to_coord_[site_i];
        }

        Coordinate bond_list_to_coord(int bond_i) const
        {
            if (name_.find("torus")==std::string::npos && (bond_i<0 || bond_i>n_bonds_total_))
            {
                return Coordinate{bond_i,bond_i,bond_i};
            }

            bond_i=(bond_i%n_bonds_total_+n_bonds_total_)%n_bonds_total_;
            return bond_list_to_coord_[bond_i];
        }

        inline const int &site_coord_to_list(Coordinate site_coord) const
        {
            site_coord[0]=(site_coord[0]%n_uc_[0]+n_uc_[0])%n_uc_[0];
            site_coord[1]=(site_coord[1]%n_uc_[1]+n_uc_[1])%n_uc_[1];
            return site_coord_to_list_[site_coord[0]][site_coord[1]][site_coord[2]];
        }
        
        const int &site_coord_to_list(int x, int y, int subsite) const
        {
            return site_coord_to_list(Coordinate{x,y,subsite});
        }

        inline const int &bond_coord_to_list(Coordinate bond_coord) const
        {
            bond_coord[0]=(bond_coord[0]%n_uc_[0]+n_uc_[0])%n_uc_[0];
            bond_coord[1]=(bond_coord[1]%n_uc_[1]+n_uc_[1])%n_uc_[1];
            return bond_coord_to_list_[bond_coord[0]][bond_coord[1]][bond_coord[2]];
        }

        const int &bond_coord_to_list(int x, int y, int subbond) const
        {
            return bond_coord_to_list(Coordinate{x,y,subbond});
        }

        const std::string &name() const { return name_; }

        //get common bond of site A and site B
        //if no common bond, return -1
        int comm_bond(int site_A, int site_B) const
        {
            for (int i=0; i<n_bonds_to_one_site_; i++)
            {
                auto bond_A=site_neighbour_bonds_[site_A][i];

                if (bond_neighbour_sites_[bond_A][0]==site_A && bond_neighbour_sites_[bond_A][1]==site_B)
                    return bond_A;

                if (bond_neighbour_sites_[bond_A][1]==site_A && bond_neighbour_sites_[bond_A][0]==site_B)
                    return bond_A;
            }
            return -1;
        }

        void print_lattice_inf() const;

        //read/write from/to file
        //void read(std::istream &s);
        //void write(std::ostream &s) const;

    protected:
        //
        //Data member
        //
        int n_sites_uc_, n_bonds_uc_, n_bonds_to_one_site_, n_sites_to_one_bond_, n_sites_total_, n_bonds_total_, n_boundary_legs_;
        std::array<int,2> n_uc_;

        //site_neighbour_sites_[i] stores the list of neighbouring sites to site i;
        //std::vector<std::vector<int>> site_neighbour_sites_;

        //site_neighbour_bonds_[i] stores the list of neibouring bonds to site i;
        std::vector<std::vector<int>> site_neighbour_bonds_;
        //bond_neighbour_sites_[i] stores the coordinate of two end sites to bond i 
        //this also defines the flow of quantum number: from bond_neighbour_sites_[0] to bond_neighbour_sites_[1]
        std::vector<std::vector<int>> bond_neighbour_sites_;
        //set boundary_neighbour_site_[boundaryi] to minus number if boundaryi connects to bond rather than site.
        //Similar for boundary_neighbour_bond_
        std::vector<int> boundary_neighbour_site_, boundary_neighbour_bond_;
        std::vector<std::vector<int>> site_neighbour_boundaries_, bond_neighbour_boundaries_;

        //translate between coordinates and list(number)
        //site/bond_list_to_coord_[i]=(x,y,i')
        std::vector<Coordinate> site_list_to_coord_, bond_list_to_coord_;
        //site/bond_coord_to_list_[x][y][i']=i
        std::vector<std::vector<std::vector<int>>> site_coord_to_list_, bond_coord_to_list_;

        std::string name_;

};

//
//square lattice on periodic BC
//
//     0|   
//     -a-1-
//      | 
//
class Square_Lattice_Torus : public Lattice_Base
{
    public:
        //
        //Constructor
        //
        //Square_Lattice_Torus() {}
        Square_Lattice_Torus(const std::array<int,2> &n_uc);

        enum Neighbour {Left=0, Up=1, Right=2, Down=3};

        //virtual void print_lattice_inf() const;

    private:
        //
        //Data Members
        //
};

//
//square lattice on cylinder (open boundary in x direction and periodic boundary on y direction)
//
class Square_Lattice_Cylinder : public Lattice_Base
{
    public:
        //Square_Lattice_Cylinder() {}
        Square_Lattice_Cylinder(const std::array<int,2> &n_uc);

        enum Neighbour {Left=0, Up=1, Right=2, Down=3};
};


//
//Square lattice with OBC
//
class Square_Lattice_Open : public Lattice_Base
{
    public:
        //Square_Lattice_Open() {}
        Square_Lattice_Open(const std::array<int,2> &n_uc);

        enum Neighbour {Left=0, Up=1, Right=2, Down=3};
};


//
//Honeycomb lattice on torus
//The convention follows fig on honeycomb_peps_notes.pdf
//three bonds (x,y,i) are connected with (x,y,0)
//
class Honeycomb_Lattice_Torus : public Lattice_Base
{
    public:
        //Honeycomb_Lattice_Torus() {}
        Honeycomb_Lattice_Torus(const std::array<int,2> &n_uc);

        enum Neighbour {NeighA=0, NeighB=1, NeighC=2};
};

//
//Kagome lattice with cirac geometry on torus
//The convention follows fig on kagome_lattice_cirac_notes.pdf
//Here, "bonds" are actually plaquette center
//
class Kagome_Cirac_Lattice_Torus : public Lattice_Base
{
    public:
        //Kagome_Cirac_Lattice_Torus() {}
        Kagome_Cirac_Lattice_Torus(const std::array<int,2> &n_uc);
};


#endif
