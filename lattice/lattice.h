
#ifndef _LATTICE_H_
#define _LATTICE_H_

#include "utilities.h"

//
//Lattice_Torus_Base
//

class Lattice_Torus_Base
{
    public:
        //
        //Constructors
        //
        Lattice_Torus_Base() {}

        Lattice_Torus_Base(const int &n_sites_uc, const int &n_bonds_uc, const std::array<int,2> &n_uc);

        //
        //Accessor Methods
        //
        inline const int &n_sites_uc() const
        {
            return n_sites_uc_;
        };

        inline const int &n_bonds_uc() const
        {
            return n_bonds_uc_;
        };

        inline const int &n_bonds_to_one_site() const
        {
            return n_bonds_to_one_site_;
        };

        inline const int &n_sites_total() const
        {
            return n_sites_total_;
        }

        inline const int &n_bonds_total() const
        {
            return n_bonds_total_;
        }

        inline const std::array<int,2> &n_uc() const
        {
            return n_uc_;
        };

        //get jth neighbour site of site_i
        inline const int &site_neighbour_sites(const int &site_i, const int &j) const
        {
            return site_neighbour_sites_[site_i][j];
        }

        //get jth bond connected to site_i
        inline const int &site_neighbour_bonds(const int &site_i, const int &j) const
        {
            return site_neighbour_bonds_[site_i][j];
        }

        inline const std::vector<int> &site_neighbour_bonds(const int &site_i) const
        {
            return site_neighbour_bonds_[site_i];
        }

        //get jth end site of bond_i
        inline const int &bond_end_sites(const int &bond_i, const int &j) const
        {
            return bond_end_sites_[bond_i][j];
        }

        inline const std::array<int,2> &bond_end_sites(const int &bond_i) const
        {
            return bond_end_sites_[bond_i];
        }

        inline const Coordinate &site_list_to_coord(int site_i) const
        {
            //ensure 0<=site_i<n_sites_total_
            site_i=(site_i%n_sites_total_+n_sites_total_)%n_sites_total_;
            return site_list_to_coord_[site_i];
        }

        inline const Coordinate &bond_list_to_coord(int bond_i) const
        {
            bond_i=(bond_i%n_bonds_total_+n_bonds_total_)%n_bonds_total_;
            return bond_list_to_coord_[bond_i];
        }

        inline const int &site_coord_to_list(Coordinate site_coord) const
        {
            site_coord[0]=(site_coord[0]%n_uc_[0]+n_uc_[0])%n_uc_[0];
            site_coord[1]=(site_coord[1]%n_uc_[1]+n_uc_[1])%n_uc_[1];
            return site_coord_to_list_[site_coord[0]][site_coord[1]][site_coord[2]];
        }

        inline const int &bond_coord_to_list(Coordinate bond_coord) const
        {
            bond_coord[0]=(bond_coord[0]%n_uc_[0]+n_uc_[0])%n_uc_[0];
            bond_coord[1]=(bond_coord[1]%n_uc_[1]+n_uc_[1])%n_uc_[1];
            return bond_coord_to_list_[bond_coord[0]][bond_coord[1]][bond_coord[2]];
        }


        virtual void print_lattice_inf()=0;

    protected:
        //
        //Data member
        //
        int n_sites_uc_, n_bonds_uc_, n_bonds_to_one_site_, n_sites_total_, n_bonds_total_;
        std::array<int,2> n_uc_;

        //site_neighbour_sites_[i] stores the list of neighbouring sites to site i;
        //site_neighbour_bonds_[i] stores the list of neibouring bonds to site i;
        std::vector<std::vector<int>> site_neighbour_sites_, site_neighbour_bonds_;
        //bond_end_sites_[i] stores the coordinate of two end sites to bond i 
        //this also defines the flow of quantum number: from bond_end_sites_[0] to bond_end_sites_[1]
        std::vector<std::array<int,2>> bond_end_sites_;

        //translate between coordinates and list(number)
        //site/bond_list_to_coord_[i]=(x,y,i')
        std::vector<Coordinate> site_list_to_coord_, bond_list_to_coord_;
        //site/bond_coord_to_list_[x][y][i']=i
        std::vector<std::vector<std::vector<int>>> site_coord_to_list_, bond_coord_to_list_;

};

//
//square lattice on periodic BC
//
//     0|   
//     -a-1-
//      | 
//
class Square_Lattice_Torus : public Lattice_Torus_Base
{
    public:
        //
        //Constructor
        //
        Square_Lattice_Torus(const std::array<int,2> &n_uc);

        enum Neighbour {Left=0, Up=1, Right=2, Down=3};

        virtual void print_lattice_inf();

    private:
        //
        //Data Members
        //
};



#endif
