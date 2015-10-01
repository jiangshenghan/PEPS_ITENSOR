
#include "lattice.h"

Lattice_Base::Lattice_Base(const int &n_sites_uc, const int &n_bonds_uc, const std::array<int,2> &n_uc):
    n_sites_uc_(n_sites_uc),
    n_bonds_uc_(n_bonds_uc),
    n_bonds_to_one_site_(2*n_bonds_uc/n_sites_uc),
    n_sites_total_(n_sites_uc*n_uc[0]*n_uc[1]),
    n_bonds_total_(n_bonds_uc*n_uc[0]*n_uc[1]),
    n_uc_(n_uc),
    site_neighbour_sites_(n_sites_total_,std::vector<int>(n_bonds_to_one_site_)),
    site_neighbour_bonds_(n_sites_total_,std::vector<int>(n_bonds_to_one_site_)),
    bond_end_sites_(n_bonds_total_),
    site_list_to_coord_(n_sites_total_),
    bond_list_to_coord_(n_bonds_total_),
    site_coord_to_list_(n_uc[0],std::vector<std::vector<int>>(n_uc[1],std::vector<int>(n_sites_uc))),
    bond_coord_to_list_(n_uc[0],std::vector<std::vector<int>>(n_uc[1],std::vector<int>(n_bonds_uc)))
{
    //
    //Initialize list and coordinate translation
    //
    for (int site_i=0; site_i<n_sites_total_; site_i++)
    {
        int lx_i=site_i/n_sites_uc_%n_uc[0],
            ly_i=site_i/n_sites_uc_/n_uc[0],
            subsite_i=site_i%n_sites_uc_;

        site_list_to_coord_[site_i]=Coordinate{lx_i,ly_i,subsite_i};
        site_coord_to_list_[lx_i][ly_i][subsite_i]=site_i;
    }

    for (int bond_i=0; bond_i<n_bonds_total_; bond_i++)
    {
        int lx_i=bond_i/n_bonds_uc_%n_uc[0],
            ly_i=bond_i/n_bonds_uc_/n_uc[0],
            subbond_i=bond_i%n_bonds_uc_;

        bond_list_to_coord_[bond_i]=Coordinate{lx_i,ly_i,subbond_i};
        bond_coord_to_list_[lx_i][ly_i][subbond_i]=bond_i;
    }
}


Square_Lattice_Torus::Square_Lattice_Torus(const std::array<int,2> &n_uc): Lattice_Base(1,2,n_uc)
{

    name_="square lattice on torus";

    for (int site_i=0; site_i<n_sites_total_; site_i++)
    {
        Coordinate site_i_coord=site_list_to_coord(site_i);
        std::vector<Coordinate> neigh_sites_coord(n_bonds_to_one_site_,site_i_coord);
        std::vector<Coordinate> neigh_bonds_coord(n_bonds_to_one_site_);

        neigh_sites_coord[Left][0]=site_i_coord[0]-1;
        neigh_sites_coord[Up][1]=site_i_coord[1]+1;
        neigh_sites_coord[Right][0]=site_i_coord[0]+1;
        neigh_sites_coord[Down][1]=site_i_coord[1]-1;

        neigh_bonds_coord[Left]=Coordinate{site_i_coord[0]-1,site_i_coord[1],1};
        neigh_bonds_coord[Up]=Coordinate{site_i_coord[0],site_i_coord[1],0};
        neigh_bonds_coord[Right]=Coordinate{site_i_coord[0],site_i_coord[1],1};
        neigh_bonds_coord[Down]=Coordinate{site_i_coord[0],site_i_coord[1]-1,0};

        for (int j=0; j<n_bonds_to_one_site_; j++)
        {
            site_neighbour_sites_[site_i][j]=site_coord_to_list(neigh_sites_coord[j]);
            site_neighbour_bonds_[site_i][j]=bond_coord_to_list(neigh_bonds_coord[j]);
        }
    }

    for (int bond_i=0; bond_i<n_bonds_total_; bond_i++)
    {
        Coordinate bond_i_coord=bond_list_to_coord(bond_i);
        std::array<Coordinate,2> end_sites_coord;

        if (bond_i_coord[2]==0) //vertical bonds, direction from down to up
        {
            end_sites_coord[0]=Coordinate{bond_i_coord[0],bond_i_coord[1],0};
            end_sites_coord[1]=Coordinate{bond_i_coord[0],bond_i_coord[1]+1,0};
        }
        else //horizontal bonds, direction from left to right
        {
            end_sites_coord[0]=Coordinate{bond_i_coord[0],bond_i_coord[1],0};
            end_sites_coord[1]=Coordinate{bond_i_coord[0]+1,bond_i_coord[1],0};
        }

        bond_end_sites_[bond_i][0]=site_coord_to_list(end_sites_coord[0]);
        bond_end_sites_[bond_i][1]=site_coord_to_list(end_sites_coord[1]);
    }
}

void Square_Lattice_Torus::print_lattice_inf()
{
    cout << name_ << endl << "lattice size " << n_uc_[0] << "x" << n_uc_[1] << " and " << n_sites_total_ << " sites, " << n_bonds_total_ << " bonds." << endl;

    cout << "Check neighbouring sites: " << endl << endl;
    for (int site_i=0; site_i<n_sites_total_; site_i++)
    {
        cout << "Neighbour sites of site " << site_list_to_coord(site_i) << " are:" << endl;
        for (auto j : site_neighbour_sites_[site_i])
        {
            cout << site_list_to_coord(j) << endl;
        }
        cout << endl;
    }
    cout << endl;

    cout << "Check neighbouring bonds: " << endl << endl;
    for (int site_i=0; site_i<n_sites_total_; site_i++)
    {
        //cout << "Bonds connected to site " << site_list_to_coord(site_i) << " are:" << endl;
        cout << "Bonds connected to site " << site_i << " are:" << endl;
        for (auto j : site_neighbour_bonds_[site_i])
        {
            cout << j << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "Check end sites of bonds: " << endl << endl;
    for (int bond_i=0; bond_i<n_bonds_total_; bond_i++)
    {
        //cout << "The two end sites of bond " << bond_list_to_coord(bond_i) << " are site " << site_list_to_coord(bond_end_sites(bond_i,0)) << " and site " << site_list_to_coord(bond_end_sites(bond_i,1)) <<endl;
        cout << "The two end sites of bond " << bond_i << " are site " << bond_end_sites(bond_i,0) << " and site " << bond_end_sites(bond_i,1) <<endl;
    }
    cout << endl;
}



Square_Lattice_Cylinder::Square_Lattice_Cylinder(const std::array<int,2> &n_uc): Lattice_Base(1,2,n_uc)
{
    name_="square lattice on cylinder";

    for (int site_i=0; site_i<n_sites_total_; site_i++)
    {
        Coordinate site_i_coord=site_list_to_coord(site_i);
        std::vector<Coordinate> neigh_sites_coord(n_bonds_to_one_site_,site_i_coord);
        std::vector<Coordinate> neigh_bonds_coord(n_bonds_to_one_site_);

        neigh_sites_coord[Left][0]=site_i_coord[0]-1;
        neigh_sites_coord[Up][1]=site_i_coord[1]+1;
        neigh_sites_coord[Right][0]=site_i_coord[0]+1;
        neigh_sites_coord[Down][1]=site_i_coord[1]-1;

        neigh_bonds_coord[Left]=Coordinate{site_i_coord[0]-1,site_i_coord[1],1};
        neigh_bonds_coord[Up]=Coordinate{site_i_coord[0],site_i_coord[1],0};
        neigh_bonds_coord[Right]=Coordinate{site_i_coord[0],site_i_coord[1],1};
        neigh_bonds_coord[Down]=Coordinate{site_i_coord[0],site_i_coord[1]-1,0};

        for (int j=0; j<n_bonds_to_one_site_; j++)
        {
            //site_i is at leftmost, then, we set the left neigh site and bond to be -1
            if (neigh_sites_coord[j][0]<0) 
            {
                site_neighbour_sites_[site_i][j]=-1;
                site_neighbour_bonds_[site_i][j]=-1;
                continue;
            }

            //site_i is at rightmost, then, we set the right neigh site to be -2, but right neigh bond still exist
            if (neigh_sites_coord[j][0]>=n_uc_[0])
            {
                site_neighbour_sites_[site_i][j]=-2;
                site_neighbour_bonds_[site_i][j]=bond_coord_to_list(neigh_bonds_coord[j]);
                continue;
            }

            site_neighbour_sites_[site_i][j]=site_coord_to_list(neigh_sites_coord[j]);
            site_neighbour_bonds_[site_i][j]=bond_coord_to_list(neigh_bonds_coord[j]);
        }
        //cout << "site " << site_i << endl << "neigh: " << site_neighbour_sites_[site_i] << endl;
    }

    for (int bond_i=0; bond_i<n_bonds_total_; bond_i++)
    {
        Coordinate bond_i_coord=bond_list_to_coord(bond_i);
        std::array<Coordinate,2> end_sites_coord;

        if (bond_i_coord[2]==0) //vertical bonds, direction from down to up
        {
            end_sites_coord[0]=Coordinate{bond_i_coord[0],bond_i_coord[1],0};
            end_sites_coord[1]=Coordinate{bond_i_coord[0],bond_i_coord[1]+1,0};
        }
        else //horizontal bonds, direction from left to right
        {
            end_sites_coord[0]=Coordinate{bond_i_coord[0],bond_i_coord[1],0};
            end_sites_coord[1]=Coordinate{bond_i_coord[0]+1,bond_i_coord[1],0};
        }

        bond_end_sites_[bond_i][0]=site_coord_to_list(end_sites_coord[0]);
        bond_end_sites_[bond_i][1]=site_coord_to_list(end_sites_coord[1]);
        // if the bond is at rightmost, set the end_site_coord[1] to be -1
        if (end_sites_coord[1][0]>=n_uc_[0])
            bond_end_sites_[bond_i][1]=-2;
    }
}


void Square_Lattice_Cylinder::print_lattice_inf()
{
    cout << name_ << endl << "lattice size " << n_uc_[0] << "x" << n_uc_[1] << " and " << n_sites_total_ << " sites, " << n_bonds_total_ << " bonds, " << n_boundary_legs() << " boudnary legs." << endl << endl;

    cout << "Check neighbouring sites: " << endl << endl;
    for (int site_i=0; site_i<n_sites_total_; site_i++)
    {
        cout << "Neighbour sites of site " << site_list_to_coord(site_i) << " are:" << endl;
        cout << site_neighbour_sites_[site_i] << endl;
        //for (auto j : site_neighbour_sites_[site_i])
        //{
        //    cout << site_list_to_coord(j) << endl;
        //}
        cout << endl;
    }
    cout << endl;

    cout << "Check neighbouring bonds: " << endl << endl;
    for (int site_i=0; site_i<n_sites_total_; site_i++)
    {
        //cout << "Bonds connected to site " << site_list_to_coord(site_i) << " are:" << endl;
        cout << "Bonds connected to site " << site_i << " are:" << endl;
        for (auto j : site_neighbour_bonds_[site_i])
        {
            cout << j << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "Check end sites of bonds: " << endl << endl;
    for (int bond_i=0; bond_i<n_bonds_total_; bond_i++)
    {
        //cout << "The two end sites of bond " << bond_list_to_coord(bond_i) << " are site " << site_list_to_coord(bond_end_sites(bond_i,0)) << " and site " << site_list_to_coord(bond_end_sites(bond_i,1)) <<endl;
        cout << "The two end sites of bond " << bond_i << " are site " << bond_end_sites(bond_i,0) << " and site " << bond_end_sites(bond_i,1) <<endl;
    }
    cout << endl;
}
