
#include "lattice.h"

Lattice_Base::Lattice_Base(int n_sites_uc, int n_bonds_uc, const std::array<int,2> &n_uc, int n_boundary_legs, int n_sites_to_one_bond):
    n_sites_uc_(n_sites_uc),
    n_bonds_uc_(n_bonds_uc),
    n_bonds_to_one_site_(n_sites_to_one_bond*n_bonds_uc/n_sites_uc),
    n_sites_to_one_bond_(n_sites_to_one_bond),
    n_sites_total_(n_sites_uc*n_uc[0]*n_uc[1]),
    n_bonds_total_(n_bonds_uc*n_uc[0]*n_uc[1]),
    n_boundary_legs_(n_boundary_legs),
    n_uc_(n_uc),
    //site_neighbour_sites_(n_sites_total_,std::vector<int>(n_bonds_to_one_site_)),
    site_neighbour_bonds_(n_sites_total_,std::vector<int>(n_bonds_to_one_site_)),
    bond_neighbour_sites_(n_bonds_total_,std::vector<int>(n_sites_to_one_bond_)),
    boundary_neighbour_site_(n_boundary_legs_),
    boundary_neighbour_bond_(n_boundary_legs_),
    site_neighbour_boundaries_(n_sites_total_,std::vector<int>(n_bonds_to_one_site_,-1)),
    bond_neighbour_boundaries_(n_bonds_total_,std::vector<int>(n_sites_to_one_bond_,-1)),
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

void Lattice_Base::print_lattice_inf() const
{
    cout << "\n==============================================\n" << endl;
    cout << name_ << endl << "lattice size " << n_uc_[0] << "x" << n_uc_[1] << " and " << n_sites_total_ << " sites, " << n_bonds_total_ << " bonds, " << n_boundary_legs() << " boudnary legs." << endl << endl;

    cout << "Check neighbouring bonds of a site: " << endl << endl;
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

    cout << "Check neighbouring sites of bonds: " << endl << endl;
    for (int bond_i=0; bond_i<n_bonds_total_; bond_i++)
    {
        //cout << "The two end sites of bond " << bond_i << " are site " << bond_neighbour_sites(bond_i,0) << " and site " << bond_neighbour_sites(bond_i,1) <<endl;
        cout << "The neighbouring sites of bond " << bond_i << " are ";
        for (auto site_i : bond_neighbour_sites_[bond_i])
            cout << "site " << site_i << " ";
        cout << endl;
    }
    cout << endl;

    if (name_.find("cylinder")!=std::string::npos || name_.find("open")!=std::string::npos)
    {
        cout << "Check neighbouring boundary legs of a site: " << endl << endl;
        for (int site_i=0; site_i<n_sites_total_; site_i++)
        {
            for (auto boundary_i : site_neighbour_boundaries_[site_i])
            {
                if (boundary_i<0) continue;
                cout << "Site " << site_i << " connects to boundary " << boundary_i << endl;
            }
        }
        cout << endl;

        cout << "Check neighbouring boundary legs of a bond: " << endl;
        for (int bond_i=0; bond_i<n_bonds_total_; bond_i++)
        {
            for (auto boundary_i : bond_neighbour_boundaries_[bond_i])
            {
                if (boundary_i<0) continue;
                cout << "Bond " << bond_i << " connects to boundary " << boundary_i << endl;
            }
        }
        cout << endl;

        cout << "Check neighbouring site and bond of boundary: " << endl;
        for (int boundary_i=0; boundary_i<n_boundary_legs_; boundary_i++)
        {
            cout << "Boundary " << boundary_i << " connects to site " << boundary_neighbour_site_[boundary_i] << endl;
            cout << "Boundary " << boundary_i << " connects to bond " << boundary_neighbour_bond_[boundary_i] << endl;
        }
        cout << endl;
    }

    cout << "Finish printing lattice information!" << endl;
    cout << "\n==============================================\n" << endl;
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
            //site_neighbour_sites_[site_i][j]=site_coord_to_list(neigh_sites_coord[j]);
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

        bond_neighbour_sites_[bond_i][0]=site_coord_to_list(end_sites_coord[0]);
        bond_neighbour_sites_[bond_i][1]=site_coord_to_list(end_sites_coord[1]);
    }
}

//void Square_Lattice_Torus::print_lattice_inf() const
//{
//    cout << name_ << endl << "lattice size " << n_uc_[0] << "x" << n_uc_[1] << " and " << n_sites_total_ << " sites, " << n_bonds_total_ << " bonds." << endl;
//
//    cout << "Check neighbouring sites: " << endl << endl;
//    for (int site_i=0; site_i<n_sites_total_; site_i++)
//    {
//        cout << "Neighbour sites of site " << site_list_to_coord(site_i) << " are:" << endl;
//        //for (auto j : site_neighbour_sites_[site_i])
//        //{
//        //    cout << site_list_to_coord(j) << endl;
//        //}
//        cout << endl;
//    }
//    cout << endl;
//
//    cout << "Check neighbouring bonds: " << endl << endl;
//    for (int site_i=0; site_i<n_sites_total_; site_i++)
//    {
//        //cout << "Bonds connected to site " << site_list_to_coord(site_i) << " are:" << endl;
//        cout << "Bonds connected to site " << site_i << " are:" << endl;
//        for (auto j : site_neighbour_bonds_[site_i])
//        {
//            cout << j << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;
//
//    cout << "Check end sites of bonds: " << endl << endl;
//    for (int bond_i=0; bond_i<n_bonds_total_; bond_i++)
//    {
//        //cout << "The two end sites of bond " << bond_list_to_coord(bond_i) << " are site " << site_list_to_coord(bond_neighbour_sites(bond_i,0)) << " and site " << site_list_to_coord(bond_neighbour_sites(bond_i,1)) <<endl;
//        cout << "The two end sites of bond " << bond_i << " are site " << bond_neighbour_sites(bond_i,0) << " and site " << bond_neighbour_sites(bond_i,1) <<endl;
//    }
//    cout << endl;
//}



Square_Lattice_Cylinder::Square_Lattice_Cylinder(const std::array<int,2> &n_uc): Lattice_Base(1,2,n_uc,n_uc[1]*2)
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
            //site_i is at leftmost, then, we set the left neighbour site and bond to be -1
            if (neigh_sites_coord[j][0]<0) 
            {
                //site_neighbour_sites_[site_i][j]=-1;
                site_neighbour_bonds_[site_i][j]=-1;
                boundary_neighbour_site_[site_i_coord[1]]=site_i;
                boundary_neighbour_bond_[site_i_coord[1]]=-1;
                site_neighbour_boundaries_[site_i][j]=site_i_coord[1];
                continue;
            }

            //site_i is at rightmost, but right neigh bond still exist
            if (neigh_sites_coord[j][0]>=n_uc_[0])
            {
                //site_neighbour_sites_[site_i][j]=-2;
                int bond_no=bond_coord_to_list(neigh_bonds_coord[j]);
                site_neighbour_bonds_[site_i][j]=bond_no;
                boundary_neighbour_site_[n_uc_[1]+site_i_coord[1]]=-2;
                boundary_neighbour_bond_[n_uc_[1]+site_i_coord[1]]=bond_no;
                bond_neighbour_boundaries_[bond_no][1]=n_uc_[1]+site_i_coord[1];
                continue;
            }

            //site_neighbour_sites_[site_i][j]=site_coord_to_list(neigh_sites_coord[j]);
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

        bond_neighbour_sites_[bond_i][0]=site_coord_to_list(end_sites_coord[0]);
        bond_neighbour_sites_[bond_i][1]=site_coord_to_list(end_sites_coord[1]);
        // if the bond is at rightmost, set the bond_neighbour_sites_[1] to be -2
        if (end_sites_coord[1][0]>=n_uc_[0])
            bond_neighbour_sites_[bond_i][1]=-2;
    }
}


//
//Square_Lattice_Open
//
Square_Lattice_Open::Square_Lattice_Open(const std::array<int,2> &n_uc): Lattice_Base(1,2,n_uc,(n_uc[0]+n_uc[1])*2)
{
    name_="square lattice with open boundary condition";

    for (int sitei=0; sitei<n_sites_total_; sitei++)
    {
        Coordinate sitei_coord=site_list_to_coord(sitei);
        std::vector<Coordinate> neigh_sites_coord(n_bonds_to_one_site_,sitei_coord);
        std::vector<Coordinate> neigh_bonds_coord(n_bonds_to_one_site_);

        neigh_sites_coord[Left][0]=sitei_coord[0]-1;
        neigh_sites_coord[Up][1]=sitei_coord[1]+1;
        neigh_sites_coord[Right][0]=sitei_coord[0]+1;
        neigh_sites_coord[Down][1]=sitei_coord[1]-1;

        neigh_bonds_coord[Left]=Coordinate{sitei_coord[0]-1,sitei_coord[1],1};
        neigh_bonds_coord[Up]=Coordinate{sitei_coord[0],sitei_coord[1],0};
        neigh_bonds_coord[Right]=Coordinate{sitei_coord[0],sitei_coord[1],1};
        neigh_bonds_coord[Down]=Coordinate{sitei_coord[0],sitei_coord[1]-1,0};

        for (int j=0; j<n_bonds_to_one_site_; j++)
        {
            //sitei is at leftmost. Then we set left boundary connected to site
            if (neigh_sites_coord[j][0]<0)
            {
                //site_neighbour_sites_[sitei][j]=-1;
                site_neighbour_bonds_[sitei][j]=-1;
                boundary_neighbour_site_[sitei_coord[1]]=sitei;
                boundary_neighbour_bond_[sitei_coord[1]]=-1;
                site_neighbour_boundaries_[sitei][j]=sitei_coord[1];
                continue;
            }
            //sitei is at upmost. Then we set up boundary connected to bond
            if (neigh_sites_coord[j][1]>=n_uc_[1])
            {
                //site_neighbour_sites_[sitei][j]=-2;
                int bond_no=bond_coord_to_list(neigh_bonds_coord[j]);
                site_neighbour_bonds_[sitei][j]=bond_no;
                boundary_neighbour_site_[n_uc_[1]+sitei_coord[0]]=-2;
                boundary_neighbour_bond_[n_uc_[1]+sitei_coord[0]]=bond_no;
                bond_neighbour_boundaries_[bond_no][1]=n_uc_[1]+sitei_coord[0];
                continue;
            }
            //sitei is at rightmost. Then, we set right boundary connected to bond
            if (neigh_sites_coord[j][0]>=n_uc[0])
            {
                //site_neighbour_sites_[sitei][j]=-3;
                int bond_no=bond_coord_to_list(neigh_bonds_coord[j]);
                site_neighbour_bonds_[sitei][j]=bond_no;
                boundary_neighbour_site_[n_uc_[0]+n_uc_[1]+sitei_coord[1]]=-3;
                boundary_neighbour_bond_[n_uc_[0]+n_uc_[1]+sitei_coord[1]]=bond_no;
                bond_neighbour_boundaries_[bond_no][1]=n_uc_[0]+n_uc_[1]+sitei_coord[1];
                continue;
            }
            //sitei is at downmost. Then, we set down boundary connected to site
            if (neigh_sites_coord[j][1]<0)
            {
                site_neighbour_bonds_[sitei][j]=-4;
                boundary_neighbour_site_[n_uc_[0]+2*n_uc_[1]+sitei_coord[0]]=sitei;
                boundary_neighbour_bond_[n_uc_[0]+2*n_uc_[1]+sitei_coord[0]]=-4;
                site_neighbour_boundaries_[sitei][j]=n_uc_[0]+2*n_uc_[1]+sitei_coord[0];
                continue;
            }

            //site_neighbour_sites_[sitei][j]=site_coord_to_list(neigh_sites_coord[j]);
            site_neighbour_bonds_[sitei][j]=bond_coord_to_list(neigh_bonds_coord[j]);
        }

        for (int bondi=0; bondi<n_bonds_total_; bondi++)
        {
            Coordinate bondi_coord=bond_list_to_coord(bondi);
            std::array<Coordinate,2> end_sites_coord;

            if (bondi_coord[2]==0) //vertical bonds, direction from down to up
            {
                end_sites_coord[0]=Coordinate{bondi_coord[0],bondi_coord[1],0};
                end_sites_coord[1]=Coordinate{bondi_coord[0],bondi_coord[1]+1,0};
            }
            if (bondi_coord[2]==1) //horizontal bonds, direction from left to right
            {
                end_sites_coord[0]=Coordinate{bondi_coord[0],bondi_coord[1],0};
                end_sites_coord[1]=Coordinate{bondi_coord[0]+1,bondi_coord[1],0};
            }

            bond_neighbour_sites_[bondi][0]=site_coord_to_list(end_sites_coord[0]);
            bond_neighbour_sites_[bondi][1]=site_coord_to_list(end_sites_coord[1]);

            //if the bond is at upmost (rightmost), set the end_sites_coord[1] to be -2 (-3)
            if (end_sites_coord[1][1]>=n_uc_[1]) bond_neighbour_sites_[bondi][1]=-2;
            if (end_sites_coord[1][0]>=n_uc_[0]) bond_neighbour_sites_[bondi][1]=-3;

        }
    }
}



//
//Honeycomb_Lattice_Torus
//
Honeycomb_Lattice_Torus::Honeycomb_Lattice_Torus(const std::array<int,2> &n_uc): Lattice_Base(2,3,n_uc)
{
    name_="honeycomb lattice on torus";

    std::vector<std::vector<int>> site_neighbour_sites(n_sites_total_,std::vector<int>(n_bonds_to_one_site_));

    for (int sitei=0; sitei<n_sites_total_; sitei++)
    {
        Coordinate sitei_coord=site_list_to_coord(sitei);
        std::vector<Coordinate> neigh_sites_coord(n_bonds_to_one_site_);
        int sublattice_pm=2*sitei_coord[2]-1;

        neigh_sites_coord[NeighA]=Coordinate{sitei_coord[0]+sublattice_pm,sitei_coord[1]+sublattice_pm,1-sitei_coord[2]};
        neigh_sites_coord[NeighB]=Coordinate{sitei_coord[0],sitei_coord[1]+sublattice_pm,1-sitei_coord[2]};
        neigh_sites_coord[NeighC]=Coordinate{sitei_coord[0]+sublattice_pm,sitei_coord[1],1-sitei_coord[2]};

        std::vector<Coordinate> neigh_bonds_coord(n_bonds_to_one_site_);
        if (sitei_coord[2]==0)
        {
            neigh_bonds_coord[NeighA]=Coordinate{sitei_coord[0],sitei_coord[1],0};
            neigh_bonds_coord[NeighB]=Coordinate{sitei_coord[0],sitei_coord[1],1};
            neigh_bonds_coord[NeighC]=Coordinate{sitei_coord[0],sitei_coord[1],2};
        }
        else //sitei_coord[2]==1
        {
            neigh_bonds_coord[NeighA]=Coordinate{sitei_coord[0]+1,sitei_coord[1]+1,0};
            neigh_bonds_coord[NeighB]=Coordinate{sitei_coord[0],sitei_coord[1]+1,1};
            neigh_bonds_coord[NeighC]=Coordinate{sitei_coord[0]+1,sitei_coord[1],2};
        }

        for (int j=0; j<n_bonds_to_one_site_; j++)
        {
            site_neighbour_sites[sitei][j]=site_coord_to_list(neigh_sites_coord[j]);
            site_neighbour_bonds_[sitei][j]=bond_coord_to_list(neigh_bonds_coord[j]);
        }
    }

    //we assume bond is from subsite 0 to subsite 1
    for (int bondi=0; bondi<n_bonds_total_; bondi++)
    {
        Coordinate bondi_coord=bond_list_to_coord(bondi);

        bond_neighbour_sites_[bondi][0]=site_coord_to_list(Coordinate{bondi_coord[0],bondi_coord[1],0});
        bond_neighbour_sites_[bondi][1]=site_neighbour_sites[bond_neighbour_sites_[bondi][0]][bondi_coord[2]];
    }
}

