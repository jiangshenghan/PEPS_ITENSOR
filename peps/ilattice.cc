
#include "ilattice.h"

iLattice_Base::iLattice_Base(int& n_sites_uc, int& n_links_uc):
    n_sites_uc_(n_sites_uc),
    n_links_uc_(n_links_uc),
    n_site_to_links_(2*n_links_uc/n_sites_uc)
    neigbour_sites_(n_sites_uc,vector<Coordinate>(n_site_to_links)),
    site_to_links_(n_sites_uc,vector<Coordinate>(n_site_to_links)),
    link_to_sites_(n_links_uc)
{
}



iSquare_Lattice::iSquare_Lattice():
    iLattice_Base(2,4),
{
    neighbour_sites_[0][Left]=Coordinate{ {-1,-1,1} };
    neighbour_sites_[0][Up]=Coordinate{ {-1,0,1} };
    neighbour_sites_[0][Right]=Coordinate{ {0,0,1} };
    neighbour_sites_[0][Down]=Coordinate{ {0,-1,1} };

    neighbour_sites_[1][Left]=Coordinate{ {0,0,0} };
    neighbour_sites_[1][Up]=Coordinate{ {0,1,0} };
    neighbour_sites_[1][Right]=Coordinate{ {1,1,0} };
    neighbour_sites_[1][Down]=Coordinate{ {1,0,0} };


    site_to_links_[0][Left]=Coordinate{ {0,0,0} };
    site_to_links_[0][Up]=Coordinate{ {0,0,1} };
    site_to_links_[0][Right]=Coordinate{ {0,0,2} };
    site_to_links_[0][Down]=Coordinate{ {0,0,3} };

    site_to_links_[1][Left]=Coordinate{ {0,0,2} };
    site_to_links_[1][Up]=Coordinate{ {0,1,3} };
    site_to_links_[1][Right]=Coordinate{ {1,1,0} };
    site_to_links_[1][Down]=Coordinate{ {1,0,1} };


    link_to_sites_[Left][0]=Coordinate{ {0,0,0} };
    link_to_sites_[Left][1]=Coordinate{ {-1,-1,1} };

    link_to_sites_[Up][0]=Coordinate{ {0,0,0} };
    link_to_sites_[Up][1]=Coordinate{ {-1,0,1} };

    link_to_sites_[Right][0]=Coordinate{ {0,0,0} };
    link_to_sites_[Right][1]=Coordinate{ {0,0,1} };

    link_to_sites_[Down][0]=Coordinate{ {0,0,0} };
    link_to_sites_[Down][1]=Coordinate{ {0,-1,1} };
}
