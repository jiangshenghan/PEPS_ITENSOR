
#include "lattice.h"

Lattice_Base::Lattice_Base(int& n_site_to_links, int &n_sites, int& n_bulk_links, int& n_boundary_links):
    n_site_to_links_(n_site_to_links),
    n_sites_(n_sites),
    n_bulk_links_(n_bulk_links),
    n_boundary_links_(n_boundary_links),
    neighbour_sites_(n_sites,vector<int>(n_site_to_links)),
    site_to_bulk_links_(n_sites,vector<int>(n_site_to_links)),
    site_to_boundary_links_(n_sites,vector<int>(n_site_to_links)),
    bulk_link_to_sites_(n_bulk_links_),
    boundary_links_to_sites_(n_boundary_links_)
{
}

Square_Lattice::Square_Lattice(int& lx, int& ly): 
    Lattice_Base(4,lx*ly,(lx-1)*ly+lx*(ly-1),2*lx+2*ly),
    lx_(lx),
    ly_(ly),
    site_to_coord_list_(lx*ly),
    coord_to_site_list_(lx,vector<int>(ly))
{
    for (int site_i=0; site_i<n_sites_; site_i++)
    {
        site_to_coord_list_[site_i][0]=i/lx;
        site_to_coord_list_[site_i][1]=i%lx;
    }

    for (int row_i=0; row_i<lx; row_i++)
        for(int col_i=0; col_i<ly; col_i++)
            coord_to_site_list[row_i][col_i]=row_i*lx+col_i;

    for (int site_i=0; site_i<n_sites_; site_i++)
    {
        neighbour_sites_[0]=
    }

}
