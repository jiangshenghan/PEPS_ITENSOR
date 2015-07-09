
#include "lattice.h"

Square_Lattice::Square_Lattice(int& lx, int& ly, std::array<Boundary,2>& boundary_condition={Open,Open}): n_site_to_links_(4), lx_(lx), ly_(ly), boundary_condtion_(boundary_condition) 
{
    n_bulk_links_=2*lx_*ly_-lx_-ly_;

    n_boundary_links_=0;
    if (boundary_condition_[0]==Open)
        n_boundary_links_+=2*ly_;
    if (boundary_condtions_[1]==Open)
        n_boundary_links+=2*lx_;
}
