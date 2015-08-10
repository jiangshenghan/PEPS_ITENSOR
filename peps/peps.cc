
#include "peps.h"

PEPS::PEPS(int d, int D, const Lattice_Base& lattice):
    d_(d),
    D_(D),
    lattice_(lattice),
    site_legs_(lattice.n_sites()),
    bulk_link_legs_(lattice.n_bulk_links()),
    boundary_link_legs_(lattice.n_boundary_links()),
    site_indexsets_(lattice.n_sites_()),
    site_tensors_(lattice.n_sites_())
{
    //Initialize legs
    for (int i=0; i<lattice_.n_sites_(); i++)
    {
        site_legs_[i]=Index(nameint("Site_leg",i),d_,Site);
    }
    for (int i=0; i<lattice_.n_bulk_links(); i++)
    {
        bulk_link_legs_[i]=Index(nameint("Bulk_link_leg",i),D_,Link);
    }
    for (int i=0; i<lattice_.n_boundary_links(); i++)
    {
        boundary_link_legs_[i]=Index(nameint("Boundary_link_leg",i),D_,Link);
    }

    //Initialze tensors
    for (int i=0; i<lattice_.n_sites(); i++)
    {
        //Adding site legs and link legs (bulk and boundary to ith siteindex sets)
        site_indexsets_[i].addIndex(site_legs_[i]);
        for (int j=0; j<lattice_.n_site_to_links(); j++)
        {
            int add_bulk_leg_j=lattice_.site_to_bulk_links(i,j),
                add_boundary_leg_j=lattice_.site_to_boundary_links(i,j);
            if (add_bulk_leg_j!=-1)
            {
                site_indexsets_[i].addIndex(bulk_link_legs_[add_bulk_leg_j]);
            }
            if (add_boundary_leg_j!=-1)
            {
                site_indexsets_[i].addIndex(boundary_link_legs_[add_boundary_legs_j]);
            }
        }

        site_tensors_[i]=ITensor(site_indexsets_[i]);
    }
}
