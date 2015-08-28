
#include "ipeps.h"

iPEPS(int& d, int& D, const Lattice_uc_Base& lattice_uc):
    d_(d),
    D_(D),
    lattice_uc_(lattice_uc),
    site_legs_(lattice_uc.n_sites_uc()),
    link_legs_(lattice_uc.n_links_uc()),
    site_tensors_(lattice_uc.n_sites_uc()),
    link_tensors_(lattice_uc.n_links_uc())
{
    //
    //Initialize site and link legs 
    //
    for (int i=0; i<lattice_uc_.n_sites_uc(); i++)
    {
        site_legs_[i]=Index(nameint("Site_leg",i),d_,Site);
    }

    for (int i=0; i<lattice_uc_.n_links_uc(); i++)
    {
        link_legs_[i][0]=Index(nameint("Link_leg",2*i),D_,Link);
        link_legs_[i][1]=Index(nameint("Link_leg",2*i+1),D_,Link);
    }

    //
    //Initialize site tensors
    //
    for (int i=0; i<lattice_uc_.n_sites_uc(); i++)
    {
        IndexSet<Index> site_indexset;
        site_indexset.addIndex(site_legs_[i]);

        for (int j=0; j<lattice_uc_.n_links_to_one_site(); j++)
        {
            //link_ind stores the sublattice index of jth link which connects to site i
            Coordinate link_ind=lattice_uc_.site_neighbour_links(i,j);
            //link_end_ind tells which end of the above link_ind is connected to site i
            int link_end_ind;

            if (lattice_uc_.link_end_sites(link_ind[2],0)==i)
                link_end_ind=0;
            else 
                link_end_ind=1;

            site_indexset.addIndex(link_legs_[link_ind][link_end_ind]);
        }

        site_tensors_[i]=ITensor(site_indexset);
    }

    //
    //Initialize link tensors
    //
    for (int i=0; i<lattice_uc_.n_links_uc(); i++)
    {
        link_tensors_[i]=ITensor(link_tensors_[i][0],link_tensors_[i][1]);
    }
}

