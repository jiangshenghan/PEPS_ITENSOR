
iPEPS_Simple_Update(int& d, int& D, const iLattice_Base& ilattice):
    d_(d),
    D_(D),
    ilattice_(ilattice),
    site_physical_legs_(ilattice.n_sites_uc()),
    site_virtual_legs_(ilattice.n_sites_uc(),vector<Index>(ilattice.n_site_to_links())),
    link_legs_(ilattice.n_links_uc()),
    site_tensors(ilattice.n_sites_uc()),
    link_tensors(ilattice.n_links_uc())
{
    //
    //Initialize legs which form site tensors in one unit cell
    //
    for (int i=0; i<ilattice_.n_sites_uc(); i++)
    {
        site_physical_legs_[i]=Index(nameint("Site_physical_leg",i),d_,Site);
        for(int j=0; j<ilattice_.n_site_to_links(); j++)
        {
            site_virtual_legs_[i][j]=Index(nameint("Site_virtual_leg",i*ilattice_.n_sites_uc()+j),D_,Link);
        }
    }

    //
    //Initialize legs forming link tensors in one unit cell
    for (int )
}

