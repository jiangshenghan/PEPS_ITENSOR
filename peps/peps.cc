
#include "peps.h"

template <class TensorT>
PEPSt<TensorT>::PEPSt(const Lattice_Base &lattice):
    lattice_ptr_(new Lattice_Base(lattice)),
    indexset_ptr_(new PEPSt_IndexSet_Base<IndexT>),
    site_tensors_(lattice.n_sites_total()),
    bond_tensors_(lattice.n_bonds_total()),
    boundary_tensors_(lattice.n_boundary_legs()),
    name_(lattice.name()+' ')
{}
template
PEPSt<ITensor>::PEPSt(const Lattice_Base &lattice);
template
PEPSt<IQTensor>::PEPSt(const Lattice_Base &lattice);

template <class TensorT>
PEPSt<TensorT>::PEPSt(const Lattice_Base &lattice, const PEPSt_IndexSet_Base<IndexT> &index_set):
    d_(index_set.d()),
    D_(index_set.D()),
    lattice_ptr_(new Lattice_Base(lattice)),
    indexset_ptr_(new PEPSt_IndexSet_Base<IndexT>(index_set)),
    site_tensors_(lattice.n_sites_total()),
    bond_tensors_(lattice.n_bonds_total()),
    boundary_tensors_(lattice.n_boundary_legs()),
    name_(lattice.name()+' '+index_set.name())
{
    //set two col of boundary tensors if cylinder geometry (left and right boundary)
    //each col contains n_uc_[1] tensors

    new_site_tensors();
    new_bond_tensors();
    new_boundary_tensors();
    //random_site_tensors();
}
template
PEPSt<ITensor>::PEPSt(const Lattice_Base &lattice, const PEPSt_IndexSet_Base<ITensor::IndexT> &index_set);
template
PEPSt<IQTensor>::PEPSt(const Lattice_Base &lattice, const PEPSt_IndexSet_Base<IQTensor::IndexT> &index_set);

template <class TensorT>
PEPSt<TensorT>::PEPSt(const Lattice_Base &lattice, const PEPSt_IndexSet_Base<IndexT> &index_set, std::vector<TensorT> site_tensors_uc, std::vector<TensorT> bond_tensors_uc, double mu_12):
    d_(index_set.d()),
    D_(index_set.D()),
    lattice_ptr_(new Lattice_Base(lattice)),
    indexset_ptr_(new PEPSt_IndexSet_Base<IndexT>(index_set)),
    site_tensors_(lattice.n_sites_total()),
    bond_tensors_(lattice.n_bonds_total()),
    boundary_tensors_(lattice.n_boundary_legs()),
    name_(lattice.name()+index_set.name())
{
    //set two col of boundary tensors if cylinder geometry (left and right boundary)
    //each col contains n_uc_[1] tensors

    new_site_tensors();
    new_bond_tensors();
    new_boundary_tensors();

    //for (const auto &tensor : site_tensors_) cout << tensor;
    //for (const auto &tensor : bond_tensors_) cout << tensor;
    //for (const auto &tensor : boundary_tensors_) cout << tensor;

    //for (const auto &tensor : site_tensors_uc) PrintDat(tensor);
    //for (const auto &tensor : bond_tensors_uc) PrintDat(tensor);

    generate_site_tensors(site_tensors_uc);
    generate_bond_tensors(bond_tensors_uc,mu_12);


}
template
PEPSt<ITensor>::PEPSt(const Lattice_Base &lattice, const PEPSt_IndexSet_Base<ITensor::IndexT> &index_set, std::vector<ITensor> site_tensors_uc, std::vector<ITensor> bond_tensors_uc, double mu_12);
template
PEPSt<IQTensor>::PEPSt(const Lattice_Base &lattice, const PEPSt_IndexSet_Base<IQTensor::IndexT> &index_set, std::vector<IQTensor> site_tensors_uc, std::vector<IQTensor> bond_tensors_uc, double mu_12);


template<class TensorT>
void PEPSt<TensorT>::generate_site_tensors(std::vector<TensorT> site_tensors_uc)
{
    for (int site_i=0; site_i<n_sites_total(); site_i++)
    {
        auto sublattice_i=lattice_ptr_->site_list_to_coord(site_i).at(2);
        tensor_assignment(site_tensors_[site_i],site_tensors_uc[sublattice_i]);

        //PrintDat(site_tensors_[site_i]);
    }
}
template
void PEPSt<ITensor>::generate_site_tensors(std::vector<ITensor> site_tensors_uc);
template
void PEPSt<IQTensor>::generate_site_tensors(std::vector<IQTensor> site_tensors_uc);


template<class TensorT>
void PEPSt<TensorT>::generate_bond_tensors(std::vector<TensorT> bond_tensors_uc, double mu_12)
{
    //for symmetric peps with half spin per site, we have 
    //B_{(x,y,i)}=\eta_{12}^{x}_{\alpha\alpha'}. (B_i)_{\alpha'...}, if leg a of B_i connect site with different y coordinate
    //So we need to double the unit cell in x direction if \mu_12==-1
    if (name_.find("half spin")!=std::string::npos && std::abs(mu_12+1)<EPSILON)
    {
        for (int bond_i=0; bond_i<n_bonds_total(); bond_i++)
        {
            auto sublattice_i=lattice_ptr_->bond_list_to_coord(bond_i)[2];
            tensor_assignment(bond_tensors_[bond_i],bond_tensors_uc[sublattice_i]);
            //if x is even, then \eta_12^x=I, so we leave bond tensors invariant
            if (lattice_ptr_->bond_list_to_coord(bond_i)[0]%2==0) continue;
            auto end_sites=lattice_ptr_->bond_neighbour_sites(bond_i);
            for (auto end_site : end_sites)
            {
                if (lattice_ptr_->site_list_to_coord(end_site)[1]%2!=0)
                {
                    obtain_tensor_after_eta_action(mu_12,bond_tensors_[bond_i],commonIndex(bond_tensors(bond_i),site_tensors(end_site)));
                }
            }

        }

        //if (name_.find("square")!=std::string::npos )
        //{
        //    for (int bond_i=0; bond_i<n_bonds_total(); bond_i++)
        //    {
        //        auto sublattice_i=lattice_ptr_->bond_list_to_coord(bond_i)[2];

        //        tensor_assignment(bond_tensors_[bond_i],bond_tensors_uc[sublattice_i]);

        //        auto end_sites=lattice_ptr_->bond_neighbour_sites(bond_i);
        //        std::array<Coordinate,2> end_sites_coord={lattice_ptr_->site_list_to_coord(end_sites[0]), lattice_ptr_->site_list_to_coord(end_sites[1])};

        //        if (std::abs(end_sites_coord[0][1]-end_sites_coord[1][1])%2==1 && std::abs(end_sites_coord[0][0]%2==1))
        //        {
        //            auto eta_12=eta_from_mu(mu_12,dag(bond_tensors_[bond_i].indices()[0]));
        //            bond_tensors_[bond_i]*=eta_12;
        //            bond_tensors_[bond_i].noprime();
        //        }
        //    }
        //}

        return;
    }

    for (int bond_i=0; bond_i<n_bonds_total(); bond_i++)
    {
        auto sublattice_i=lattice_ptr_->bond_list_to_coord(bond_i).at(2);
        tensor_assignment(bond_tensors_[bond_i],bond_tensors_uc[sublattice_i]);
    }
}
template
void PEPSt<ITensor>::generate_bond_tensors(std::vector<ITensor> bond_tensors_uc, double mu_12);
template
void PEPSt<IQTensor>::generate_bond_tensors(std::vector<IQTensor> bond_tensors_uc, double mu_12);


//According to the library, ITensor is constructed by IndexSet<Index>(indexset), while IQTensor is constructed by indexset directly
template <>
void PEPSt<ITensor>::construct_tensor(ITensor &tensor, std::vector<ITensor::IndexT> &indexset)
{
    tensor=ITensor(IndexSet<Index>(indexset));
    return;
}
template <>
void PEPSt<IQTensor>::construct_tensor(IQTensor &tensor, std::vector<IQTensor::IndexT> &indexset)
{
    tensor=IQTensor(indexset);
    return;
}


template <class TensorT>
void PEPSt<TensorT>::new_site_tensors()
{
    int site_i=0;
    for(auto &tensor : site_tensors_)
    {
        //tensor_indices stores both physical legs and virtual legs for site_tensors_[site_i]
        std::vector<IndexT> tensor_indices;
        tensor_indices.push_back(indexset_ptr_->phys_legs(site_i));

        for (int j=0; j<n_bonds_to_one_site(); j++)
        {
            IndexT neigh_virt_ind=indexset_ptr_->virt_legs(j+site_i*n_bonds_to_one_site());
            tensor_indices.push_back(neigh_virt_ind);
        }

        construct_tensor(tensor,tensor_indices);

        site_i++;
    }

    return;
}
template
void PEPSt<ITensor>::new_site_tensors();
template
void PEPSt<IQTensor>::new_site_tensors();


template <class TensorT>
void PEPSt<TensorT>::new_bond_tensors()
{
    int bond_i=0;
    int boundary_legs_no=0;

    for (auto &tensor : bond_tensors_)
    {
        std::vector<IndexT> tensor_indices;

        for (int endi=0; endi<n_sites_to_one_bond(); endi++)
        {
            int endi_site=lattice_ptr_->bond_neighbour_sites(bond_i,endi);

            //consider the boundary bonds case separetely
            if (endi_site<0)
            {
                IndexT endi_index=indexset_ptr_->virt_legs(n_sites_total()*n_bonds_to_one_site()+boundary_legs_no);
                endi_index.dag();
                tensor_indices.push_back(endi_index);
                boundary_legs_no++;
                continue;
            }

            std::vector<int> endi_site_neigh=lattice_ptr_->site_neighbour_bonds(endi_site);
            int legi=std::find(endi_site_neigh.begin(),endi_site_neigh.end(),bond_i)-endi_site_neigh.begin();
            assert(legi<n_bonds_to_one_site());

            IndexT endi_index=indexset_ptr_->virt_legs(endi_site*n_bonds_to_one_site()+legi);
            endi_index.dag();


            //for (const auto &i : endi_site_neigh)
            //    cout << i << " ";
            //cout << endl;
            //cout << bond_i << " " << endi << ": site" << endi_site << " leg" << legi << endl;
            
            tensor_indices.push_back(endi_index);
        }

        construct_tensor(tensor,tensor_indices);

        //cout << tensor << endl;

        bond_i++;
    }
}
template
void PEPSt<ITensor>::new_bond_tensors();
template
void PEPSt<IQTensor>::new_bond_tensors();

template <class TensorT>
void PEPSt<TensorT>::new_boundary_tensors()
{
    if (name_.find("torus")!=std::string::npos) return;

    if (name_.find("cylinder")!=std::string::npos || name_.find("open")!=std::string::npos)
    {
        for (int boundaryi=0; boundaryi<lattice_ptr_->n_boundary_legs(); boundaryi++)
        {
            int site_no=lattice_ptr_->boundary_neighbour_site(boundaryi);
            int bond_no=lattice_ptr_->boundary_neighbour_bond(boundaryi);
            //boundaryi is connected to site_no
            if (bond_no<0 && site_no>=0)
            {
                auto &neighbour_boundaries=lattice_ptr_->site_neighbour_boundaries(site_no);
                int boundary_ind_no=std::find(neighbour_boundaries.begin(),neighbour_boundaries.end(),boundaryi)-neighbour_boundaries.begin();
                assert(boundary_ind_no<n_bonds_to_one_site());
                boundary_ind_no+=site_no*n_bonds_to_one_site();
                auto boundary_leg=dag(indexset_ptr_->virt_legs(boundary_ind_no));
                boundary_tensors_[boundaryi]=TensorT(boundary_leg);
                continue;
            }
            //boundaryi is connected to bond_no
            if (site_no<0 && bond_no>=0)
            {
                auto &neighbour_boundaries=lattice_ptr_->bond_neighbour_boundaries(bond_no);
                int boundary_ind_no=std::find(neighbour_boundaries.begin(),neighbour_boundaries.end(),boundaryi)-neighbour_boundaries.begin();
                assert(boundary_ind_no<n_sites_to_one_bond());
                auto boundary_leg=dag(bond_tensors_[bond_no].indices()[boundary_ind_no]);
                boundary_tensors_[boundaryi]=TensorT(boundary_leg);
                continue;
            }

            cout << "Incorrect boundary setting!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    //if (name_.find("cylinder")!=std::string::npos || name_.find("open")!=std::string::npos)
    //{
    //    std::vector<bool> boundary_tensor_created(lattice_ptr_->n_boundary_legs(),false);
    //    for (int boundary_i=0; boundary_i<lattice_ptr_->n_boundary_legs(); boundary_i++)
    //    {
    //        if (boundary_tensor_created[boundary_i]) continue;
    //        int site_i=lattice_ptr_->boundary_end_site(boundary_i);
    //        std::vector<int> neighbour_sites=lattice_ptr_->site_neighbour_sites(site_i);
    //        //begin_iter is the position start to seach the "invalid" neighbour sites
    //        //should be updated to find next "invalid" one
    //        auto begin_iter=neighbour_sites.begin();
    //        for (const auto &boundary_id : lattice_ptr_->site_neighbour_boundary(site_i))
    //        {
    //            auto neighbour_boundary_iter=std::find_if(begin_iter,neighbour_sites.end(), [](int i){ return (i<0); });
    //            int neighbour_boundary_no=neighbour_boundary_iter-neighbour_sites.begin();

    //            //For cylinder geometry, we get left and right boundaries
    //            if (name_.find("cylinder")!=std::string::npos)
    //            {
    //                //left boundary legs are associated with sites, labeled as -1
    //                if (*neighbour_boundary_iter==-1)
    //                {
    //                    auto boundary_leg=dag(indexset_ptr_->virt_legs(site_i*n_bonds_to_one_site()+neighbour_boundary_no));
    //                    boundary_tensors_[boundary_id]=TensorT(boundary_leg);
    //                }
    //                //right boundary legs are associated with bonds, labeled as -2
    //                if (*neighbour_boundary_iter==-2)
    //                {
    //                    int bond_i=lattice_ptr_->site_neighbour_bonds(site_i,neighbour_boundary_no);
    //                    int bond_boundary_i=(lattice_ptr_->bond_neighbour_sites(bond_i,0)==-2 ? 0 : 1);
    //                    //if (lattice_ptr_->bond_neighbour_sites(bond_i,0)==-2) bond_boundary_i=0;
    //                    //if (lattice_ptr_->bond_neighbour_sites(bond_i,1)==-2) bond_boundary_i=1;
    //                    auto boundary_leg=dag(bond_tensors_[bond_i].indices()[bond_boundary_i]);
    //                    boundary_tensors_[boundary_id]=TensorT(boundary_leg);
    //                }
    //            }

    //            //For system with OBC, we get left, up, right and down boundaries
    //            if (name_.find("open")!=std::string::npos)
    //            {
    //                //left boundaries are associated with sites, labeled as -1
    //                //down boundaries are associated with sites, labeled as -4
    //                if (*neighbour_boundary_iter==-1 || *neighbour_boundary_iter==-4)
    //                {
    //                    auto boundary_leg=dag(indexset_ptr_->virt_legs(site_i*n_bonds_to_one_site()+neighbour_boundary_no));
    //                    boundary_tensors_[boundary_id]=TensorT(boundary_leg);
    //                }
    //                //up boundaries are associated with bonds, labeled as -2
    //                //right boundaries are associated with bonds, labeled as -3
    //                if (*neighbour_boundary_iter==-2 || *neighbour_boundary_iter==-3)
    //                {
    //                    int bond_i=lattice_ptr_->site_neighbour_bonds(site_i,neighbour_boundary_no);
    //                    int bond_boundary_i=(lattice_ptr_->bond_neighbour_sites(bond_i,0)<0 ? 0 : 1);
    //                    auto boundary_leg=dag(bond_tensors_[bond_i].indices()[bond_boundary_i]);
    //                    boundary_tensors_[boundary_id]=TensorT(boundary_leg);
    //                }
    //            }

    //            boundary_tensor_created[boundary_id]=true;
    //            begin_iter=neighbour_boundary_iter;
    //            begin_iter++;
    //        }
    //    }

    //    //int boundary_no=0;
    //    //for (int sitei=0; sitei<n_sites_total(); sitei++)
    //    //{
    //    //    std::vector<int> neighbour_bonds=lattice_ptr_->site_neighbour_bonds(sitei);
    //    //    auto begin_iter=neighbour_bonds.begin();
    //    //    //sitei may have many boundary legs
    //    //    while (begin_iter!=neighbour_bonds.end())
    //    //    {
    //    //        auto found_iter=std::find(begin_iter,neighbour_bonds.end(),-1);
    //    //        if (found_iter==neighbour_bonds.end()) break;

    //    //        int neigh_bond_no=found_iter-begin_iter;
    //    //        boundary_tensors_[boundary_no]=TensorT(dag(indexset_ptr_->virt_legs(sitei*n_bonds_to_one_site()+neigh_bond_no)));

    //    //        //cout << boundary_tensors_[0][boundary_no];

    //    //        found_iter++;
    //    //        begin_iter=found_iter;
    //    //        boundary_no++;
    //    //    }
    //    //}


    //    //for (int bondi=0; bondi<n_bonds_total(); bondi++)
    //    //{
    //    //    for (int legi=0; legi<2; legi++)
    //    //    {
    //    //        if (lattice_ptr_->bond_neighbour_sites(bondi,legi)==-2)
    //    //        {
    //    //            boundary_tensors_[boundary_no]=TensorT(dag(bond_tensors_[bondi].indices()[legi]));

    //    //            //cout << boundary_tensors_[1][boundary_no];

    //    //            boundary_no++;
    //    //        }
    //    //    }
    //    //}
    //}
}
template
void PEPSt<ITensor>::new_boundary_tensors();
template 
void PEPSt<IQTensor>::new_boundary_tensors();


template <class TensorT>
std::vector<TensorT> PEPSt<TensorT>::combined_site_tensors() const
{
    std::vector<TensorT> combined_site_tensors;
    for (int sitei=0; sitei<this->n_sites_total(); sitei++)
    {
        auto site_tensor=this->site_tensors(sitei);
        //cout << "Original tensor:" << endl << sitei << endl << site_tensor << endl;

        //multiply neighbouring bond tensors 
        for (int neighi=0; neighi<this->lattice().n_bonds_to_one_site(); neighi++)
        {
            int bond_no=this->lattice().site_neighbour_bonds(sitei,neighi);
            if (bond_no<0) continue;
            //For bulk bond, we only multiply those start from the site to avoid double counting. Namely bond_neighbour_sites(bond_no,0)==sitei
            //For boundary bond, we will always absorb the bond
            //TODO: consider boundary bonds connect multiple sites
            if (this->lattice().bond_neighbour_sites(bond_no,0)==sitei || this->lattice().bond_neighbour_sites(bond_no,0)<0)
            {
                site_tensor*=this->bond_tensors(bond_no);
            }
        }
        //multiply neighbouring boundary tensors
        for (const auto &boundary_no : this->lattice().site_neighbour_boundaries(sitei))
        {
            if (boundary_no<0) continue;
            site_tensor*=this->boundary_tensors(boundary_no);
        }
        clean(site_tensor);

        combined_site_tensors.push_back(site_tensor);

        //cout << "Single layer tensor " << sitei << endl << site_tensor << endl;
    }

    return combined_site_tensors;
}
template
std::vector<ITensor> PEPSt<ITensor>::combined_site_tensors() const;
template
std::vector<IQTensor> PEPSt<IQTensor>::combined_site_tensors() const;

template <typename TensorT>
void PEPSt<TensorT>::obtain_combined_site_tensors()
{
    combined_site_tensors_=combined_site_tensors();
}
template 
void PEPSt<ITensor>::obtain_combined_site_tensors();
template
void PEPSt<IQTensor>::obtain_combined_site_tensors();

//before reading, we should init PEPSt use lattice
template <class TensorT>
void PEPSt<TensorT>::read(std::istream &s)
{
    //if (lattice_ptr_->is_empty()) lattice_ptr_->read(s);
    indexset_ptr_->read(s);

    for (auto &tensor : site_tensors_) tensor.read(s);
    for (auto &tensor : bond_tensors_) tensor.read(s);
    for (auto &tensor : boundary_tensors_) tensor.read(s);

    int nlength;
    s.read((char*)&nlength,sizeof(nlength));
    auto newname=std::unique_ptr<char[]>(new char[nlength+1]);
    s.read(newname.get(),nlength+1);
    name_=std::string(newname.get());

    d_=indexset_ptr_->d();
    D_=indexset_ptr_->D();
}
template
void PEPSt<ITensor>::read(std::istream &s);
template
void PEPSt<IQTensor>::read(std::istream &s);

template <class TensorT>
void PEPSt<TensorT>::write(std::ostream &s) const
{
    //lattice_ptr_->write(s);
    indexset_ptr_->write(s);

    for (const auto &tensor : site_tensors_) tensor.write(s);
    for (const auto &tensor : bond_tensors_) tensor.write(s);
    for (const auto &tensor : boundary_tensors_) tensor.write(s);

    const int nlength=name_.length();
    s.write((char*)&nlength,sizeof(nlength));
    s.write(name_.c_str(),nlength+1);
}
template
void PEPSt<ITensor>::write(std::ostream &s) const;
template
void PEPSt<IQTensor>::write(std::ostream &s) const;


template <class TensorT>
Tnetwork_Storage<TensorT> peps_to_tnetwork_storage(const PEPSt<TensorT> &peps)
{
    Tnetwork_Storage<TensorT> tnetwork_storage;
    
    if (peps.name().find("torus")!=std::string::npos) tnetwork_storage._boundary_condition=1;

    //Translate from square lattice
    if (peps.name().find("square")!=std::string::npos) 
    {
        tnetwork_storage._tnetwork_type=1;
        tnetwork_storage._n_subl=1;

        int Lx=peps.n_uc()[0], Ly=peps.n_uc()[1];
        tnetwork_storage._Lx=Lx;
        tnetwork_storage._Ly=Ly;
        //Print(tnetwork_storage._Lx);
        //Print(tnetwork_storage._Ly);

        tnetwork_storage._tensor_list.set_size(peps.n_sites_total());
        auto tensor_list=peps.combined_site_tensors();
        for (int sitei=0; sitei<peps.n_sites_total(); sitei++) 
        {
            tnetwork_storage._tensor_list(sitei)=tensor_list[sitei];
            //Print(sitei);
            //Print(tnetwork_storage._tensor_list(sitei));
        }

        tnetwork_storage._coor_to_siteind.set_size(Lx,Ly);
        for (int x=0; x<Lx; x++)
        {
            for (int y=0; y<Ly; y++)
            {
                tnetwork_storage._coor_to_siteind(x,y).set_size(tnetwork_storage._n_subl);
                for (int subli=0; subli<tnetwork_storage._n_subl; subli++)
                    tnetwork_storage._coor_to_siteind(x,y)(subli)=peps.lattice().site_coord_to_list(x,y,subli);
            }
        }
    }

    //translate from kagome cirac lattice
    if (peps.name().find("kagome cirac")!=std::string::npos)
    {
        tnetwork_storage._tnetwork_type=5;
        tnetwork_storage._n_subl=5;

        int Lx=peps.n_uc()[1],
            Ly=peps.n_uc()[0];
        tnetwork_storage._Lx=Lx;
        tnetwork_storage._Ly=Ly;

        tnetwork_storage._tensor_list.set_size(peps.n_sites_total()+peps.n_bonds_total());
        tnetwork_storage._coor_to_siteind.set_size(Lx,Ly);
        for (int x=0; x<Lx; x++)
        {
            for (int y=0; y<Ly; y++)
            {
                tnetwork_storage._coor_to_siteind(x,y).set_size(tnetwork_storage._n_subl);
                for (int subli=0; subli<tnetwork_storage._n_subl; subli++)
                {
                    int siteind=(x+y*Lx)*tnetwork_storage._n_subl+subli;
                    tnetwork_storage._coor_to_siteind(x,y)(subli)=siteind;
                    if (subli==0)
                    {
                        int original_siteind=peps.lattice().site_coord_to_list(y,-x,1);
                        tnetwork_storage._tensor_list(siteind)=peps.site_tensors(original_siteind);
                    }
                    if (subli==1)
                    {
                        int original_siteind=peps.lattice().site_coord_to_list(y,-x,2);
                        tnetwork_storage._tensor_list(siteind)=peps.site_tensors(original_siteind);
                    }
                    if (subli==2)
                    {
                        int original_siteind=peps.lattice().site_coord_to_list(y,-x,0);
                        tnetwork_storage._tensor_list(siteind)=peps.site_tensors(original_siteind);
                    }
                    if (subli==3)
                    {
                        int original_siteind=peps.lattice().bond_coord_to_list(y,-x,0);
                        tnetwork_storage._tensor_list(siteind)=peps.bond_tensors(original_siteind);
                    }
                    if (subli==4)
                    {
                        int original_siteind=peps.lattice().bond_coord_to_list(y-1,-x-1,1);
                        tnetwork_storage._tensor_list(siteind)=peps.bond_tensors(original_siteind);
                    }
                }
            }
        }
    }

    //translate from kagome normal lattice
    if (peps.name().find("kagome normal")!=std::string::npos)
    {
        tnetwork_storage._tnetwork_type=8;
        tnetwork_storage._n_subl=3;

        int Lx=peps.n_uc()[1],
            Ly=peps.n_uc()[0];
        tnetwork_storage._Lx=Lx;
        tnetwork_storage._Ly=Ly;

        tnetwork_storage._tensor_list.set_size(peps.n_sites_total());
        auto tensor_list=peps.combined_site_tensors();
        for (int sitei=0; sitei<peps.n_sites_total(); sitei++)
        {
            tnetwork_storage._tensor_list(sitei)=tensor_list[sitei];
        }

        tnetwork_storage._coor_to_siteind.set_size(Lx,Ly);
        for (int x=0; x<Lx; x++)
        {
            for (int y=0; y<Ly; y++)
            {
                tnetwork_storage._coor_to_siteind(x,y).set_size(tnetwork_storage._n_subl);
                for (int subli=0; subli<tnetwork_storage._n_subl; subli++)
                {
                    tnetwork_storage._coor_to_siteind(x,y)(subli)=peps.lattice().site_coord_to_list(x,y,subli);
                }
            }
        }
    }

    return tnetwork_storage;
}
template
Tnetwork_Storage<ITensor> peps_to_tnetwork_storage(const PEPSt<ITensor> &peps);
template
Tnetwork_Storage<IQTensor> peps_to_tnetwork_storage(const PEPSt<IQTensor> &peps);



