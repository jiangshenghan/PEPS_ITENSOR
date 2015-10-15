
#include "peps.h"

template <class TensorT>
PEPSt<TensorT>::PEPSt(const Lattice_Base &lattice):
    lattice_(lattice),
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
PEPSt<TensorT>::PEPSt(const Lattice_Base &lattice, PEPSt_IndexSet_Base<IndexT> &index_set):
    d_(index_set.d()),
    D_(index_set.D()),
    lattice_(lattice),
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
PEPSt<ITensor>::PEPSt(const Lattice_Base &lattice, PEPSt_IndexSet_Base<ITensor::IndexT> &index_set);
template
PEPSt<IQTensor>::PEPSt(const Lattice_Base &lattice, PEPSt_IndexSet_Base<IQTensor::IndexT> &index_set);

template <class TensorT>
PEPSt<TensorT>::PEPSt(const Lattice_Base &lattice, PEPSt_IndexSet_Base<IndexT> &index_set, std::vector<TensorT> site_tensors_uc, std::vector<TensorT> bond_tensors_uc):
    d_(index_set.d()),
    D_(index_set.D()),
    lattice_(lattice),
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
    generate_bond_tensors(bond_tensors_uc);


}
template
PEPSt<ITensor>::PEPSt(const Lattice_Base &lattice, PEPSt_IndexSet_Base<ITensor::IndexT> &index_set, std::vector<ITensor> site_tensors_uc, std::vector<ITensor> bond_tensors_uc);
template
PEPSt<IQTensor>::PEPSt(const Lattice_Base &lattice, PEPSt_IndexSet_Base<IQTensor::IndexT> &index_set, std::vector<IQTensor> site_tensors_uc, std::vector<IQTensor> bond_tensors_uc);


template<class TensorT>
void PEPSt<TensorT>::generate_site_tensors(std::vector<TensorT> site_tensors_uc)
{
    for (int site_i=0; site_i<n_sites_total(); site_i++)
    {
        auto sublattice_i=lattice_.site_list_to_coord(site_i).at(2);
        tensor_assignment(site_tensors_[site_i],site_tensors_uc[sublattice_i]);

        //PrintDat(site_tensors_[site_i]);
    }
}
template
void PEPSt<ITensor>::generate_site_tensors(std::vector<ITensor> site_tensors_uc);
template
void PEPSt<IQTensor>::generate_site_tensors(std::vector<IQTensor> site_tensors_uc);


template<class TensorT>
void PEPSt<TensorT>::generate_bond_tensors(std::vector<TensorT> bond_tensors_uc)
{
    //for symmetric peps with half spin per site, we have 
    //B_{(x,y,i)}=\eta_{12}^{x}B_i, if B_i connect site with different y coordinate
    //So we need to double the unit cell in x direction if \mu_12==-1
    if (name_.find("half spin")!=std::string::npos && std::abs(mu_12+1)<EPSILON)
    {
        for (int bond_i=0; bond_i<n_bonds_total(); bond_i++)
        {
            auto sublattice_i=lattice_.bond_list_to_coord(bond_i)[2];

            tensor_assignment(bond_tensors_[bond_i],bond_tensors_uc[sublattice_i]);

            auto end_sites=lattice_.bond_end_sites(bond_i);
            std::array<Coordinate,2> end_sites_coord={lattice_.site_list_to_coord(end_sites[0]), lattice_.site_list_to_coord(end_sites[1])};

            if (std::abs(end_sites_coord[0][1]-end_sites_coord[1][1])%2==1)
            {
                auto eta_12=eta_from_mu(mu_12,dag(bond_tensors_[bond_i].indices()[0]));
                bond_tensors_[bond_i]*=eta_12;
                bond_tensors_[bond_i].noprime();
            }
        }

        return;
    }

    for (int bond_i=0; bond_i<n_bonds_total(); bond_i++)
    {
        auto sublattice_i=lattice_.bond_list_to_coord(bond_i).at(2);
        tensor_assignment(bond_tensors_[bond_i],bond_tensors_uc[sublattice_i]);
    }
}


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

        for (int endi=0; endi<2; endi++)
        {
            int endi_site=lattice_.bond_end_sites(bond_i,endi);

            //consider the boundary bonds case separetely
            if (endi_site<0)
            {
                IndexT endi_index=indexset_ptr_->virt_legs(n_sites_total()*n_bonds_to_one_site()+boundary_legs_no);
                endi_index.dag();
                tensor_indices.push_back(endi_index);
                boundary_legs_no++;
                continue;
            }
            
            std::vector<int> endi_site_neigh=lattice_.site_neighbour_bonds(endi_site);
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

    if (name_.find("cylinder")!=std::string::npos)
    {
        //left boundary legs are associated with sites, labeled as -1
        //right boundary legs are associated with bonds, labeled as -2
        std::vector<bool> boundary_tensor_created(lattice_.n_boundary_legs(),false);
        for (int boundary_i=0; boundary_i<lattice_.n_boundary_legs(); boundary_i++)
        {
            if (boundary_tensor_created[boundary_i]) continue;
            int site_i=lattice_.boundary_end_site(boundary_i);
            std::vector<int> neighbour_sites=lattice_.site_neighbour_sites(site_i);
            //begin_iter is the position start to seach the "invalid" neighbour sites
            //should be updated to find next "invalid" one
            auto begin_iter=neighbour_sites.begin();
            for (const auto &boundary_id : lattice_.site_neighbour_boundary(site_i))
            {
                auto neighbour_boundary_iter=std::find_if(begin_iter,neighbour_sites.end(), [](int i){ return (i<0); });
                int neighbour_boundary_no=neighbour_boundary_iter-begin_iter;

                //left boundary, which connects to site tensor
                if (*neighbour_boundary_iter==-1)
                {
                    auto boundary_leg=dag(indexset_ptr_->virt_legs(site_i*n_bonds_to_one_site()+neighbour_boundary_no));
                    boundary_tensors_[boundary_id]=TensorT(boundary_leg);
                }

                //right boundary, which connects to bond tensor
                if (*neighbour_boundary_iter==-2)
                {
                    int bond_i=lattice_.site_neighbour_bonds(site_i,neighbour_boundary_no);
                    int bond_boundary_i=(lattice_.bond_end_sites(bond_i,0)==-2 ? 0 : 1);
                    //if (lattice_.bond_end_sites(bond_i,0)==-2) bond_boundary_i=0;
                    //if (lattice_.bond_end_sites(bond_i,1)==-2) bond_boundary_i=1;
                    auto boundary_leg=dag(bond_tensors_[bond_i].indices()[bond_boundary_i]);
                    boundary_tensors_[boundary_id]=TensorT(boundary_leg);
                }

                boundary_tensor_created[boundary_id]=true;
                begin_iter=neighbour_boundary_iter;
                begin_iter++;
            }
        }

        //int boundary_no=0;
        //for (int sitei=0; sitei<n_sites_total(); sitei++)
        //{
        //    std::vector<int> neighbour_bonds=lattice_.site_neighbour_bonds(sitei);
        //    auto begin_iter=neighbour_bonds.begin();
        //    //sitei may have many boundary legs
        //    while (begin_iter!=neighbour_bonds.end())
        //    {
        //        auto found_iter=std::find(begin_iter,neighbour_bonds.end(),-1);
        //        if (found_iter==neighbour_bonds.end()) break;

        //        int neigh_bond_no=found_iter-begin_iter;
        //        boundary_tensors_[boundary_no]=TensorT(dag(indexset_ptr_->virt_legs(sitei*n_bonds_to_one_site()+neigh_bond_no)));

        //        //cout << boundary_tensors_[0][boundary_no];

        //        found_iter++;
        //        begin_iter=found_iter;
        //        boundary_no++;
        //    }
        //}


        //for (int bondi=0; bondi<n_bonds_total(); bondi++)
        //{
        //    for (int legi=0; legi<2; legi++)
        //    {
        //        if (lattice_.bond_end_sites(bondi,legi)==-2)
        //        {
        //            boundary_tensors_[boundary_no]=TensorT(dag(bond_tensors_[bondi].indices()[legi]));

        //            //cout << boundary_tensors_[1][boundary_no];

        //            boundary_no++;
        //        }
        //    }
        //}
    }
}
template
void PEPSt<ITensor>::new_boundary_tensors();
template 
void PEPSt<IQTensor>::new_boundary_tensors();


//before reading, we should init PEPSt use lattice
template <class TensorT>
void PEPSt<TensorT>::read(std::istream &s)
{
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


void randomize_spin_sym_square_peps(IQPEPS &spin_peps)
{
    //generate site tensors by random parameters
    //random number generator
    std::default_random_engine generator(std::time(0));
    std::uniform_real_distribution<double> distribution(-1.0,1.0);
    auto rand_param=std::bind(distribution,generator);

    //Init a symmetric site0, then generating all other sites
    Singlet_Tensor_Basis site0_tensor_basis(spin_peps.site_tensors(0).indices());
    std::vector<Complex> site0_tensor_params(site0_tensor_basis.dim(),0.);

    //cout << "Basis No. equals to " << site0_tensor_basis.dim() << endl << endl;
    //for (const auto &tensor : site0_tensor_basis)
    //{
    //    PrintDat(tensor);
    //}

    //site0_basis_visited==true means params of the basis has been assigned
    std::vector<bool> site0_basis_visited(site0_tensor_basis.dim(),false); 

    IQIndex site0_phys_indice;
    IndexSet<IQIndex> site0_virt_indices;
    std::vector<int> virt_dims;
    for (const auto &indice : spin_peps.site_tensors(0).indices())
    {
        if (indice.type()==Site) 
        {
            site0_phys_indice=indice;
            continue;
        }
        site0_virt_indices.addindex(indice);
        virt_dims.push_back(indice.m());
    }

    for (int i=0; i<site0_tensor_basis.dim(); i++)
    {
        if (site0_basis_visited[i]) continue;
        site0_basis_visited[i]=true;

        auto spin_list=site0_tensor_basis.spin_configs(i);
        auto deg_list=site0_tensor_basis.deg_configs(i);
        int fusion_channel=site0_tensor_basis.fusion_channel(i);
        site0_tensor_params[i]=rand_param();

        cout << "Spins: " << spin_list << endl
             << "Colors: " << deg_list << endl
             << "Params: " << site0_tensor_params[i] << endl << endl;
        //cout   << "Basis:" << site0_tensor_basis[i] << endl;


        //find the first nonzero elem for site0_tensor_basis[i] for later convient (compare the sign of different basis related by rotation symmetry)
        int val_num=0;
        double basis_elem;
        std::vector<int> val_list;

        //cout << virt_dims << endl;
        do
        {
            val_list=list_from_num(val_num,virt_dims);

            basis_elem=site0_tensor_basis[i](site0_phys_indice(1),site0_virt_indices[0](val_list[0]+1),site0_virt_indices[1](val_list[1]+1),site0_virt_indices[2](val_list[2]+1),site0_virt_indices[3](val_list[3]+1));
            //cout << val_list << endl;
            //cout << basis_elem << endl;
            val_num++;
        }
        while (std::abs(basis_elem)<EPSILON);

        //TODO: Debug the following
        //spin_oddness[i]==1 means spin_list[i] stores half-int spin
        std::vector<int> spin_oddness;
        for (const auto &S : spin_list)
            spin_oddness.push_back(S%2);

        //params may be imaginary when Theta_c4 is imaginary as well as the following condition holds
        if (spin_oddness[1]!=spin_oddness[3])
            site0_tensor_params[i]*=Theta_c4;

        //generate rotation symmetry related params
        for (int j=1; j<4; j++)
        {

            std::rotate(spin_list.begin(),spin_list.begin()+1,spin_list.end());
            std::rotate(deg_list.begin(),deg_list.begin()+1,deg_list.end());
            std::rotate(val_list.begin(),val_list.begin()+1,val_list.end());
            cout << "Spins: " << spin_list << endl
                << "Colors: " << deg_list << endl;

            auto rotate_basis_no=site0_tensor_basis.spin_deg_list_to_basis_no(spin_list,deg_list,fusion_channel);
            //mark the rotated basis as visited.
            if (site0_basis_visited[rotate_basis_no]) continue;
            site0_basis_visited[rotate_basis_no]=true;

            //when generating rotated params, be cautious that the singlet basis may differ by -1, which we saved in rotate_basis_ratio
            auto rotate_basis_ratio=site0_tensor_basis[rotate_basis_no](site0_phys_indice(1),site0_virt_indices[0](val_list[0]+1),site0_virt_indices[1](val_list[1]+1),site0_virt_indices[2](val_list[2]+1),site0_virt_indices[3](val_list[3]+1))/basis_elem;

            site0_tensor_params[rotate_basis_no]=site0_tensor_params[i]*std::pow(chi_c4*Theta_c4,j)/rotate_basis_ratio;

            cout << "Basis Ratios: " << rotate_basis_ratio << endl
                << "Params: " << site0_tensor_params[i] << endl << endl;
        }//generate rotated params
    }//generate all params

    //generate all translational related site tensors
    spin_peps.generate_site_tensors({singlet_tensor_from_basis_params(site0_tensor_basis,site0_tensor_params)});


    //generate bond tensors
    //init bond0, and generated all other bonds
    Singlet_Tensor_Basis bond0_tensor_basis(spin_peps.bond_tensors(0).indices());
    std::vector<Complex> bond0_tensor_params(bond0_tensor_basis.dim(),0.);

    for (int i=0; i<bond0_tensor_basis.dim(); i++)
    {
        auto spin_list=bond0_tensor_basis.spin_configs(i);
        auto deg_list=bond0_tensor_basis.deg_configs(i);
        
        //we choose gauge such that bond tensor in extra deg space is either diagonal (symmetric) or ~i\sigma^y\otimes\mathrm{I}_{n/2}
        int deg_dim=bond0_tensor_basis.spin_degs(0)[spin_list[0]];
        if (deg_list[0]!=deg_list[1]) //not diagonal
        {
            //not the other case
            if (deg_dim%2==1 || std::abs(deg_list[0]-deg_list[1])!=deg_dim/2) continue;
        }

        //singlet form by two integer spins
        //chi_c4=1: real sym
        //chi_c4=-1: real antisym
        if (spin_list[0]%2==0)
        {
            if (std::abs(chi_c4-1)<EPSILON)
            {
                if (deg_list[0]!=deg_list[1]) continue;

                //TODO: consider the case with -1
                bond0_tensor_params[i]=1.;
            }
            if (std::abs(chi_c4+1)<EPSILON)
            {
                assert(deg_dim%2==0);

                if (deg_list[0]-deg_list[1]==deg_dim/2) 
                    bond0_tensor_params[i]=-1.;
                if (deg_list[1]-deg_list[0]==deg_dim/2)
                    bond0_tensor_params[i]=1.;
            }
            continue;
        }

        //singlet form by two half integer spins
        //mu_t1T=1: real
        //mu_t1T=-1: imag
        //mu_t2c4.chi_c4=-1: sym
        //mu_t2c4.chi_c4=+1: antisym
        if (spin_list[0]%2==1)
        {
            if (std::abs(mu_t2c4*chi_c4+1)<EPSILON)
            {
                if (deg_list[0]!=deg_list[1]) continue;

                //TODO: consider the case with -1
                bond0_tensor_params[i]=1.;
            }
            if (std::abs(mu_t2c4*chi_c4-1)<EPSILON)
            {
                assert(deg_dim%2==0);

                if (deg_list[0]-deg_list[1]==deg_dim/2) 
                    bond0_tensor_params[i]=-1.;
                if (deg_list[1]-deg_list[0]==deg_dim/2)
                    bond0_tensor_params[i]=1.;
            }
            if (std::abs(mu_t1T+1)<EPSILON)
            {
                bond0_tensor_params[i]*=Complex_i;
            }
        }
    }
    spin_peps.bond_tensors(0)=singlet_tensor_from_basis_params(bond0_tensor_basis,bond0_tensor_params);

    //generate the horizontal bond B_1 by vertical bond B_0
    tensor_assignment(spin_peps.bond_tensors(1),spin_peps.bond_tensors(0));
    spin_peps.bond_tensors(1)*=chi_c4;

    //generate all translational related bond tensors
    spin_peps.generate_bond_tensors({spin_peps.bond_tensors(0),spin_peps.bond_tensors(1)});
}
