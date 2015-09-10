
#include "peps.h"

//
//class PEPSt_Torus
//

//
//Constructors
//
template <class TensorT>
PEPSt_Torus<TensorT>::PEPSt_Torus(const Lattice_Torus_Base &lattice, const PEPSt_IndexSet_Base<IndexT> &index_set):
    d_(index_set.d()),
    D_(index_set.D()),
    lattice_(lattice),
    index_set_(index_set),
    site_tensors_(lattice.n_sites_total()),
    bond_tensors_(lattice.n_bonds_total())
{
    new_site_tensors();
    new_bond_tensors();
    //random_site_tensors();
}
template
PEPSt_Torus<ITensor>::PEPSt_Torus(const Lattice_Torus_Base &lattice, const PEPSt_IndexSet_Base<ITensor::IndexT> &index_set);
template
PEPSt_Torus<IQTensor>::PEPSt_Torus(const Lattice_Torus_Base &lattice, const PEPSt_IndexSet_Base<IQTensor::IndexT> &index_set);

template <class TensorT>
PEPSt_Torus<TensorT>::PEPSt_Torus(const Lattice_Torus_Base &lattice, const PEPSt_IndexSet_Base<IndexT> &index_set, std::vector<TensorT> &site_tensors_uc, std::vector<TensorT> &bond_tensors_uc):
    d_(index_set.d()),
    D_(index_set.D()),
    lattice_(lattice),
    index_set_(index_set),
    site_tensors_(lattice.n_sites_total()),
    bond_tensors_(lattice.n_bonds_total())
{
    new_site_tensors();
    new_bond_tensors();

    for (int site_i=0; site_i<n_sites_total(); site_i++)
    {
        auto sublattice_i=lattice.site_list_to_coord(site_i).at(2);
        tensor_assignment(site_tensors_[site_i],site_tensors_uc[sublattice_i]);
    }

    for (int bond_i=0; bond_i<n_bonds_total(); bond_i++)
    {
        auto sublattice_i=lattice.bond_list_to_coord(bond_i).at(2);
        tensor_assignment(bond_tensors_[bond_i],bond_tensors_uc[sublattice_i]);
    }
}
template
PEPSt_Torus<ITensor>::PEPSt_Torus(const Lattice_Torus_Base &lattice, const PEPSt_IndexSet_Base<ITensor::IndexT> &index_set, std::vector<ITensor> &site_tensors_uc, std::vector<ITensor> &bond_tensors_uc);
template
PEPSt_Torus<IQTensor>::PEPSt_Torus(const Lattice_Torus_Base &lattice, const PEPSt_IndexSet_Base<IQTensor::IndexT> &index_set, std::vector<IQTensor> &site_tensors_uc, std::vector<IQTensor> &bond_tensors_uc);



template<class TensorT>
void PEPSt_Torus<TensorT>::tensor_assignment(TensorT &TA, TensorT &TB)
{
    TensorT tensor_tmp=TB;

    for (int leg_i=0; leg_i<TA.r(); leg_i++)
    {
        tensor_tmp.replaceIndex(tensor_tmp.indices()[leg_i],TA.indices()[leg_i]);
    }
    TA=tensor_tmp;
    return;
}
template
void PEPSt_Torus<ITensor>::tensor_assignment(ITensor &TA, ITensor &TB);
template
void PEPSt_Torus<IQTensor>::tensor_assignment(IQTensor &TA, IQTensor &TB);

//According to the library, ITensor is constructed by IndexSet<Index>(indexset), while IQTensor is constructed by indexset directly
template <>
void PEPSt_Torus<ITensor>::construct_tensor(ITensor &tensor, std::vector<ITensor::IndexT> &indexset)
{
    tensor=ITensor(IndexSet<Index>(indexset));
    return;
}
template<>
void PEPSt_Torus<IQTensor>::construct_tensor(IQTensor &tensor, std::vector<IQTensor::IndexT> &indexset)
{
    tensor=IQTensor(indexset);
    return;
}


template <class TensorT>
void PEPSt_Torus<TensorT>::new_site_tensors()
{
    int site_i=0;
    for(auto &tensor : site_tensors_)
    {
        //tensor_indices stores both physical legs and virtual legs for site_tensors_[site_i]
        std::vector<IndexT> tensor_indices;
        tensor_indices.push_back(index_set_.phys_legs(site_i));

        for (int j=0; j<n_bonds_to_one_site(); j++)
        {
            IndexT neigh_virt_ind=index_set_.virt_legs(j+site_i*n_bonds_to_one_site());
            tensor_indices.push_back(neigh_virt_ind);
        }

        construct_tensor(tensor,tensor_indices);

        site_i++;
    }

    return;
}
template
void PEPSt_Torus<ITensor>::new_site_tensors();
template
void PEPSt_Torus<IQTensor>::new_site_tensors();


template <class TensorT>
void PEPSt_Torus<TensorT>::new_bond_tensors()
{
    int bond_i=0;

    for (auto &tensor : bond_tensors_)
    {
        std::vector<IndexT> tensor_indices;

        for (int endi=0; endi<2; endi++)
        {
            auto endi_site=lattice_.bond_end_sites(bond_i,endi);
            auto endi_site_neigh=lattice_.site_neighbour_bonds(endi_site);
            int legi=std::find(endi_site_neigh.begin(),endi_site_neigh.end(),bond_i)-endi_site_neigh.begin();
            assert(legi<n_bonds_to_one_site());

            IndexT endi_index=index_set_.virt_legs(endi_site*n_bonds_to_one_site()+legi);
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
void PEPSt_Torus<ITensor>::new_bond_tensors();
template
void PEPSt_Torus<IQTensor>::new_bond_tensors();

//template <class TensorT>
//void PEPSt_Torus<TensorT>::random_site_tensors()
//{
//    new_site_tensors();
//    for(auto& tensor : site_tensors_)
//    {
//        tensor.randomize();
//    }
//}
//template
//void PEPSt_Torus<ITensor>::random_site_tensors();
//template
//void PEPSt_Torus<IQTensor>::random_site_tensors();
