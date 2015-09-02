
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
    tensors_div_(lattice.n_sites_total())
{
    new_site_tensors();
    //random_site_tensors();
}
template
PEPSt_Torus<ITensor>::PEPSt_Torus(const Lattice_Torus_Base &lattice, const PEPSt_IndexSet_Base<ITensor::IndexT> &index_set);
template
PEPSt_Torus<IQTensor>::PEPSt_Torus(const Lattice_Torus_Base &lattice, const PEPSt_IndexSet_Base<IQTensor::IndexT> &index_set);

template <class TensorT>
PEPSt_Torus<TensorT>::PEPSt_Torus(const Lattice_Torus_Base &lattice, const PEPSt_IndexSet_Base<IndexT> &index_set, std::vector<TensorT> &site_tensors_uc):
    d_(index_set.d()),
    D_(index_set.D()),
    lattice_(lattice),
    index_set_(index_set),
    site_tensors_(lattice.n_sites_total()),
    tensors_div_(lattice.n_sites_total())
{
    new_site_tensors();

    for (int site_i=0; site_i<n_sites_total(); site_i++)
    {
        auto sublattice_i=lattice.site_list_to_coord(site_i).at(2);
        tensor_assignment(site_tensors_[site_i],site_tensors_uc[sublattice_i]);
    }
}
template
PEPSt_Torus<ITensor>::PEPSt_Torus(const Lattice_Torus_Base &lattice, const PEPSt_IndexSet_Base<ITensor::IndexT> &index_set, std::vector<ITensor> &site_tensors_uc);
template
PEPSt_Torus<IQTensor>::PEPSt_Torus(const Lattice_Torus_Base &lattice, const PEPSt_IndexSet_Base<IQTensor::IndexT> &index_set, std::vector<IQTensor> &site_tensors_uc);



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
    for (int site_i=0; site_i<n_sites_total(); site_i++)
    {
        //sitei_indexset stores both physical legs and virtual legs for site_tensors_[site_i]
        std::vector<IndexT> sitei_indexset;
        sitei_indexset.push_back(index_set_.phys_legs(site_i));

        for (int j=0; j<n_bonds_to_one_site(); j++)
        {
            int neigh_bond=lattice_.site_neighbour_bonds(site_i,j);
            IndexT neigh_virt_ind=index_set_.virt_legs(neigh_bond);
            
            //set the direction of neigh_virt_ind according to the bond direction
            if (lattice_.bond_end_sites(neigh_bond,0)==site_i)
            {
                sitei_indexset.push_back(neigh_virt_ind);
            }
            else
            {
                assert(lattice_.bond_end_sites(neigh_bond,1)==site_i);
                neigh_virt_ind.dag();
                sitei_indexset.push_back(neigh_virt_ind);
            }

        }
        construct_tensor(site_tensors_[site_i],sitei_indexset);
    }
    return;
}
template
void PEPSt_Torus<ITensor>::new_site_tensors();
template
void PEPSt_Torus<IQTensor>::new_site_tensors();


template <class TensorT>
void PEPSt_Torus<TensorT>::random_site_tensors()
{
    new_site_tensors();
    for(auto& tensor : site_tensors_)
    {
        tensor.randomize();
    }
}
template
void PEPSt_Torus<ITensor>::random_site_tensors();
template
void PEPSt_Torus<IQTensor>::random_site_tensors();
