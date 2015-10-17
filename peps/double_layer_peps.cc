
#include "double_layer_peps.h"

template <class TensorT>
Double_Layer_PEPSt<TensorT>::Double_Layer_PEPSt(const Lattice_Base &lattice):
    single_layer_peps_(lattice),
    single_layer_tensors_(lattice.n_sites_total()),
    double_layer_tensors_(lattice.n_sites_total()),
    virt_leg_combiners_(lattice.n_sites_total())
{}
template
Double_Layer_PEPSt<ITensor>::Double_Layer_PEPSt(const Lattice_Base &lattice);
template
Double_Layer_PEPSt<IQTensor>::Double_Layer_PEPSt(const Lattice_Base &lattice);

template <class TensorT>
Double_Layer_PEPSt<TensorT>::Double_Layer_PEPSt(const PEPSt<TensorT> &peps):
    single_layer_peps_(peps),
    single_layer_tensors_(peps.n_sites_total()),
    double_layer_tensors_(peps.n_sites_total()),
    virt_leg_combiners_(peps.n_sites_total())
{
    obtain_single_layer_tensors();
    
    obtain_layered_tensors_with_combined_legs();
}
template
Double_Layer_PEPSt<ITensor>::Double_Layer_PEPSt(const PEPSt<ITensor> &peps);
template
Double_Layer_PEPSt<IQTensor>::Double_Layer_PEPSt(const PEPSt<IQTensor> &peps);


template <class TensorT>
typename TensorT::CombinerT Double_Layer_PEPSt<TensorT>::decombine_double_virt_indice(int sitei, const IndexT &double_virt_ind, IndexT &lower_ind)
{
    assert(hasindex(double_layer_tensors_[sitei],double_virt_ind));

    for (const auto &combiner : virt_leg_combiners_[sitei])
    {
        if (combiner.right()==double_virt_ind)
        {
            if (left_leg_of_combiners(combiner,0).primeLevel()==0) lower_ind=left_leg_of_combiners(combiner,0);
            if (left_leg_of_combiners(combiner,1).primeLevel()==0) lower_ind=left_leg_of_combiners(combiner,1);

            return combiner;
        }
    }

    cout << "Fail to find lower tensor on site " << sitei << endl;
    exit(EXIT_FAILURE);
}
template
typename ITensor::CombinerT Double_Layer_PEPSt<ITensor>::decombine_double_virt_indice(int sitei, const ITensor::IndexT &double_virt_ind, ITensor::IndexT &lower_ind);
template
typename IQTensor::CombinerT Double_Layer_PEPSt<IQTensor>::decombine_double_virt_indice(int sitei, const IQTensor::IndexT &double_virt_ind, IQTensor::IndexT &lower_ind);


template <class TensorT>
void Double_Layer_PEPSt<TensorT>::obtain_peps_sandwich_single_site_operators(std::vector<TensorT> single_site_operators, const std::vector<int> &acting_sites, std::vector<TensorT>& sandwiched_tensors)
{
    //auto sandwiched_tensors=double_layer_tensors_;

    for (auto site_no : acting_sites)
    {
        //get operator act on site_no
        IndexT phys_leg;
        for (const auto &leg : single_layer_tensors_[site_no].indices())
        {
            if (leg.type()==Site)
            {
                phys_leg=leg;
                break;
            }
        }
        for (int legi=0; legi<2; legi++)
        {
            auto old_leg=single_site_operators[site_no].indices()[legi];
            if (old_leg.primeLevel()==0) single_site_operators[site_no].replaceIndex(old_leg,dag(phys_leg));
            if (old_leg.primeLevel()==1) single_site_operators[site_no].replaceIndex(old_leg,prime(phys_leg));
        }

        //acting operator and replace the corresponding sandwiched_tensors
        TensorT lower_tensor=single_layer_tensors_[site_no],
                upper_tensor=prime(dag(lower_tensor));
        sandwiched_tensors[site_no]=lower_tensor*single_site_operators[site_no]*upper_tensor;
        for (const auto &combiner : virt_leg_combiners_[site_no])
        {
            sandwiched_tensors[site_no]=sandwiched_tensors[site_no]*combiner;
        }
    }
}
template
void Double_Layer_PEPSt<ITensor>::obtain_peps_sandwich_single_site_operators(std::vector<ITensor> single_site_operators, const std::vector<int> &acting_sites, std::vector<ITensor>& sandwiched_tensors);
template
void Double_Layer_PEPSt<IQTensor>::obtain_peps_sandwich_single_site_operators(std::vector<IQTensor> single_site_operators, const std::vector<int> &acting_sites, std::vector<IQTensor>& sandwiched_tensors);


template <class TensorT>
void Double_Layer_PEPSt<TensorT>::obtain_single_layer_tensors()
{
    for (int sitei=0; sitei<single_layer_peps_.n_sites_total(); sitei++)
    {
        auto site_tensor=single_layer_peps_.site_tensors(sitei);

        //multiply neighbouring bond tensors 
        for (int neighi=0; neighi<this->lattice().n_bonds_to_one_site(); neighi++)
        {
            int bond_no=this->lattice().site_neighbour_bonds(sitei,neighi);
            if (bond_no<0) continue;
            //For bulk bond, we only multiply those start from the site to avoid double counting. Namely bond_end_sites(bond_no,0)==sitei
            //For boundary bond, we will always absorb to the bond
            if (this->lattice().bond_end_sites(bond_no,0)==sitei ||
                this->lattice().bond_end_sites(bond_no,0)<0)
            {
                site_tensor*=single_layer_peps_.bond_tensors(bond_no);
            }
        }
        //multiply neighbouring boundary tensors
        for (const auto &boundary_no : this->lattice().site_neighbour_boundary(sitei))
        {
            site_tensor*=single_layer_peps_.boundary_tensors(boundary_no);
        }
        clean(site_tensor);

        single_layer_tensors_[sitei]=site_tensor;

        //cout << "Single layer tensor " << sitei << endl << site_tensor << endl;
    }
}
template
void Double_Layer_PEPSt<ITensor>::obtain_single_layer_tensors();
template
void Double_Layer_PEPSt<IQTensor>::obtain_single_layer_tensors();


template <class TensorT>
void Double_Layer_PEPSt<TensorT>::obtain_layered_tensors_with_combined_legs()
{
    //create combiners for virt_leg and prime(dag(virt_leg)), where virt_leg connects sitei and neighbour_sites[neighbour_i]
    //every combiner should appear twice, so we stores combined indiceswhich only appear once, and delete it if it already appear twice
    std::vector<IndexT> combined_indices;
    std::vector<CombinerT> indices_combiners;

    int combine_i=0;
    for (int sitei=0; sitei<this->lattice().n_sites_total(); sitei++)
    {
        for (const auto &virt_leg : single_layer_tensors_[sitei].indices())
        {
            //for physical leg, we do not combine
            if (virt_leg.type()==Site) continue;

            auto virt_leg_iter=std::find(combined_indices.begin(),combined_indices.end(),virt_leg);

            //if the leg has already been combined, add the combiner to virt_leg_combiners_[sitei], and delete the legs and combiners in combined_indices(_combiners) since each leg only appear twice
            if (virt_leg_iter!=combined_indices.end())
            {
                std::swap(*virt_leg_iter,*(combined_indices.end()-1));
                std::swap(indices_combiners[virt_leg_iter-combined_indices.begin()],*(indices_combiners.end()-1));

                virt_leg_combiners_[sitei].push_back(dag(*(indices_combiners.end()-1)));

                combined_indices.pop_back();
                indices_combiners.pop_back();

                continue;
            }

            //for the case where the virtual leg appears the first time
            auto leg_combiner=CombinerT(virt_leg,dag(virt_leg).prime());
            leg_combiner.init(nameint("leg_combiner ",combine_i),Link,virt_leg.dir());
            virt_leg_combiners_[sitei].push_back(leg_combiner);
            combine_i++;

            combined_indices.push_back(virt_leg);
            indices_combiners.push_back(leg_combiner);
        }
    }

    //for (int sitei=0; sitei<this->lattice().n_sites_total(); sitei++)
    //{
    //    cout << "Combiners for site " << sitei << endl;
    //    for (const auto &leg_combiner : virt_leg_combiners_[sitei])
    //    {
    //        cout << leg_combiner;
    //    }
    //    cout << endl;
    //}


    //create layered tensors with virtual legs combined
    std::vector<TensorT> upper_tensors(single_layer_tensors_);
    for (auto &tensor : upper_tensors)
    {
        tensor.dag();
        tensor.prime(Link);
    }

    for (int sitei=0; sitei<this->lattice().n_sites_total(); sitei++)
    {
        double_layer_tensors_[sitei]=single_layer_tensors_[sitei]*upper_tensors[sitei];
        for (const auto &combiner : virt_leg_combiners_[sitei])
        {
            double_layer_tensors_[sitei]=double_layer_tensors_[sitei]*combiner;
        }
        clean(double_layer_tensors_[sitei]);
    }

}
template
void Double_Layer_PEPSt<ITensor>::obtain_layered_tensors_with_combined_legs();
template
void Double_Layer_PEPSt<IQTensor>::obtain_layered_tensors_with_combined_legs();


template <class TensorT>
void Double_Layer_PEPSt<TensorT>::read(std::istream &s)
{
    single_layer_peps_.read(s);
    
    for (auto &tensor : single_layer_tensors_) tensor.read(s);

    //reconstruct double_layer_tensors_ and virt_leg_combiners_
    obtain_layered_tensors_with_combined_legs();
}
template
void Double_Layer_PEPSt<ITensor>::read(std::istream &s);
template
void Double_Layer_PEPSt<IQTensor>::read(std::istream &s);

template <class TensorT>
void Double_Layer_PEPSt<TensorT>::write(std::ostream &s) const
{
    single_layer_peps_.write(s);

    for (const auto &tensor : single_layer_tensors_) tensor.write(s);
}
template
void Double_Layer_PEPSt<ITensor>::write(std::ostream &s) const;
template
void Double_Layer_PEPSt<IQTensor>::write(std::ostream &s) const;


