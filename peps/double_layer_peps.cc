
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
            if (left_leg_of_combiners(combiner,0).primeLevel()==0) lower_ind=dag(left_leg_of_combiners(combiner,0));
            if (left_leg_of_combiners(combiner,1).primeLevel()==0) lower_ind=dag(left_leg_of_combiners(combiner,1));

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
typename TensorT::CombinerT Double_Layer_PEPSt<TensorT>::decombine_double_virt_indice(std::array<int,2> sites_no, IndexT &lower_ind)
{
    for (const auto &combineri : virt_leg_combiners_[sites_no[0]])
    {
        for (const auto &combinerj : virt_leg_combiners_[sites_no[1]])
        {
            if (combineri.right()==combinerj.right())
            {
                if (left_leg_of_combiners(combineri,0).primeLevel()==0) lower_ind=dag(left_leg_of_combiners(combineri,0));
                if (left_leg_of_combiners(combineri,1).primeLevel()==0) lower_ind=dag(left_leg_of_combiners(combineri,1));

                return combineri;
            }
        }
    }

    cout << "Fail to find lower tensor on site " << sites_no[0] << endl;
    exit(EXIT_FAILURE);
}
template
typename ITensor::CombinerT Double_Layer_PEPSt<ITensor>::decombine_double_virt_indice(std::array<int,2> sites_no, ITensor::IndexT &lower_ind);
template
typename IQTensor::CombinerT Double_Layer_PEPSt<IQTensor>::decombine_double_virt_indice(std::array<int,2> sites_no, IQTensor::IndexT &lower_ind);



template <class TensorT>
void Double_Layer_PEPSt<TensorT>::obtain_peps_sandwich_single_site_operators(std::vector<TensorT> direct_prod_operators, const std::vector<int> &acting_sites_list, std::vector<TensorT>& sandwiched_tensors)
{
    if (sandwiched_tensors.empty()) sandwiched_tensors=double_layer_tensors_;

    for (int i=0; i<direct_prod_operators.size(); i++)
    {
        //get operator act on site_no
        IndexT phys_leg;
        int site_no=acting_sites_list[i];
        for (const auto &leg : single_layer_tensors_[site_no].indices())
        {
            if (leg.type()==Site)
            {
                phys_leg=leg;
                //Print(phys_leg);
                break;
            }
        }
        for (int legi=0; legi<2; legi++)
        {
            auto old_leg=direct_prod_operators[i].indices()[legi];
            if (old_leg.primeLevel()==0) direct_prod_operators[i].replaceIndex(old_leg,dag(phys_leg));
            if (old_leg.primeLevel()==1) direct_prod_operators[i].replaceIndex(old_leg,prime(phys_leg));

            //Print(old_leg);
        }

        //acting operator and replace the corresponding sandwiched_tensors
        TensorT lower_tensor=single_layer_tensors_[site_no],
                upper_tensor=prime(dag(lower_tensor));
        sandwiched_tensors[site_no]=lower_tensor*direct_prod_operators[i]*upper_tensor;
        for (const auto &combiner : virt_leg_combiners_[site_no])
        {
            sandwiched_tensors[site_no]=sandwiched_tensors[site_no]*combiner;
        }
        //PrintDat(sandwiched_tensors[site_no]);
    }
}
template
void Double_Layer_PEPSt<ITensor>::obtain_peps_sandwich_single_site_operators(std::vector<ITensor> direct_prod_operators, const std::vector<int> &acting_sites_list, std::vector<ITensor>& sandwiched_tensors);
template
void Double_Layer_PEPSt<IQTensor>::obtain_peps_sandwich_single_site_operators(std::vector<IQTensor> direct_prod_operators, const std::vector<int> &acting_sites_list, std::vector<IQTensor>& sandwiched_tensors);

template <class TensorT>
void Double_Layer_PEPSt<TensorT>::obtain_peps_sandwich_tensor_prod_operators(std::vector<TPOt<TensorT>> tensor_prod_operators, const std::vector<std::vector<int>> &acting_sites_list, std::vector<TensorT> &sandwiched_tensors)
{
    if (sandwiched_tensors.empty()) sandwiched_tensors=double_layer_tensors_;

    for (int opi=0; opi<tensor_prod_operators.size(); opi++)
    {
        std::vector<int> acting_sites=acting_sites_list[opi];
        for (int j=0; j<acting_sites.size(); j++)
        {
            IndexT phys_leg=findtype(single_layer_tensors_[acting_sites[j]],Site);

            auto old_leg=tensor_prod_operators[opi].phys_legs(j);
            tensor_prod_operators[opi].site_tensors(j).replaceIndex(old_leg,dag(phys_leg));
            tensor_prod_operators[opi].site_tensors(j).replaceIndex(dag(prime(old_leg)),prime(phys_leg));

            PrintDat(tensor_prod_operators[opi].site_tensors(j));

            //sandwich tensor_prod_operators
            TensorT lower_tensor=single_layer_tensors_[acting_sites[j]],
                    upper_tensor=prime(dag(lower_tensor));

            //Print(double_layer_tensor_from_lower_upper_tensors(lower_tensor,upper_tensor,virt_leg_combiners_[acting_sites[j]]));

            sandwiched_tensors[acting_sites[j]]=double_layer_tensor_from_lower_upper_tensors(lower_tensor,upper_tensor,virt_leg_combiners_[acting_sites[j]])*tensor_prod_operators[opi].site_tensors(j);
            //Print(sandwiched_tensors[acting_sites[j]]);
        }
    }

}
template
void Double_Layer_PEPSt<ITensor>::obtain_peps_sandwich_tensor_prod_operators(std::vector<TPOt<ITensor>> tensor_prod_operators, const std::vector<std::vector<int>> &acting_sites_list, std::vector<ITensor> &sandwiched_tensors);
template
void Double_Layer_PEPSt<IQTensor>::obtain_peps_sandwich_tensor_prod_operators(std::vector<TPOt<IQTensor>> tensor_prod_operators, const std::vector<std::vector<int>> &acting_sites_list, std::vector<IQTensor> &sandwiched_tensors);


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
            //For bulk bond, we only multiply those start from the site to avoid double counting. Namely bond_neighbour_sites(bond_no,0)==sitei
            //For boundary bond, we will always absorb to the bond
            if (this->lattice().bond_neighbour_sites(bond_no,0)==sitei ||
                this->lattice().bond_neighbour_sites(bond_no,0)<0)
            {
                site_tensor*=single_layer_peps_.bond_tensors(bond_no);
            }
        }
        //multiply neighbouring boundary tensors
        for (const auto &boundary_no : this->lattice().site_neighbour_boundaries(sitei))
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

            //if the leg has already been combined, add the combiner to virt_leg_combiners_[sitei]. We also delete the legs and combiners in combined_indices(_combiners) since each leg only appear twice
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
TensorT Double_Layer_PEPSt<TensorT>::double_layer_tensor_from_lower_upper_tensors(TensorT lower_tensor, TensorT upper_tensor, const std::vector<CombinerT> &pair_combiners)
{
    //combine virtual legs of lower(upper) tensor
    std::vector<CombinerT> leg_combiners;
    auto indset=lower_tensor.indices();
    for (const auto &leg : indset)
    {
        if (leg.type()==Site) continue;
        if (leg_combiners.empty()) 
        {
            leg_combiners.push_back(CombinerT(leg));
        }
        else
        {
            leg_combiners.push_back(CombinerT(leg,leg_combiners.back().right()));
        }

        //Print(lower_tensor);
        //Print(upper_tensor);
        //Print(leg_combiners.back());

        lower_tensor=lower_tensor*leg_combiners.back();
        upper_tensor=upper_tensor*prime(dag(leg_combiners.back()));
    }

    //decombine vitual legs one by one and recombine pairs of them
    auto combined_tensor=lower_tensor*upper_tensor;
    while (leg_combiners.empty()==false)
    {
        auto combiner=leg_combiners.back();
        combined_tensor=combined_tensor*dag(combiner)*prime(combiner);
        for (const auto pair_combiner : pair_combiners)
        {
            if (hasindex(combined_tensor.indices(),left_leg_of_combiners(pair_combiner,0)))
            {
                combined_tensor=combined_tensor*pair_combiner;
                break;
            }
        }
        leg_combiners.pop_back();
    }
    return combined_tensor;
}
template 
ITensor Double_Layer_PEPSt<ITensor>::double_layer_tensor_from_lower_upper_tensors(ITensor lower_tensor, ITensor upper_tensor, const std::vector<ITensor::CombinerT> &pair_combiners);
template 
IQTensor Double_Layer_PEPSt<IQTensor>::double_layer_tensor_from_lower_upper_tensors(IQTensor lower_tensor, IQTensor upper_tensor, const std::vector<IQTensor::CombinerT> &pair_combiners);



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


