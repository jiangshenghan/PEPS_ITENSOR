
#include "double_layer_peps.h"

template <class TensorT>
Double_Layer_PEPSt<TensorT>::Double_Layer_PEPSt(const PEPSt<TensorT> &peps):
    DD_(peps.D()*peps.D()),
    lattice_(peps.lattice()),
    layered_site_tensors_(peps.n_sites_total()),
    virt_leg_combiners_(peps.n_sites_total())
{
    std::vector<TensorT> combined_site_tensors;
    obtain_combined_site_tensors(peps,combined_site_tensors);
    obtain_layered_tensors_with_combined_legs(combined_site_tensors);
}
template
Double_Layer_PEPSt<ITensor>::Double_Layer_PEPSt(const PEPSt<ITensor> &peps);
template
Double_Layer_PEPSt<IQTensor>::Double_Layer_PEPSt(const PEPSt<IQTensor> &peps);


//template <class TensorT>
//Double_Layer_PEPSt<TensorT>::Double_Layer_PEPSt(const Lattice_Base &square_lattice_cylinder, const TensorT &square_site_tensor, const TensorT &boundary_tensor):
//    lattice_(square_lattice_cylinder)
//{
//    for (const auto &leg : square_site_tensor.indices())
//    {
//        if (leg.type()==Link)
//        {
//            DD_=leg.m()*leg.m();
//            break;
//        }
//    }
//
//    std::vector<IndexT> virt_legs(lattice_.n_sites_total()-2*lattice_.n_uc()[1]);
//
//    TensorT layered_bulk_tensor, layered_left_boundary_tensor, layered_right_boundary_tensor;
//
//    for (const auto &leg : square_site_tensor.indices())
//    {
//    }
//
//    for (int sitei=0; sitei<lattice_.n_sites_total(); sitei++)
//}


template <class TensorT>
void Double_Layer_PEPSt<TensorT>::obtain_combined_site_tensors(const PEPSt<TensorT> &peps, std::vector<TensorT> &combined_site_tensors)
{
    for (int sitei=0; sitei<peps.n_sites_total(); sitei++)
    {
        auto site_tensor=peps.site_tensors(sitei);

        //multiply neighbouring bond tensors 
        for (int neighi=0; neighi<lattice_.n_bonds_to_one_site(); neighi++)
        {
            int bond_no=lattice_.site_neighbour_bonds(sitei,neighi);
            if (bond_no<0) continue;
            //For bulk bond, we only multiply those start from the site to avoid double counting. Namely bond_end_sites(bond_no,0)==sitei
            //For boundary bond, we will always absorb to the bond
            if (lattice_.bond_end_sites(bond_no,0)==sitei ||
                lattice_.bond_end_sites(bond_no,0)<0)
            {
                site_tensor*=peps.bond_tensors(bond_no);
            }
        }
        //multiply neighbouring boundary tensors
        for (const auto &boundary_no : lattice_.site_neighbour_boundary(sitei))
        {
            site_tensor*=peps.boundary_tensors(boundary_no);
        }
        clean(site_tensor);

        combined_site_tensors.push_back(site_tensor);

        //cout << "Combined tensor " << sitei << endl << site_tensor << endl;
    }
}
template
void Double_Layer_PEPSt<ITensor>::obtain_combined_site_tensors(const PEPSt<ITensor> &peps, std::vector<ITensor> &combined_site_tensors);
template
void Double_Layer_PEPSt<IQTensor>::obtain_combined_site_tensors(const PEPSt<IQTensor> &peps, std::vector<IQTensor> &combined_site_tensors);


template <class TensorT>
void Double_Layer_PEPSt<TensorT>::obtain_layered_tensors_with_combined_legs(const std::vector<TensorT> &lower_tensors)
{
    //create combiners for virt_leg and prime(dag(virt_leg)), where virt_leg connects sitei and neighbour_sites[neighbour_i]
    //every combiner should appear twice, so we stores combined indiceswhich only appear once, and delete it if it already appear twice
    std::vector<IndexT> combined_indices;
    std::vector<CombinerT> indices_combiners;

    int combine_i=0;
    for (int sitei=0; sitei<lattice_.n_sites_total(); sitei++)
    {
        for (const auto &virt_leg : lower_tensors[sitei].indices())
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

    //for (int sitei=0; sitei<lattice_.n_sites_total(); sitei++)
    //{
    //    cout << "Combiners for site " << sitei << endl;
    //    for (const auto &leg_combiner : virt_leg_combiners_[sitei])
    //    {
    //        cout << leg_combiner;
    //    }
    //    cout << endl;
    //}


    //create layered tensors with virtual legs combined
    std::vector<TensorT> upper_tensors(lower_tensors);
    for (auto &tensor : upper_tensors)
    {
        tensor.dag();
        tensor.prime(Link);
    }

    for (int sitei=0; sitei<lattice_.n_sites_total(); sitei++)
    {
        layered_site_tensors_[sitei]=lower_tensors[sitei]*upper_tensors[sitei];
        for (const auto &combiner : virt_leg_combiners_[sitei])
        {
            layered_site_tensors_[sitei]=layered_site_tensors_[sitei]*combiner;
        }
        clean(layered_site_tensors_[sitei]);
    }

}
template
void Double_Layer_PEPSt<ITensor>::obtain_layered_tensors_with_combined_legs(const std::vector<ITensor> &lower_tensors);
template
void Double_Layer_PEPSt<IQTensor>::obtain_layered_tensors_with_combined_legs(const std::vector<IQTensor> &lower_tensors);

