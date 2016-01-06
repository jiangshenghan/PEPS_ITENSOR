
#include "simple_update_patch_general.h"

template <class TensorT>
General_Patch_RDM<TensorT>::General_Patch_RDM(const std::string &name, PEPSt<TensorT> &peps, const TensorT &env_tens, const std::vector<int> &patch_sites, const std::vector<int> &cutting_sites):
    name_(peps.name()+name),
    patch_sites_(patch_sites),
    cutting_sites_(cutting_sites),
    env_tens_(env_tens),
    peps_(peps)
{
    init_patch_dim();
    init_env_tens();
    init_legs_combiners();
    init_patch_double_layer_tensors();
    obtain_RDM_and_patch_norm();

    //Print(name_);
    //Print(patch_sites_);
    //Print(cutting_sites_);
    //Print(patch_dim_);
    //PrintDat(env_tens_);
    //for (int i=0; i<patch_sites_.size(); i++)
    //{
    //    Print(patch_sites_[i]);
    //    Print(virt_legs_combiners_[i]);
    //}
    //Print(phys_legs_combiners_);
    //Print(patch_double_layer_tensors_);
    //Print(RDM_);
    //Print(patch_norm_);

}
template
General_Patch_RDM<ITensor>::General_Patch_RDM(const std::string &name, PEPSt<ITensor> &peps, const ITensor &env_tens, const std::vector<int> &patch_sites, const std::vector<int> &cutting_sites);
template
General_Patch_RDM<IQTensor>::General_Patch_RDM(const std::string &name, PEPSt<IQTensor> &peps, const IQTensor &env_tens, const std::vector<int> &patch_sites, const std::vector<int> &cutting_sites);


template <class TensorT>
void General_Patch_RDM<TensorT>::obtain_RDM_and_patch_norm()
{
    RDM_=tensor_from_contract_patch(patch_double_layer_tensors_);
    RDM_=RDM_*dag(phys_legs_combiners_[0])*dag(phys_legs_combiners_[1]);
    //get patch_norm_
    auto RDM_trace=RDM_;
    for (int cuti=0; cuti<cutting_sites_.size(); cuti++)
    {
        const auto &phys_leg=peps_.phys_legs(cutting_sites_[cuti]);
        RDM_trace.trace(phys_leg,dag(prime(phys_leg)));
    }
    patch_norm_=RDM_trace.toComplex().real();
}
template
void General_Patch_RDM<ITensor>::obtain_RDM_and_patch_norm();
template
void General_Patch_RDM<IQTensor>::obtain_RDM_and_patch_norm();


template <class TensorT>
std::vector<TensorT> General_Patch_RDM<TensorT>::double_layer_tensors_from_replaced_cutting_sites_tensors(std::vector<std::array<TensorT,2>> replaced_tensors_ket_bra)
{
    auto double_layer_tensors=patch_double_layer_tensors_;
    int comm_bond_no=peps_.lattice().comm_bond(cutting_sites_);
    for (int cuti=0; cuti<cutting_sites_.size(); cuti++)
    {
        absorb_env_tens(cutting_sites_[cuti],replaced_tensors_ket_bra[cuti]);
        int cut_coord=std::find(patch_sites_.begin(),patch_sites_.end(),cutting_sites_[cuti])-patch_sites_.begin();
        double_layer_tensors[cut_coord]=replaced_tensors_ket_bra[cuti][0]*replaced_tensors_ket_bra[cuti][1];
        absorb_neigh_bonds(cutting_sites_[cuti],double_layer_tensors[cut_coord],comm_bond_no);
    }
    return double_layer_tensors;
}
template
std::vector<ITensor> General_Patch_RDM<ITensor>::double_layer_tensors_from_replaced_cutting_sites_tensors(std::vector<std::array<ITensor,2>> replaced_tensors_ket_bra);
template
std::vector<IQTensor> General_Patch_RDM<IQTensor>::double_layer_tensors_from_replaced_cutting_sites_tensors(std::vector<std::array<IQTensor,2>> replaced_tensors_ket_bra);


template <class TensorT>
TensorT General_Patch_RDM<TensorT>::tensor_from_contract_patch(const std::vector<TensorT> &double_layer_tensors)
{
    //For square lattice, we separate the patch to left and right part
    if (name_.find("square")!=std::string::npos)
    {
        return tensor_from_contraction_part_patch(double_layer_tensors,0)*tensor_from_contraction_part_patch(double_layer_tensors,1);
    }
    //For Cirac kagome lattice

    return TensorT();
}
template
ITensor General_Patch_RDM<ITensor>::tensor_from_contract_patch(const std::vector<ITensor> &double_layer_tensors);
template
IQTensor General_Patch_RDM<IQTensor>::tensor_from_contract_patch(const std::vector<IQTensor> &double_layer_tensors);


template <class TensorT>
TensorT General_Patch_RDM<TensorT>::tensor_from_contraction_part_patch(const std::vector<TensorT> &double_layer_tensors, int parti)
{
    TensorT half_patch_tensor;
    //RDM for square lattice
    if (name_.find("square")!=std::string::npos)
    {
        if (name_.find("regular shape")!=std::string::npos)
        {
            //for parti=0/1, we start from first/last col, and end at col of cutting_sites_[parti]
            int start_col=parti*(patch_dim_[0]-1),
                final_col=(std::find(patch_sites_.begin(),patch_sites_.end(),cutting_sites_[parti])-patch_sites_.begin())%(patch_dim_[0]);
            for (int coli=start_col; ; coli+=1-2*parti)
            {
                for (int rowi=0; rowi<patch_dim_.size(); rowi++)
                {
                    if (!half_patch_tensor.valid())
                    {
                        half_patch_tensor=double_layer_tensors[rowi*patch_dim_[0]+coli];
                    }
                    else
                    {
                        half_patch_tensor*=double_layer_tensors[rowi*patch_dim_[0]+coli];
                    }
                }
                if (coli==final_col) break;
            }
        }
        if (name_.find("special shape I")!=std::string::npos)
        {
            if (parti==0)
            {
                half_patch_tensor=double_layer_tensors[0]*(double_layer_tensors[6]*double_layer_tensors[2])*double_layer_tensors[4];
            }
            if (parti==1)
            {
                half_patch_tensor=double_layer_tensors[1]*(double_layer_tensors[7]*double_layer_tensors[3])*double_layer_tensors[5];
            }
        }
    }

    //RDM for cirac kagome lattice

    return half_patch_tensor;
}
template
ITensor General_Patch_RDM<ITensor>::tensor_from_contraction_part_patch(const std::vector<ITensor> &double_layer_tensors, int parti);
template
IQTensor General_Patch_RDM<IQTensor>::tensor_from_contraction_part_patch(const std::vector<IQTensor> &double_layer_tensors, int parti);


template <class TensorT>
void General_Patch_RDM<TensorT>::init_patch_dim()
{
    if (name_.find("square")!=std::string::npos && name_.find("regular shape")!=std::string::npos)
    {
        int row_no=patch_sites_[0]/peps_.n_uc()[0],
            site_no=0;
        for (auto site : patch_sites_)
        {
            if (site/peps_.n_uc()[0]==row_no)
            {
                site_no++;
            }
            else
            {
                patch_dim_.push_back(site_no);
                row_no=site/peps_.n_uc()[0];
                site_no=1;
            }
        }
        patch_dim_.push_back(site_no);
    }
}
template
void General_Patch_RDM<ITensor>::init_patch_dim();
template
void General_Patch_RDM<IQTensor>::init_patch_dim();


template<>
void General_Patch_RDM<ITensor>::init_env_tens()
{
    for (int i=0; i<2; i++)
    {
        auto env_leg=env_tens_.indices()[i];
        Index new_env_leg("env_leg",env_leg.m(),env_leg.type(),env_leg.primeLevel());
        env_tens_.replaceIndex(env_leg,new_env_leg);
    }
}

template <>
void General_Patch_RDM<IQTensor>::init_env_tens()
{
    for (int i=0; i<2; i++)
    {
        auto env_leg=env_tens_.indices()[i];
        std::vector<IndexQN> new_indqns;
        for (const auto &indqn : env_leg.indices())
        {
            new_indqns.push_back(IndexQN(Index("env_indqn",indqn.m(),indqn.type(),indqn.primeLevel()),indqn.qn));
        }
        IQIndex new_env_leg("env_leg",new_indqns,env_leg.dir(),env_leg.primeLevel());
        env_tens_.replaceIndex(env_leg,new_env_leg);
    }
}


template <class TensorT>
void General_Patch_RDM<TensorT>::init_legs_combiners()
{
    //init virt_legs_combiners
    for (int sitei=0; sitei<patch_sites_.size(); sitei++)
    {
        std::vector<CombinerT> virt_legs_combiners_one_site;
        for (int legi=0; legi<peps_.lattice().n_bonds_to_one_site(); legi++)
        {
            auto curr_leg=peps_.virt_legs(patch_sites_[sitei],legi);
            virt_legs_combiners_one_site.push_back(CombinerT(curr_leg,dag(prime(curr_leg))));
        }
        virt_legs_combiners_.push_back(virt_legs_combiners_one_site);
    }

    //init phys_legs_combiners
    for (int i=0; i<cutting_sites_.size(); i++)
    {
        auto phys_leg=peps_.phys_legs(cutting_sites_[i]);
        phys_legs_combiners_.push_back(CombinerT(phys_leg,dag(prime(phys_leg))));
    }
}
template
void General_Patch_RDM<ITensor>::init_legs_combiners();
template
void General_Patch_RDM<IQTensor>::init_legs_combiners();


template <class TensorT>
void General_Patch_RDM<TensorT>::init_patch_double_layer_tensors()
{
    //get double layer site tensors
    for (int i=0; i<patch_sites_.size(); i++)
    {
        int sitei=patch_sites_[i];
        std::array<TensorT,2> ket_bra_tensors={peps_.site_tensors(sitei),prime(dag(peps_.site_tensors(sitei))).noprime(Site)};

        //obtain ket_bra_tensors with env_tens_ absorbed
        absorb_env_tens(sitei,ket_bra_tensors);

        std::vector<CombinerT> temp_combiner;
        auto cut_iter=std::find(cutting_sites_.begin(),cutting_sites_.end(),sitei);
        //if sitei is at cutting_sites_, then we make phys_legs of ket bra tensor different
        if (cut_iter!=cutting_sites_.end()) 
        {
            ket_bra_tensors[1].prime(Site);
            //take care of legs_num after contraction to prevent overflow.
            //we assume only tensors at BULK cutting_sites_ are possible to overflow
            int legs_num=legs_num_after_contraction<TensorT>({ket_bra_tensors[0],ket_bra_tensors[1]});
            if (legs_num>NMAX)
            {
                int n_legs=peps_.lattice().n_bonds_to_one_site();
                temp_combiner.push_back(CombinerT(peps_.virt_legs(sitei,n_legs-2),peps_.virt_legs(sitei,n_legs-1)));
                ket_bra_tensors[0]=ket_bra_tensors[0]*temp_combiner[0];
                ket_bra_tensors[1]=ket_bra_tensors[1]*dag(prime(temp_combiner[0]));
            }
            //obtain double layer site tensors with all top and bottom legs combined
            int cuti=cut_iter-cutting_sites_.begin();
            patch_double_layer_tensors_.push_back(ket_bra_tensors[0]*ket_bra_tensors[1]*phys_legs_combiners_[cuti]);
            for (const auto &virt_leg_combiner : virt_legs_combiners_[i])
            {
                const auto &indices=patch_double_layer_tensors_[i].indices();
                if (std::find(indices.begin(),indices.end(),*(virt_leg_combiner.left().begin()))!=indices.end())
                {
                    patch_double_layer_tensors_[i]=patch_double_layer_tensors_[i]*virt_leg_combiner;
                    if (!temp_combiner.empty())
                    {
                        patch_double_layer_tensors_[i]=patch_double_layer_tensors_[i]*prime(*(temp_combiner.end()-1))*dag(*(temp_combiner.end()-1));
                        temp_combiner.pop_back();
                    }
                    //Print(i);
                    //Print(sitei);
                    //Print(virt_leg_combiner);
                    //Print(patch_double_layer_tensors_[i]);
                }
            }
        }
        //for tensors not at cutting_sites_
        else
        {
            patch_double_layer_tensors_.push_back(ket_bra_tensors[0]*ket_bra_tensors[1]);
            for (const auto &virt_leg_combiner : virt_legs_combiners_[i])
            {
                const auto &indices=patch_double_layer_tensors_[i].indices();
                if (std::find(indices.begin(),indices.end(),*(virt_leg_combiner.left().begin()))!=indices.end())
                {
                    patch_double_layer_tensors_[i]=patch_double_layer_tensors_[i]*virt_leg_combiner;
                }
            }
        }
        //absorb the neighbouring bonds
        absorb_neigh_bonds(sitei,patch_double_layer_tensors_[i]);

        //Print(i);
        //Print(sitei);
        //Print(ket_bra_tensors[0]);
        //Print(ket_bra_tensors[1]);
        //Print(patch_double_layer_tensors_[i]);
            
    }

}
template
void General_Patch_RDM<ITensor>::init_patch_double_layer_tensors();
template
void General_Patch_RDM<IQTensor>::init_patch_double_layer_tensors();


template <class TensorT>
void General_Patch_RDM<TensorT>::absorb_env_tens(int sitei, std::array<TensorT,2> &ket_bra_tensors)
{
    int legi=0;
    for (int neigh_bond : peps_.lattice().site_neighbour_bonds(sitei))
    {
        auto curr_leg=peps_.virt_legs(sitei,legi);
        for (int neigh_site : peps_.lattice().bond_neighbour_sites(neigh_bond))
        {
            if (neigh_site==sitei) continue;
            auto site_iter=std::find(patch_sites_.begin(),patch_sites_.end(),neigh_site); 
            //if neigh_site is outside patch then absorb env_tens on boundary leg
            if (site_iter==patch_sites_.end())
            {
                obtain_env_dressed_tensor(ket_bra_tensors[0],env_tens_,curr_leg);
                ket_bra_tensors[1].noprime(dag(prime(curr_leg)));
                obtain_env_dressed_tensor(ket_bra_tensors[1],dag(env_tens_),dag(curr_leg));
                break;
            }
        }

        //Print(sitei);
        //Print(neigh_bond);
        //Print(ket_bra_tensors[0]);
        //Print(ket_bra_tensors[1]);

        legi++;
    }
}
template
void General_Patch_RDM<ITensor>::absorb_env_tens(int sitei, std::array<ITensor,2> &ket_bra_tensors);
template
void General_Patch_RDM<IQTensor>::absorb_env_tens(int sitei, std::array<IQTensor,2> &ket_bra_tensors);


template <class TensorT>
void General_Patch_RDM<TensorT>::absorb_neigh_bonds(int sitei, TensorT &double_layer_tensor, int forbid_bond)
{
    for (auto bondi : peps_.lattice().site_neighbour_bonds(sitei))
    {
        if (bondi==forbid_bond) continue;
        //make sure every bond only be absorbed once
        if (sitei!=peps_.lattice().bond_neighbour_sites(bondi,0)) continue;
        //make sure only absorb bonds in the bulk
        bool bond_in_bulk=true;
        auto double_layer_bondi_tensor=peps_.bond_tensors(bondi)*dag(prime(peps_.bond_tensors(bondi)));
        for (auto bondi_neigh_site : peps_.lattice().bond_neighbour_sites(bondi))
        {
            auto bondi_neigh_site_iter=std::find(patch_sites_.begin(),patch_sites_.end(),bondi_neigh_site);
            if (bondi_neigh_site_iter==patch_sites_.end()) 
            {
                bond_in_bulk=false;
                break;
            }
            int patch_site_no=bondi_neigh_site_iter-patch_sites_.begin();
            const auto &legs=peps_.lattice().site_neighbour_bonds(bondi_neigh_site);
            int leg_no=std::find(legs.begin(),legs.end(),bondi)-legs.begin();
            double_layer_bondi_tensor=double_layer_bondi_tensor*dag(virt_legs_combiners_[patch_site_no][leg_no]);
        }
        if (bond_in_bulk)
            double_layer_tensor*=double_layer_bondi_tensor;
    }
}
template
void General_Patch_RDM<ITensor>::absorb_neigh_bonds(int sitei, ITensor &double_layer_tensor, int forbid_bond);
template
void General_Patch_RDM<IQTensor>::absorb_neigh_bonds(int sitei, IQTensor &double_layer_tensor, int forbid_bond);
