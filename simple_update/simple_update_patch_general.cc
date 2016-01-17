
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
    //PrintDat(RDM_);
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
    for (int cuti=0; cuti<cutting_sites_.size(); cuti++)
    {
        RDM_=RDM_*dag(phys_legs_combiners_[cuti]);
    }
    //RDM_=RDM_*dag(phys_legs_combiners_[0])*dag(phys_legs_combiners_[1]);

    //get patch_norm_
    auto RDM_trace=RDM_;
    for (int cuti=0; cuti<cutting_sites_.size(); cuti++)
    {
        const auto &phys_leg=peps_.phys_legs(cutting_sites_[cuti]);
        RDM_trace.trace(phys_leg,dag(prime(phys_leg)));
    }
    patch_norm_=RDM_trace.toComplex().real();
    patch_norm_=sqrt(patch_norm_);
}
template
void General_Patch_RDM<ITensor>::obtain_RDM_and_patch_norm();
template
void General_Patch_RDM<IQTensor>::obtain_RDM_and_patch_norm();


template <class TensorT>
std::vector<TensorT> General_Patch_RDM<TensorT>::double_layer_tensors_from_replaced_cutting_sites_tensors(std::vector<std::array<TensorT,2>> replaced_tensors_ket_bra)
{
    //init replaced_tensors_bra
    for (auto &tensor_ket_bra : replaced_tensors_ket_bra)
    {
        std::array<IndexT,2> phys_leg_ket_bra={findtype(tensor_ket_bra[0],Site),findtype(tensor_ket_bra[1],Site)};
        if (phys_leg_ket_bra[0].dir()==phys_leg_ket_bra[1].dir())
        {
            tensor_ket_bra[1].dag();
            tensor_ket_bra[1].prime(Link);
        }
    }

    auto double_layer_tensors=patch_double_layer_tensors_;
    int comm_bond_no=peps_.lattice().comm_bond(cutting_sites_);
    auto comm_bond_tensor=peps_.bond_tensors(comm_bond_no);
    for (int cuti=0; cuti<cutting_sites_.size(); cuti++)
    {
        absorb_env_tens(cutting_sites_[cuti],replaced_tensors_ket_bra[cuti]);
        int cut_coord=std::find(patch_sites_.begin(),patch_sites_.end(),cutting_sites_[cuti])-patch_sites_.begin();
        double_layer_tensors[cut_coord]=replaced_tensors_ket_bra[cuti][0]*replaced_tensors_ket_bra[cuti][1];

        //combine top and bottom virtual legs. 
        //Notice, we do not touch legs of comm_bond
        for (const auto &virt_leg_combiner : virt_legs_combiners_[cut_coord])
        {
            const auto &indices=double_layer_tensors[cut_coord].indices();
            if (std::find(indices.begin(),indices.end(),*(virt_leg_combiner.left().begin()))!=indices.end() && std::find(comm_bond_tensor.indices().begin(),comm_bond_tensor.indices().end(),*(virt_leg_combiner.left().begin()))==comm_bond_tensor.indices().end())
            {
                double_layer_tensors[cut_coord]=double_layer_tensors[cut_coord]*virt_leg_combiner;
            }
        }

        //absorb neighbouring bonds (exclude comm bond)
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
        return tensor_from_contract_part_patch(double_layer_tensors,0)*tensor_from_contract_part_patch(double_layer_tensors,1);
    }
    //For Cirac kagome lattice, we separate the patch to three parts (u,v,w)
    if (name_.find("kagome cirac")!=std::string::npos)
    {
        if (name_.find("tree shape I")!=std::string::npos)
        {
            return double_layer_tensors[0]*double_layer_tensors[1]*double_layer_tensors[2];
        }
    }

    return TensorT();
}
template
ITensor General_Patch_RDM<ITensor>::tensor_from_contract_patch(const std::vector<ITensor> &double_layer_tensors);
template
IQTensor General_Patch_RDM<IQTensor>::tensor_from_contract_patch(const std::vector<IQTensor> &double_layer_tensors);


template <class TensorT>
TensorT General_Patch_RDM<TensorT>::tensor_from_contract_part_patch(const std::vector<TensorT> &double_layer_tensors, int parti)
{
    TensorT part_patch_tensor;
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
                    if (!part_patch_tensor.valid())
                    {
                        part_patch_tensor=double_layer_tensors[rowi*patch_dim_[0]+coli];
                    }
                    else
                    {
                        part_patch_tensor*=double_layer_tensors[rowi*patch_dim_[0]+coli];
                    }
                }
                if (coli==final_col) break;
            }
        }
        if (name_.find("special shape I")!=std::string::npos)
        {
            if (parti==0)
            {
                auto temp_tens=double_layer_tensors[6]*double_layer_tensors[2];
                //Print(double_layer_tensors[6]);
                //Print(temp_tens);
                part_patch_tensor=double_layer_tensors[0]*temp_tens*double_layer_tensors[4];
            }
            if (parti==1)
            {
                auto temp_tens=double_layer_tensors[7]*double_layer_tensors[3];
                part_patch_tensor=double_layer_tensors[1]*temp_tens*double_layer_tensors[5];
            }
        }
    }

    //RDM for cirac kagome lattice
    if (name_.find("kagome cirac")!=std::string::npos)
    {
        if (name_.find("tree shape I")!=std::string::npos)
        {
            return double_layer_tensors[parti];
        }
    }

    return part_patch_tensor;
}
template
ITensor General_Patch_RDM<ITensor>::tensor_from_contract_part_patch(const std::vector<ITensor> &double_layer_tensors, int parti);
template
IQTensor General_Patch_RDM<IQTensor>::tensor_from_contract_part_patch(const std::vector<IQTensor> &double_layer_tensors, int parti);


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
        //PrintDat(patch_double_layer_tensors_[i]);
            
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




void spin_square_peps_patch_simple_update(IQPEPS &square_peps, const Evolution_Params &square_su_params, std::vector<int> patch_sites, std::array<int,2> evolved_sites, std::string patch_name)
{
    //Initialize trotter gate
    NN_Heisenberg_Trotter_Gate evolve_gate({square_peps.phys_legs(evolved_sites[0]),square_peps.phys_legs(evolved_sites[1])});
    //Print(evolve_gate.site_tensors(0));

    //Initialize env_tens
    //env_tens[i] stores environment of site i, which is direct product of matrices for the legs of site i 
    //there is no env matrix between A and B
    std::array<std::vector<IQTensor>,2> env_tens;

    int comm_bond=square_peps.lattice().comm_bond(evolved_sites[0],evolved_sites[1]);
    auto comm_bond_tensor=square_peps.bond_tensors(comm_bond);
    //Print(comm_bond_tensor);

    //leg gates is used to approx evolve_gate, which formed by three indices, two in one out (primed). 
    //two in legs are contracting to the site tensor which has applied by trotter gate
    //one out leg is to contract comm bond tensor
    std::array<IndexSet<IQIndex>,2> leg_gates_indices;
    std::array<Singlet_Tensor_Basis,2> leg_gates_basis;
    for (int i=0; i<2; i++)
    {
        leg_gates_indices[i].addindex(commonIndex(dag(square_peps.site_tensors(evolved_sites[i])),comm_bond_tensor));
        leg_gates_indices[i].addindex(commonIndex(dag(evolve_gate.site_tensors(i)),evolve_gate.bond_tensors(0)));
        leg_gates_indices[i].addindex(commonIndex(square_peps.site_tensors(evolved_sites[i]),dag(comm_bond_tensor)).prime());

        leg_gates_basis[i]=Singlet_Tensor_Basis(leg_gates_indices[i]);

        //PrintDat(leg_gates_basis[i]);
    }

    //leg_gates_for_one_site is used to update site_tensors[evolved_sites[0]]
    std::vector<IQTensor> leg_gates_for_one_site;
    auto indice_from_evolve_gate=commonIndex(dag(evolve_gate.site_tensors(0)),evolve_gate.bond_tensors(0));
    for (const auto &leg_indice : square_peps.site_tensors(evolved_sites[0]).indices())
    {
        if (leg_indice.type()==Site) continue;

        std::vector<IQIndex> leg_gates_indices_for_one_site{dag(leg_indice),indice_from_evolve_gate,prime(leg_indice)};

        leg_gates_for_one_site.push_back(IQTensor(leg_gates_indices_for_one_site));
    }
    //Print(leg_gates_for_one_site);

    //leg_gate_params stores tunable parameters for the leg gate.
    std::vector<double> leg_gate_params(leg_gates_basis[0].dim());
    //random generate leg_gate_params
    for (auto &param : leg_gate_params) param=rand_gen();
    //Print(leg_gate_params);
    
    for (int iter=0; iter<square_su_params.iter_nums; iter++)
    {
        evolve_gate.change_time(square_su_params.ts[iter]);

        for (int step=0; step<square_su_params.steps_nums[iter]; step++)
        {
            Print(iter);
            Print(step);
            Print(square_su_params.ts[iter]);

            get_env_tensor_minimization(square_peps.site_tensors(evolved_sites[0])*comm_bond_tensor,square_peps.site_tensors(evolved_sites[1]),env_tens);

            //using general patch RDM
            General_Patch_RDM<IQTensor> square_RDM(patch_name,square_peps,env_tens[0][0],patch_sites,{evolved_sites[0],evolved_sites[1]});
            //measure energy by RDM
            Print(heisenberg_energy_from_RDM(square_RDM.RDM()));



            //if we cannot obtain a reasonable leg_gate, we try smaller time separation
            double cutoff=square_su_params.ts[iter]/10.;
            if (cutoff>1e-5) cutoff=1e-5;
            //double cutoff=1e-5;
            if (!obtain_spin_sym_leg_gates_params_minimization_from_RDM(square_RDM,evolve_gate,leg_gates_basis,leg_gate_params,cutoff)) break;

            //using leg_gate_params to generate all leg gates
            auto leg_gate_sample=singlet_tensor_from_basis_params(leg_gates_basis[0],leg_gate_params);

            for (auto &gate : leg_gates_for_one_site) 
            {
                tensor_assignment(gate,leg_gate_sample); 
                //Print(gate.indices());
            }


            //updated site_tensors[evolved_sites[0]]
            auto updated_site_tens=square_peps.site_tensors(evolved_sites[0]);

            for (const auto &leg_gate : leg_gates_for_one_site)
            {
                updated_site_tens*=evolve_gate.site_tensors(0)*leg_gate;
                updated_site_tens.noprime();
            }
            //we should never change order of indices of site tensor
            auto updated_site_tens_ordered_ind=square_peps.site_tensors(evolved_sites[0]);
            tensor_assignment_diff_order(updated_site_tens_ordered_ind,updated_site_tens);

            //check symmetry of updated_site_tens
            //auto sym_updated_site_tens=updated_site_tens_ordered_ind;
            //rotation_symmetrize_square_rvb_site_tensor(sym_updated_site_tens);
            //sym_updated_site_tens*=0.25;
            //Print(updated_site_tens_ordered_ind.norm());
            //Print(sym_updated_site_tens.norm());
            //Print((updated_site_tens_ordered_ind-sym_updated_site_tens).norm());
            //Print(diff_tensor_and_singlet_projection(updated_site_tens_ordered_ind));
            //Print(diff_tensor_and_singlet_projection(sym_updated_site_tens));

            //symmetrize updated_site_tens_ordered_ind, and using updated_site_tens_ordered_ind to generate all site tensors of peps
            rotation_symmetrize_square_rvb_site_tensor(updated_site_tens_ordered_ind);
            //updated_site_tens_ordered_ind*=0.25;
            //we always keep the same norm as original state
            updated_site_tens_ordered_ind*=(square_peps.site_tensors(evolved_sites[0]).norm())/(updated_site_tens_ordered_ind.norm());
            updated_site_tens_ordered_ind.clean();
            //Print(updated_site_tens_ordered_ind.norm());
            //Print(square_peps.site_tensors(evolved_sites[0]).norm());
            //Print((updated_site_tens_ordered_ind-square_peps.site_tensors(evolved_sites[0])).norm());

            square_peps.generate_site_tensors({updated_site_tens_ordered_ind});

            //stores as PEPS, which is used for further evolution
            if ((step+1)*10%square_su_params.steps_nums[iter]==0)
            {
                std::stringstream ss;

                //zero flux state
                if (std::abs(square_psg::mu_12-1)<EPSILON)
                {
                    ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_rvb_D=" << square_peps.D() << "_Lx=" << square_peps.n_uc()[0] << "_Ly=" << square_peps.n_uc()[1] << "_iter=" << iter << "_step=" << step << "_" << square_RDM.patch_name();
                }
                //pi flux state
                if (std::abs(square_psg::mu_12+1)<EPSILON)
                {
                    ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_pi_rvb_D=" << square_peps.D() << "_Lx=" << square_peps.n_uc()[0] << "_Ly=" << square_peps.n_uc()[1] << "_iter=" << iter << "_step=" << step << square_RDM.patch_name();
                }

                std::string file_name=ss.str();
                writeToFile(file_name,square_peps);
                
                //reinit leg_gate_params to get rid of local minimal
                //for (auto &param : leg_gate_params) param=rand_gen();
            }

        }//trotter steps

    }//trotter iters
}


void spin_kagome_cirac_peps_patch_simple_update(IQPEPS &kagome_rvb, const Evolution_Params &su_params, std::vector<int> patch_sites, std:vector<int> evolved_sites, std::string patch_name)
{
    //Initialize env_tens
    std::array<std::vector<IQTensor>,2> env_tens;

    int comm_bond_no=kagome_rvb.lattice().comm_bond(evolved_sites);

    std::array<std::vector<double>,2> leg_gates_params;

    double site_norm=kagome_rvb.site_tensors(evolved_sites[0]).norm(),
           bond_norm=kagome_rvb.bond_tensors(comm_bond_no).norm();

    for (int iter=0; iter<su_params.iter_nums; iter++)
    {
        //Init evolve gate
        IQTPO evolve_gate=trotter_gate_kagome_cirac({kagome_rvb.phys_legs(evolved_sites[0]),kagome_rvb.phys_legs(evolved_sites[1]),kagome_rvb.phys_legs(evolved_sites[2])},su_params.ts[iter]);

        //init leg gates for sites and bonds, which are used to approx evolve_gate
        //every leg gate is formed by two in legs and one out leg
        std::array<std::vector<Singlet_Tensor_Basis>,2> leg_gates_basis;
        for (int evolvei=0; evolvei<evolve_sites.size(); evolvei++)
        {
            IndexSet<IQIndex> leg_gate_indices;
            leg_gate_indices.addindex(commonIndex(dag(kagome_rvb.site_tensors(evolved_sites[evolvei])),kagome_rvb.bond_tensors(comm_bond_no)));
            leg_gate_indices.addindex(commonIndex(dag(evolve_gate.site_tensors(evolvei)),evolve_gate.bond_tensors(0)));
            leg_gate_indices.addindex(commonIndex(kagome_rvb.site_tensors(evolved_sites[evolvei]),dag(kagome_rvb.bond_tensors(comm_bond_no))).prime());
            leg_gates_basis[0].push_back(Singlet_Tensor_Basis(leg_gate_indices));
            leg_gates_basis[1].push_back(Singlet_Tensor_Basis(leg_gate_indices.dag()));
        }

        //init leg_gates_params
        for (int i=0; i<2; i++)
        {
            if (leg_gates_params[i].empty())
            {
                for (int parami=0; parami<leg_gates_basis[i][0].size(); parami++) leg_gates_params[i].push_back(rand_gen());
            }
        }

        //leg_gates_for_one_tensor[0/1] is for site/plaquette
        std::array<std::vector<IQTensor>,2> leg_gates_for_one_tensor;
        //init for site legs
        auto indice_from_evolve_gate=commonIndex(dag(evolve_gate.site_tensors(0)),evolve_gate.bond_tensors(0));
        for (const auto &virt_leg : kagome_rvb.site_tensors(evolved_sites[0]).indices())
        {
            if (virt_leg.type()==Site) continue;
            std::vector<IQIndex> leg_gate_indices{dag(virt_leg),indice_from_evolve_gate,prime(virt_leg)};
            leg_gates_for_one_tensor[0].push_back(IQTensor(leg_gate_indices));
        }
        //init for plaquette legs
        for (int legi=0; legi<evolved_sites.size(); legi++)
        {
            auto peps_virt_leg=commonIndex(kagome_rvb.bond_tensors(comm_bond_no),kagome_rvb.site_tensors(evolved_sites[legi]));
            auto evolve_gate_virt_leg=commonIndex(evolve_gate.bond_tensors(0),evolve_gate.site_tensors(legi));
           std::vector<IQIndex> leg_gate_indices{dag(peps_virt_leg),dag(evolve_gate_virt_leg),prime(peps_virt_leg)};
           leg_gates_for_one_tensor[1].push_back(IQTensor(leg_gate_indices));
        }

        for (int step=0; step<su_params.steps_nums[iter]; step++)
        {
            Print(iter);
            Print(step);
            Print(su_params.ts[iter]);

            //obtain env tensors, and measure energy by patch
            IQTensor site_tensA=kagome_rvb.site_tensors(0)*kagome_rvb.bond_tensors(0)*kagome_rvb.site_tensors(1)*kagome_rvb.site_tensors(2);
            IQTensor site_tensB=kagome_rvb.site_tensors({1,0,0})*kagome_rvb.bond_tensors({0,-1,1})*kagome_rvb.site_tensors({1,-1,1});
            get_env_tensor_minimization(site_tensA,site_tensB,env_tens);

            General_Patch_RDM<IQTensor> kagome_patch_RDM(patch_name,kagome_rvb,env_tens[0][0],patch_sites,evolved_sites);
            Print(heisenberg_energy_from_RDM(kagome_patch_RDM));

            //cutoff for leg gate 
            double cutoff=1e-5;
            if (!obtain_kagome_cirac_leg_gates_params_minimization(kagome_patch_RDM,evolve_gate,leg_gates_basis,leg_gates_params,cutoff)) break;

            //using leg_gates_params to generate all leg_gates for one site and one plaquette
            for (int i=0; i<2; i++)
            {
                auto leg_gate_sample=singlet_tensor_from_basis_params(leg_gates_basis[i][0],leg_gates_params[i]);

                for (const auto &gate : leg_gates_for_one_tensor[i])
                {
                    tensor_assignment(gate,leg_gate_sample);
                }

            }

            //updte site_tensors
            IQTensor updated_site_tens_unordered=kagome_rvb.site_tensors(evolved_sites[0]);
            for (const auto &site_leg_gate : leg_gates_for_one_tensor[0]) 
            {
                updated_site_tens_unordered*=evolve_gate.site_tensors(0)*site_leg_gate; 
            }
            updated_site_tens_unordered.noprime();
            //we should never change order of indices of site tensors
            auto updated_site_tens=kagome_rvb.site_tensors(evolved_sites[0]);
            tensor_assignment_diff_order(updated_site_tens,updated_site_tens_unordered);
            rotation_symmetrize_kagome_site_tensor(updated_site_tens);
            //keep the same norm
            updated_site_tens*=site_norm/(updated_site_tens.norm());
            kagome_rvb.generate_site_tensors(updated_site_tens,updated_site_tens,tensor_permutation({0,2,1},updated_site_tens));

            //updated bond tensors
            IQTensor updated_bond_tens_unordered=kagome_rvb.bond_tensors(comm_bond_no);
            for (const auto &bond_leg_gate : leg_gates_for_one_tensor[1])
            {
                updated_bond_tens_unordered*=evolve_gate.bond_tensors(0)*bond_leg_gate;
            }
            updated_bond_tens_unordered.noprime();
            //we should never change order of inds
            auto updated_bond_tens=kagome_rvb.bond_tensors(comm_bond_no);
            tensor_assignment_diff_order(updated_bond_tens,updated_bond_tens_unordered);
            rotation_symmetrize_kagome_bond_tensor(updated_bond_tens);
            //keep the same norm
            updated_bond_tens*=bond_norm/(updated_bond_tens.norm());
            kagome_rvb.generate_bond_tensors({bond_tensor,tensor_permutation({2,0,1},bond_tensor)},kagome_psg::mu_12);

            //stores as PEPS
            if ((step+1)*10%su_params.step_nums[iter]==0)
            {
                std::stringstream ss;

                //zero flux state
                ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/kagome_rvb_D=" << kagome_rvb.D() << "_Lx=" << kagome_rvb.n_uc()[0] << "_Ly=" << kagome_rvb.n_uc()[1] << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_iter=" << iter << "_step=" << step << "_" << kagome_patch_RDM.patch_name();

                std::string file_name=ss.str();
                writeToFile(file_name,kagome_rvb);

            }

        }//end of trotter steps

    }// end of trotter iters
}



bool obtain_spin_sym_leg_gates_params_minimization_from_RDM(General_Patch_RDM<IQTensor> &square_RDM, const Trotter_Gate &trotter_gate, const std::array<Singlet_Tensor_Basis,2> &leg_gates_basis, std::vector<double> &leg_gate_params, double cutoff)
{
    //init leg_gate_params
    if (leg_gate_params.empty())
    {
        for (int i=0; i<leg_gates_basis[0].dim(); i++)
        {
            leg_gate_params.push_back(rand_gen());
        }
        Print(leg_gate_params);
    }

    //prepare evolved site tensors and bond tensors
    std::array<IQTensor,2> site_tensors_evolved={square_RDM.cutting_site_tensors(0),square_RDM.cutting_site_tensors(1)};
    IQTensor bond_tensor_evolved(square_RDM.cutting_bond_tensor());

    for (int i=0; i<2; i++)
    {
        site_tensors_evolved[i]*=trotter_gate.site_tensors(i);
        site_tensors_evolved[i].noprime();
    }
    bond_tensor_evolved*=trotter_gate.bond_tensors(0);

    //evolve_legs_combiners are used to combine legs of peps site_tensors and trotter gate site tensors
    std::array<IQCombiner,2> evolve_legs_combiners={IQCombiner(square_RDM.cutting_virt_legs(0),trotter_gate.virt_legs(0)[0]),IQCombiner(square_RDM.cutting_virt_legs(1),trotter_gate.virt_legs(1)[0])};

    //two method to obtain evolved_wf_norm
    //1. directly use two_sites_RDM
    double evolved_wf_norm=((square_RDM.RDM()*dag(trotter_gate.site_tensors(0)*trotter_gate.bond_tensors(0)*trotter_gate.site_tensors(1)).prime()).mapprime(2,1)*(trotter_gate.site_tensors(0)*trotter_gate.bond_tensors(0)*trotter_gate.site_tensors(1))).toComplex().real();
    evolved_wf_norm=sqrt(evolved_wf_norm);
    //Print(evolved_wf_norm);
    //2. replace cutting_tensors with tensors_evolved, and obtain expectation value
    //evolved_wf_norm=sqrt(square_RDM.expect_val_from_replaced_tensors({{site_tensors_evolved[0]*bond_tensor_evolved*dag(evolve_legs_combiners[1]),site_tensors_evolved[0]*bond_tensor_evolved*dag(evolve_legs_combiners[1])},{site_tensors_evolved[1]*evolve_legs_combiners[1],site_tensors_evolved[1]*evolve_legs_combiners[1]}}).real());
    //Print(evolved_wf_norm);


    //get M_{ijkl}=\langle\phi_{ij}|\phi_{kl}\rangle, store in a vector updated_wf_basis_overlap
    std::array<std::vector<IQTensor>,2> site_tensors_updated_basis;
    for (int i=0; i<2; i++)
    {
        for (const auto &leg_base : leg_gates_basis[i])
        {
            site_tensors_updated_basis[i].push_back(noprime(site_tensors_evolved[i]*leg_base));
        }
        //Print(i);
        //Print(site_tensors_updated_basis[i]);
    }
    int N_leg_basis=site_tensors_updated_basis[0].size();

    //get left half and right half respectively
    std::vector<IQTensor> left_half_basis, right_half_basis;
    for (int base_num=0; base_num<N_leg_basis*N_leg_basis; base_num++)
    {
        std::array<int,2> base_list={base_num/N_leg_basis,base_num%N_leg_basis};

        std::vector<std::array<IQTensor,2>> left_replaced_tensors_ket_bra={{site_tensors_updated_basis[0][base_list[1]]*square_RDM.cutting_bond_tensor(),site_tensors_updated_basis[0][base_list[0]]*square_RDM.cutting_bond_tensor()},{site_tensors_updated_basis[1][0],site_tensors_updated_basis[1][0]}},
                                            right_replaced_tensors_ket_bra={{site_tensors_updated_basis[0][0]*square_RDM.cutting_bond_tensor(),site_tensors_updated_basis[0][0]*square_RDM.cutting_bond_tensor()},{site_tensors_updated_basis[1][base_list[1]],site_tensors_updated_basis[1][base_list[0]]}};

        auto double_layer_tensors_for_left_half_basis=square_RDM.double_layer_tensors_from_replaced_cutting_sites_tensors(left_replaced_tensors_ket_bra);
        auto double_layer_tensors_for_right_half_basis=square_RDM.double_layer_tensors_from_replaced_cutting_sites_tensors(right_replaced_tensors_ket_bra);

        left_half_basis.push_back(square_RDM.tensor_from_contract_part_patch(double_layer_tensors_for_left_half_basis,0));
        right_half_basis.push_back(square_RDM.tensor_from_contract_part_patch(double_layer_tensors_for_right_half_basis,1));

        //TODO:debug the speed
        //Print(base_num);
        //Print(double_layer_tensors_for_left_half_basis);
        //Print(double_layer_tensors_for_right_half_basis);
        //Print(*(left_half_basis.end()-1));
        //Print(*(right_half_basis.end()-1));
    }

    int total_base_num=1;
    std::vector<int> max_base_list;
    for (int i=0; i<4; i++)
    {
        total_base_num*=N_leg_basis;
        max_base_list.push_back(N_leg_basis);
    }
    std::vector<double> updated_wf_basis_overlap;
    //int nonzero_M=0;
    for (int base_num=0; base_num<total_base_num; base_num++)
    {
        auto base_list=list_from_num(base_num,max_base_list);
        //only get nonzero result when spin of contraction indice matches
        if (leg_gates_basis[0].spin_configs(base_list[0])[2]!=leg_gates_basis[1].spin_configs(base_list[1])[2])
        {
            updated_wf_basis_overlap.push_back(0);
            continue;
        }
        if (leg_gates_basis[0].spin_configs(base_list[2])[2]!=leg_gates_basis[1].spin_configs(base_list[3])[2])
        {
            updated_wf_basis_overlap.push_back(0);
            continue;
        }

        auto M_ijkl=(left_half_basis[base_list[0]*N_leg_basis+base_list[2]]*right_half_basis[base_list[1]*N_leg_basis+base_list[3]]).toComplex();
        //Print(base_num);
        //Print(M_ijkl);
        updated_wf_basis_overlap.push_back(M_ijkl.real());

    }

    //get w_ij=\langle\phi_{ij}|\psi\rangle+\langle\psi|\phi_{ij}\rangle,
    //stored in updated_wf_basis_evolved_wf_overlap
    std::vector<double> updated_wf_basis_evolved_wf_overlap;
    for (int base_num=0; base_num<N_leg_basis*N_leg_basis; base_num++)
    {
        std::array<int,2> base_list={base_num/N_leg_basis,base_num%N_leg_basis};
        if (leg_gates_basis[0].spin_configs(base_list[0])[2]!=leg_gates_basis[1].spin_configs(base_list[1])[2]) 
        {
            updated_wf_basis_evolved_wf_overlap.push_back(0);
            continue;
        }
        auto w_ij=square_RDM.expect_val_from_replaced_tensors({{site_tensors_evolved[0]*bond_tensor_evolved*dag(evolve_legs_combiners[1]),site_tensors_updated_basis[0][base_list[0]]*square_RDM.cutting_bond_tensor()},{site_tensors_evolved[1]*evolve_legs_combiners[1],site_tensors_updated_basis[1][base_list[1]]}});

        //Print(leg_gates_basis[0].spin_configs(base_list[0]));
        //Print(leg_gates_basis[1].spin_configs(base_list[1]));
        //Print(w_ij);

        updated_wf_basis_evolved_wf_overlap.push_back(2*w_ij.real());
    }
    //Print(updated_wf_basis_evolved_wf_overlap);


    //using conjugate gradient minimization to minimize distance square between updated wf and time evolved wf
    int find_min_status;
    int iter=0, max_iter=1e4;

    const gsl_multimin_fdfminimizer_type *minimize_T;
    gsl_multimin_fdfminimizer *s;

    Wf_Distance_Params *wf_distance_params=new Wf_Distance_Params(N_leg_basis,evolved_wf_norm,updated_wf_basis_overlap,updated_wf_basis_evolved_wf_overlap);

    //x stores coefficient for leg gates
    gsl_vector *x;
    gsl_multimin_function_fdf wf_distance_func;

    wf_distance_func.n=leg_gate_params.size();
    wf_distance_func.f=wf_distance_f;
    wf_distance_func.df=wf_distance_df;
    wf_distance_func.fdf=wf_distance_fdf;
    wf_distance_func.params=wf_distance_params;

    x=gsl_vector_alloc(leg_gate_params.size());
    for (int i=0; i<leg_gate_params.size(); i++) gsl_vector_set(x,i,leg_gate_params[i]);

    minimize_T=gsl_multimin_fdfminimizer_conjugate_fr;
    s=gsl_multimin_fdfminimizer_alloc(minimize_T,leg_gate_params.size());

    gsl_multimin_fdfminimizer_set(s,&wf_distance_func,x,0.1,0.1);

    do
    {
        iter++;
        find_min_status=gsl_multimin_fdfminimizer_iterate(s);
        if (find_min_status) break;
        find_min_status=gsl_multimin_test_gradient(s->gradient,cutoff);

        //Print(iter);
        //Print(s->f);
    }
    while (find_min_status==GSL_CONTINUE && iter<max_iter);

    //Print(iter);
    //Print(s->f);
    //normalized distance
    double wf_norm=square_RDM.wf_norm(),
           wf_evolved_wf_overlap=2*square_RDM.expect_val_from_replaced_tensors({{site_tensors_evolved[0]*bond_tensor_evolved*dag(evolve_legs_combiners[1]),square_RDM.cutting_site_tensors(0)*square_RDM.cutting_bond_tensor()},{site_tensors_evolved[1]*evolve_legs_combiners[1],square_RDM.cutting_site_tensors(1)}}).real(),
           //we get distance with original wf resize to the norm of evolved_wf
           wf_evolved_wf_distance=std::sqrt(2*evolved_wf_norm*evolved_wf_norm-evolved_wf_norm/wf_norm*wf_evolved_wf_overlap);
    double updated_wf_evolved_wf_distance=std::sqrt(wf_distance_f(s->x,wf_distance_params));
    Print(wf_evolved_wf_distance/evolved_wf_norm);
    Print(updated_wf_evolved_wf_distance/evolved_wf_norm);

    if (iter==max_iter)
    {
        cout << "Leg gate is not good enough, may be trapped in local minima!" << endl << "try smaller time step!" << endl;
        gsl_multimin_fdfminimizer_free(s);
        gsl_vector_free(x);
        delete wf_distance_params;

        //leg_gate_params.clear();
        //obtain_spin_sym_leg_gates_params_minimization(site_tensors,bond_tensor,trotter_gate,leg_gates_basis,leg_gate_params);
        return false;
    }

    for (int i=0; i<leg_gate_params.size(); i++)
        leg_gate_params[i]=gsl_vector_get(s->x,i);
    Print(leg_gate_params);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
    delete wf_distance_params;

    return true;
}

bool obtain_kagome_cirac_leg_gates_params_minimization(General_Patch_RDM<IQTensor> &kagome_patch_RDM, const IQTPO &evolve_gate, const std::array<std::vector<Singlet_Tensor_Basis>,2> &leg_gates_basis, std::array<std::vector<double>,2> &leg_gates_params, double cutoff=1E-5)
{
    //construct site tensors and bond tensors after applied by evolve_gate
    std::vector<IQTensor> site_tensors_evolved;
    //evolve_legs_combiners are used to combine legs of peps site tensors and evolve gate site tensors
    std::vector<IQCombiner> evolve_legs_combiners;
    for (int cuti=0; cuti<kagome_patch_RDM.cutting_sites_no(); cuti++)
    {
        site_tensors_evolved.push_back(kagome_patch_RDM.cutting_site_tensors(cuti)*evolve_gate.site_tensors(cuti));
        site_tensors_evolved[cuti].noprime();
        auto evolve_gate_virt_leg=commonIndex(evolve_gate.site_tensors(cuti),evolve_gate.bond_tensors(cuti));
        evolve_legs_combiners.push_back(IQCombiner(kagome_patch_RDM.cutting_virt_legs(cuti),evolve_gate_virt_leg));
    }
    auto bond_tensor_evolved=kagome_patch_RDM.cutting_bond_tensor()*evolve_gate.bond_tensors(0);

    //obtain evolved_wf_norm
    IQTensor evolve_gate_tensor=evolve_gate.bond_tensors(0);
    for (const auto &tens : evolve_gate.site_tensors()) evolve_gate_tensor*=tens;
    double evolved_wf_norm=sqrt((kagome_patch_RDM.RDM()*(evolve_gate_tensor*dag(swapPrime(evolve_gate_tensor,0,2))).mapprime(2,1)).real());

    //using conjugate gradient methods to minimize distance square between updated_wf and evolved_wf
    int find_min_status;
    int iter=0, max_iter=1E4;

    const gsl_multimin_fdfminimizer_type *minimize_T;
    gsl_multimin_fdfminimizer *s;

    //params to do minimization
    Kagome_Cirac_Wf_Distance_Params kagome_cirac_wf_distance_params(evolved_wf_norm,kagome_patch_RDM,site_tensors_evolved,bond_tensor_evolved,evolve_legs_combiners,leg_gates_basis);
    
    //x stores coefficient for site leg gates and plaquette leg gates
    gsl_vector *x;
    gsl_multimin_function_fdf wf_distance_sq_func;

    wf_distance_sq_func.n=leg_gates_params[0].size()+leg_gates_params[1].size();
    wf_distance_sq_func.f=kagome_cirac_wf_distance_sq_f;
    wf_distance_sq_func.df=kagome_cirac_wf_distance_sq_df;
    wf_distance_sq_func.fdf=kagome_cirac_wf_distance_sq_fdf;

    x=gsl_vector_alloc(wf_distance_sq_func.n);
    for (int i=0; i<leg_gates_params[0].size(); i++) gsl_vector_set(x,i,leg_gates_params[0][i]);
    for (int i=0; i<leg_gates_params[1].size(); i++) gsl_vector_set(x,leg_gates_params[0].size()+i,leg_gates_params[1][i]);

    minimize_T=gsl_multimin_fdfminimizer_conjugate_fr;
    s=gsl_multimin_fdfminimizer_alloc(minimize_T,wf_distance_sq_func.n);
    gsl_multimin_fdfminimizer_set(s,&wf_distance_sq_func,x,0.1,0.1);

    do
    {
        iter++;
        find_min_status=gsl_multimin_fdfminimizer_iterate(s);
        if (find_min_status) break;
        find_min_status=gsl_multimin_test_gradient(s->gradient,cutoff);
    }
    while (find_min_status==GSL_CONTINUE && iter<max_iter);

    double wf_norm=kagome_patch_RDM.wf_norm(),
           wf_evolved_wf_overlap=2*(kagome_patch_RDM.expect_val_from_RDM(evolve_gate_tensor)).real(),
           wf_evolved_wf_distance=std::sqrt(2*pow(evolved_wf_norm,2.)-evolved_wf_norm/wf_norm*wf_evolved_wf_overlap),
           updated_wf_evolved_wf_distance=sqrt(wf_distance_f(s->x,kagome_cirac_wf_distance_sq_params));

    Print(wf_evolved_wf_distance/evolved_wf_norm);
    Print(updated_wf_basis_evolved_wf_overlap/evolved_wf_norm);

    if (iter==max_iter)
    {
        cout << "Leg gate is not good enough, may be trapped in local minima!" << endl << "try smaller time step!" << endl;
        gsl_multimin_fdfminimizer_free(s);
        gsl_vector_free(x);
        delete kagome_cirac_wf_distance_params;

        return false;
    }

    for (int i=0; i<leg_gates_params[0].size(); i++)
        leg_gates_params[0][i]=gsl_vector_get(s->x,i);
    for (int i=0; i<leg_gates_params[1].size(); i++)
        leg_gates_params[1][i]=gsl_vector_get(s->x,leg_gates_params[0].size()+i);

    Print(leg_gates_params[0]);
    Print(leg_gates_params[1]);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
    delete kagome_cirac_wf_distance_params;

    return true;
}


double kagome_cirac_wf_distance_f(const gsl_vector *x, void *params)
{
    std::array<std::vector<double>,2> leg_gates_params;
    Kagome_Cirac_Wf_Distance_Params *kagome_cirac_wf_distance_params=(Kagome_Cirac_Wf_Distance_Params *)params;
    std::array<int,2> N_legs_basis={kagome_cirac_wf_distance_params->leg_gates_basis[0][0].size(),kagome_cirac_wf_distance_params->leg_gates_basis[1][0].size()};
   
    for (int i=0; i<N_legs_basis[0]; i++)
        leg_gates_params[0].push_back(gsl_vector_get(x,i));
    for (int i=0; i<N_legs_basis[1]; i++)
        leg_gates_params[1].push_back(gsl_vector_get(x,i+N_legs_basis[0]));

    //construct leg gates
    std::array<std::vector<IQTensor>,2> leg_gates;
    for (int typei=0; typei<2; typei++)
    {
        for (const auto &singlet_basis : kagome_cirac_wf_distance_params->leg_gates_basis[typei])
        {
            leg_gates[typei].push_back(singlet_tensor_from_basis_params(singlet_basis,leg_gates_params[typei]));
        }
    }

    //updated sites and bonds are obtained by multiplication of evolved_sites (bonds) with leg gates
    auto evolved_site_tensors=kagome_cirac_wf_distance_params->evolved_site_tensors;
    auto evolved_bond_tensor=kagome_cirac_wf_distance_params->evolved_bond_tensor;
    std::vector<IQTensor> updated_site_tensors;
    IQTensor updated_bond_tensor=evolved_bond_tensor;
    const auto &evolve_legs_combiners=kagome_cirac_wf_distance_params->evolve_legs_combiners;
    for (int cuti=0; cuti<3; cuti++)
    {
        updated_site_tensors.push_back(evolved_site_tensors[cuti]*leg_gates[0][cuti]);
        updated_site_tensors[cuti].noprime();
        updated_bond_tensor*=leg_gates[1][cuti];
        updated_bond_tensor.noprime();

        //combine virtual legs of evolved_tensors
        evolved_site_tensors[cuti]=evolved_site_tensors[cuti]*dag(evolve_legs_combiners[cuti]);
        evolved_bond_tensor=evolved_bond_tensor*dag(evolved_legs_combiner[cuti]);
    }


    double updated_wf_norm_sq=(kagome_cirac_wf_distance_params->kagome_patch_RDM.expect_val_from_replaced_tensors({updated_site_tensors[0]*updated_bond_tensor,updated_site_tensors[1],updated_site_tensors[2]})).real();

    double updated_evolved_wf_overlap=(2*kagome_cirac_wf_distance_params->kagome_patch_RDM.expect_val_from_replaced_tensors({{updated_site_tensors[0]*updated_bond_tensor,evolved_site_tensors[0]*evolved_bond_tensor},{updated_site_tensors[1],evolved_site_tensors[1]},{updated_site_tensors[2],evolved_site_tensors[2]}})).real();

    double distance_sq=updated_wf_norm_sq+pow(kagome_cirac_wf_distance_params->evolved_wf_norm,2.)-updated_evolved_wf_overlap;

    return distance_sq;
}

void kagome_cirac_wf_distance_df(const gsl_vector *x, void *params, gsl_vector *df)
{
    int x_size=x->size;

    gsl_vector *x_plus_dx;
    x_plux_dx=gsl_vector_alloc(x_size);
    gsl_vector_memcpy(x_plux_dx,x);

    double f_x=kagome_cirac_wf_distance_f(x,params);
    for (int i=0; i<x_size; i++)
    {
        double dxi=1E-10;
        gsl_vector_set(x_plus_dx,i,gsl_vector_get(x,i),dxi);
        double f_x_plus_dxi=kagome_cirac_wf_distance_f(x_plus_dx,params);
        gsl_vector_set(df,i,(f_x_plus_dxi-f_x)/dxi);
        gsl_vector_set(x_plus_dx,i,gsl_vector_get(x,i));
    }

    gsl_vector_free(x_plus_dx);
}

void kagome_cirac_wf_distance_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
    *f=kagome_cirac_wf_distance_f(x,params);
    wf_distance_df(x,params,df);
}

double heisenberg_energy_from_RDM(const General_Patch_RDM<IQTensor> &patch_rdm)
{
    std::vector<IQIndex> phys_legs;
    for (int cuti=0; cuti<patch_rdm.cutting_sites().size(); cuti++)
    {
        phys_legs.push_back(patch_rdm.cutting_phys_legs(cuti));
    }

    double energy=0;

    if (patch_rdm.peps_name().find("square")!=std::string::npos)
    {
        NN_Heisenberg_Hamiltonian heisenberg_gate({phys_legs[0],phys_legs[1]});
        energy=(patch_rdm.RDM()*(heisenberg_gate.site_tensors(0)*heisenberg_gate.bond_tensor()*heisenberg_gate.site_tensors(1))).toComplex().real();
    }

    //Print(patch_rdm.patch_name());
    if (patch_rdm.peps_name().find("kagome cirac")!=std::string::npos)
    {
        IQTPO heisenberg_gate=SpinSpin_kagome_cirac(phys_legs);
        energy=(patch_rdm.RDM()*heisenberg_gate.site_tensors(0)*heisenberg_gate.site_tensors(1)*heisenberg_gate.bond_tensors(0)*heisenberg_gate.site_tensors(2)).toComplex().real();
    }

    return energy/std::pow(patch_rdm.wf_norm(),2.);
}

