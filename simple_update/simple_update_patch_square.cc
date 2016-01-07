
#include "simple_update_patch_square.h"

//class Square_Patch_RDM
Square_Patch_RDM::Square_Patch_RDM(const IQPEPS &square_peps, const IQTensor &env_tens, const std::vector<std::vector<int>> &patch_sites, const std::array<int,2> cutting_sites):
    patch_dim_{patch_sites.size(),patch_sites[0].size()},
    patch_sites_(patch_sites),
    cutting_sites_(cutting_sites),
    env_tens_(env_tens),
    square_peps_(square_peps)
{
    init_patch_cutting_coords();
    init_env_tens();
    init_patch_tensors();
    init_legs_combiners();
    init_patch_double_layer_tensors();
    obtain_two_sites_RDM();
    wf_norm_=std::sqrt(trace(two_sites_RDM_,cutting_phys_legs(0),prime(dag(cutting_phys_legs(0)))).trace(cutting_phys_legs(1),prime(dag(cutting_phys_legs(1)))).toComplex().real());
}

void Square_Patch_RDM::obtain_two_sites_RDM()
{
    if (patch_dim_[0]>=patch_dim_[1])
    {
        for (int rowi=0; rowi<patch_dim_[0]; rowi++)
        {
            for (int coli=0; coli<patch_dim_[1]; coli++)
            {
                if (rowi==0 && coli==0) 
                    two_sites_RDM_=patch_double_layer_tensors_[rowi][coli];
                else
                    two_sites_RDM_*=patch_double_layer_tensors_[rowi][coli];

                //Print(rowi);
                //Print(coli);
                //Print(patch_double_layer_tensors_[rowi][coli]);
                //Print(two_sites_RDM_);
            }
        }
    }
    else
    {
        for (int coli=0; coli<patch_dim_[1]; coli++)
        {
            for (int rowi=0; rowi<patch_dim_[0]; rowi++)
            {
                if (rowi==0 && coli==0) 
                    two_sites_RDM_=patch_double_layer_tensors_[rowi][coli];
                else
                    two_sites_RDM_*=patch_double_layer_tensors_[rowi][coli];

                //Print(rowi);
                //Print(coli);
                //Print(patch_double_layer_tensors_[rowi][coli]);
                //Print(two_sites_RDM_);
            }
        }
    }

    //decombine physical legs
    two_sites_RDM_=two_sites_RDM_*dag(phys_legs_combiners_[0])*dag(phys_legs_combiners_[1]);
}

std::vector<std::vector<IQTensor>> Square_Patch_RDM::double_layer_tensors_from_replaced_tensors(std::array<IQTensor,2> replaced_tensors, std::array<IQTensor,2> replaced_tensors_dag) const
{
    for (auto &tensor : replaced_tensors_dag)
    {
        tensor=dag(tensor).prime().noprime(Site);
    }

    //combine common legs of replaced_tensors(_dag) to avoid leg overall flow
    combine_comm_legs(replaced_tensors);
    combine_comm_legs(replaced_tensors_dag);

    //Print(replaced_tensors[0]);
    //Print(replaced_tensors[1]);
    //Print(replaced_tensors_dag[0]);
    //Print(replaced_tensors_dag[1]);

    std::vector<std::vector<IQTensor>> replaced_double_layer_tensors(patch_double_layer_tensors_);

    for (int i=0; i<2; i++)
    {
        std::vector<IQIndex> boundary_legs;
        //we only combine top and bottom legs on combine_dir=true direction
        std::vector<bool> combine_dir(4,true);
        //down edge case
        if (patch_cutting_coords_[i][0]==0)
        {
            boundary_legs.push_back(square_peps_.virt_legs(cutting_sites_[i],3));
            combine_dir[3]=false;
        }
        //up edge
        if (patch_cutting_coords_[i][0]==patch_dim_[0]-1)
        {
            boundary_legs.push_back(square_peps_.virt_legs(cutting_sites_[i],1));
            combine_dir[1]=false;
        }
        //left edge
        if (patch_cutting_coords_[i][1]==0)
        {
            boundary_legs.push_back(square_peps_.virt_legs(cutting_sites_[i],0));
            combine_dir[0]=false;
        }
        //right edge
        if (patch_cutting_coords_[i][1]==patch_dim_[1]-1)
        {
            boundary_legs.push_back(square_peps_.virt_legs(cutting_sites_[i],2));
            combine_dir[2]=false;
        }

        //Print(combine_dir);
        //Print(boundary_legs);

        //modify to boundary tensors
        for (const auto &boundary_leg : boundary_legs)
        {
            obtain_env_dressed_tensor(replaced_tensors[i],env_tens_,boundary_leg);
            replaced_tensors_dag[i].noprime(dag(prime(boundary_leg)));
            obtain_env_dressed_tensor(replaced_tensors_dag[i],dag(env_tens_),dag(boundary_leg));
            //Print(replaced_tensors[i]);
            //Print(replaced_tensors_dag[i]);
        }


        //we already count the bond between two cutting sites in the replaced tensor, so we do not absorb that bond
        int diff_row=patch_cutting_coords_[1-i][0]-patch_cutting_coords_[i][0],
            diff_col=patch_cutting_coords_[1-i][1]-patch_cutting_coords_[i][1];
        if (diff_row==1 && diff_col==0) combine_dir[1]=false;
        if (diff_row==-1 && diff_col==0) combine_dir[3]=false;
        if (diff_row==0 && diff_col==1) combine_dir[2]=false;
        if (diff_row==0 && diff_col==-1) combine_dir[0]=false;

        //Print(diff_row);
        //Print(diff_col);
        //Print(combine_dir);
        //Print(replaced_tensors[i]);
        //Print(replaced_tensors_dag[i]);

        //absorb bond legs and combine top and bottom legs
        auto &curr_row=patch_cutting_coords_[i][0],
             &curr_col=patch_cutting_coords_[i][1];
        auto &curr_tensor=replaced_double_layer_tensors[curr_row][curr_col];
        curr_tensor=replaced_tensors[i]*replaced_tensors_dag[i];

        if (combine_dir[0])
        {
            curr_tensor=curr_tensor*left_legs_combiners_[curr_row][curr_col];
        }
        if (combine_dir[1])
        {
            int ubond=square_peps_.lattice().comm_bond(patch_sites_[curr_row][curr_col],patch_sites_[curr_row+1][curr_col]);
            curr_tensor=curr_tensor*square_peps_.bond_tensors(ubond)*dag(prime(square_peps_.bond_tensors(ubond)))*dag(down_legs_combiners_[curr_row+1][curr_col]);
        }
        if (combine_dir[2])
        {
            int rbond=square_peps_.lattice().comm_bond(patch_sites_[curr_row][curr_col],patch_sites_[curr_row][curr_col+1]);
            curr_tensor=curr_tensor*square_peps_.bond_tensors(rbond)*dag(prime(square_peps_.bond_tensors(rbond)))*dag(left_legs_combiners_[curr_row][curr_col+1]);
        }
        if (combine_dir[3])
        {
            curr_tensor=curr_tensor*down_legs_combiners_[curr_row][curr_col];
        }

        //Print(curr_tensor);

    }

    return replaced_double_layer_tensors;

}

Complex Square_Patch_RDM::expect_val_from_replaced_tensors(std::array<IQTensor,2> replaced_tensors, std::array<IQTensor,2> replaced_tensors_dag) const
{
    auto replaced_double_layer_tensors=double_layer_tensors_from_replaced_tensors(replaced_tensors,replaced_tensors_dag);
    return expect_val_from_double_layer_tensors(replaced_double_layer_tensors);
}

Complex Square_Patch_RDM::expect_val_from_double_layer_tensors(const std::vector<std::vector<IQTensor>> &double_layer_tensors) const
{
    IQTensor network_tensor;
    if (patch_dim_[0]>=patch_dim_[1])
    {
        for (int rowi=0; rowi<patch_dim_[0]; rowi++)
        {
            for (int coli=0; coli<patch_dim_[1]; coli++)
            {
                if (rowi==0 && coli==0)
                    network_tensor=double_layer_tensors[0][0];
                else
                    network_tensor*=double_layer_tensors[rowi][coli];
            }
        }
    }
    else
    {
        for (int coli=0; coli<patch_dim_[1]; coli++)
        {
            for (int rowi=0; rowi<patch_dim_[0]; rowi++)
            {
                if (rowi==0 && coli==0)
                    network_tensor=double_layer_tensors[0][0];
                else
                    network_tensor*=double_layer_tensors[rowi][coli];
            }
        }
    }

    //PrintDat(network_tensor);

    return network_tensor.toComplex();
}

IQTensor Square_Patch_RDM::half_patch_tensor_from_double_layer_tensors(const std::vector<std::vector<IQTensor>> &double_layer_tensors, int dir) const
{
    int coli=(patch_dim_[1]-1)*(dir+1)/2;
    IQTensor half_patch_tensor;
    while (coli!=patch_cutting_coords_[(1-dir)/2][1])
    {
        for (int rowi=0; rowi<patch_dim_[0]; rowi++)
        {
            if (!half_patch_tensor.valid())
            {
                half_patch_tensor=double_layer_tensors[rowi][coli];
            }
            else
            {
                half_patch_tensor*=double_layer_tensors[rowi][coli];
            }
        }
        coli=coli-dir;
    }
    return half_patch_tensor;
}

//void Square_Patch_RDM::modify_env_tens(const IQTensor &env_tens)
//{
//    env_tens_=env_tens;
//
//    init_env_tens();
//    modify_boundary_patch_tensors();
//    init_patch_double_layer_tensors();
//    obtain_two_sites_RDM();
//    wf_norm_=std::sqrt(trace(two_sites_RDM_,cutting_phys_legs(0),prime(dag(cutting_phys_legs(0)))).trace(cutting_phys_legs(1),prime(dag(cutting_phys_legs(1)))).toComplex().real());
//}

void Square_Patch_RDM::init_patch_cutting_coords()
{
    for (int rowi=0; rowi<patch_dim_[0]; rowi++)
    {
        for (int coli=0; coli<patch_dim_[1]; coli++)
        {
            if (patch_sites_[rowi][coli]==cutting_sites_[0])
            {
                patch_cutting_coords_[0][0]=rowi;
                patch_cutting_coords_[0][1]=coli;
            }
            if (patch_sites_[rowi][coli]==cutting_sites_[1])
            {
                patch_cutting_coords_[1][0]=rowi;
                patch_cutting_coords_[1][1]=coli;
            }
        }
    }

    //cout << "patch_cutting_coords: " << endl << patch_cutting_coords_[0][0] << " " << patch_cutting_coords_[0][1] << endl << patch_cutting_coords_[1][0] << " " << patch_cutting_coords_[1][1] << endl;

}

void Square_Patch_RDM::init_env_tens()
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

void Square_Patch_RDM::init_patch_tensors()
{
    //init patch_tensors(_dag) with site tensors multiplying right and up bond tensors
    //However, there are no bond tensors on boundary
    for (int rowi=0; rowi<patch_dim_[0]; rowi++)
    {
        std::vector<IQTensor> tensors_one_row, tensors_one_row_dag;
        for (int coli=0; coli<patch_dim_[1]; coli++)
        {
            auto curr_tens=square_peps_.site_tensors(patch_sites_[rowi][coli]);
            //absorb up bond
            if (rowi<patch_dim_[0]-1)
            {
                int ubond=square_peps_.lattice().comm_bond(patch_sites_[rowi][coli],patch_sites_[rowi+1][coli]);
                curr_tens*=square_peps_.bond_tensors(ubond);
            }
            //absorb right bond
            if (coli<patch_dim_[1]-1)
            {
                int rbond=square_peps_.lattice().comm_bond(patch_sites_[rowi][coli],patch_sites_[rowi][coli+1]);
                curr_tens*=square_peps_.bond_tensors(rbond);
            }
            tensors_one_row.push_back(curr_tens);
            if (patch_sites_[rowi][coli]==cutting_sites_[0] || patch_sites_[rowi][coli]==cutting_sites_[1])
            {
                tensors_one_row_dag.push_back(dag(curr_tens).prime());
            }
            else
            {
                tensors_one_row_dag.push_back(dag(curr_tens).prime().noprime(Site));
            }
        }
        patch_tensors_.push_back(tensors_one_row);
        patch_tensors_dag_.push_back(tensors_one_row_dag);
    }

    modify_boundary_patch_tensors();

    //for (int rowi=0; rowi<patch_dim_[0]; rowi++)
    //{
    //    for (int coli=0; coli<patch_dim_[1]; coli++)
    //    {
    //        Print(rowi);
    //        Print(coli);
    //        Print(patch_tensors_[rowi][coli]);
    //        Print(patch_tensors_dag_[rowi][coli]);
    //    }
    //}

}

void Square_Patch_RDM::modify_boundary_patch_tensors()
{
    IQIndex boundary_leg;

    for (int rowi=0; rowi<patch_dim_[0]; rowi++)
    {
        //multiply left leg of first col with env tens
        boundary_leg=square_peps_.virt_legs(patch_sites_[rowi][0],0);
        obtain_env_dressed_tensor(patch_tensors_[rowi][0],env_tens_,boundary_leg);
        patch_tensors_dag_[rowi][0].noprime(prime(dag(boundary_leg)));
        obtain_env_dressed_tensor(patch_tensors_dag_[rowi][0],dag(env_tens_),dag(boundary_leg));
        //multiply right leg of last col with env tens
        boundary_leg=square_peps_.virt_legs(patch_sites_[rowi][patch_dim_[1]-1],2);
        obtain_env_dressed_tensor(patch_tensors_[rowi][patch_dim_[1]-1],env_tens_,boundary_leg);
        patch_tensors_dag_[rowi][patch_dim_[1]-1].noprime(prime(dag(boundary_leg)));
        obtain_env_dressed_tensor(patch_tensors_dag_[rowi][patch_dim_[1]-1],dag(env_tens_),dag(boundary_leg));
    }
    for (int coli=0; coli<patch_dim_[1]; coli++)
    {
        //multiply down leg of first row with env tens
        boundary_leg=square_peps_.virt_legs(patch_sites_[0][coli],3);
        obtain_env_dressed_tensor(patch_tensors_[0][coli],env_tens_,boundary_leg);
        patch_tensors_dag_[0][coli].noprime(prime(dag(boundary_leg)));
        obtain_env_dressed_tensor(patch_tensors_dag_[0][coli],dag(env_tens_),dag(boundary_leg));
        //multiply up leg of last row with env tens
        boundary_leg=square_peps_.virt_legs(patch_sites_[patch_dim_[0]-1][coli],1);
        obtain_env_dressed_tensor(patch_tensors_[patch_dim_[0]-1][coli],env_tens_,boundary_leg);
        patch_tensors_dag_[patch_dim_[0]-1][coli].noprime(prime(dag(boundary_leg)));
        obtain_env_dressed_tensor(patch_tensors_dag_[patch_dim_[0]-1][coli],dag(env_tens_),dag(boundary_leg));
    }

}

void Square_Patch_RDM::init_legs_combiners()
{
    //init left and down virt legs combiners
    for (int rowi=0; rowi<patch_dim_[0]; rowi++)
    {
        std::vector<IQCombiner> left_combiners_one_row, down_combiners_one_row;
        for (int coli=0; coli<patch_dim_[1]; coli++)
        {
            auto left_leg=square_peps_.virt_legs(patch_sites_[rowi][coli],0),
                 down_leg=square_peps_.virt_legs(patch_sites_[rowi][coli],3);
            left_combiners_one_row.push_back(IQCombiner(left_leg,dag(prime(left_leg))));
            down_combiners_one_row.push_back(IQCombiner(down_leg,dag(prime(down_leg))));
        }
        left_legs_combiners_.push_back(left_combiners_one_row);
        down_legs_combiners_.push_back(down_combiners_one_row);
    }

    //init phys legs combiners
    for (int i=0; i<2; i++)
    {
        auto phys_leg=square_peps_.phys_legs(cutting_sites_[i]);
        phys_legs_combiners_.push_back(IQCombiner(phys_leg,prime(dag(phys_leg))));
    }

}

void Square_Patch_RDM::init_patch_double_layer_tensors()
{
    patch_double_layer_tensors_=patch_tensors_;
    
    for (int rowi=0; rowi<patch_dim_[0]; rowi++)
    {
        for (int coli=0; coli<patch_dim_[1]; coli++)
        {
            auto &curr_tensor=patch_double_layer_tensors_[rowi][coli];
            //if no comm indice, the the tensor is at cutting sites in patch bulk. We should treat carefully about indices number overflow
            if (!commonIndex(patch_tensors_[rowi][coli],patch_tensors_dag_[rowi][coli]).valid())
            {
                IQCombiner lr_combiner(square_peps_.virt_legs(patch_sites_[rowi][coli],0),dag(square_peps_.virt_legs(patch_sites_[rowi][coli+1],0)));
                curr_tensor=(patch_tensors_[rowi][coli]*lr_combiner)*(patch_tensors_dag_[rowi][coli]*dag(prime(lr_combiner)));

                curr_tensor=curr_tensor*down_legs_combiners_[rowi][coli]*dag(down_legs_combiners_[rowi+1][coli]);
                curr_tensor=curr_tensor*dag(lr_combiner)*prime(lr_combiner);
                curr_tensor=curr_tensor*left_legs_combiners_[rowi][coli]*dag(left_legs_combiners_[rowi][coli+1]);

                //Print(rowi);
                //Print(coli);
                //Print(curr_tensor);

                continue;
            }

            curr_tensor*=patch_tensors_dag_[rowi][coli];
            //combine top and bottom virt legs
            if (rowi!=0)
                curr_tensor=curr_tensor*down_legs_combiners_[rowi][coli];
            if (rowi!=patch_dim_[0]-1)
                curr_tensor=curr_tensor*dag(down_legs_combiners_[rowi+1][coli]);
            if (coli!=0)
                curr_tensor=curr_tensor*left_legs_combiners_[rowi][coli];
            if (coli!=patch_dim_[1]-1)
                curr_tensor=curr_tensor*dag(left_legs_combiners_[rowi][coli+1]);

            //Print(rowi);
            //Print(coli);
            //Print(curr_tensor);

        }
    }

    //combine phys legs
    for (int i=0; i<2; i++)
    {
        auto &curr_tensor=patch_double_layer_tensors_[patch_cutting_coords_[i][0]][patch_cutting_coords_[i][1]];
        curr_tensor=curr_tensor*phys_legs_combiners_[i];
    }

}


template <class TensorT>
void obtain_env_dressed_tensor(TensorT &dressed_tens, const TensorT &env_tens, const typename TensorT::IndexT &boundary_leg)
{
    //std::vector<typename TensorT::IndexT> env_tens_leg={env_tens.indices()[0],env_tens.indices()[1]};
    //if (env_tens_leg[0].primeLevel()==1) 
    //{
    //    env_tens_leg[0]=env_tens.indices()[1];
    //    env_tens_leg[1]=env_tens.indices()[0];
    //}
    dressed_tens=tensor_contraction<TensorT,TensorT::IndexT>(dressed_tens,env_tens,{boundary_leg},{env_tens.indices()[0]});
    dressed_tens.replaceIndex(env_tens.indices()[1],boundary_leg);
}
template
void obtain_env_dressed_tensor(ITensor &dressed_tens, const ITensor &env_tens, const ITensor::IndexT &boundary_leg);
template
void obtain_env_dressed_tensor(IQTensor &dressed_tens, const IQTensor &env_tens, const IQTensor::IndexT &boundary_leg);

template <class TensorT>
void combine_comm_legs(std::array<TensorT,2> &tensors)
{
    std::vector<typename TensorT::IndexT> comm_legs;
    for (const auto &leg: tensors[0].indices())
    {
        if (hasindex(tensors[1],leg)) 
            comm_legs.push_back(leg);
    }
    if (comm_legs.size()==0)
    {
        cout << "no comm legs!" << endl;
        exit(EXIT_FAILURE);
    }
    if (comm_legs.size()==1) return;
    typename TensorT::CombinerT comm_legs_combiner(comm_legs[0]);
    for (int legi=1; legi<comm_legs.size(); legi++) 
        comm_legs_combiner.addleft(comm_legs[legi]);
    tensors[0]=tensors[0]*comm_legs_combiner;
    tensors[1]=tensors[1]*dag(comm_legs_combiner);
}
template
void combine_comm_legs(std::array<ITensor,2> &tensors);
template
void combine_comm_legs(std::array<IQTensor,2> &tensors);

void spin_square_peps_patch_simple_update(IQPEPS &square_peps, const Evolution_Params &square_su_params, std::vector<std::vector<int>> patch_sites, std::array<int,2> evolved_sites)
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

            //obtain square patch RDM
            Square_Patch_RDM square_RDM(square_peps,env_tens[0][0],patch_sites,evolved_sites);
            //using general patch RDM
            //General_Patch_RDM<IQTensor> square_RDM(patch_name,square_peps,env_tens[0][0],patch_sites,{evolved_sites[0],evolved_sites[1]});
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
            Print(updated_site_tens_ordered_ind.norm());
            Print(square_peps.site_tensors(evolved_sites[0]).norm());
            Print((updated_site_tens_ordered_ind-square_peps.site_tensors(evolved_sites[0])).norm());

            square_peps.generate_site_tensors({updated_site_tens_ordered_ind});

            //stores as PEPS, which is used for further evolution
            if ((step+1)*10%square_su_params.steps_nums[iter]==0)
            {
                std::stringstream ss;

                //zero flux state
                if (std::abs(square_psg::mu_12-1)<EPSILON)
                {
                    ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_rvb_D=" << square_peps.D() << "_Lx=" << square_peps.n_uc()[0] << "_Ly=" << square_peps.n_uc()[1] << "_iter=" << iter << "_step=" << step << "_patch=" << patch_sites.size() << "x" << patch_sites[0].size();
                }
                //pi flux state
                if (std::abs(square_psg::mu_12+1)<EPSILON)
                {
                    ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_pi_rvb_D=" << square_peps.D() << "_Lx=" << square_peps.n_uc()[0] << "_Ly=" << square_peps.n_uc()[1] << "_iter=" << iter << "_step=" << step << "_patch=" << patch_sites.size() << "x" << patch_sites[0].size();
                }

                std::string file_name=ss.str();
                writeToFile(file_name,square_peps);
                
                //reinit leg_gate_params to get rid of local minimal
                //for (auto &param : leg_gate_params) param=rand_gen();
            }

        }//trotter steps

    }//trotter iters
}


bool obtain_spin_sym_leg_gates_params_minimization_from_RDM(const Square_Patch_RDM &square_RDM, const Trotter_Gate &trotter_gate, const std::array<Singlet_Tensor_Basis,2> &leg_gates_basis, std::vector<double> &leg_gate_params, double cutoff)
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
    std::array<IQTensor,2> site_tensors_evolved(square_RDM.cutting_site_tensors());
    IQTensor bond_tensor_evolved(square_RDM.cutting_bond_tensor());

    for (int i=0; i<2; i++)
    {
        site_tensors_evolved[i]*=trotter_gate.site_tensors(i);
        site_tensors_evolved[i].noprime();
    }
    bond_tensor_evolved*=trotter_gate.bond_tensors(0);

    //two method to obtain evolved_wf_norm
    //1. directly use two_sites_RDM
    double evolved_wf_norm=((square_RDM.two_sites_RDM()*dag(trotter_gate.site_tensors(0)*trotter_gate.bond_tensors(0)*trotter_gate.site_tensors(1)).prime()).mapprime(2,1)*(trotter_gate.site_tensors(0)*trotter_gate.bond_tensors(0)*trotter_gate.site_tensors(1))).toComplex().real();
    evolved_wf_norm=sqrt(evolved_wf_norm);
    //Print(evolved_wf_norm);
    //2. replace cutting_tensors with tensors_evolved, and obtain expectation value
    //evolved_wf_norm=sqrt(square_RDM.expect_val_from_replaced_tensors({site_tensors_evolved[0]*bond_tensor_evolved,site_tensors_evolved[1]},{site_tensors_evolved[0]*bond_tensor_evolved,site_tensors_evolved[1]}).real());
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

        std::array<IQTensor,2> left_replaced_tensors={site_tensors_updated_basis[0][base_list[1]]*square_RDM.cutting_bond_tensor(),site_tensors_updated_basis[1][0]},
                               left_replaced_tensors_dag={site_tensors_updated_basis[0][base_list[0]]*square_RDM.cutting_bond_tensor(),site_tensors_updated_basis[1][0]},
                               right_replaced_tensors={site_tensors_updated_basis[0][0]*square_RDM.cutting_bond_tensor(),site_tensors_updated_basis[1][base_list[1]]},
                               right_replaced_tensors_dag={site_tensors_updated_basis[0][0]*square_RDM.cutting_bond_tensor(),site_tensors_updated_basis[1][base_list[0]]};

        auto double_layer_tensors_for_left_half_basis=square_RDM.double_layer_tensors_from_replaced_tensors(left_replaced_tensors,left_replaced_tensors_dag);
        auto double_layer_tensors_for_right_half_basis=square_RDM.double_layer_tensors_from_replaced_tensors(right_replaced_tensors,right_replaced_tensors_dag);

        left_half_basis.push_back(square_RDM.half_patch_tensor_from_double_layer_tensors(double_layer_tensors_for_left_half_basis,-1));
        right_half_basis.push_back(square_RDM.half_patch_tensor_from_double_layer_tensors(double_layer_tensors_for_right_half_basis,1));
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

        //auto M_ijkl=square_RDM.expect_val_from_replaced_tensors({site_tensors_updated_basis[0][base_list[2]]*square_RDM.cutting_bond_tensor(),site_tensors_updated_basis[1][base_list[3]]},{site_tensors_updated_basis[0][base_list[0]]*square_RDM.cutting_bond_tensor(),site_tensors_updated_basis[1][base_list[1]]});
        auto M_ijkl=(left_half_basis[base_list[0]*N_leg_basis+base_list[2]]*right_half_basis[base_list[1]*N_leg_basis+base_list[3]]).toComplex();
        updated_wf_basis_overlap.push_back(M_ijkl.real());

        //if (std::abs(M_ijkl)>0) 
        //{
        //    Print(base_num);
        //    //Print(leg_gates_basis[0].spin_configs(base_list[0]));
        //    //Print(leg_gates_basis[1].spin_configs(base_list[1]));
        //    //Print(leg_gates_basis[0].spin_configs(base_list[2]));
        //    //Print(leg_gates_basis[1].spin_configs(base_list[3]));
        //    Print(M_ijkl);
        //    nonzero_M++;
        //}
    }
    //Print(nonzero_M);
    //Print(total_base_num);
    //Print(updated_wf_basis_overlap);

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
        auto w_ij=square_RDM.expect_val_from_replaced_tensors({site_tensors_evolved[0]*bond_tensor_evolved,site_tensors_evolved[1]},{site_tensors_updated_basis[0][base_list[0]]*square_RDM.cutting_bond_tensor(),site_tensors_updated_basis[1][base_list[1]]});

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
           wf_evolved_wf_overlap=2*square_RDM.expect_val_from_replaced_tensors({site_tensors_evolved[0]*bond_tensor_evolved,site_tensors_evolved[1]},{square_RDM.cutting_site_tensors(0)*square_RDM.cutting_bond_tensor(),square_RDM.cutting_site_tensors(1)}).real(),
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


double heisenberg_energy_from_RDM(const IQTensor &two_sites_RDM)
{
    std::vector<IQIndex> phys_legs;
    for (const auto &indice : two_sites_RDM.indices())
    {
        if (indice.primeLevel()==0) phys_legs.push_back(indice); 
    }

    auto norm_sq=trace(two_sites_RDM,phys_legs[0],prime(dag(phys_legs[0]))).trace(phys_legs[1],prime(dag(phys_legs[1]))).toComplex();
    NN_Heisenberg_Hamiltonian heisenberg_gate({phys_legs[0],phys_legs[1]});
    auto energy=(two_sites_RDM*(heisenberg_gate.site_tensors(0)*heisenberg_gate.bond_tensor()*heisenberg_gate.site_tensors(1))).toComplex();

    //Print(norm_sq);
    //Print(energy);

    return (energy/norm_sq).real();
}

double heisenberg_energy_from_RDM(const Square_Patch_RDM &square_RDM)
{
    NN_Heisenberg_Hamiltonian heisenberg_gate(square_RDM.cutting_phys_legs());
    auto energy=(square_RDM.two_sites_RDM()*(heisenberg_gate.site_tensors(0)*heisenberg_gate.bond_tensor()*heisenberg_gate.site_tensors(1))).toComplex().real();
    return energy/std::pow(square_RDM.wf_norm(),2.);

    //NN_Heisenberg_Hamiltonian heisenberg_gate(square_RDM.cutting_phys_legs());
    //std::array<IQTensor,2> combined_tensors={square_RDM.cutting_site_tensors(0)*square_RDM.cutting_bond_tensor(),square_RDM.cutting_site_tensors(1)};
    //std::array<IQTensor,2> evolved_tensors={(combined_tensors[0]*heisenberg_gate.site_tensors(0)*heisenberg_gate.bond_tensor()).noprime(),(combined_tensors[1]*heisenberg_gate.site_tensors(1)).noprime()};
    //auto norm_sq=square_RDM.expect_val_from_replaced_tensors(combined_tensors,combined_tensors);
    //auto energy=square_RDM.expect_val_from_replaced_tensors(evolved_tensors,combined_tensors);

    //Print(norm_sq);
    //Print(energy);

    //return (energy/norm_sq).real();
}
