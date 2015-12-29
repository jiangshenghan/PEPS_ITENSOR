
#include "simple_update_patch.h"

//class Square_Patch_RDM
Square_Patch_RDM::Square_Patch_RDM(const IQPEPS &square_peps, const IQTensor &env_tens, const std::vector<std::vector<int>> &patch_sites, const std::array<int,2> cutting_sites):
    patch_dim_{patch_sites.size(),patch_sites[0].size()},
    patch_sites_(patch_sites),
    cutting_sites_(cutting_sites),
    env_tens_(env_tens),
    square_peps_(square_peps)
{
    init_patch_cutting_coords();
    init_patch_tensors();
    init_legs_combiners();
    init_patch_double_layer_tensors();
    obtain_two_sites_RDM();
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

Complex Square_Patch_RDM::expect_val_from_replaced_tensors(std::array<IQTensor,2> replaced_tensors, std::array<IQTensor,2> replaced_tensors_dag)
{
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

        //modify to boundary tensors
        for (const auto &boundary_leg : boundary_legs)
        {
            obtain_env_dressed_tensor(replaced_tensors[i],env_tens_,boundary_leg);
            obtain_env_dressed_tensor(replaced_tensors_dag[i],env_tens_,dag(boundary_leg));
        }

        //we already count the bond between two cutting sites in the replaced tensor
        int neigh_cutting_site_dir=(patch_cutting_coords_[1-i][1]-patch_cutting_coords_[i][1]+1)+(-patch_cutting_coords_[1-i][0]+patch_cutting_coords_[i][0]+2);
        combine_dir[neigh_cutting_site_dir]=false;

        //absorb bond legs and combine top and bottom legs
        auto &curr_row=patch_cutting_coords_[i][0],
             &curr_col=patch_cutting_coords_[i][1];
        auto &curr_tensor=replaced_double_layer_tensors[curr_row][curr_col];
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
    }

    //calculate expectation value using replaced_double_layer_tensors
    IQTensor network_tensor;
    if (patch_dim_[0]>=patch_dim_[1])
    {
        for (int rowi=0; rowi<patch_dim_[0]; rowi++)
        {
            for (int coli=0; coli<patch_dim_[1]; coli++)
            {
                if (rowi==0 && coli==0)
                    network_tensor=replaced_double_layer_tensors[0][0];
                else
                    network_tensor*=replaced_double_layer_tensors[rowi][coli];
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
                    network_tensor=replaced_double_layer_tensors[0][0];
                else
                    network_tensor*=replaced_double_layer_tensors[rowi][coli];
            }
        }
    }

    return network_tensor.toComplex();

}

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


void obtain_env_dressed_tensor(IQTensor &dressed_tens, const IQTensor &env_tens, const IQIndex &boundary_leg)
{
    std::vector<IQIndex> env_tens_leg={env_tens.indices()[0],env_tens.indices()[1]};
    if (env_tens_leg[0].primeLevel()==1) 
    {
        env_tens_leg[0]=env_tens.indices()[1];
        env_tens_leg[1]=env_tens.indices()[0];
    }
    dressed_tens=tensor_contraction<IQTensor,IQIndex>(dressed_tens,env_tens,{boundary_leg},{env_tens_leg[0]});
    dressed_tens.replaceIndex(env_tens_leg[1],boundary_leg);
}


/*
void spin_square_peps_patch_simple_update(IQPEPS &square_peps, const Evolution_Params &square_su_params, std::vector<std::vector<int>> patch_sites, std::array<int,2> evolved_sites)
{
    //Initialize trotter gate
    std::array<IQIndex,2> gate_legs={square_peps.phys_leg(evolved_sites[0]),square_peps.phys_leg(evolved_sites[1])};
    NN_Heisenberg_Trotter_Gate evolve_gate(gate_legs);

    //Initialize env_tens
    //env_tens[i] stores environment of site i, which is direct product of matrices for the legs of site i 
    //there is no env matrix between A and B
    std::array<std::vector<IQTensor>,2> env_tens;

    int comm_bond=square_peps.lattice().comm_bond(evolved_sites[0],evolved_sites[1]);
    auto comm_bond_tensor=square_peps.bond_tensors(comm_bond);

    //leg gates is used to approx evolve_gate, which formed by three indices, two in one out (primed). 
    //two in legs are contracting to the site tensor which has applied by trotter gate
    //one out leg is to contract comm bond tensor
    std::array<IndexSet<IQIndex>,2> leg_gates_indices;
    std::array<Singlet_Tensor_Basis,2> leg_gates_basis;
    for (int i=0; i<2; i++)
    {
        leg_gates_indices[i].addindex(commonIndex(dag(square_peps.site_tensors(evolved_sites[i])),comm_bond_tensor));
        leg_gates_indices[i].addindex(commonIndex(dag(evolve_gate.site_tensors(evolved_sites[i])),evolve_gate.bond_tensors(0)));
        leg_gates_indices[i].addindex(commonIndex(square_peps.site_tensors(evolved_sites[i]),dag(comm_bond_tensor)).prime());

        leg_gates_basis[i]=Singlet_Tensor_Basis(leg_gates_indices[i]);

        //PrintDat(leg_gates_basis[i]);
    }

    //leg_gates_for_one_site is used to update site_tensors[evolved_sites[0]]
    std::vector<IQTensor> leg_gates_for_one_site;
    auto indice_from_evolve_gate=commonIndex(dag(evolve_gate.site_tensors(evolved_sites[0])),evolve_gate.bond_tensors(0));
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

            //obtain two sites RDM
            IQTensor two_sites_RDM=square_peps_two_sites_RDM_simple_update(square_peps,env_tens[0][0],patch_sites,evolved_sites);

            //check wf_distance_f, wf_distance_df
            //wf_distance_func_check(site_env_tens,comm_bond_tensor,evolve_gate,leg_gates_basis,leg_gate_params);

            //measure energy by site_env_tens
            Print(heisenberg_energy_from_RDM(two_sites_RDM));

            //if we cannot obtain a reasonable leg_gate, we try smaller time separation
            double cutoff=square_su_params.ts[iter]/10.;
            if (cutoff>1e-5) cutoff=1e-5;
            //double cutoff=1e-5;
            if (!obtain_spin_sym_leg_gates_params_minimization_from_RDM(two_sites_RDM,evolve_gate,leg_gates_basis,leg_gate_params,cutoff)) break;

            //using leg_gate_params to generate all leg gates
            auto leg_gate_sample=singlet_tensor_from_basis_params(leg_gates_basis[0],leg_gate_params);

            for (auto &gate : leg_gates_for_one_site) 
            {
                tensor_assignment(gate,leg_gate_sample); 
                //Print(gate.indices());
            }


            //updated site_tensors[0]
            auto updated_site_tens=square_peps.site_tensors(evolved_sites[0]);

            for (const auto &leg_gate : leg_gates_for_one_site)
            {
                updated_site_tens*=evolve_gate.site_tensors(0)*leg_gate;
                updated_site_tens.noprime();
                //PrintDat(leg_gate);
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
                //ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_rvb_D=" << square_peps.D() << "_Lx=" << square_peps.n_uc()[0] << "_Ly=" << square_peps.n_uc()[1] << "_iter=" << iter << "_step=" << step;
                //pi flux state
                ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_pi_rvb_D=" << square_peps.D() << "_Lx=" << square_peps.n_uc()[0] << "_Ly=" << square_peps.n_uc()[1] << "_iter=" << iter << "_step=" << step;

                std::string file_name=ss.str();
                writeToFile(file_name,square_peps);
                
                //reinit leg_gate_params to get rid of local minimal
                //for (auto &param : leg_gate_params) param=rand_gen();
            }


        }//trotter steps

    }//trotter iters
}


bool obtain_spin_sym_leg_gates_params_minimization_from_RDM(const IQTensor &two_sites_RDM, const Trotter_Gate &trotter_gate, const std::array<Singlet_Tensor_Basis,2> &leg_gates_basis, std::vector<double> &leg_gate_params, double cutoff=1E-5)
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

    int N_basis=leg_gate_params.size();

    //prepare evolved/updated site tensors, which is used to obtain minimization params
    std::array<IQTensor,2> site_tensors_evolved();

    //get time evloved tensors: evolved_wf_norm
    //first method: directly using two-sites RDM
    //get combined trotter gate for bra/ket physical indices
    //we set the out physical leg of trotter gate with primeLevel=2
    IQTensor combined_trotter_gate_bra=trotter_gate.site_tensors(0)*trotter_gate.bond_tensors(0)*trotter_gate.site_tensors(1),
             combined_trotter_gate_ket=combined_trotter_gate_bra;
    combined_trotter_gate_ket.mapprime(1,2);
    combined_trotter_gate_bra.dag().prime();
    double evolved_wf_norm=two_sites_RDM*combined_trotter_gate_bra*combined_trotter_gate_ket;

    //second method, using RDM variation

    //get M_{ijkl}=\langle\phi_{ij}|\phi_{kl}\rangle, store in a vector updated_wf_basis_overlap
    int total_base_num=1;
    std::vector<int> max_base_list;
    for (int i=0; i<4; i++)
    {
        total_base_num*=N_leg_basis;
        max_base_list.push_back(N_leg_basis);
    }
    std::vector<double> updated_wf_basis_overlap;
    for (int base_num=0; base_num<total_base_num; base_num++)
    {
        auto base_list=list_from_num(base_num,max_base_list);
        auto M_ijkl=two_sites_RDM*;
        updated_wf_overlap.push_back(M_ijkl.toComplex().real());
    }

    //get w_ij=\langle\phi_{ij}|\psi\rangle+\langle\psi|\phi_{ij}\rangle


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
    double wf_distance=std::sqrt(wf_distance_f(s->x,wf_distance_params))/wf_distance_params->evolved_wf_norm;
    Print(wf_distance);

    //if the leg_gate is not a good approx of trotter gate, we may be trapped in a local minima, thus, we retry to find leg gate
    //if (wf_distance>cutoff*30)
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
*/


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

    Print(norm_sq);
    Print(energy);

    return (energy/norm_sq).real();
}
