
#include "simple_update.h"


void spin_square_peps_simple_update(IQPEPS &square_peps, const Evolution_Params &square_su_params)
{
    //Initialize trotter gate
    std::array<IQIndex,2> site01_legs{{square_peps.phys_legs(0),square_peps.phys_legs(1)}};
    NN_Heisenberg_Hamiltonian hamiltonian_gate(site01_legs);
    NN_Heisenberg_Trotter_Gate evolve_gate(site01_legs);

    //Initialize env_tens
    //env_tens[i] stores environment of site i, which is direct product of matrices for the legs of site i 
    //there is no env matrix between A and B
    std::array<std::vector<IQTensor>,2> env_tens;

    int comm_bond=square_peps.lattice().comm_bond(0,1);
    auto comm_bond_tensor=square_peps.bond_tensors(comm_bond);

    //leg gates is used to approx evolve_gate, which formed by three indices, two in one out (primed). 
    //two in legs are contracting to the site tensor which has applied by trotter gate
    //one out leg is to contract comm bond tensor
    std::array<IndexSet<IQIndex>,2> leg_gates_indices;
    std::array<Singlet_Tensor_Basis,2> leg_gates_basis;
    for (int i=0; i<2; i++)
    {
        leg_gates_indices[i].addindex(commonIndex(dag(square_peps.site_tensors(i)),comm_bond_tensor));
        leg_gates_indices[i].addindex(commonIndex(dag(evolve_gate.site_tensors(i)),evolve_gate.bond_tensors(0)));
        leg_gates_indices[i].addindex(commonIndex(square_peps.site_tensors(i),dag(comm_bond_tensor)).prime());

        leg_gates_basis[i]=Singlet_Tensor_Basis(leg_gates_indices[i]);

        //PrintDat(leg_gates_basis[i]);
    }

    //leg_gates_for_site0 is used to update site_tensors[0]
    std::vector<IQTensor> leg_gates_site0;
    auto indice_from_evolve_gate=commonIndex(dag(evolve_gate.site_tensors(0)),evolve_gate.bond_tensors(0));
    for (const auto &leg_indice : square_peps.site_tensors(0).indices())
    {
        if (leg_indice.type()==Site) continue;

        std::vector<IQIndex> leg_gates_site0_indices{dag(leg_indice),indice_from_evolve_gate,prime(leg_indice)};

        leg_gates_site0.push_back(IQTensor(leg_gates_site0_indices));
    }
    //Print(leg_gates_site0);

    //leg_gate_params stores tunable parameters for the leg gate.
    std::vector<double> leg_gate_params(leg_gates_basis[0].dim());
    for (auto &param : leg_gate_params) param=rand_gen();
    Print(leg_gate_params);
    
    for (int iter=0; iter<square_su_params.iter_nums; iter++)
    {
        evolve_gate.change_time(square_su_params.ts[iter]);

        for (int step=0; step<square_su_params.steps_nums[iter]; step++)
        {
            Print(iter);
            Print(step);
            Print(square_su_params.ts[iter]);

            //get_env_tensor_iterative(square_peps.site_tensors(0)*comm_bond_tensor,square_peps.site_tensors(1),env_tens);
            //get_env_tensor_minimization(square_peps.site_tensors(0)*comm_bond_tensor,square_peps.site_tensors(1),env_tens);
            auto combined_site_tens0=square_peps.site_tensors(0);
            for (int neighi=0; neighi<square_peps.n_bonds_to_one_site(); neighi++)
            {
                int bondi=square_peps.lattice().site_neighbour_bonds(0,neighi);
                if (bondi==comm_bond) continue;
                combined_site_tens0*=square_peps.bond_tensors(bondi);
            }
            get_env_tensor_minimization(combined_site_tens0*comm_bond_tensor,square_peps.site_tensors(1),env_tens);

            std::array<IQTensor,2> site_env_tens{{combined_site_tens0,square_peps.site_tensors(1)}};
            
            for (int sitei=0; sitei<2; sitei++)
            {
                for (const auto &env_leg_tensor : env_tens[sitei])
                {
                    site_env_tens[sitei]*=env_leg_tensor;
                }
                site_env_tens[sitei].noprime();
            }

            //check wf_distance_f, wf_distance_df
            //wf_distance_func_check(site_env_tens,comm_bond_tensor,evolve_gate,leg_gates_basis,leg_gate_params);

            //measure energy by site_env_tens
            Print(heisenberg_energy_from_site_env_tensors(site_env_tens,comm_bond_tensor,hamiltonian_gate));

            //if we cannot obtain a reasonable leg_gate, we try smaller time separation
            if (!obtain_spin_sym_leg_gates_params_minimization(site_env_tens,comm_bond_tensor,evolve_gate,leg_gates_basis,leg_gate_params)) break;

            //using leg_gate_params to generate all leg gates
            auto leg_gate_sample=singlet_tensor_from_basis_params(leg_gates_basis[0],leg_gate_params);

            for (auto &gate : leg_gates_site0) 
            {
                tensor_assignment(gate,leg_gate_sample); 
                //Print(gate.indices());
            }


            //updated site_tensors[0]
            auto site_tens0=square_peps.site_tensors(0);

            int legi=0;
            for (const auto &leg_gate : leg_gates_site0)
            {
                site_tens0*=evolve_gate.site_tensors(0)*leg_gate;
                site_tens0.noprime();
                legi++;
                //PrintDat(leg_gate);
            }
            //site_tens0*=square_peps.site_tensors(0).norm()/site_tens0.norm();
            //we should never change order of indices of site tensor
            auto site_tens0_ordered_ind=square_peps.site_tensors(0);
            tensor_assignment_diff_order(site_tens0_ordered_ind,site_tens0);

            //check symmetry of site_tens0
            //auto sym_site_tens0=site_tens0_ordered_ind;
            //rotation_symmetrize_square_rvb_site_tensor(sym_site_tens0);
            //sym_site_tens0*=0.25;
            //Print(site_tens0_ordered_ind.norm());
            //Print(sym_site_tens0.norm());
            //Print((site_tens0_ordered_ind-sym_site_tens0).norm());
            //Print(diff_tensor_and_singlet_projection(site_tens0_ordered_ind));
            //Print(diff_tensor_and_singlet_projection(sym_site_tens0));

            //symmetrize site_tens0_ordered_ind, and using site_tens0_ordered_ind to generate all site tensors of peps
            rotation_symmetrize_square_rvb_site_tensor(site_tens0_ordered_ind);
            //site_tens0_ordered_ind*=0.25;
            //we always keep the same norm as original state
            site_tens0_ordered_ind*=(square_peps.site_tensors(0).norm())/(site_tens0_ordered_ind.norm());
            site_tens0_ordered_ind.clean();
            Print(site_tens0_ordered_ind.norm());
            Print(square_peps.site_tensors(0).norm());
            Print((site_tens0_ordered_ind-square_peps.site_tensors(0)).norm());

            square_peps.generate_site_tensors({site_tens0_ordered_ind});

            //store as Tnetwork Storage, which is used for energy measurement
            //if (step*10%square_su_params.steps_nums[iter]==0)
            //{
            //    std::stringstream ss;
            //    if (std::abs(square_psg::mu_12-1)<EPSILON) //zero-flux
            //    {
            //        ss << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_rvb_D=" << square_peps.D() << "_Lx=" << square_peps.n_uc()[0] << "_Ly=" << square_peps.n_uc()[1] << "_iter=" << iter << "_step=" << step  << ".txt";
            //    }
            //    if (std::abs(square_psg::mu_12+1)<EPSILON) //pi-flux
            //    {
            //        ss << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_pi_rvb_D=" << square_peps.D() << "_Lx=" << square_peps.n_uc()[0] << "_Ly=" << square_peps.n_uc()[1] << "_iter=" << iter << "_step=" << step  << ".txt";
            //    }
            //    std::string file_name=ss.str();
            //    Tnetwork_Storage<IQTensor> square_rvb_storage=peps_to_tnetwork_storage(square_peps);
            //    writeToFile(file_name,square_rvb_storage);
            //}

        }//trotter steps

        //stores as PEPS, which is used for further evolution
        //std::stringstream ss;
        //if (std::abs(square_psg::mu_12-1)<EPSILON) //zero-flux state
        //{
        //    ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_rvb_D=" << square_peps.D() << "_Lx=" << square_peps.n_uc()[0] << "_Ly=" << square_peps.n_uc()[1] << "_iter=" << iter << ".txt";
        //}
        //if (std::abs(square_psg::mu_12+1)<EPSILON) //pi-flux state
        //{
        //    ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_pi_rvb_D=" << square_peps.D() << "_Lx=" << square_peps.n_uc()[0] << "_Ly=" << square_peps.n_uc()[1] << "_iter=" << iter << ".txt";
        //}
        //std::string file_name=ss.str();
        //writeToFile(file_name,square_peps);

    }//finish all simple update
}


void obtain_spin_sym_leg_gates_params_iterative(const std::array<IQTensor,2> &site_tensors, const IQTensor &bond_tensor, const Trotter_Gate &evolve_gate, const std::array<Singlet_Tensor_Basis,2> &leg_gates_basis, std::vector<double> &leg_gate_params)
{
    //get time evloved tensors
    std::array<IQTensor,2> site_tensors_evolved(site_tensors);
    IQTensor bond_tensor_evolved(bond_tensor);

    for (int i=0; i<2; i++)
    {
        site_tensors_evolved[i]*=evolve_gate.site_tensors(i);
        site_tensors_evolved[i].noprime();
        site_tensors_evolved[i].clean();
    }
    bond_tensor_evolved*=evolve_gate.bond_tensors(0);
    bond_tensor_evolved.clean();

    
    //Print(diff_tensor_and_singlet_projection(site_tensors_evolved[0]));
    //Print(diff_tensor_and_singlet_projection(site_tensors_evolved[1]));
    //Print(diff_tensor_and_singlet_projection(bond_tensor_evolved));

    double evolved_tensor_norm=(site_tensors_evolved[0]*bond_tensor_evolved*site_tensors_evolved[1]).norm();
    //Print(evolved_tensor_norm);


    //site0_updated_basis[i] stores result of multiplication of site_tensors[0] with each linear basis
    //site0_updated_basis_dag stores dag of site0_updated_basis where the virt_leg connecting bond tensor is primed
    std::vector<IQTensor> site0_updated_basis,site0_updated_basis_dag;
    //we always set the inner leg of ket no prime and inner leg of bra primed
    IQTensor bond_tensor_dag=prime(dag(bond_tensor));
    for (const auto &tensor_base : leg_gates_basis[0])
    {
        IQTensor temp_tens=site_tensors_evolved[0]*tensor_base;

        site0_updated_basis.push_back(noprime(temp_tens));
        site0_updated_basis_dag.push_back(dag(temp_tens));
        //Print(site0_updated_basis.back());
        //Print(diff_tensor_and_singlet_projection(site0_updated_basis.back()));
    }


    //fix another gate (gate v), and get updated site1
    int N_leg_basis=leg_gates_basis[0].dim();
    arma::Col<Complex> vec_params(leg_gate_params.size());
    for (int i=0; i<leg_gate_params.size(); i++)
        vec_params(i)=leg_gate_params[i];
    //matN_{ij}=\langle T_i|T_j\rangle
    //vecb_i=\langle T_i|\psi_0\rangle
    //where T_i=site0_updated_basis[i]*bond_tensor*site1_updated
    //\psi_0 is the wavefunction after time evolving
    //Then we get 
    //(abs(|\psi\rangle-|\psi_0\rangle))^2=\sum_{i,j}p_ip_j.matN_{ij}-\sum_ip_i.(vecb_i+vecb_i^*)+\langle\psi_0|\psi_0\rangle
    arma::Mat<Complex> matN(N_leg_basis,N_leg_basis);
    arma::Col<Complex> vecb(N_leg_basis);

    int iter=0, max_iter=30;
    while (iter<max_iter)
    {
        IQTensor gate_v=singlet_tensor_from_basis_params(leg_gates_basis[1],arma::real(vec_params));
        IQTensor site1_updated=site_tensors_evolved[1]*gate_v;
        IQTensor site1_updated_dag=dag(site1_updated);
        site1_updated.noprime();


        for (int i=0; i<N_leg_basis; i++)
        {
            auto b_i=(site_tensors_evolved[0]*bond_tensor_evolved)*(site0_updated_basis_dag[i]*bond_tensor_dag)*(site_tensors_evolved[1]*site1_updated_dag);
            vecb(i)=b_i.toComplex();
            for (int j=0; j<N_leg_basis; j++)
            {
                auto N_ij=(site0_updated_basis_dag[i]*bond_tensor_dag)*(site0_updated_basis[j]*bond_tensor)*(site1_updated_dag*site1_updated);
                matN(i,j)=N_ij.toComplex();
            }
        }

        auto updated_tensor_norm=std::sqrt(arma::dot(vec_params,matN*vec_params).real());
        //vec_params=vec_params/updated_tensor_norm*evolved_tensor_norm;
        //updated_tensor_norm=std::sqrt(arma::dot(vec_params,matN*vec_params).real());
        auto wf_overlap=arma::dot(vec_params,vecb);
        auto diff_sq=updated_tensor_norm*updated_tensor_norm+evolved_tensor_norm*evolved_tensor_norm-wf_overlap-std::conj(wf_overlap);
        //auto diff_sq=arma::dot(vec_params,matN*vec_params)-arma::dot(vec_params,vecb)-arma::dot(arma::conj(vecb),vec_params)+evolved_tensor_norm*evolved_tensor_norm;

        vec_params=arma::inv(matN)*vecb;

        //Print(iter);
        //Print(arma::norm(arma::imag(matN),2));
        //Print(arma::real(matN));
        //Print(arma::real(vecb).t());
        //Print(arma::real(vec_params).t());
        //Print(evolved_tensor_norm);
        //Print(updated_tensor_norm);
        Print(wf_overlap/(evolved_tensor_norm*updated_tensor_norm));
        Print(diff_sq);

        iter++;
    }

    for (int i=0; i<N_leg_basis; i++)
    {
        leg_gate_params[i]=vec_params(i).real();
    }

}


bool obtain_spin_sym_leg_gates_params_minimization(const std::array<IQTensor,2> &site_tensors, const IQTensor &bond_tensor, const Trotter_Gate &trotter_gate, const std::array<Singlet_Tensor_Basis,2> &leg_gates_basis, std::vector<double> &leg_gate_params)
{
    //init leg_gate_params
    if (leg_gate_params.empty())
    {
        for (int i=0; i<leg_gates_basis[0].dim(); i++)
            leg_gate_params.push_back(rand_gen());
        Print(leg_gate_params);
    }

    //get time evloved tensors
    std::array<IQTensor,2> site_tensors_evolved(site_tensors);
    IQTensor bond_tensor_evolved(bond_tensor);

    for (int i=0; i<2; i++)
    {
        site_tensors_evolved[i]*=trotter_gate.site_tensors(i);
        site_tensors_evolved[i].noprime();
        site_tensors_evolved[i].clean();
    }
    bond_tensor_evolved*=trotter_gate.bond_tensors(0);
    bond_tensor_evolved.clean();

    double evolved_wf_norm=(site_tensors_evolved[0]*bond_tensor_evolved*site_tensors_evolved[1]).norm();

    //Print(site_tensors_evolved[0]);
    //Print(site_tensors_evolved[1]);
    //Print(bond_tensor_evolved);

    //site_tensors_updated_basis stores result of applying leg_gates_basis to site_tensors
    //we always set the inner leg of ket no prime and inner leg of bra primed
    std::array<std::vector<IQTensor>,2> site_tensors_updated_basis, site_tensors_updated_basis_dag;
    for (int i=0; i<2; i++)
    {
        for (const auto &leg_base : leg_gates_basis[i])
        {
            IQTensor temp_tens=site_tensors_evolved[i]*leg_base;
            site_tensors_updated_basis[i].push_back(noprime(temp_tens));
            site_tensors_updated_basis_dag[i].push_back(dag(temp_tens));
        }
        //Print(site_tensors_updated_basis[i]);
    }
    IQTensor bond_tensor_dag=prime(dag(bond_tensor));

    int N_leg_basis=leg_gate_params.size();

    //get M_{ijkl}=\langle\phi_{ij}|\phi_{kl}\rangle, store in a vector
    int total_base_num=1;
    std::vector<int> max_base_list;
    std::vector<double> updated_wf_basis_overlap;
    for (int i=0; i<4; i++)
    {
        total_base_num*=N_leg_basis;
        max_base_list.push_back(N_leg_basis);
    }
    for (int base_num=0; base_num<total_base_num; base_num++)
    {
        auto base_list=list_from_num(base_num,max_base_list);
        auto M_ijkl=(site_tensors_updated_basis_dag[0][base_list[0]]*site_tensors_updated_basis[0][base_list[2]])*bond_tensor_dag*bond_tensor*(site_tensors_updated_basis_dag[1][base_list[1]]*site_tensors_updated_basis[1][base_list[3]]);
        updated_wf_basis_overlap.push_back(M_ijkl.toComplex().real());
        //PrintDat(M_ijkl);
    }
    //Print(updated_wf_basis_overlap);

    //get w_ij=\langle\phi_{ij}|\psi\rangle+\langle\psi|\phi_{ij}\rangle
    std::vector<double> updated_wf_basis_evolved_wf_overlap;
    for (int base_num=0; base_num<N_leg_basis*N_leg_basis; base_num++)
    {
        std::array<int,2> base_list={base_num/N_leg_basis,base_num%N_leg_basis};
        auto w_ij=(site_tensors_updated_basis_dag[0][base_list[0]]*site_tensors_evolved[0])*bond_tensor_dag*bond_tensor_evolved*(site_tensors_updated_basis_dag[1][base_list[1]]*site_tensors_evolved[1]);
        w_ij+=dag(w_ij);
        updated_wf_basis_evolved_wf_overlap.push_back(w_ij.toComplex().real());
    }
    //Print(updated_wf_basis_evolved_wf_overlap);

    //using conjugate gradient minimization to minimize distance square between updated wf and time evolved wf
    int find_min_status;
    int iter=0, max_iter=100;

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
        find_min_status=gsl_multimin_test_gradient(s->gradient,0.1);

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
    if (wf_distance>0.1)
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


double wf_distance_f(const gsl_vector *x, void *params)
{
    std::vector<double> leg_gate_params;
    Wf_Distance_Params *wf_distance_params=(Wf_Distance_Params *)params;
    int N_leg_basis=wf_distance_params->N_leg_basis;

    //Print(N_leg_basis);
    //Print(wf_distance_params->M);
    //Print(wf_distance_params->w);

    for (int i=0; i<N_leg_basis; i++)
        leg_gate_params.push_back(gsl_vector_get(x,i));

    //distance_sq=x_i.x_j.x_k.x_l.M_{ijkl}-x_i.x_j.w_{ij}+evolved_wf_norm^2
    //=updated_wf_norm^2+evolved_wf_norm^2-wf_overlap
    int total_base_num=1;
    std::vector<int> max_base_list;
    for (int i=0; i<4; i++)
    {
        total_base_num*=N_leg_basis;
        max_base_list.push_back(N_leg_basis);
    }
    double updated_wf_norm=0;
    for (int base_num=0; base_num<total_base_num; base_num++)
    {
        auto base_list=list_from_num(base_num,max_base_list);
        updated_wf_norm+=leg_gate_params[base_list[0]]*leg_gate_params[base_list[1]]*leg_gate_params[base_list[2]]*leg_gate_params[base_list[3]]*wf_distance_params->M[base_num];
        //distance_sq+=leg_gate_params[base_list[0]]*leg_gate_params[base_list[1]]*leg_gate_params[base_list[2]]*leg_gate_params[base_list[3]]*wf_distance_params->M[base_num];
    }
    updated_wf_norm=std::sqrt(updated_wf_norm);
    
    double updated_evolved_wf_overlap=0;
    for (int base_num=0; base_num<N_leg_basis*N_leg_basis; base_num++)
    {
        std::array<int,2> base_list={base_num/N_leg_basis,base_num%N_leg_basis};
        updated_evolved_wf_overlap+=leg_gate_params[base_list[0]]*leg_gate_params[base_list[1]]*wf_distance_params->w[base_num];
        //distance_sq-=leg_gate_params[base_list[0]]*leg_gate_params[base_list[1]]*wf_distance_params->w[base_num];
    }

    double distance_sq=std::pow(updated_wf_norm,2)+std::pow(wf_distance_params->evolved_wf_norm,2)-updated_evolved_wf_overlap;

    //Print(wf_distance_params->evolved_wf_norm);
    //Print(updated_wf_norm);
    //Print(updated_evolved_wf_overlap);
    //Print(updated_evolved_wf_overlap/(2.*updated_wf_norm*wf_distance_params->evolved_wf_norm));
    //Print(std::sqrt(distance_sq));

    return distance_sq;
}

void wf_distance_df(const gsl_vector *x, void *params, gsl_vector *df)
{
    std::vector<double> leg_gate_params;
    Wf_Distance_Params *wf_distance_params=(Wf_Distance_Params *)params;
    int N_leg_basis=wf_distance_params->N_leg_basis;

    for (int i=0; i<N_leg_basis; i++)
        leg_gate_params.push_back(gsl_vector_get(x,i));

    //distance_sq_df[i]=\sum_{jkl}c_j*c_k*c_l*(M_{ijkl}+M_{jikl}+M_{jkil}+M_{jkli})-\sum_jc_j*(w_ij+w_ji)
    std::vector<double> distance_sq_df(N_leg_basis,0);
    int total_base_num=1;
    std::vector<int> max_base_list;
    for (int i=0; i<4; i++)
    {
        total_base_num*=N_leg_basis;
        max_base_list.push_back(N_leg_basis);
    }
    for (int base_num=0; base_num<total_base_num; base_num++)
    {
        auto base_list=list_from_num(base_num,max_base_list);
        auto base_list_rotate=base_list;
        auto base_num_rotate=base_num;
        double M_rotate_sum=0;
        for (int rotatei=0; rotatei<4; rotatei++)
        {
            M_rotate_sum+=wf_distance_params->M[base_num_rotate];
            std::rotate(base_list_rotate.begin(),base_list_rotate.begin()+1,base_list_rotate.end());
            base_num_rotate=num_from_list(base_list_rotate,max_base_list);
        }
        distance_sq_df[base_list[0]]+=leg_gate_params[base_list[1]]*leg_gate_params[base_list[2]]*leg_gate_params[base_list[3]]*M_rotate_sum;
    }
    for (int base_num=0; base_num<N_leg_basis*N_leg_basis; base_num++)
    {
        std::array<int,2> base_list={base_num/N_leg_basis,base_num%N_leg_basis};
        int base_num_rotate=base_list[0]+base_list[1]*N_leg_basis;
        distance_sq_df[base_list[0]]-=leg_gate_params[base_list[1]]*(wf_distance_params->w[base_num]+wf_distance_params->w[base_num_rotate]);
    }

    for (int i=0; i<N_leg_basis; i++)
        gsl_vector_set(df,i,distance_sq_df[i]);

}

void wf_distance_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
    *f=wf_distance_f(x,params);
    wf_distance_df(x,params,df);
}

void wf_distance_func_check(const std::array<IQTensor,2> &site_tensors, const IQTensor &bond_tensor, const Trotter_Gate &trotter_gate, const std::array<Singlet_Tensor_Basis,2> &leg_gates_basis, std::vector<double> &leg_gate_params)
{
    //get evolved_wf information
    std::array<IQTensor,2> site_tensors_evolved(site_tensors);
    IQTensor bond_tensor_evolved(bond_tensor);
    for (int i=0; i<2; i++)
    {
        site_tensors_evolved[i]*=trotter_gate.site_tensors(i);
        site_tensors_evolved[i].noprime();
        site_tensors_evolved[i].clean();
    }
    bond_tensor_evolved*=trotter_gate.bond_tensors(0);
    bond_tensor_evolved.clean();
    double evolved_wf_norm=(site_tensors_evolved[0]*bond_tensor_evolved*site_tensors_evolved[1]).norm();

    //get updated_wf information
    std::array<IQTensor,2> leg_gates={singlet_tensor_from_basis_params(leg_gates_basis[0],leg_gate_params),singlet_tensor_from_basis_params(leg_gates_basis[1],leg_gate_params)};
    std::array<IQTensor,2> site_tensors_updated={(site_tensors_evolved[0]*leg_gates[0]).noprime(),(site_tensors_evolved[1]*leg_gates[1]).noprime()};
    std::array<IQTensor,2> site_tensors_updated_dag={(site_tensors_evolved[0]*leg_gates[0]).dag(),(site_tensors_evolved[1]*leg_gates[1]).dag()};
    IQTensor bond_tensor_dag=dag(bond_tensor).prime();
    double updated_wf_norm=(site_tensors_updated[0]*bond_tensor*site_tensors_updated[1]).norm();

    //get overlap between updated_wf and evolved_wf: langle\phi|\psi\rangle+\langle\psi|\phi|ranlge
    double wf_overlap=2*((site_tensors_evolved[0]*site_tensors_updated_dag[0])*bond_tensor_evolved*bond_tensor_dag*(site_tensors_evolved[1]*site_tensors_updated_dag[1])).toComplex().real();


    //get information from wf_distance_f
    int N_leg_basis=leg_gate_params.size();
    //get M_{ijkl}=\langle\phi_{ij}|\phi_{kl}\rangle, store in a vector
    int total_base_num=1;
    std::vector<int> max_base_list;
    std::vector<double> updated_wf_basis_overlap;
    for (int i=0; i<4; i++)
    {
        total_base_num*=N_leg_basis;
        max_base_list.push_back(N_leg_basis);
    }
    //site_tensors_updated_basis stores result of applying leg_gates_basis to site_tensors
    //we always set the inner leg of ket no prime and inner leg of bra primed
    std::array<std::vector<IQTensor>,2> site_tensors_updated_basis, site_tensors_updated_basis_dag;
    for (int i=0; i<2; i++)
    {
        for (const auto &leg_base : leg_gates_basis[i])
        {
            IQTensor temp_tens=site_tensors_evolved[i]*leg_base;
            site_tensors_updated_basis[i].push_back(noprime(temp_tens));
            site_tensors_updated_basis_dag[i].push_back(dag(temp_tens));
        }
        //Print(site_tensors_updated_basis[i]);
    }
    for (int base_num=0; base_num<total_base_num; base_num++)
    {
        auto base_list=list_from_num(base_num,max_base_list);
        auto M_ijkl=(site_tensors_updated_basis_dag[0][base_list[0]]*site_tensors_updated_basis[0][base_list[2]])*bond_tensor_dag*bond_tensor*(site_tensors_updated_basis_dag[1][base_list[1]]*site_tensors_updated_basis[1][base_list[3]]);
        updated_wf_basis_overlap.push_back(M_ijkl.toComplex().real());
    }
    //get w_ij=\langle\phi_{ij}|\psi\rangle+\langle\psi|\phi_{ij}\rangle
    std::vector<double> updated_wf_basis_evolved_wf_overlap;
    for (int base_num=0; base_num<N_leg_basis*N_leg_basis; base_num++)
    {
        std::array<int,2> base_list={base_num/N_leg_basis,base_num%N_leg_basis};
        auto w_ij=(site_tensors_updated_basis_dag[0][base_list[0]]*site_tensors_evolved[0])*bond_tensor_dag*bond_tensor_evolved*(site_tensors_updated_basis_dag[1][base_list[1]]*site_tensors_evolved[1]);
        w_ij+=dag(w_ij);
        updated_wf_basis_evolved_wf_overlap.push_back(w_ij.toComplex().real());
    }

    //params for wf_distance_f, wf_distance_df
    gsl_vector *x;
    x=gsl_vector_alloc(leg_gate_params.size());
    for (int i=0; i<leg_gate_params.size(); i++) 
        gsl_vector_set(x,i,leg_gate_params[i]);
    Wf_Distance_Params *wf_distance_params=new Wf_Distance_Params(N_leg_basis,evolved_wf_norm,updated_wf_basis_overlap,updated_wf_basis_evolved_wf_overlap);

    wf_distance_f(x,wf_distance_params);

    //Compare result 
    cout << "Compare result with direct tensor contraction" << endl;
    Print(evolved_wf_norm);
    Print(updated_wf_norm);
    Print(wf_overlap);

    //Check for wf_distance_df
    cout << "Compare the reuslt for derivative" << endl;
    gsl_vector *df;
    df=gsl_vector_alloc(N_leg_basis);
    wf_distance_df(x,wf_distance_params,df);
    for (int i=0; i<N_leg_basis; i++)
    {
        Print(i);
        Print(gsl_vector_get(df,i));

        double dxi=1E-6,
               x_plus_dxi=gsl_vector_get(x,i)+dxi,
               origin_x=gsl_vector_get(x,i);

        double f_x=wf_distance_f(x,wf_distance_params);
        gsl_vector_set(x,i,x_plus_dxi);
        double f_x_plus_dxi=wf_distance_f(x,wf_distance_params);
        Print((f_x_plus_dxi-f_x)/dxi);
        gsl_vector_set(x,i,origin_x);
    }

    delete wf_distance_params;
    gsl_vector_free(x);
    gsl_vector_free(df);
}


double heisenberg_energy_from_site_env_tensors(const std::array<IQTensor,2> &site_env_tens, const IQTensor &comm_bond_tensor, const NN_Heisenberg_Hamiltonian &hamiltonian_gate)
{
    //PrintDat(hamiltonian_gate.site_tensors(0)*hamiltonian_gate.bond_tensor()*hamiltonian_gate.site_tensors(1));
    std::array<IQTensor,2> site_env_tens_dag={dag(site_env_tens[0]),dag(site_env_tens[1])};
    auto comm_bond_tensor_dag=dag(comm_bond_tensor).prime();
    site_env_tens_dag[0].prime(commonIndex(site_env_tens_dag[0],comm_bond_tensor));
    site_env_tens_dag[1].prime(commonIndex(site_env_tens_dag[1],comm_bond_tensor));

    double wf_norm_sq=((site_env_tens_dag[0]*site_env_tens[0])*comm_bond_tensor*comm_bond_tensor_dag*(site_env_tens_dag[1]*site_env_tens[1])).toComplex().real();

    site_env_tens_dag[0].prime(Site);
    site_env_tens_dag[1].prime(Site);

    double energy=((site_env_tens[0]*hamiltonian_gate.site_tensors(0)*site_env_tens_dag[0])*hamiltonian_gate.bond_tensor()*comm_bond_tensor*comm_bond_tensor_dag*(site_env_tens[1]*hamiltonian_gate.site_tensors(1)*site_env_tens_dag[1])).toComplex().real();

    return energy/wf_norm_sq;
}

double heisenberg_energy_from_site_env_tensors(const std::array<IQTensor,2> &site_env_tens, const NN_Heisenberg_Hamiltonian &hamiltonian_gate)
{
    std::array<IQTensor,2> site_env_tens_dag={dag(site_env_tens[0]),dag(site_env_tens[1])};
    auto comm_ind=commonIndex(site_env_tens_dag[0],dag(site_env_tens_dag[1]));
    site_env_tens_dag[0].prime(comm_ind);
    site_env_tens_dag[1].prime(dag(comm_ind));

    double wf_norm_sq=((site_env_tens_dag[0]*site_env_tens[0])*(site_env_tens_dag[1]*site_env_tens[1])).toComplex().real();

    site_env_tens_dag[0].prime(Site);
    site_env_tens_dag[1].prime(Site);

    double energy=((site_env_tens[0]*hamiltonian_gate.site_tensors(0)*site_env_tens_dag[0])*hamiltonian_gate.bond_tensor()*(site_env_tens[1]*hamiltonian_gate.site_tensors(1)*site_env_tens_dag[1])).toComplex().real();

    return energy/wf_norm_sq;
}
