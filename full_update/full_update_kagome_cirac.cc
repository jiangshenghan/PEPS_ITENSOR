
#include "tnetwork.h"
#include "tnetwork_ctm.h"
#include "full_update_ctm.h"
#include "full_update.h"

void kagome_cirac_rvb_fast_full_update(IQPEPS &kagome_rvb, const Evolution_Params &su_params, const std::array<double,2> &bond_param_norms)
{
    //init environment
    //options: 0:mode. If mode=1: chi only, if mode=2: cutoff only
    //1: if mode=1: this is chi. if mode=2: cut off is 1.E-(this number). eg. this number=14
    //2: number of update steps, 
    //[3,4]: starting_corr (e.g.: could be [1,1])
    //[5,6]: unit cell size (e.g. could be [1,2] for pi_flux)
    //7: svd_method, 1 if using Itensor svd, 2 if using lapack svd.
    int Lx=kagome_rvb.lattice().n_uc()[0], Ly=kagome_rvb.lattice().n_uc()[1];
    iVec env_options="1,30,5,1,1,1,1,2";
    if (abs(kagome_psg::mu_12+1.)<EPSILON) env_options(6)=2;

    //if (!kagome_rvb.site_tensors(0).valid())
    //{
    //    random_init_kagome_rvb_cirac_peps(kagome_rvb,bond_param_norms);
    //}
    Tnetwork_Storage<IQTensor> kagome_rvb_tnetwork_storage=peps_to_tnetwork_storage(kagome_rvb);
    full_update_ctm<IQTensor> full_update_env_tensors(kagome_rvb_tnetwork_storage,env_options);
    auto env_tensors_temp=full_update_env_tensors.output_env_MPOs();

    
    //get a good init rvb state
    while (SzSz_measure(full_update_env_tensors)(0).real()>-0.04)
    {
        random_init_kagome_rvb_cirac_peps(kagome_rvb,bond_param_norms);
        kagome_rvb_tnetwork_storage=peps_to_tnetwork_storage(kagome_rvb);
        full_update_env_tensors=full_update_ctm<IQTensor>(kagome_rvb_tnetwork_storage,env_options);
        full_update_env_tensors.output_env_MPOs();
    }

    //Print(SzSz_measure(full_update_env_tensors));
    //std::vector<IQTensor> env_tensors={full_update_env_tensors.env_tensors(0),full_update_env_tensors.env_tensors(1),full_update_env_tensors.env_tensors(2)};
    std::vector<IQTensor> env_tensors={env_tensors_temp(0),env_tensors_temp(1),env_tensors_temp(2)};
    env_options(2)=1;

    //init params
    std::vector<int> cutting_sites={3*((Ly-1)*Lx+1),3*((Ly-1)*Lx+1)+1,3*((Ly-1)*Lx+1)+2}, cutting_bonds={2*((Ly-1)*Lx+1)};
    std::vector<double> leg_gate_params, bond_params;
    Singlet_Tensor_Basis cutting_bond_basis(kagome_rvb.bond_tensors(cutting_bonds[0]).indices()),
                         cutting_site_basis(kagome_rvb.site_tensors(cutting_sites[0]).indices());
    //bond_params_index stores information of symmetry related bond params. For example, if bond_params_index[j]=a, then the corresponding params for singlet j equals sgn(a)*bond_params[abs(a)-1]
    std::vector<int> bond_params_index(cutting_bond_basis.dim(),0);
    //init bond_params, where we only store parameters not related by symmetry
    //bond singlet params are all params regardless of rotation symmetry
    std::vector<double> bond_singlet_params;
    for (int basei=0; basei<cutting_bond_basis.dim(); basei++)
    {
        double param=(kagome_rvb.bond_tensors(cutting_bonds[0])*dag(cutting_bond_basis[basei])).toReal();
        bond_singlet_params.push_back(param);
    }
    for (int basei=0; basei<cutting_bond_basis.dim(); basei++)
    {
        int bond_params_size=bond_params.size();
        bool params_appeared=false;
        for (int parami=0; parami<bond_params_size; parami++)
        {
            //Print(basei);
            //Print(parami);
            //Print(std::abs(bond_singlet_params[basei]/bond_params[parami]-1));
            //Print(std::abs(bond_singlet_params[basei]/bond_params[parami]+1));
            if (std::abs(bond_singlet_params[basei]/bond_params[parami]-1)<1E-5)
            {
                bond_params_index[basei]=parami+1;
                params_appeared=true;
                break;
            }
            if (std::abs(bond_singlet_params[basei]/bond_params[parami]+1)<1E-5)
            {
                bond_params_index[basei]=-(parami+1);
                params_appeared=true;
                break;
            }
        }
        if (!params_appeared)
        {
            bond_params.push_back(bond_singlet_params[basei]);
            bond_params_index[basei]=bond_params_size+1;
        }
    }
    Print(bond_params);
    Print(bond_params_index);
    Print(bond_singlet_params);

    //init kagome RDM
    PEPSt_RDM<IQTensor> kagome_cirac_rdm("tree shape I",cutting_sites,cutting_bonds,env_tensors,kagome_rvb);
    //Print(kagome_cirac_rdm.RDM());
    //Print(kagome_cirac_rdm.wf_norm());


    for (int iter=0; iter<su_params.iter_nums; iter++)
    {
        //Init evolve gates
        IQTPO evolve_gate=trotter_gate_kagome_cirac({kagome_rvb.phys_legs(cutting_sites[0]),kagome_rvb.phys_legs(cutting_sites[1]),kagome_rvb.phys_legs(cutting_sites[2])},su_params.ts[iter]);

        //reset env_options
        //env_options(2)=(su_params.iter_nums-1-iter)/2+1;
        //env_options(2)=su_params.iter_nums-iter;
        
        //init leg gates for sites, which are used to approx evolve_gate
        //every leg gate is formed by two in legs and one out leg
        std::vector<Singlet_Tensor_Basis> leg_gates_basis;
        for (int evolvei=0; evolvei<cutting_sites.size(); evolvei++)
        {
            IndexSet<IQIndex> leg_gate_indices;
            leg_gate_indices.addindex(commonIndex(dag(kagome_rvb.site_tensors(cutting_sites[evolvei])),kagome_rvb.bond_tensors(cutting_bonds[0])));
            leg_gate_indices.addindex(commonIndex(dag(evolve_gate.site_tensors(evolvei)),evolve_gate.bond_tensors(0)));
            leg_gate_indices.addindex(commonIndex(kagome_rvb.site_tensors(cutting_sites[evolvei]),dag(kagome_rvb.bond_tensors(cutting_bonds[0]))).prime());
            leg_gates_basis.push_back(Singlet_Tensor_Basis(leg_gate_indices));
        }
        //Print(leg_gates_basis);
        
        //init leg_gate_params 
        if (leg_gate_params.empty())
        {
            for (int parami=0; parami<leg_gates_basis[0].dim(); parami++) 
            {
                //leg_gate_params.push_back(rand_gen());
                auto spin_config=leg_gates_basis[0].spin_configs(parami);
                //init leg_gate_params to be identity
                if (spin_config[1]==0)
                {
                    leg_gate_params.push_back(sqrt((spin_config[0]+1)*1.));
                }
                else
                {
                    leg_gate_params.push_back(0);
                }
            }
        }

        //init leg gates for one site tensor
        std::vector<IQTensor> leg_gates_for_one_site;
        auto indice_from_evolve_gate=commonIndex(dag(evolve_gate.site_tensors(0)),evolve_gate.bond_tensors(0));
        for (const auto &virt_leg : kagome_rvb.site_tensors(cutting_sites[0]).indices())
        {
            if (virt_leg.type()==Site) continue;
            std::vector<IQIndex> leg_gate_indices{dag(virt_leg),indice_from_evolve_gate,prime(virt_leg)};
            leg_gates_for_one_site.push_back(IQTensor(leg_gate_indices));
        }
        //Print(leg_gates_for_one_site);

        for (int step=0; step<su_params.steps_nums[iter]; step++)
        {
            Print(iter);
            Print(step);
            Print(su_params.ts[iter]);

            //update env tensors, and measure energy by RDM
            kagome_rvb_tnetwork_storage=peps_to_tnetwork_storage(kagome_rvb);
            full_update_env_tensors.update_env(kagome_rvb_tnetwork_storage,env_options);
            env_tensors_temp=full_update_env_tensors.output_env_MPOs();
            //full_update_env_tensors=full_update_ctm<IQTensor>(kagome_rvb_tnetwork_storage,env_options);
            //env_tensors={full_update_env_tensors.env_tensors(0),full_update_env_tensors.env_tensors(1),full_update_env_tensors.env_tensors(2)};
            env_tensors={env_tensors_temp(0),env_tensors_temp(1),env_tensors_temp(2)};
            kagome_cirac_rdm.update_peps_rdm(env_tensors,kagome_rvb);
            Print(SzSz_measure(full_update_env_tensors));
            Print(heisenberg_energy_from_RDM(kagome_cirac_rdm));

            //get leg_gate_params and bond_params
            double cutoff=0.1*su_params.ts[iter];
            obtain_kagome_cirac_site_leg_gate_bond_params_minimization(kagome_cirac_rdm,evolve_gate,leg_gates_basis,cutting_bond_basis,leg_gate_params,bond_params,bond_params_index,cutoff);

            //using leg_gate_params to generate all leg_gates for one site
            auto leg_gate_sample=singlet_tensor_from_basis_params(leg_gates_basis[0],leg_gate_params);
            for (auto &gate : leg_gates_for_one_site)
            {
                tensor_assignment(gate,leg_gate_sample);
            }

            //update site_tensors
            IQTensor updated_site_tens_unordered=kagome_rvb.site_tensors(cutting_sites[0]);
            for (const auto &site_leg_gate : leg_gates_for_one_site) 
            {
                updated_site_tens_unordered*=evolve_gate.site_tensors(0)*site_leg_gate; 
                updated_site_tens_unordered.noprime();
            }
            //we should never change order of indices of site tensors
            auto updated_site_tens=kagome_rvb.site_tensors(cutting_sites[0]);
            tensor_assignment_diff_order(updated_site_tens,updated_site_tens_unordered);
            rotation_symmetrize_kagome_rvb_cirac_site_tensor(updated_site_tens);
            //keep the same norm
            double site_norm=kagome_rvb.site_tensors(cutting_sites[0]).norm();
            updated_site_tens*=site_norm/(updated_site_tens.norm());
            Print(site_norm);
            Print((updated_site_tens-kagome_rvb.site_tensors(cutting_sites[0])).norm());
            kagome_rvb.generate_site_tensors({updated_site_tens,updated_site_tens,tensor_permutation({0,2,1},updated_site_tens)});

            std::vector<Complex> site_params;
            for (const auto &basis : cutting_site_basis) site_params.push_back((dag(basis)*updated_site_tens).toComplex());
            Print(site_params);

            //update_bond_tensors
            std::vector<double> bond_singlet_params;
            for (int i=0; i<bond_params_index.size(); i++)
            {
                int parami=abs(bond_params_index[i])-1,
                    param_sgn=(parami+1)/bond_params_index[i];
                bond_singlet_params.push_back(param_sgn*bond_params[parami]);
            }
            IQTensor updated_bond_tensor=singlet_tensor_from_basis_params(cutting_bond_basis,bond_singlet_params);
            //rotation_symmetrize_kagome_rvb_cirac_site_tensor(updated_bond_tensor);
            //fix_ratio_kagome_rvb_bond_tensor(updated_bond_tensor,comm_bond_basis,bond_param_norms);
            kagome_rvb.generate_bond_tensors({updated_bond_tensor,tensor_permutation({2,0,1},updated_bond_tensor)},kagome_psg::mu_12);

            //measure energy accurately
            //if ((step*su_params.ts[iter]-floor(step*su_params.ts[iter])<su_params.ts[iter]*0.01 && step%10==0) || step%100==0)
            if (step%10==0)
            {
                cout << "Measure energy accurately!" << endl;
                int env_step=env_options(2);
                env_options(2)=6;
                kagome_rvb_tnetwork_storage=peps_to_tnetwork_storage(kagome_rvb);
                full_update_env_tensors.update_env(kagome_rvb_tnetwork_storage,env_options);
                full_update_env_tensors.output_whole_env();
                Print(SzSz_measure(full_update_env_tensors));
                env_options(2)=env_step;
            }

            //stores as PEPS 
            if ((step+1)*10%su_params.steps_nums[iter]==0)
            {
                std::stringstream ss;
                ss << "/home/jiangsb/tn_ying/tensor_network/result/peps_storage/kagome_rvb_D=" << kagome_rvb.D() << "_Lx=" << kagome_rvb.n_uc()[0] << "_Ly=" << kagome_rvb.n_uc()[1] << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_iter=" << iter << "_step=" << step; 
                std::string file_name=ss.str();
                writeToFile(file_name,kagome_rvb);
            }

        }//end of trotter steps

    }//end of trotter iters

}


//std::vector<IQTensor> cluster_env_tensors_for_kagome_cirac_geometry(const std::vector<int> &cutting_sites, const std::vector<int> &cutting_bonds, const PEPSt<IQTensor> &kagome_rvb, std::array<std::vector<IQTensor>,2> &env_mats)
//{
//    //obtain env matrix
//    IQTensor site_tensA=kagome_rvb.site_tensors(0)*kagome_rvb.bond_tensors(0)*kagome_rvb.site_tensors(1)*kagome_rvb.site_tensors(2);
//    IQTensor site_tensB=kagome_rvb.site_tensors({1,0,0})*kagome_rvb.bond_tensors({0,-1,1})*kagome_rvb.site_tensors(1)*kagome_rvb.site_tensors(2);
//    get_env_tensor_minimization(sit_tensA,site_tensB,env_mats);
//
//    //construct MPO from env_mats
//    Coordinate cutting_bond_coord=kagome_rvb.lattice().bond_list_to_coord(cutting_bonds[0]),
//               cutting_site_coord=kagome_rvb.lattice().site_list_to_coord(cutting_sites[0]);
//    //multiply env mat on boundary leg of bond tensor
//    auto env_mat_double_layer=(env_mat[0][0]*(dag(env_mat[0][0]).mapprime(0,2))).mapprime(2,1);
//
//}



bool obtain_kagome_cirac_site_leg_gate_bond_params_minimization(PEPSt_RDM<IQTensor> &kagome_cirac_rdm, const IQTPO &evolve_gate, const std::vector<Singlet_Tensor_Basis> &leg_gates_basis, const Singlet_Tensor_Basis &cutting_bond_basis, std::vector<double> &leg_gate_params, std::vector<double> &bond_params, const std::vector<int> &bond_params_index, double cutoff)
{
    //construct site and bond tensors after applied by evolve gate
    std::vector<IQTensor> site_tensors_evolved;
    //evolve_legs_combiners are used to combine legs of peps site tensors and evolve gate site tensors
    std::vector<IQCombiner> evolve_legs_combiners;
    for (int cuti=0; cuti<kagome_cirac_rdm.cutting_sites_no(); cuti++)
    {
        site_tensors_evolved.push_back(kagome_cirac_rdm.cutting_site_tensors(cuti)*evolve_gate.site_tensors(cuti));
        site_tensors_evolved[cuti].noprime();
        auto evolve_gate_virt_leg=commonIndex(evolve_gate.site_tensors(cuti),evolve_gate.bond_tensors(0));
        IQIndex cutting_virt_leg=commonIndex(kagome_cirac_rdm.cutting_site_tensors(cuti),kagome_cirac_rdm.cutting_bond_tensors(0));
        evolve_legs_combiners.push_back(IQCombiner(cutting_virt_leg,evolve_gate_virt_leg));
    }
    auto bond_tensor_evolved=kagome_cirac_rdm.cutting_bond_tensors(0)*evolve_gate.bond_tensors(0);
    //Print(site_tensors_evolved);
    //Print(bond_tensor_evolved);

    //obtain evolved_wf_norm
    IQTensor evolve_gate_tensor=evolve_gate.bond_tensors(0);
    for (const auto &tens : evolve_gate.site_tensors()) evolve_gate_tensor*=tens;
    double evolved_wf_norm=sqrt((kagome_cirac_rdm.RDM()*(evolve_gate_tensor*dag(swapPrime(evolve_gate_tensor,0,2))).mapprime(2,1)).toComplex().real());
    //Print(evolved_wf_norm);

    //obtain distance between wf and wf_evolved
    double wf_norm=kagome_cirac_rdm.wf_norm(),
           wf_evolved_wf_overlap=2.*(kagome_cirac_rdm.RDM()*evolve_gate_tensor).toComplex().real(),
           wf_evolved_wf_distance=std::sqrt(2.*pow(evolved_wf_norm,2.)-evolved_wf_norm/wf_norm*wf_evolved_wf_overlap);
    //Print(wf_norm);
    //Print(wf_evolved_wf_overlap/(2.*wf_norm*wf_norm));
    //Print(wf_evolved_wf_distance);

    //using conjugate gradient methods to minimize distance square between updated_wf and evolved_wf
    int find_min_status;
    int iter=0, max_iter=200;

    const gsl_multimin_fdfminimizer_type *minimize_T;
    gsl_multimin_fdfminimizer *s;

    //params to do minimization
    Kagome_Cirac_Wf_Distance_Params_FU *kagome_cirac_wf_distance_params_fu=new Kagome_Cirac_Wf_Distance_Params_FU(evolved_wf_norm,kagome_cirac_rdm,site_tensors_evolved,bond_tensor_evolved,evolve_legs_combiners,leg_gates_basis,cutting_bond_basis,bond_params_index);

    //x stores coefficient for site leg gates and plaquette leg gates
    gsl_vector *x;
    gsl_multimin_function_fdf wf_distance_sq_func;

    wf_distance_sq_func.n=leg_gate_params.size()+bond_params.size();
    wf_distance_sq_func.f=kagome_cirac_wf_distance_sq_fu_f;
    wf_distance_sq_func.df=kagome_cirac_wf_distance_sq_fu_df;
    wf_distance_sq_func.fdf=kagome_cirac_wf_distance_sq_fu_fdf;
    wf_distance_sq_func.params=kagome_cirac_wf_distance_params_fu;

    x=gsl_vector_alloc(wf_distance_sq_func.n);
    for (int i=0; i<leg_gate_params.size(); i++) gsl_vector_set(x,i,leg_gate_params[i]);
    for (int i=0; i<bond_params.size(); i++) gsl_vector_set(x,leg_gate_params.size()+i,bond_params[i]);

    minimize_T=gsl_multimin_fdfminimizer_conjugate_fr;
    s=gsl_multimin_fdfminimizer_alloc(minimize_T,wf_distance_sq_func.n);
    gsl_multimin_fdfminimizer_set(s,&wf_distance_sq_func,x,0.1,0.1);
    do
    {
        iter++;
        find_min_status=gsl_multimin_fdfminimizer_iterate(s);
        if (find_min_status) break;
        find_min_status=gsl_multimin_test_gradient(s->gradient,pow(kagome_cirac_rdm.wf_norm(),2.)*cutoff);

        Print(iter);
        Print(s->f);
        Print(pow(wf_evolved_wf_distance,2.));
    }
    while ((find_min_status==GSL_CONTINUE || sqrt(s->f)>0.5*wf_evolved_wf_distance) && iter<max_iter);
    //while (find_min_status==GSL_CONTINUE && iter<max_iter);


    Print(iter);
    Print(s->f);

    double updated_wf_evolved_wf_distance=sqrt(kagome_cirac_wf_distance_sq_fu_f(s->x,kagome_cirac_wf_distance_params_fu));

    Print(wf_evolved_wf_distance/evolved_wf_norm);
    Print(updated_wf_evolved_wf_distance/evolved_wf_norm);

    //if (iter==max_iter)
    //{
    //    cout << "Leg gate is not good enough, may be trapped in local minima!" << endl << "try smaller time step!" << endl;
    //    gsl_multimin_fdfminimizer_free(s);
    //    gsl_vector_free(x);
    //    delete kagome_cirac_wf_distance_params_fu;

    //    return false;
    //}

    for (int i=0; i<leg_gate_params.size(); i++)
        leg_gate_params[i]=gsl_vector_get(s->x,i);
    for (int i=0; i<bond_params.size(); i++)
        bond_params[i]=gsl_vector_get(s->x,leg_gate_params.size()+i);

    Print(leg_gate_params);
    Print(bond_params);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
    delete kagome_cirac_wf_distance_params_fu;

    return true;
}

double kagome_cirac_wf_distance_sq_fu_f(const gsl_vector *x, void *params)
{
    std::vector<double> leg_gate_params, bond_params;
    Kagome_Cirac_Wf_Distance_Params_FU *kagome_cirac_wf_distance_params_fu=(Kagome_Cirac_Wf_Distance_Params_FU *)params;
    int N_leg_params=kagome_cirac_wf_distance_params_fu->leg_gates_basis[0].dim();

    for (int i=0; i<N_leg_params; i++) 
        leg_gate_params.push_back(gsl_vector_get(x,i));
    for (int i=N_leg_params; i<x->size; i++)
        bond_params.push_back(gsl_vector_get(x,i));
    //Print(leg_gate_params);
    //Print(bond_params);

    //construct leg gates
    std::vector<IQTensor> leg_gates;
    for (const auto &leg_gate_basis : kagome_cirac_wf_distance_params_fu->leg_gates_basis)
    {
        leg_gates.push_back(singlet_tensor_from_basis_params(leg_gate_basis,leg_gate_params));
    }

    auto evolved_site_tensors=kagome_cirac_wf_distance_params_fu->evolved_site_tensors;
    auto evolved_bond_tensor=kagome_cirac_wf_distance_params_fu->evolved_bond_tensor;

    //combine virtual legs of evolved tensors
    const auto &evolve_legs_combiners=kagome_cirac_wf_distance_params_fu->evolve_legs_combiners;
    for (int cuti=0; cuti<3; cuti++)
    {
        evolved_site_tensors[cuti]=evolved_site_tensors[cuti]*evolve_legs_combiners[cuti];
        evolved_bond_tensor=evolved_bond_tensor*dag(evolve_legs_combiners[cuti]);
    }
    //Print(evolved_site_tensors);
    //Print(evolved_bond_tensor);

    //obtain updated_site_tensors and updated_bond_tensors
    std::vector<IQTensor> updated_site_tensors;
    for (int cuti=0; cuti<3; cuti++)
    {
        updated_site_tensors.push_back(kagome_cirac_wf_distance_params_fu->evolved_site_tensors[cuti]*leg_gates[cuti]);
        updated_site_tensors[cuti].noprime();
    }

    //init bond_singlet_params
    std::vector<double> bond_singlet_params;
    const auto &bond_params_index=kagome_cirac_wf_distance_params_fu->bond_params_index;
    for (int i=0; i<bond_params_index.size(); i++)
    {
        int parami=abs(bond_params_index[i])-1,
            param_sgn=(parami+1)/bond_params_index[i];
        bond_singlet_params.push_back(param_sgn*bond_params[parami]);
    }
    IQTensor updated_bond_tensor=singlet_tensor_from_basis_params(kagome_cirac_wf_distance_params_fu->bond_tensor_basis,bond_singlet_params);
    //Print(updated_site_tensors);
    //Print(updated_bond_tensor);
    //auto test_updated_bond_tensor=updated_bond_tensor;
    //rotation_symmetrize_kagome_rvb_cirac_bond_tensor(test_updated_bond_tensor);
    //Print((updated_bond_tensor-test_updated_bond_tensor).norm()/updated_bond_tensor.norm());

    double updated_wf_norm_sq=(kagome_cirac_wf_distance_params_fu->kagome_cirac_rdm.expect_val_from_replaced_tensors(updated_site_tensors,{updated_bond_tensor})).real();
    double updated_evolved_wf_overlap=(2.*kagome_cirac_wf_distance_params_fu->kagome_cirac_rdm.expect_val_from_replaced_tensors({updated_site_tensors,evolved_site_tensors},{std::vector<IQTensor>{updated_bond_tensor},std::vector<IQTensor>{evolved_bond_tensor}})).real();
    double distance_sq=updated_wf_norm_sq+pow(kagome_cirac_wf_distance_params_fu->evolved_wf_norm,2.)-updated_evolved_wf_overlap;
    //Print(kagome_cirac_wf_distance_params_fu->evolved_wf_norm);
    //Print(sqrt(updated_wf_norm_sq));
    //Print(updated_evolved_wf_overlap);

    //gsl_vector *df;
    //df=gsl_vector_alloc(x->size);
    //kagome_cirac_wf_distance_sq_fu_df(x,params,df);
    //for (int i=0; i<x->size; i++)
    //    Print(gsl_vector_get(df,i));
    //kagome_cirac_wf_distance_sq_fu_df_check(x,params,df);
    //for (int i=0; i<x->size; i++)
    //    Print(gsl_vector_get(df,i));
    //gsl_vector_free(df);

    return distance_sq;
}

void kagome_cirac_wf_distance_sq_fu_df(const gsl_vector *x, void *params, gsl_vector *df)
{
    std::vector<double> leg_gate_params, bond_params;
    Kagome_Cirac_Wf_Distance_Params_FU *kagome_cirac_wf_distance_params_fu=(Kagome_Cirac_Wf_Distance_Params_FU *)params;
    int N_leg_params=kagome_cirac_wf_distance_params_fu->leg_gates_basis[0].dim();

    for (int i=0; i<N_leg_params; i++)
        leg_gate_params.push_back(gsl_vector_get(x,i));
    for (int i=N_leg_params; i<x->size; i++)
        bond_params.push_back(gsl_vector_get(x,i));

    //obtain evolved tensors with combined virtual legs
    auto evolved_site_tensors=kagome_cirac_wf_distance_params_fu->evolved_site_tensors;
    auto evolved_bond_tensor=kagome_cirac_wf_distance_params_fu->evolved_bond_tensor;
    const auto &evolve_legs_combiners=kagome_cirac_wf_distance_params_fu->evolve_legs_combiners;
    for (int cuti=0; cuti<3; cuti++)
    {
        evolved_site_tensors[cuti]=evolved_site_tensors[cuti]*evolve_legs_combiners[cuti];
        evolved_bond_tensor=evolved_bond_tensor*dag(evolve_legs_combiners[cuti]);
    }

    //init bond_singlet_params
    std::vector<double> bond_singlet_params;
    const auto &bond_params_index=kagome_cirac_wf_distance_params_fu->bond_params_index;
    for (int i=0; i<bond_params_index.size(); i++)
    {
        int parami=abs(bond_params_index[i])-1,
            param_sgn=(parami+1)/bond_params_index[i];
        bond_singlet_params.push_back(param_sgn*bond_params[parami]);
    }

    //init leg gates
    std::vector<IQTensor> leg_gates;
    for (const auto &leg_gate_basis : kagome_cirac_wf_distance_params_fu->leg_gates_basis)
    {
        leg_gates.push_back(singlet_tensor_from_basis_params(leg_gate_basis,leg_gate_params));
    }

    //init updated_site_tensors and updated_bond_tensor
    std::vector<IQTensor> updated_site_tensors;
    for (int cuti=0; cuti<3; cuti++)
    {
        updated_site_tensors.push_back(kagome_cirac_wf_distance_params_fu->evolved_site_tensors[cuti]*leg_gates[cuti]);
        updated_site_tensors[cuti].noprime();
    }
    IQTensor updated_bond_tensor=singlet_tensor_from_basis_params(kagome_cirac_wf_distance_params_fu->bond_tensor_basis,bond_singlet_params);

    //init updated_wf_norm_sq and updated_evolved_wf_overlap 
    double updated_wf_norm_sq=(kagome_cirac_wf_distance_params_fu->kagome_cirac_rdm.expect_val_from_replaced_tensors(updated_site_tensors,{updated_bond_tensor})).real();
    //here updated_evolved_wf_overlap=\langle\psi|\phi\rangle
    Complex updated_evolved_wf_overlap=kagome_cirac_wf_distance_params_fu->kagome_cirac_rdm.expect_val_from_replaced_tensors({updated_site_tensors,evolved_site_tensors},{std::vector<IQTensor>{updated_bond_tensor},std::vector<IQTensor>{evolved_bond_tensor}});

    //Print(updated_wf_norm_sq);
    //Print(updated_evolved_wf_overlap);

    double dxi=1E-3;
    //get derivative for leg_gate_params
    const auto &leg_gates_basis=kagome_cirac_wf_distance_params_fu->leg_gates_basis;
    for (int i=0; i<N_leg_params; i++)
    {
        double dfi=0;
        std::vector<double> leg_gate_params_dx=leg_gate_params;
        leg_gate_params_dx[i]=leg_gate_params[i]+dxi;
        auto leg_gates_dx=leg_gates;
        auto updated_site_tensors_dx=updated_site_tensors;

        for (int cuti=0; cuti<3; cuti++)
        {
            leg_gates_dx[cuti]=singlet_tensor_from_basis_params(leg_gates_basis[cuti],leg_gate_params_dx);
            updated_site_tensors_dx[cuti]=kagome_cirac_wf_distance_params_fu->evolved_site_tensors[cuti]*leg_gates_dx[cuti];
            updated_site_tensors_dx[cuti].noprime();

            auto updated_wf_norm_sq_dx=kagome_cirac_wf_distance_params_fu->kagome_cirac_rdm.expect_val_from_replaced_tensors({updated_site_tensors_dx,updated_site_tensors},{std::vector<IQTensor>{updated_bond_tensor},std::vector<IQTensor>{updated_bond_tensor}});
            auto updated_evolved_wf_overlap_dx=kagome_cirac_wf_distance_params_fu->kagome_cirac_rdm.expect_val_from_replaced_tensors({updated_site_tensors_dx,evolved_site_tensors},{std::vector<IQTensor>{updated_bond_tensor},std::vector<IQTensor>{evolved_bond_tensor}});

            dfi+=2.*((updated_wf_norm_sq_dx-updated_wf_norm_sq)/dxi-(updated_evolved_wf_overlap_dx-updated_evolved_wf_overlap)/dxi).real();

            //recover leg_gates_dx and updated_site_tensors_dx
            leg_gates_dx[cuti]=leg_gates[cuti];
            updated_site_tensors_dx[cuti]=updated_site_tensors[cuti];

            //Print((updated_wf_norm_sq_dx-updated_wf_norm_sq)/dxi);
            //Print((updated_evolved_wf_overlap_dx-updated_evolved_wf_overlap)/dxi);
        }

        //Print(i);
        //Print(dfi);
        //gsl_vector *x_plus_dxi;
        //x_plus_dxi=gsl_vector_alloc(x->size);
        //gsl_vector_memcpy(x_plus_dxi,x);
        //gsl_vector_set(x_plus_dxi,i,gsl_vector_get(x,i)+1E-10);
        //Print((kagome_cirac_wf_distance_sq_fu_f(x_plus_dxi,params)-kagome_cirac_wf_distance_sq_fu_f(x,params))/1E-10);
        //gsl_vector_set(x_plus_dxi,i,gsl_vector_get(x,i)+1E-12);
        //Print((kagome_cirac_wf_distance_sq_fu_f(x_plus_dxi,params)-kagome_cirac_wf_distance_sq_fu_f(x,params))/1E-12);
        //gsl_vector_free(x_plus_dxi);

        gsl_vector_set(df,i,dfi);
    }

    //get derivative for bond_params
    for (int i=0; i<x->size-N_leg_params; i++)
    {
        double dfi=0;
        std::vector<double> bond_params_dx=bond_params;
        bond_params_dx[i]=bond_params[i]+dxi;

        for (int basei=0; basei<bond_params_index.size(); basei++)
        {
            if (abs(bond_params_index[basei])!=(i+1)) continue;
            std::vector<double> bond_singlet_params_dx=bond_singlet_params;
            bond_singlet_params_dx[basei]+=dxi*(bond_params_index[basei]/(i+1));
            auto updated_bond_tensor_dx=singlet_tensor_from_basis_params(kagome_cirac_wf_distance_params_fu->bond_tensor_basis,bond_singlet_params_dx);

            auto updated_wf_norm_sq_dx=kagome_cirac_wf_distance_params_fu->kagome_cirac_rdm.expect_val_from_replaced_tensors({updated_site_tensors,updated_site_tensors},{std::vector<IQTensor>{updated_bond_tensor_dx},std::vector<IQTensor>{updated_bond_tensor}});
            auto updated_evolved_wf_overlap_dx=kagome_cirac_wf_distance_params_fu->kagome_cirac_rdm.expect_val_from_replaced_tensors({updated_site_tensors,evolved_site_tensors},{std::vector<IQTensor>{updated_bond_tensor_dx},std::vector<IQTensor>{evolved_bond_tensor}});

            //Print((updated_wf_norm_sq_dx-updated_wf_norm_sq)/dxi);
            //Print((updated_evolved_wf_overlap_dx-updated_evolved_wf_overlap)/dxi);

            dfi+=2.*((updated_wf_norm_sq_dx-updated_wf_norm_sq)/dxi-(updated_evolved_wf_overlap_dx-updated_evolved_wf_overlap)/dxi).real();

        }

        //Print(i+N_leg_params);
        //Print(dfi);
        //gsl_vector *x_plus_dxi;
        //x_plus_dxi=gsl_vector_alloc(x->size);
        //gsl_vector_memcpy(x_plus_dxi,x);
        //gsl_vector_set(x_plus_dxi,i+N_leg_params,gsl_vector_get(x,i+N_leg_params)+1E-10);
        //Print((kagome_cirac_wf_distance_sq_fu_f(x_plus_dxi,params)-kagome_cirac_wf_distance_sq_fu_f(x,params))/1E-10);
        //gsl_vector_set(x_plus_dxi,i+N_leg_params,gsl_vector_get(x,i+N_leg_params)+1E-12);
        //Print((kagome_cirac_wf_distance_sq_fu_f(x_plus_dxi,params)-kagome_cirac_wf_distance_sq_fu_f(x,params))/1E-12);
        //gsl_vector_free(x_plus_dxi);

        gsl_vector_set(df,i+N_leg_params,dfi);
    }
}

void kagome_cirac_wf_distance_sq_fu_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
    *f=kagome_cirac_wf_distance_sq_fu_f(x,params);
    kagome_cirac_wf_distance_sq_fu_df(x,params,df);
}

void kagome_cirac_wf_distance_sq_fu_df_check(const gsl_vector *x, void *params, gsl_vector *df)
{
    int x_size=x->size;

    gsl_vector *x_plus_dx;
    x_plus_dx=gsl_vector_alloc(x_size);
    gsl_vector_memcpy(x_plus_dx,x);

    double f_x=kagome_cirac_wf_distance_sq_fu_f(x,params);
    for (int i=0; i<x_size; i++)
    {
        double dxi=1E-10;
        gsl_vector_set(x_plus_dx,i,gsl_vector_get(x,i)+dxi);
        double f_x_plus_dxi=kagome_cirac_wf_distance_sq_fu_f(x_plus_dx,params);
        gsl_vector_set(df,i,(f_x_plus_dxi-f_x)/dxi);
        gsl_vector_set(x_plus_dx,i,gsl_vector_get(x,i));
        //Print(gsl_vector_get(df,i));
    }

    gsl_vector_free(x_plus_dx);
}
