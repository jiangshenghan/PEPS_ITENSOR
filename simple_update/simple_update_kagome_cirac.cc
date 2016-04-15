
#include "simple_update_kagome_cirac.h"

void spin_kagome_cirac_peps_patch_simple_update(IQPEPS &kagome_rvb, const Evolution_Params &su_params, std::vector<int> patch_sites, std::vector<int> evolved_sites, std::string patch_name, std::array<double,2> bond_param_norms)
{
    //Initialize env_tens
    std::array<std::vector<IQTensor>,2> env_tens;

    int comm_bond_no=kagome_rvb.lattice().comm_bond(evolved_sites);

    std::array<std::vector<double>,2> leg_gates_params;

    double site_norm=kagome_rvb.site_tensors(evolved_sites[0]).norm(),
           bond_norm=kagome_rvb.bond_tensors(comm_bond_no).norm();
    Singlet_Tensor_Basis comm_bond_basis(kagome_rvb.bond_tensors(comm_bond_no).indices());

    for (int iter=0; iter<su_params.iter_nums; iter++)
    {
        //Init evolve gate
        IQTPO evolve_gate=trotter_gate_kagome_cirac({kagome_rvb.phys_legs(evolved_sites[0]),kagome_rvb.phys_legs(evolved_sites[1]),kagome_rvb.phys_legs(evolved_sites[2])},su_params.ts[iter]);

        //init leg gates for sites and bonds, which are used to approx evolve_gate
        //every leg gate is formed by two in legs and one out leg
        std::array<std::vector<Singlet_Tensor_Basis>,2> leg_gates_basis;
        for (int evolvei=0; evolvei<evolved_sites.size(); evolvei++)
        {
            IndexSet<IQIndex> leg_gate_indices;
            leg_gate_indices.addindex(commonIndex(dag(kagome_rvb.site_tensors(evolved_sites[evolvei])),kagome_rvb.bond_tensors(comm_bond_no)));
            leg_gate_indices.addindex(commonIndex(dag(evolve_gate.site_tensors(evolvei)),evolve_gate.bond_tensors(0)));
            leg_gate_indices.addindex(commonIndex(kagome_rvb.site_tensors(evolved_sites[evolvei]),dag(kagome_rvb.bond_tensors(comm_bond_no))).prime());
            leg_gates_basis[0].push_back(Singlet_Tensor_Basis(leg_gate_indices));
            leg_gate_indices.dag();
            leg_gates_basis[1].push_back(Singlet_Tensor_Basis(leg_gate_indices));
        }
        //Print(leg_gates_basis[0]);
        //Print(leg_gates_basis[1]);

        //init leg_gates_params
        for (int i=0; i<2; i++)
        {
            if (leg_gates_params[i].empty())
            {
                for (int parami=0; parami<leg_gates_basis[i][0].dim(); parami++) leg_gates_params[i].push_back(rand_gen());
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
        //Print(leg_gates_for_one_tensor[0]);
        //Print(leg_gates_for_one_tensor[1]);

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

                for (auto &gate : leg_gates_for_one_tensor[i])
                {
                    tensor_assignment(gate,leg_gate_sample);
                }

            }

            //updte site_tensors
            IQTensor updated_site_tens_unordered=kagome_rvb.site_tensors(evolved_sites[0]);
            for (const auto &site_leg_gate : leg_gates_for_one_tensor[0]) 
            {
                updated_site_tens_unordered*=evolve_gate.site_tensors(0)*site_leg_gate; 
                updated_site_tens_unordered.noprime();
            }
            //we should never change order of indices of site tensors
            auto updated_site_tens=kagome_rvb.site_tensors(evolved_sites[0]);
            tensor_assignment_diff_order(updated_site_tens,updated_site_tens_unordered);
            rotation_symmetrize_kagome_rvb_cirac_site_tensor(updated_site_tens);
            //keep the same norm
            updated_site_tens*=site_norm/(updated_site_tens.norm());
            kagome_rvb.generate_site_tensors({updated_site_tens,updated_site_tens,tensor_permutation({0,2,1},updated_site_tens)});

            //update bond tensors
            IQTensor updated_bond_tens_unordered=kagome_rvb.bond_tensors(comm_bond_no)*evolve_gate.bond_tensors(0);
            for (const auto &bond_leg_gate : leg_gates_for_one_tensor[1])
            {
                updated_bond_tens_unordered*=bond_leg_gate;
            }
            updated_bond_tens_unordered.noprime();
            //we should never change order of inds
            auto updated_bond_tens=kagome_rvb.bond_tensors(comm_bond_no);
            tensor_assignment_diff_order(updated_bond_tens,updated_bond_tens_unordered);
            rotation_symmetrize_kagome_rvb_cirac_bond_tensor(updated_bond_tens);
            //keep the same norm
            //updated_bond_tens*=bond_norm/(updated_bond_tens.norm());
            fix_ratio_kagome_rvb_bond_tensor(updated_bond_tens,comm_bond_basis,bond_param_norms);
            kagome_rvb.generate_bond_tensors({updated_bond_tens,tensor_permutation({2,0,1},updated_bond_tens)},kagome_psg::mu_12);

            //stores as PEPS
            if ((step+1)*10%su_params.steps_nums[iter]==0)
            {
                std::stringstream ss;

                ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/kagome_rvb_D=" << kagome_rvb.D() << "_Lx=" << kagome_rvb.n_uc()[0] << "_Ly=" << kagome_rvb.n_uc()[1] << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_iter=" << iter << "_step=" << step << "_" << kagome_patch_RDM.patch_name();

                std::string file_name=ss.str();
                writeToFile(file_name,kagome_rvb);

            }

        }//end of trotter steps

    }// end of trotter iters
}


void spin_kagome_cirac_peps_patch_simple_update_no_bond_leg_gates(IQPEPS &kagome_rvb, const Evolution_Params &su_params, std::vector<int> patch_sites, std::vector<int> evolved_sites, std::string patch_name, std::array<double,2> bond_param_norms)
{
    std::array<std::vector<IQTensor>,2> env_tens;
    int comm_bond_no=kagome_rvb.lattice().comm_bond(evolved_sites);
    double site_norm=kagome_rvb.site_tensors(evolved_sites[0]).norm();
    std::vector<double> leg_gate_params, bond_params;
    Singlet_Tensor_Basis comm_bond_basis(kagome_rvb.bond_tensors(comm_bond_no).indices());
    IQTensor &comm_bond_tensor=kagome_rvb.bond_tensors(comm_bond_no);
    //bond_params_index stores information of symmetry related bond params. For example, if bond_params_index[j]=a, then the corresponding params for singlet j equals sgn(a)*bond_params[abs(a)-1]
    std::vector<int> bond_params_index(comm_bond_basis.dim(),0);
    //init bond_params, where we only store parameters not related by symmetry
    std::vector<double> bond_singlet_params;
    for (int basei=0; basei<comm_bond_basis.dim(); basei++)
    {
        double param=(comm_bond_tensor*dag(comm_bond_basis[basei])).toReal();
        bond_singlet_params.push_back(param);
    }
    for (int basei=0; basei<comm_bond_basis.dim(); basei++)
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


    //time evolutions
    for (int iter=0; iter<su_params.iter_nums; iter++)
    {
        //Init evolve gates
        IQTPO evolve_gate=trotter_gate_kagome_cirac({kagome_rvb.phys_legs(evolved_sites[0]),kagome_rvb.phys_legs(evolved_sites[1]),kagome_rvb.phys_legs(evolved_sites[2])},su_params.ts[iter]);
        
        //init leg gates for sites, which are used to approx evolve_gate
        //every leg gate is formed by two in legs and one out leg
        std::vector<Singlet_Tensor_Basis> leg_gates_basis;
        for (int evolvei=0; evolvei<evolved_sites.size(); evolvei++)
        {
            IndexSet<IQIndex> leg_gate_indices;
            leg_gate_indices.addindex(commonIndex(dag(kagome_rvb.site_tensors(evolved_sites[evolvei])),kagome_rvb.bond_tensors(comm_bond_no)));
            leg_gate_indices.addindex(commonIndex(dag(evolve_gate.site_tensors(evolvei)),evolve_gate.bond_tensors(0)));
            leg_gate_indices.addindex(commonIndex(kagome_rvb.site_tensors(evolved_sites[evolvei]),dag(kagome_rvb.bond_tensors(comm_bond_no))).prime());
            leg_gates_basis.push_back(Singlet_Tensor_Basis(leg_gate_indices));
        }
        //Print(leg_gates_basis);
        
        //init leg_gates_params
        if (leg_gate_params.empty())
        {
            for (int parami=0; parami<leg_gates_basis[0].dim(); parami++) leg_gate_params.push_back(rand_gen());
        }

        //init leg gates for one site tensor
        std::vector<IQTensor> leg_gates_for_one_site;
        auto indice_from_evolve_gate=commonIndex(dag(evolve_gate.site_tensors(0)),evolve_gate.bond_tensors(0));
        for (const auto &virt_leg : kagome_rvb.site_tensors(evolved_sites[0]).indices())
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

            //obtain env tensors, and measure energy by patch
            IQTensor site_tensA=kagome_rvb.site_tensors(0)*kagome_rvb.bond_tensors(0)*kagome_rvb.site_tensors(1)*kagome_rvb.site_tensors(2);
            IQTensor site_tensB=kagome_rvb.site_tensors({1,0,0})*kagome_rvb.bond_tensors({0,-1,1})*kagome_rvb.site_tensors({1,-1,1});
            get_env_tensor_minimization(site_tensA,site_tensB,env_tens);

            General_Patch_RDM<IQTensor> kagome_patch_RDM(patch_name,kagome_rvb,env_tens[0][0],patch_sites,evolved_sites);
            Print(heisenberg_energy_from_RDM(kagome_patch_RDM));

            //cutoff for leg gate 
            double cutoff=1e-5;
            if (!obtain_kagome_cirac_site_leg_gate_bond_params_minimization(kagome_patch_RDM,evolve_gate,leg_gates_basis,comm_bond_basis,leg_gate_params,bond_params,bond_params_index,cutoff)) break;

            //using leg_gate_params to generate all leg_gates for one site
            auto leg_gate_sample=singlet_tensor_from_basis_params(leg_gates_basis[0],leg_gate_params);
            for (auto &gate : leg_gates_for_one_site)
            {
                tensor_assignment(gate,leg_gate_sample);
            }

            //update site_tensors
            IQTensor updated_site_tens_unordered=kagome_rvb.site_tensors(evolved_sites[0]);
            for (const auto &site_leg_gate : leg_gates_for_one_site) 
            {
                updated_site_tens_unordered*=evolve_gate.site_tensors(0)*site_leg_gate; 
                updated_site_tens_unordered.noprime();
            }
            //we should never change order of indices of site tensors
            auto updated_site_tens=kagome_rvb.site_tensors(evolved_sites[0]);
            tensor_assignment_diff_order(updated_site_tens,updated_site_tens_unordered);
            rotation_symmetrize_kagome_rvb_cirac_site_tensor(updated_site_tens);
            //keep the same norm
            updated_site_tens*=site_norm/(updated_site_tens.norm());
            kagome_rvb.generate_site_tensors({updated_site_tens,updated_site_tens,tensor_permutation({0,2,1},updated_site_tens)});

            //update_bond_tensors
            std::vector<double> bond_singlet_params;
            for (int i=0; i<bond_params_index.size(); i++)
            {
                int parami=abs(bond_params_index[i])-1,
                    param_sgn=(parami+1)/bond_params_index[i];
                bond_singlet_params.push_back(param_sgn*bond_params[parami]);
            }
            IQTensor updated_bond_tensor=singlet_tensor_from_basis_params(comm_bond_basis,bond_singlet_params);
            //rotation_symmetrize_kagome_rvb_cirac_site_tensor(updated_bond_tensor);
            //fix_ratio_kagome_rvb_bond_tensor(updated_bond_tensor,comm_bond_basis,bond_param_norms);
            kagome_rvb.generate_bond_tensors({updated_bond_tensor,tensor_permutation({2,0,1},updated_bond_tensor)},kagome_psg::mu_12);

            //stores as PEPS
            if ((step+1)*10%su_params.steps_nums[iter]==0)
            {
                std::stringstream ss;

                ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/kagome_rvb_D=" << kagome_rvb.D() << "_Lx=" << kagome_rvb.n_uc()[0] << "_Ly=" << kagome_rvb.n_uc()[1] << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_iter=" << iter << "_step=" << step << "_" << kagome_patch_RDM.patch_name();

                std::string file_name=ss.str();
                writeToFile(file_name,kagome_rvb);

            }
        }//end of trotter steps

    }//end of trotter iters
}


bool obtain_kagome_cirac_leg_gates_params_minimization(General_Patch_RDM<IQTensor> &kagome_patch_RDM, const IQTPO &evolve_gate, const std::array<std::vector<Singlet_Tensor_Basis>,2> &leg_gates_basis, std::array<std::vector<double>,2> &leg_gates_params, double cutoff)
{
    //construct site tensors and bond tensors after applied by evolve_gate
    std::vector<IQTensor> site_tensors_evolved;
    //evolve_legs_combiners are used to combine legs of peps site tensors and evolve gate site tensors
    std::vector<IQCombiner> evolve_legs_combiners;
    for (int cuti=0; cuti<kagome_patch_RDM.cutting_sites_no(); cuti++)
    {
        site_tensors_evolved.push_back(kagome_patch_RDM.cutting_site_tensors(cuti)*evolve_gate.site_tensors(cuti));
        site_tensors_evolved[cuti].noprime();
        auto evolve_gate_virt_leg=commonIndex(evolve_gate.site_tensors(cuti),evolve_gate.bond_tensors(0));
        evolve_legs_combiners.push_back(IQCombiner(kagome_patch_RDM.cutting_virt_legs(cuti),evolve_gate_virt_leg));
    }
    auto bond_tensor_evolved=kagome_patch_RDM.cutting_bond_tensor()*evolve_gate.bond_tensors(0);
    //PrintDat(site_tensors_evolved);
    //PrintDat(bond_tensor_evolved);
    //Print(evolve_legs_combiners);

    //obtain evolved_wf_norm
    IQTensor evolve_gate_tensor=evolve_gate.bond_tensors(0);
    for (const auto &tens : evolve_gate.site_tensors()) evolve_gate_tensor*=tens;
    double evolved_wf_norm=sqrt((kagome_patch_RDM.RDM()*(evolve_gate_tensor*dag(swapPrime(evolve_gate_tensor,0,2))).mapprime(2,1)).real());
    //Print(evolved_wf_norm);

    //using conjugate gradient methods to minimize distance square between updated_wf and evolved_wf
    int find_min_status;
    int iter=0, max_iter=1E4;

    const gsl_multimin_fdfminimizer_type *minimize_T;
    gsl_multimin_fdfminimizer *s;

    //params to do minimization
    Kagome_Cirac_Wf_Distance_Params *kagome_cirac_wf_distance_params=new Kagome_Cirac_Wf_Distance_Params(evolved_wf_norm,kagome_patch_RDM,site_tensors_evolved,bond_tensor_evolved,evolve_legs_combiners,leg_gates_basis);

    //x stores coefficient for site leg gates and plaquette leg gates
    gsl_vector *x;
    gsl_multimin_function_fdf wf_distance_sq_func;

    wf_distance_sq_func.n=leg_gates_params[0].size()+leg_gates_params[1].size();
    wf_distance_sq_func.f=kagome_cirac_wf_distance_sq_f;
    wf_distance_sq_func.df=kagome_cirac_wf_distance_sq_df;
    wf_distance_sq_func.fdf=kagome_cirac_wf_distance_sq_fdf;
    wf_distance_sq_func.params=kagome_cirac_wf_distance_params;

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
        find_min_status=gsl_multimin_test_gradient(s->gradient,kagome_patch_RDM.wf_norm()*cutoff);

        //Print(iter);
        //Print(s->f);
    }
    while (find_min_status==GSL_CONTINUE && iter<max_iter);


    //using minimization methods without derivative to minimize distance square between updated_wf and evolved_wf
    //gsl_vector *ss;
    //ss=gsl_vector_alloc(x->size);
    //minimize_T=gsl_multimin_fminimizer_nmsimplex2;
    //gsl_multimin_fminimizer_set(s,&wf_distance_sq_func,x,0.1,0.1);


    Print(iter);
    Print(s->f);

    double wf_norm=kagome_patch_RDM.wf_norm(),
           wf_evolved_wf_overlap=2*(kagome_patch_RDM.expect_val_from_RDM(evolve_gate_tensor)).real(),
           wf_evolved_wf_distance=std::sqrt(2*pow(evolved_wf_norm,2.)-evolved_wf_norm/wf_norm*wf_evolved_wf_overlap),
           updated_wf_evolved_wf_distance=sqrt(kagome_cirac_wf_distance_sq_f(s->x,kagome_cirac_wf_distance_params));

    Print(wf_evolved_wf_distance/evolved_wf_norm);
    Print(updated_wf_evolved_wf_distance/evolved_wf_norm);

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

bool obtain_kagome_cirac_site_leg_gate_bond_params_minimization(General_Patch_RDM<IQTensor> &kagome_patch_RDM, const IQTPO &evolve_gate, const std::vector<Singlet_Tensor_Basis> &leg_gates_basis, const Singlet_Tensor_Basis &comm_bond_basis, std::vector<double> &leg_gate_params, std::vector<double> &bond_params, const std::vector<int> &bond_params_index, double cutoff)
{
    //construct site and bond tensors after applied by evolve gate
    std::vector<IQTensor> site_tensors_evolved;
    //evolve_legs_combiners are used to combine legs of peps site tensors and evolve gate site tensors
    std::vector<IQCombiner> evolve_legs_combiners;
    for (int cuti=0; cuti<kagome_patch_RDM.cutting_sites_no(); cuti++)
    {
        site_tensors_evolved.push_back(kagome_patch_RDM.cutting_site_tensors(cuti)*evolve_gate.site_tensors(cuti));
        site_tensors_evolved[cuti].noprime();
        auto evolve_gate_virt_leg=commonIndex(evolve_gate.site_tensors(cuti),evolve_gate.bond_tensors(0));
        evolve_legs_combiners.push_back(IQCombiner(kagome_patch_RDM.cutting_virt_legs(cuti),evolve_gate_virt_leg));
    }
    auto bond_tensor_evolved=kagome_patch_RDM.cutting_bond_tensor()*evolve_gate.bond_tensors(0);
    //Print(site_tensors_evolved);
    //Print(bond_tensor_evolved);

    //obtain evolved_wf_norm
    IQTensor evolve_gate_tensor=evolve_gate.bond_tensors(0);
    for (const auto &tens : evolve_gate.site_tensors()) evolve_gate_tensor*=tens;
    double evolved_wf_norm=sqrt((kagome_patch_RDM.RDM()*(evolve_gate_tensor*dag(swapPrime(evolve_gate_tensor,0,2))).mapprime(2,1)).real());
    //Print(evolved_wf_norm);

    //obtain distance between wf and wf_evolved
    double wf_norm=kagome_patch_RDM.wf_norm(),
           wf_evolved_wf_overlap=2*(kagome_patch_RDM.expect_val_from_RDM(evolve_gate_tensor)).real(),
           wf_evolved_wf_distance=std::sqrt(2*pow(evolved_wf_norm,2.)-evolved_wf_norm/wf_norm*wf_evolved_wf_overlap);

    //using conjugate gradient methods to minimize distance square between updated_wf and evolved_wf
    int find_min_status;
    int iter=0, max_iter=1E4;

    const gsl_multimin_fdfminimizer_type *minimize_T;
    gsl_multimin_fdfminimizer *s;

    //params to do minimization
    Kagome_Cirac_Wf_Distance_Params_NBLG *kagome_cirac_wf_distance_params_nblg=new Kagome_Cirac_Wf_Distance_Params_NBLG(evolved_wf_norm,kagome_patch_RDM,site_tensors_evolved,bond_tensor_evolved,evolve_legs_combiners,leg_gates_basis,comm_bond_basis,bond_params_index);

    //x stores coefficient for site leg gates and plaquette leg gates
    gsl_vector *x;
    gsl_multimin_function_fdf wf_distance_sq_nblg_func;

    wf_distance_sq_nblg_func.n=leg_gate_params.size()+bond_params.size();
    wf_distance_sq_nblg_func.f=kagome_cirac_wf_distance_sq_nblg_f;
    wf_distance_sq_nblg_func.df=kagome_cirac_wf_distance_sq_nblg_df;
    wf_distance_sq_nblg_func.fdf=kagome_cirac_wf_distance_sq_nblg_fdf;
    wf_distance_sq_nblg_func.params=kagome_cirac_wf_distance_params_nblg;

    x=gsl_vector_alloc(wf_distance_sq_nblg_func.n);
    for (int i=0; i<leg_gate_params.size(); i++) gsl_vector_set(x,i,leg_gate_params[i]);
    for (int i=0; i<bond_params.size(); i++) gsl_vector_set(x,leg_gate_params.size()+i,bond_params[i]);

    minimize_T=gsl_multimin_fdfminimizer_conjugate_fr;
    s=gsl_multimin_fdfminimizer_alloc(minimize_T,wf_distance_sq_nblg_func.n);
    gsl_multimin_fdfminimizer_set(s,&wf_distance_sq_nblg_func,x,0.1,0.1);
    do
    {
        iter++;
        find_min_status=gsl_multimin_fdfminimizer_iterate(s);
        if (find_min_status) break;
        find_min_status=gsl_multimin_test_gradient(s->gradient,kagome_patch_RDM.wf_norm()*cutoff);

        Print(iter);
        Print(s->f);
    }
    while (find_min_status==GSL_CONTINUE && iter<max_iter);


    Print(iter);
    Print(s->f);

    double updated_wf_evolved_wf_distance=sqrt(kagome_cirac_wf_distance_sq_nblg_f(s->x,kagome_cirac_wf_distance_params_nblg));

    Print(wf_evolved_wf_distance/evolved_wf_norm);
    Print(updated_wf_evolved_wf_distance/evolved_wf_norm);

    if (iter==max_iter)
    {
        cout << "Leg gate is not good enough, may be trapped in local minima!" << endl << "try smaller time step!" << endl;
        gsl_multimin_fdfminimizer_free(s);
        gsl_vector_free(x);
        delete kagome_cirac_wf_distance_params_nblg;

        return false;
    }

    for (int i=0; i<leg_gate_params.size(); i++)
        leg_gate_params[i]=gsl_vector_get(s->x,i);
    for (int i=0; i<bond_params.size(); i++)
        bond_params[i]=gsl_vector_get(s->x,leg_gate_params.size()+i);

    Print(leg_gate_params);
    Print(bond_params);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
    delete kagome_cirac_wf_distance_params_nblg;

    return true;


}

double kagome_cirac_wf_distance_sq_f(const gsl_vector *x, void *params)
{
    std::array<std::vector<double>,2> leg_gates_params;
    Kagome_Cirac_Wf_Distance_Params *kagome_cirac_wf_distance_params=(Kagome_Cirac_Wf_Distance_Params *)params;
    std::array<int,2> N_legs_basis={kagome_cirac_wf_distance_params->leg_gates_basis[0][0].dim(),kagome_cirac_wf_distance_params->leg_gates_basis[1][0].dim()};

    //Print(kagome_cirac_wf_distance_params->evolved_wf_norm);
    //Print(kagome_cirac_wf_distance_params->kagome_patch_RDM);
    //Print(kagome_cirac_wf_distance_params->evolved_site_tensors);
    //Print(kagome_cirac_wf_distance_params->evolved_bond_tensor);
    //Print(kagome_cirac_wf_distance_params->evolve_legs_combiners);
    //Print(kagome_cirac_wf_distance_params->leg_gates_basis[0]);
    //Print(kagome_cirac_wf_distance_params->leg_gates_basis[1]);
    //Print(N_legs_basis[0]);
    //Print(N_legs_basis[1]);
   
    for (int i=0; i<N_legs_basis[0]; i++)
        leg_gates_params[0].push_back(gsl_vector_get(x,i));
    for (int i=0; i<N_legs_basis[1]; i++)
        leg_gates_params[1].push_back(gsl_vector_get(x,i+N_legs_basis[0]));

    //Print(leg_gates_params[0]);
    //Print(leg_gates_params[1]);

    //construct leg gates
    std::array<std::vector<IQTensor>,2> leg_gates;
    for (int typei=0; typei<2; typei++)
    {
        for (const auto &singlet_basis : kagome_cirac_wf_distance_params->leg_gates_basis[typei])
        {
            leg_gates[typei].push_back(singlet_tensor_from_basis_params(singlet_basis,leg_gates_params[typei]));
        }
    }
    //PrintDat(leg_gates[0]);
    //PrintDat(leg_gates[1]);

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
        evolved_site_tensors[cuti]=evolved_site_tensors[cuti]*evolve_legs_combiners[cuti];
        evolved_bond_tensor=evolved_bond_tensor*dag(evolve_legs_combiners[cuti]);
    }
    //PrintDat(evolved_site_tensors);
    //PrintDat(evolved_bond_tensor);


    double updated_wf_norm_sq=(kagome_cirac_wf_distance_params->kagome_patch_RDM.expect_val_from_replaced_tensors({updated_site_tensors[0]*updated_bond_tensor,updated_site_tensors[1],updated_site_tensors[2]})).real();

    double updated_evolved_wf_overlap=(2.*kagome_cirac_wf_distance_params->kagome_patch_RDM.expect_val_from_replaced_tensors({{updated_site_tensors[0]*updated_bond_tensor,evolved_site_tensors[0]*evolved_bond_tensor},{updated_site_tensors[1],evolved_site_tensors[1]},{updated_site_tensors[2],evolved_site_tensors[2]}})).real();

    double distance_sq=updated_wf_norm_sq+pow(kagome_cirac_wf_distance_params->evolved_wf_norm,2.)-updated_evolved_wf_overlap;

    //Print(updated_evolved_wf_overlap);
    //Print(kagome_cirac_wf_distance_params->evolved_wf_norm);
    //Print(sqrt(updated_wf_norm_sq));
    //Print(sqrt(distance_sq));

    return distance_sq;
}

void kagome_cirac_wf_distance_sq_df(const gsl_vector *x, void *params, gsl_vector *df)
{
    int x_size=x->size;

    gsl_vector *x_plus_dx;
    x_plus_dx=gsl_vector_alloc(x_size);
    gsl_vector_memcpy(x_plus_dx,x);

    double f_x=kagome_cirac_wf_distance_sq_f(x,params);
    for (int i=0; i<x_size; i++)
    {
        double dxi=1E-10;
        if (dxi>f_x/10.) dxi=f_x/10.;
        gsl_vector_set(x_plus_dx,i,gsl_vector_get(x,i)+dxi);
        double f_x_plus_dxi=kagome_cirac_wf_distance_sq_f(x_plus_dx,params);
        gsl_vector_set(df,i,(f_x_plus_dxi-f_x)/dxi);
        gsl_vector_set(x_plus_dx,i,gsl_vector_get(x,i));
        //Print(gsl_vector_get(df,i));
    }

    gsl_vector_free(x_plus_dx);
}

void kagome_cirac_wf_distance_sq_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
    *f=kagome_cirac_wf_distance_sq_f(x,params);
    kagome_cirac_wf_distance_sq_df(x,params,df);
}

double kagome_cirac_wf_distance_sq_nblg_f(const gsl_vector *x, void *params)
{
    std::vector<double> leg_gate_params, bond_params, bond_singlet_params;
    Kagome_Cirac_Wf_Distance_Params_NBLG *kagome_cirac_wf_distance_params_nblg=(Kagome_Cirac_Wf_Distance_Params_NBLG *)params;
    int N_leg_params=kagome_cirac_wf_distance_params_nblg->leg_gates_basis[0].dim();

    for (int i=0; i<N_leg_params; i++) 
        leg_gate_params.push_back(gsl_vector_get(x,i));
    for (int i=N_leg_params; i<x->size; i++)
        bond_params.push_back(gsl_vector_get(x,i));
    //Print(leg_gate_params);
    //Print(bond_params);

    //construct leg gates
    std::vector<IQTensor> leg_gates;
    for (const auto &leg_gate_basis : kagome_cirac_wf_distance_params_nblg->leg_gates_basis)
    {
        leg_gates.push_back(singlet_tensor_from_basis_params(leg_gate_basis,leg_gate_params));
    }

    auto evolved_site_tensors=kagome_cirac_wf_distance_params_nblg->evolved_site_tensors;
    auto evolved_bond_tensor=kagome_cirac_wf_distance_params_nblg->evolved_bond_tensor;

    //combine virtual legs of evolved tensors
    const auto &evolve_legs_combiners=kagome_cirac_wf_distance_params_nblg->evolve_legs_combiners;
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
        updated_site_tensors.push_back(kagome_cirac_wf_distance_params_nblg->evolved_site_tensors[cuti]*leg_gates[cuti]);
        updated_site_tensors[cuti].noprime();
    }

    //init bond_singlet_params
    const auto &bond_params_index=kagome_cirac_wf_distance_params_nblg->bond_params_index;
    for (int i=0; i<bond_params_index.size(); i++)
    {
        int parami=abs(bond_params_index[i])-1,
            param_sgn=(parami+1)/bond_params_index[i];
        bond_singlet_params.push_back(param_sgn*bond_params[parami]);
    }
    IQTensor updated_bond_tensor=singlet_tensor_from_basis_params(kagome_cirac_wf_distance_params_nblg->bond_tensor_basis,bond_singlet_params);
    //Print(updated_site_tensors);
    //Print(updated_bond_tensor);
    //auto test_updated_bond_tensor=updated_bond_tensor;
    //rotation_symmetrize_kagome_rvb_cirac_bond_tensor(test_updated_bond_tensor);
    //Print((updated_bond_tensor-test_updated_bond_tensor).norm()/updated_bond_tensor.norm());

    double updated_wf_norm_sq=(kagome_cirac_wf_distance_params_nblg->kagome_patch_RDM.expect_val_from_replaced_tensors({updated_site_tensors[0]*updated_bond_tensor,updated_site_tensors[1],updated_site_tensors[2]})).real();
    double updated_evolved_wf_overlap=(2.*kagome_cirac_wf_distance_params_nblg->kagome_patch_RDM.expect_val_from_replaced_tensors({{updated_site_tensors[0]*updated_bond_tensor,evolved_site_tensors[0]*evolved_bond_tensor},{updated_site_tensors[1],evolved_site_tensors[1]},{updated_site_tensors[2],evolved_site_tensors[2]}})).real();
    double distance_sq=updated_wf_norm_sq+pow(kagome_cirac_wf_distance_params_nblg->evolved_wf_norm,2.)-updated_evolved_wf_overlap;
    //Print(kagome_cirac_wf_distance_params_nblg->evolved_wf_norm);
    //Print(sqrt(updated_wf_norm_sq));
    //Print(updated_evolved_wf_overlap);
    //Print(sqrt(distance_sq));

    //gsl_vector *df;
    //df=gsl_vector_alloc(x->size);
    //kagome_cirac_wf_distance_sq_nblg_df(x,params,df);
    //for (int i=0; i<x->size; i++)
    //    Print(gsl_vector_get(df,i));
    //kagome_cirac_wf_distance_sq_nblg_df_check(x,params,df);
    //for (int i=0; i<x->size; i++)
    //    Print(gsl_vector_get(df,i));
    //gsl_vector_free(df);

    return distance_sq;
}

void kagome_cirac_wf_distance_sq_nblg_df(const gsl_vector *x, void *params, gsl_vector *df)
{
    std::vector<double> leg_gate_params, bond_params, bond_singlet_params;
    Kagome_Cirac_Wf_Distance_Params_NBLG *kagome_cirac_wf_distance_params_nblg=(Kagome_Cirac_Wf_Distance_Params_NBLG *)params;
    int N_leg_params=kagome_cirac_wf_distance_params_nblg->leg_gates_basis[0].dim();

    for (int i=0; i<N_leg_params; i++)
        leg_gate_params.push_back(gsl_vector_get(x,i));
    for (int i=N_leg_params; i<x->size; i++)
        bond_params.push_back(gsl_vector_get(x,i));

    //obtain evolved tensors with combined virtual legs
    auto evolved_site_tensors=kagome_cirac_wf_distance_params_nblg->evolved_site_tensors;
    auto evolved_bond_tensor=kagome_cirac_wf_distance_params_nblg->evolved_bond_tensor;
    const auto &evolve_legs_combiners=kagome_cirac_wf_distance_params_nblg->evolve_legs_combiners;
    for (int cuti=0; cuti<3; cuti++)
    {
        evolved_site_tensors[cuti]=evolved_site_tensors[cuti]*evolve_legs_combiners[cuti];
        evolved_bond_tensor=evolved_bond_tensor*dag(evolve_legs_combiners[cuti]);
    }

    //init bond_singlet_params
    const auto &bond_params_index=kagome_cirac_wf_distance_params_nblg->bond_params_index;
    for (int i=0; i<bond_params_index.size(); i++)
    {
        int parami=abs(bond_params_index[i])-1,
            param_sgn=(parami+1)/bond_params_index[i];
        bond_singlet_params.push_back(param_sgn*bond_params[parami]);
    }

    //init leg gates
    std::vector<IQTensor> leg_gates;
    for (const auto &leg_gate_basis : kagome_cirac_wf_distance_params_nblg->leg_gates_basis)
    {
        leg_gates.push_back(singlet_tensor_from_basis_params(leg_gate_basis,leg_gate_params));
    }

    //init updated_site_tensors and updated_bond_tensor
    std::vector<IQTensor> updated_site_tensors;
    for (int cuti=0; cuti<3; cuti++)
    {
        updated_site_tensors.push_back(kagome_cirac_wf_distance_params_nblg->evolved_site_tensors[cuti]*leg_gates[cuti]);
        updated_site_tensors[cuti].noprime();
    }
    IQTensor updated_bond_tensor=singlet_tensor_from_basis_params(kagome_cirac_wf_distance_params_nblg->bond_tensor_basis,bond_singlet_params);

    //init updated_wf_norm_sq and updated_evolved_wf_overlap 
    double updated_wf_norm_sq=(kagome_cirac_wf_distance_params_nblg->kagome_patch_RDM.expect_val_from_replaced_tensors({updated_site_tensors[0]*updated_bond_tensor,updated_site_tensors[1],updated_site_tensors[2]})).real();
    //here updated_evolved_wf_overlap=\langle\psi|\phi\rangle
    Complex updated_evolved_wf_overlap=kagome_cirac_wf_distance_params_nblg->kagome_patch_RDM.expect_val_from_replaced_tensors({{updated_site_tensors[0]*updated_bond_tensor,evolved_site_tensors[0]*evolved_bond_tensor},{updated_site_tensors[1],evolved_site_tensors[1]},{updated_site_tensors[2],evolved_site_tensors[2]}});

    //Print(updated_wf_norm_sq);
    //Print(updated_evolved_wf_overlap);

    double dxi=1E-3;
    //get derivative for leg_gate_params
    const auto &leg_gates_basis=kagome_cirac_wf_distance_params_nblg->leg_gates_basis;
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
            updated_site_tensors_dx[cuti]=kagome_cirac_wf_distance_params_nblg->evolved_site_tensors[cuti]*leg_gates_dx[cuti];
            updated_site_tensors_dx[cuti].noprime();

            auto updated_wf_norm_sq_dx=kagome_cirac_wf_distance_params_nblg->kagome_patch_RDM.expect_val_from_replaced_tensors({{updated_site_tensors_dx[0]*updated_bond_tensor,updated_site_tensors[0]*updated_bond_tensor},{updated_site_tensors_dx[1],updated_site_tensors[1]},{updated_site_tensors_dx[2],updated_site_tensors[2]}});
            auto updated_evolved_wf_overlap_dx=kagome_cirac_wf_distance_params_nblg->kagome_patch_RDM.expect_val_from_replaced_tensors({{updated_site_tensors_dx[0]*updated_bond_tensor,evolved_site_tensors[0]*evolved_bond_tensor},{updated_site_tensors_dx[1],evolved_site_tensors[1]},{updated_site_tensors_dx[2],evolved_site_tensors[2]}});

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
        //Print((kagome_cirac_wf_distance_sq_nblg_f(x_plus_dxi,params)-kagome_cirac_wf_distance_sq_nblg_f(x,params))/1E-10);
        //gsl_vector_set(x_plus_dxi,i,gsl_vector_get(x,i)+1E-12);
        //Print((kagome_cirac_wf_distance_sq_nblg_f(x_plus_dxi,params)-kagome_cirac_wf_distance_sq_nblg_f(x,params))/1E-12);
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
            auto updated_bond_tensor_dx=singlet_tensor_from_basis_params(kagome_cirac_wf_distance_params_nblg->bond_tensor_basis,bond_singlet_params_dx);

            auto updated_wf_norm_sq_dx=kagome_cirac_wf_distance_params_nblg->kagome_patch_RDM.expect_val_from_replaced_tensors({{updated_site_tensors[0]*updated_bond_tensor_dx,updated_site_tensors[0]*updated_bond_tensor},{updated_site_tensors[1],updated_site_tensors[1]},{updated_site_tensors[2],updated_site_tensors[2]}});
            auto updated_evolved_wf_overlap_dx=kagome_cirac_wf_distance_params_nblg->kagome_patch_RDM.expect_val_from_replaced_tensors({{updated_site_tensors[0]*updated_bond_tensor_dx,evolved_site_tensors[0]*evolved_bond_tensor},{updated_site_tensors[1],evolved_site_tensors[1]},{updated_site_tensors[2],evolved_site_tensors[2]}});

            //Print((updated_wf_norm_sq_dx-updated_wf_norm_sq)/dxi);
            //Print((updated_evolved_wf_overlap_dx-updated_evolved_wf_overlap)/dxi);

            dfi+=2*((updated_wf_norm_sq_dx-updated_wf_norm_sq)/dxi-(updated_evolved_wf_overlap_dx-updated_evolved_wf_overlap)/dxi).real();

        }

        //Print(i+N_leg_params);
        //Print(dfi);
        //gsl_vector *x_plus_dxi;
        //x_plus_dxi=gsl_vector_alloc(x->size);
        //gsl_vector_memcpy(x_plus_dxi,x);
        //gsl_vector_set(x_plus_dxi,i+N_leg_params,gsl_vector_get(x,i+N_leg_params)+1E-10);
        //Print((kagome_cirac_wf_distance_sq_nblg_f(x_plus_dxi,params)-kagome_cirac_wf_distance_sq_nblg_f(x,params))/1E-10);
        //gsl_vector_set(x_plus_dxi,i+N_leg_params,gsl_vector_get(x,i+N_leg_params)+1E-12);
        //Print((kagome_cirac_wf_distance_sq_nblg_f(x_plus_dxi,params)-kagome_cirac_wf_distance_sq_nblg_f(x,params))/1E-12);
        //gsl_vector_free(x_plus_dxi);

        gsl_vector_set(df,i+N_leg_params,dfi);
    }
}

void kagome_cirac_wf_distance_sq_nblg_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
    *f=kagome_cirac_wf_distance_sq_nblg_f(x,params);
    kagome_cirac_wf_distance_sq_nblg_df(x,params,df);
}

void kagome_cirac_wf_distance_sq_nblg_df_check(const gsl_vector *x, void *params, gsl_vector *df)
{
    int x_size=x->size;

    gsl_vector *x_plus_dx;
    x_plus_dx=gsl_vector_alloc(x_size);
    gsl_vector_memcpy(x_plus_dx,x);

    double f_x=kagome_cirac_wf_distance_sq_nblg_f(x,params);
    for (int i=0; i<x_size; i++)
    {
        double dxi=1E-10;
        gsl_vector_set(x_plus_dx,i,gsl_vector_get(x,i)+dxi);
        double f_x_plus_dxi=kagome_cirac_wf_distance_sq_nblg_f(x_plus_dx,params);
        gsl_vector_set(df,i,(f_x_plus_dxi-f_x)/dxi);
        gsl_vector_set(x_plus_dx,i,gsl_vector_get(x,i));
        //Print(gsl_vector_get(df,i));
    }

    gsl_vector_free(x_plus_dx);
}
