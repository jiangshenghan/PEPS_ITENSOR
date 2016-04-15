
#include "tnetwork.h"
#include "tnetwork_ctm.h"
#include "full_update_ctm.h"
#include "full_update.h"
#include "cluster_env.h"



void kagome_normal_rvb_fast_full_update(IQPEPS &kagome_rvb, const Evolution_Params &su_params)
{
    int Lx=kagome_rvb.lattice().n_uc()[0], 
        Ly=kagome_rvb.lattice().n_uc()[1];

    //using SIMPLE UPDATE environment
    /*
    //init params
    std::vector<int> cutting_sites={3*(Lx+1),3*(Lx+1)+2}, 
                     cutting_bonds={6*(Lx+1),6*(Lx+1)+1,6*(Lx+1)+4};
    std::vector<double> leg_gate_params;
    Singlet_Tensor_Basis cutting_site_basis(kagome_rvb.site_tensors(cutting_sites[0]).indices());

    std::array<std::vector<IQTensor>,2> su_env_mats;
    get_env_tensor_minimization(kagome_rvb.site_tensors(cutting_sites[0])*kagome_rvb.bond_tensors(cutting_bonds[0])*kagome_rvb.bond_tensors(cutting_bonds[1]),kagome_rvb.site_tensors(cutting_sites[1])*kagome_rvb.bond_tensors(cutting_bonds[2]),su_env_mats);
    //obtain env_tensors
    std::vector<IQTensor> env_tensors;
    obtain_kagome_normal_env_MPO(0,cutting_sites,cutting_bonds,su_env_mats,kagome_rvb,env_tensors);

    //init kagome RDM
    PEPSt_RDM<IQTensor> kagome_normal_rdm("two sites shape",cutting_sites,cutting_bonds,env_tensors,kagome_rvb);
    //Print(kagome_normal_rdm.RDM());
    //Print(kagome_normal_rdm.wf_norm_sq());

    //get a good init rvb state
    while (heisenberg_energy_from_RDM(kagome_normal_rdm)>-0.08)
    {
        Print(heisenberg_energy_from_RDM(kagome_normal_rdm));
        random_init_kagome_rvb_normal_peps(kagome_rvb);
        get_env_tensor_minimization(kagome_rvb.site_tensors(cutting_sites[0])*kagome_rvb.bond_tensors(cutting_bonds[0])*kagome_rvb.bond_tensors(cutting_bonds[1]),kagome_rvb.site_tensors(cutting_sites[1])*kagome_rvb.bond_tensors(cutting_bonds[2]),su_env_mats);
        obtain_kagome_normal_env_MPO(0,cutting_sites,cutting_bonds,su_env_mats,kagome_rvb,env_tensors);
        kagome_normal_rdm.update_peps_rdm(env_tensors,kagome_rvb);
    }

    //su_env_option=0 : simple update env, correspond to env_contract_seq={0,0,0,0,0,0}
    //su_env_option=1 : triangle shape, correspond to env_contract_seq={0,0,0,0,1}
    //su_env_option=2 : large tree shape, correspond to env_contract_seq={0,0,1}
    int su_env_option=1;
    std::vector<int> env_contract_seq={0,0,0,0,1};
    Print(su_env_option);

    // */

    //
    //using CLUSTER UPDATE
    // /*
    //init params
    std::vector<int> cutting_sites={3*(Lx+1),3*(Lx+1)+2}, 
                     cutting_bonds={6*(Lx+1),6*(Lx+1)+1,6*(Lx+1)+4};
    std::vector<double> leg_gate_params;
    Singlet_Tensor_Basis cutting_site_basis(kagome_rvb.site_tensors(cutting_sites[0]).indices());

    //su_env_option=0 : simple update env, correspond to env_contract_seq={0,0,0,0,0,0}
    //su_env_option=1 : triangle shape, correspond to env_contract_seq={0,0,0,0,1}
    //su_env_option=2 : large tree shape, correspond to env_contract_seq={0,0,1}
    int su_env_option=1;
    std::vector<int> env_contract_seq={0,0,0,0,1};
    Print(su_env_option);

    //init cluster direct product env
    Cluster_Env kagome_cluster_env(kagome_rvb.name()+"triangle shape",{0,2},{kagome_rvb.site_tensors(0),kagome_rvb.site_tensors(1)*kagome_rvb.bond_tensors(0)*kagome_rvb.bond_tensors(2),kagome_rvb.site_tensors(2)*kagome_rvb.bond_tensors(1)});
    kagome_cluster_env.obtain_sl_env_iterative_nodeg();
    std::array<std::vector<IQTensor>,2> su_env_mats;
    init_env_tensor(kagome_rvb.site_tensors(cutting_sites[0])*kagome_rvb.bond_tensors(cutting_bonds[0])*kagome_rvb.bond_tensors(cutting_bonds[1]),kagome_rvb.site_tensors(cutting_sites[1])*kagome_rvb.bond_tensors(cutting_bonds[2]),kagome_cluster_env.sl_env_tensor(),su_env_mats);
    //obtain env_tensors
    std::vector<IQTensor> env_tensors;
    obtain_kagome_normal_env_MPO(su_env_option,cutting_sites,cutting_bonds,su_env_mats,kagome_rvb,env_tensors);

    //init kagome RDM
    PEPSt_RDM<IQTensor> kagome_normal_rdm("two sites shape",cutting_sites,cutting_bonds,env_tensors,kagome_rvb);
    //Print(kagome_normal_rdm.RDM());
    //Print(kagome_normal_rdm.wf_norm_sq());

    //get a good init rvb state
    while (heisenberg_energy_from_RDM(kagome_normal_rdm)>-0.08)
    {
        Print(heisenberg_energy_from_RDM(kagome_normal_rdm));
        random_init_kagome_rvb_normal_peps(kagome_rvb);
        kagome_cluster_env.update_site_tensors({kagome_rvb.site_tensors(0),kagome_rvb.site_tensors(1)*kagome_rvb.bond_tensors(0)*kagome_rvb.bond_tensors(2),kagome_rvb.site_tensors(2)*kagome_rvb.bond_tensors(1)});
        kagome_cluster_env.obtain_sl_env_iterative_nodeg();
        init_env_tensor(kagome_rvb.site_tensors(cutting_sites[0])*kagome_rvb.bond_tensors(cutting_bonds[0])*kagome_rvb.bond_tensors(cutting_bonds[1]),kagome_rvb.site_tensors(cutting_sites[1])*kagome_rvb.bond_tensors(cutting_bonds[2]),kagome_cluster_env.sl_env_tensor(),su_env_mats);
        obtain_kagome_normal_env_MPO(su_env_option,cutting_sites,cutting_bonds,su_env_mats,kagome_rvb,env_tensors);
        kagome_normal_rdm.update_peps_rdm(env_tensors,kagome_rvb);
    }
    // */


    //
    //using FAST FULL UPDATE
    /*
    //init params
    //make sure we choose cutting sites right
    std::vector<int> cutting_sites={3*(Lx+1),3*(Lx+1)+2}, 
                     cutting_bonds={6*(Lx+1),6*(Lx+1)+1,6*(Lx+1)+4};
    std::vector<double> leg_gate_params;
    Singlet_Tensor_Basis cutting_site_basis(kagome_rvb.site_tensors(cutting_sites[0]).indices());

    //obtain good init state
    std::array<std::vector<IQTensor>,2> su_env_mats;
    get_env_tensor_minimization(kagome_rvb.site_tensors(cutting_sites[0])*kagome_rvb.bond_tensors(cutting_bonds[0])*kagome_rvb.bond_tensors(cutting_bonds[1]),kagome_rvb.site_tensors(cutting_sites[1])*kagome_rvb.bond_tensors(cutting_bonds[2]),su_env_mats);
    //obtain env_tensors
    std::vector<IQTensor> env_tensors;
    obtain_kagome_normal_env_MPO(0,cutting_sites,cutting_bonds,su_env_mats,kagome_rvb,env_tensors);

    //init kagome RDM
    PEPSt_RDM<IQTensor> kagome_normal_rdm("two sites shape",cutting_sites,cutting_bonds,env_tensors,kagome_rvb);
    //Print(kagome_normal_rdm.RDM());
    //Print(kagome_normal_rdm.wf_norm_sq());

    //get a good init rvb state
    while (heisenberg_energy_from_RDM(kagome_normal_rdm)>-0.08)
    {
        Print(heisenberg_energy_from_RDM(kagome_normal_rdm));
        random_init_kagome_rvb_normal_peps(kagome_rvb);
        get_env_tensor_minimization(kagome_rvb.site_tensors(cutting_sites[0])*kagome_rvb.bond_tensors(cutting_bonds[0])*kagome_rvb.bond_tensors(cutting_bonds[1]),kagome_rvb.site_tensors(cutting_sites[1])*kagome_rvb.bond_tensors(cutting_bonds[2]),su_env_mats);
        obtain_kagome_normal_env_MPO(0,cutting_sites,cutting_bonds,su_env_mats,kagome_rvb,env_tensors);
        kagome_normal_rdm.update_peps_rdm(env_tensors,kagome_rvb);
    }

    //init full update environment
    //options: 0:mode. If mode=1: chi only; if mode=2: cutoff only; if mode=3, both chi and cutoff
    //1: if mode=1: this is chi. if mode=2: cut off is 1.E-(this number). eg. this number=14. if mode=3, this is chi
    //2: number of update steps, 
    //[3,4]: starting_corr (e.g.: could be [1,1])
    //[5,6]: unit cell size (e.g. could be [2,1] for pi_flux)
    //7: svd_method, 1 if using Itensor svd, 2 if using lapack svd.
    iVec env_options="3,15,10,1,1,1,1,2";
    Print(env_options[1]);
    if (abs(kagome_psg::mu_12+1.)<EPSILON) env_options(5)=2;

    Tnetwork_Storage<IQTensor> kagome_rvb_tnetwork_storage=peps_to_tnetwork_storage(kagome_rvb);
    full_update_itebd<IQTensor> full_update_env_tensors(kagome_rvb_tnetwork_storage,env_options,"double layer X");
    auto env_tensors_temp=full_update_env_tensors.output_env_MPOs();
    //Print(env_tensors_temp);
    SzSz_measure(full_update_env_tensors);
    env_tensors.clear();
    for (int i=0; i<env_tensors_temp.size(); i++) env_tensors.push_back(env_tensors_temp(i));
    env_options(2)=2;

    //reinit kagome RDM use full update env
    kagome_normal_rdm=PEPSt_RDM<IQTensor>("two sites shape",cutting_sites,cutting_bonds,{env_tensors[0]*env_tensors[1],env_tensors[5],env_tensors[2]*env_tensors[3],env_tensors[4]},kagome_rvb,{0,1,0,1});
    //Print(kagome_normal_rdm.wf_norm_sq());
    // */



    for (int iter=0; iter<su_params.iter_nums; iter++)
    {
        //Init evolve gates
        IQTPO evolve_gate=trotter_gate_NN_Heisenberg({kagome_rvb.phys_legs(cutting_sites[0]),kagome_rvb.phys_legs(cutting_sites[1])},su_params.ts[iter]);

        //init leg gates for sites, which are used to approx evolve_gate
        //every leg gate is formed by two in legs and one out leg
        std::vector<Singlet_Tensor_Basis> leg_gates_basis;
        for (int evolvei=0; evolvei<cutting_sites.size(); evolvei++)
        {
            IndexSet<IQIndex> leg_gate_indices;
            leg_gate_indices.addindex(commonIndex(dag(kagome_rvb.site_tensors(cutting_sites[evolvei])),kagome_rvb.bond_tensors(cutting_bonds[1])));
            leg_gate_indices.addindex(commonIndex(dag(evolve_gate.site_tensors(evolvei)),evolve_gate.bond_tensors(0)));
            leg_gate_indices.addindex(commonIndex(kagome_rvb.site_tensors(cutting_sites[evolvei]),dag(kagome_rvb.bond_tensors(cutting_bonds[1]))).prime());
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
                //Print(spin_config);
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
            Print(leg_gate_params);
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

            //FAST FULL UPDATE
            /*
            //update env tensors
            kagome_rvb_tnetwork_storage=peps_to_tnetwork_storage(kagome_rvb);
            full_update_env_tensors.update_env(kagome_rvb_tnetwork_storage,env_options);
            env_tensors_temp=full_update_env_tensors.output_env_MPOs();
            env_tensors.clear();
            for (int i=0; i<env_tensors_temp.size(); i++) 
            {
                env_tensors.push_back(env_tensors_temp(i));
                Print(env_tensors.back().norm());
            }
            //Print((env_tensors[0]*env_tensors[1]).norm());
            //Print((env_tensors[2]*env_tensors[3]).norm());
            kagome_normal_rdm.update_peps_rdm({env_tensors[0]*env_tensors[1],env_tensors[5],env_tensors[2]*env_tensors[3],env_tensors[4]},kagome_rvb);
            Print(SzSz_measure(full_update_env_tensors));
            // */

            //SIMPLE UPDATE
            /*
            //update env tensors
            get_env_tensor_minimization(kagome_rvb.site_tensors(cutting_sites[0])*kagome_rvb.bond_tensors(cutting_bonds[0])*kagome_rvb.bond_tensors(cutting_bonds[1]),kagome_rvb.site_tensors(cutting_sites[1])*kagome_rvb.bond_tensors(cutting_bonds[2]),su_env_mats);
            obtain_kagome_normal_env_MPO(su_env_option,cutting_sites,cutting_bonds,su_env_mats,kagome_rvb,env_tensors);
            kagome_normal_rdm.update_peps_rdm(env_tensors,kagome_rvb,env_contract_seq);
            // */

            //CLUSTER UPDATE
            // /*
            //update env tensors
            kagome_cluster_env.update_site_tensors({kagome_rvb.site_tensors(0),kagome_rvb.site_tensors(1)*kagome_rvb.bond_tensors(0)*kagome_rvb.bond_tensors(2),kagome_rvb.site_tensors(2)*kagome_rvb.bond_tensors(1)});
            kagome_cluster_env.obtain_sl_env_iterative_nodeg();
            init_env_tensor(kagome_rvb.site_tensors(cutting_sites[0])*kagome_rvb.bond_tensors(cutting_bonds[0])*kagome_rvb.bond_tensors(cutting_bonds[1]),kagome_rvb.site_tensors(cutting_sites[1])*kagome_rvb.bond_tensors(cutting_bonds[2]),kagome_cluster_env.sl_env_tensor(),su_env_mats);
            obtain_kagome_normal_env_MPO(su_env_option,cutting_sites,cutting_bonds,su_env_mats,kagome_rvb,env_tensors);
            kagome_normal_rdm.update_peps_rdm(env_tensors,kagome_rvb,env_contract_seq);
            // */

            Print(heisenberg_energy_from_RDM(kagome_normal_rdm));

            //get leg_gate_params
            double cutoff=1e-5*su_params.ts[iter];
            obtain_kagome_normal_leg_gate_params_minimization(kagome_normal_rdm,evolve_gate,leg_gates_basis,leg_gate_params,cutoff);

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
            rotation_reflection_symmetrize_kagome_rvb_normal_site_tensor(updated_site_tens);
            //keep the same norm
            double site_norm=kagome_rvb.site_tensors(cutting_sites[0]).norm();
            updated_site_tens*=site_norm/(updated_site_tens.norm());
            Print(site_norm);
            Print((updated_site_tens-kagome_rvb.site_tensors(cutting_sites[0])).norm());
            kagome_rvb.generate_site_tensors({updated_site_tens,updated_site_tens,updated_site_tens});

            //Print site tensor params
            //for (int basei=0; basei<cutting_site_basis.dim(); basei++) 
            //{
            //    Complex site_param=(dag(cutting_site_basis[basei])*updated_site_tens).toComplex();
            //    Print(cutting_site_basis.spin_configs(basei));
            //    Print(cutting_site_basis.fusion_channel(basei));
            //    Print(site_param);
            //}

            //FAST FULL UPDATE
            /*
            //measure energy accurately
            //if ((step*su_params.ts[iter]-floor(step*su_params.ts[iter])<su_params.ts[iter]*0.01 && step%10==0) || step%100==0)
            if (step%20==0)
            {
                cout << "Measure energy accurately!" << endl;
                int env_chi=env_options(1),
                    env_step=env_options(2);
                env_options(2)=10;
                kagome_rvb_tnetwork_storage=peps_to_tnetwork_storage(kagome_rvb);
                full_update_env_tensors=full_update_itebd<IQTensor>(kagome_rvb_tnetwork_storage,env_options,"double layer X");
                full_update_env_tensors.output_env_MPOs();
                Print(SzSz_measure(full_update_env_tensors));
                env_options(1)=env_chi;
                env_options(2)=env_step;

                //std::stringstream ss;
                //ss << "/home/jiangsb/tn_ying/tensor_network/result/peps_storage/kagome_rvb_normal_D=" << kagome_rvb.D() << "_Lx=" << kagome_rvb.n_uc()[0] << "_Ly=" << kagome_rvb.n_uc()[1] << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_fast_full_chi=" << env_options[1] << "_iter=" << iter << "_step=" << step; 
                //std::string file_name=ss.str();
                //writeToFile(file_name,kagome_rvb);
            }
            // */

            //stores as PEPS 
            if ((step+1)*10%su_params.steps_nums[iter]==0)
            {
                std::stringstream ss;
                //SIMPLE UPDATE
                //ss << "/home/jiangsb/tn_ying/tensor_network/result/peps_storage/kagome_rvb_normal_D=" << kagome_rvb.D() << "_Lx=" << kagome_rvb.n_uc()[0] << "_Ly=" << kagome_rvb.n_uc()[1] << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_su_env_option=" << su_env_option << "_iter=" << iter << "_step=" << step; 
                //CLUSTER UPDATE
                ss << "/home/jiangsb/tn_ying/tensor_network/result/peps_storage/kagome_rvb_normal_D=" << kagome_rvb.D() << "_Lx=" << kagome_rvb.n_uc()[0] << "_Ly=" << kagome_rvb.n_uc()[1] << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_cluster_su_env_option=" << su_env_option << "_iter=" << iter << "_step=" << step; 
                //FAST FULL UPDATE
                //ss << "/home/jiangsb/tn_ying/tensor_network/result/peps_storage/kagome_rvb_normal_D=" << kagome_rvb.D() << "_Lx=" << kagome_rvb.n_uc()[0] << "_Ly=" << kagome_rvb.n_uc()[1] << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_fast_full_chi=" << env_options[1] << "_iter=" << iter << "_step=" << step; 
                std::string file_name=ss.str();
                writeToFile(file_name,kagome_rvb);
            }

        }//end of trotter steps

    }//end of trotter iters

}


bool obtain_kagome_normal_leg_gate_params_minimization(PEPSt_RDM<IQTensor> &kagome_normal_rdm, const IQTPO &evolve_gate, const std::vector<Singlet_Tensor_Basis> &leg_gates_basis, std::vector<double> &leg_gate_params, double cutoff)
{
    //construct site and bond tensors after applied by evolve gate
    std::vector<IQTensor> site_tensors_evolved;
    //evolve_legs_combiners are used to combine legs of peps site tensors and evolve gate site tensors
    std::vector<IQCombiner> evolve_legs_combiners;
    for (int cuti=0; cuti<kagome_normal_rdm.cutting_sites_no(); cuti++)
    {
        site_tensors_evolved.push_back(kagome_normal_rdm.cutting_site_tensors(cuti)*evolve_gate.site_tensors(cuti));
        site_tensors_evolved[cuti].noprime();
        auto evolve_gate_virt_leg=commonIndex(evolve_gate.site_tensors(cuti),evolve_gate.bond_tensors(0));
        IQIndex cutting_virt_leg=commonIndex(kagome_normal_rdm.cutting_site_tensors(cuti),kagome_normal_rdm.cutting_bond_tensors(1));
        evolve_legs_combiners.push_back(IQCombiner(cutting_virt_leg,evolve_gate_virt_leg));
    }
    auto bond_tensor_evolved=kagome_normal_rdm.cutting_bond_tensors(1)*evolve_gate.bond_tensors(0);
    //Print(site_tensors_evolved);
    //Print(bond_tensor_evolved);

    //obtain evolved_wf_norm
    IQTensor evolve_gate_tensor=evolve_gate.bond_tensors(0);
    for (const auto &tens : evolve_gate.site_tensors()) evolve_gate_tensor*=tens;
    Complex evolved_wf_norm_sq=(kagome_normal_rdm.RDM()*(evolve_gate_tensor*dag(swapPrime(evolve_gate_tensor,0,2))).mapprime(2,1)).toComplex();
    //Print(evolved_wf_norm);

    //obtain distance between wf and wf_evolved
    Complex wf_norm_sq=kagome_normal_rdm.wf_norm_sq(),
           wf_evolved_wf_overlap=2.*(kagome_normal_rdm.RDM()*evolve_gate_tensor).toComplex(),
           //wf_evolved_wf_distance=std::sqrt(2.*pow(evolved_wf_norm,2.)-evolved_wf_norm/wf_norm*wf_evolved_wf_overlap);
           //TODO:normalize evolved_wf to get distance
           wf_evolved_wf_distance_sq=wf_norm_sq+evolved_wf_norm_sq-wf_evolved_wf_overlap;
    //Print(wf_norm_sq);
    //Print(wf_evolved_wf_distance_sq);

    //using conjugate gradient methods to minimize distance square between updated_wf and evolved_wf
    int find_min_status;
    int iter=0, max_iter=30;

    const gsl_multimin_fdfminimizer_type *minimize_T;
    gsl_multimin_fdfminimizer *s;

    //params to do minimization
    Kagome_Normal_Wf_Distance_Params_FU *kagome_normal_wf_distance_params_fu=new Kagome_Normal_Wf_Distance_Params_FU(evolved_wf_norm_sq,kagome_normal_rdm,site_tensors_evolved,bond_tensor_evolved,evolve_legs_combiners,leg_gates_basis);

    //x stores coefficient for site leg gates and plaquette leg gates
    gsl_vector *x;
    gsl_multimin_function_fdf wf_distance_sq_func;

    wf_distance_sq_func.n=leg_gate_params.size();
    wf_distance_sq_func.f=kagome_normal_wf_distance_sq_fu_f;
    wf_distance_sq_func.df=kagome_normal_wf_distance_sq_fu_df;
    wf_distance_sq_func.fdf=kagome_normal_wf_distance_sq_fu_fdf;
    wf_distance_sq_func.params=kagome_normal_wf_distance_params_fu;

    x=gsl_vector_alloc(wf_distance_sq_func.n);
    for (int i=0; i<leg_gate_params.size(); i++) gsl_vector_set(x,i,leg_gate_params[i]);

    //minimize_T=gsl_multimin_fdfminimizer_conjugate_fr;
    minimize_T=gsl_multimin_fdfminimizer_vector_bfgs2;
    s=gsl_multimin_fdfminimizer_alloc(minimize_T,wf_distance_sq_func.n);
    gsl_multimin_fdfminimizer_set(s,&wf_distance_sq_func,x,0.1,0.1);
    do
    {
        iter++;
        find_min_status=gsl_multimin_fdfminimizer_iterate(s);
        if (find_min_status) break;
        find_min_status=gsl_multimin_test_gradient(s->gradient,cutoff);

        Print(iter);
        Print(s->f);
        Print(wf_evolved_wf_distance_sq/wf_norm_sq);

        //Print(find_min_status==GSL_SUCCESS);
        //Print(find_min_status==GSL_CONTINUE);
        //Print(find_min_status==GSL_FAILURE);
    }
    //while ((find_min_status==GSL_CONTINUE || sqrt(s->f)>0.3*wf_evolved_wf_distance) && iter<max_iter);
    while (find_min_status==GSL_CONTINUE && iter<max_iter);
    //Print(cutoff);


    //Print(iter);
    //Print(s->f);

    //Print(kagome_normal_rdm.wf_norm_sq());
    Print(sqrt((wf_evolved_wf_distance_sq/wf_norm_sq).real()));
    Print(sqrt(kagome_normal_wf_distance_sq_fu_f(s->x,kagome_normal_wf_distance_params_fu)));

    //if (iter==max_iter)
    //{
    //    cout << "Leg gate is not good enough, may be trapped in local minima!" << endl << "try smaller time step!" << endl;
    //    gsl_multimin_fdfminimizer_free(s);
    //    gsl_vector_free(x);
    //    delete kagome_normal_wf_distance_params_fu;

    //    return false;
    //}

    for (int i=0; i<leg_gate_params.size(); i++)
        leg_gate_params[i]=gsl_vector_get(s->x,i);

    Print(leg_gate_params);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
    delete kagome_normal_wf_distance_params_fu;

    return true;
}


double kagome_normal_wf_distance_sq_fu_f(const gsl_vector *x, void *params)
{
    std::vector<double> leg_gate_params;
    Kagome_Normal_Wf_Distance_Params_FU *kagome_normal_wf_distance_params_fu=(Kagome_Normal_Wf_Distance_Params_FU *)params;
    int N_leg_params=kagome_normal_wf_distance_params_fu->leg_gates_basis[0].dim();

    for (int i=0; i<N_leg_params; i++) 
        leg_gate_params.push_back(gsl_vector_get(x,i));
    //Print(leg_gate_params);

    //construct leg gates
    std::vector<IQTensor> leg_gates;
    for (const auto &leg_gate_basis : kagome_normal_wf_distance_params_fu->leg_gates_basis)
    {
        leg_gates.push_back(singlet_tensor_from_basis_params(leg_gate_basis,leg_gate_params));
        //PrintDat(leg_gates);
    }

    auto evolved_site_tensors=kagome_normal_wf_distance_params_fu->evolved_site_tensors;
    auto evolved_bond_tensor=kagome_normal_wf_distance_params_fu->evolved_bond_tensor;

    //combine virtual legs of evolved tensors
    const auto &evolve_legs_combiners=kagome_normal_wf_distance_params_fu->evolve_legs_combiners;
    for (int cuti=0; cuti<2; cuti++)
    {
        evolved_site_tensors[cuti]=evolved_site_tensors[cuti]*evolve_legs_combiners[cuti];
        evolved_bond_tensor=evolved_bond_tensor*dag(evolve_legs_combiners[cuti]);
    }
    //Print(evolved_site_tensors);
    //Print(evolved_bond_tensor);

    //obtain updated_site_tensors
    std::vector<IQTensor> updated_site_tensors;
    for (int cuti=0; cuti<2; cuti++)
    {
        updated_site_tensors.push_back(kagome_normal_wf_distance_params_fu->evolved_site_tensors[cuti]*leg_gates[cuti]);
        updated_site_tensors[cuti].noprime();
    }

    //get cut bond tensors
    auto &kagome_normal_rdm=kagome_normal_wf_distance_params_fu->kagome_normal_rdm;
    const std::vector<IQTensor> &bond_tensors=kagome_normal_rdm.cutting_bond_tensors();

    //Print(updated_site_tensors);

    Complex wf_norm_sq=kagome_normal_wf_distance_params_fu->kagome_normal_rdm.wf_norm_sq();
    Complex updated_wf_norm_sq=kagome_normal_rdm.expect_val_from_replaced_tensors(updated_site_tensors);
    double updated_evolved_wf_overlap_rel=2.*(kagome_normal_rdm.expect_val_from_replaced_tensors({updated_site_tensors,evolved_site_tensors},{bond_tensors,{bond_tensors[0],evolved_bond_tensor,bond_tensors[2]}})/wf_norm_sq).real();
    Complex distance_sq_rel=updated_wf_norm_sq/wf_norm_sq+kagome_normal_wf_distance_params_fu->evolved_wf_norm_sq/wf_norm_sq-updated_evolved_wf_overlap_rel;
    //Print(kagome_normal_wf_distance_params_fu->evolved_wf_norm_sq/wf_norm_sq);
    //Print(updated_wf_norm_sq);
    //Print(updated_evolved_wf_overlap);

    //Print(distance_sq_rel);

    return distance_sq_rel.real();
}

void kagome_normal_wf_distance_sq_fu_df(const gsl_vector *x, void *params, gsl_vector *df)
{
    std::vector<double> leg_gate_params;
    Kagome_Normal_Wf_Distance_Params_FU *kagome_normal_wf_distance_params_fu=(Kagome_Normal_Wf_Distance_Params_FU *)params;
    int N_leg_params=kagome_normal_wf_distance_params_fu->leg_gates_basis[0].dim();

    for (int i=0; i<N_leg_params; i++)
        leg_gate_params.push_back(gsl_vector_get(x,i));
    //Print(leg_gate_params);

    //obtain evolved tensors with combined virtual legs
    auto evolved_site_tensors=kagome_normal_wf_distance_params_fu->evolved_site_tensors;
    auto evolved_bond_tensor=kagome_normal_wf_distance_params_fu->evolved_bond_tensor;
    const auto &evolve_legs_combiners=kagome_normal_wf_distance_params_fu->evolve_legs_combiners;
    for (int cuti=0; cuti<2; cuti++)
    {
        evolved_site_tensors[cuti]=evolved_site_tensors[cuti]*evolve_legs_combiners[cuti];
        evolved_bond_tensor=evolved_bond_tensor*dag(evolve_legs_combiners[cuti]);
    }


    //init leg gates
    std::vector<IQTensor> leg_gates;
    for (const auto &leg_gate_basis : kagome_normal_wf_distance_params_fu->leg_gates_basis)
    {
        leg_gates.push_back(singlet_tensor_from_basis_params(leg_gate_basis,leg_gate_params));
    }

    //init updated_site_tensors
    std::vector<IQTensor> updated_site_tensors;
    for (int cuti=0; cuti<2; cuti++)
    {
        updated_site_tensors.push_back(kagome_normal_wf_distance_params_fu->evolved_site_tensors[cuti]*leg_gates[cuti]);
        updated_site_tensors[cuti].noprime();
    }

    //get cut bond tensors
    auto &kagome_normal_rdm=kagome_normal_wf_distance_params_fu->kagome_normal_rdm;
    const std::vector<IQTensor> &bond_tensors=kagome_normal_rdm.cutting_bond_tensors();

    Complex wf_norm_sq=kagome_normal_wf_distance_params_fu->kagome_normal_rdm.wf_norm_sq();
    //init updated_wf_norm_sq and updated_evolved_wf_overlap 
    Complex updated_wf_norm_sq_rel=kagome_normal_rdm.expect_val_from_replaced_tensors(updated_site_tensors)/wf_norm_sq;
    //here updated_evolved_wf_overlap=\langle\psi|\phi\rangle
    Complex updated_evolved_wf_overlap_rel=kagome_normal_wf_distance_params_fu->kagome_normal_rdm.expect_val_from_replaced_tensors({updated_site_tensors,evolved_site_tensors},{bond_tensors,{bond_tensors[0],evolved_bond_tensor,bond_tensors[2]}})/wf_norm_sq;

    //Print(updated_wf_norm_sq);
    //Print(updated_evolved_wf_overlap);

    double dxi=1E-3;
    //get derivative for leg_gate_params
    const auto &leg_gates_basis=kagome_normal_wf_distance_params_fu->leg_gates_basis;
    double df_sq=0;
    for (int i=0; i<N_leg_params; i++)
    {
        Complex dfi=0;
        std::vector<double> leg_gate_params_dx=leg_gate_params;
        leg_gate_params_dx[i]=leg_gate_params[i]+dxi;
        auto leg_gates_dx=leg_gates;
        auto updated_site_tensors_dx=updated_site_tensors;

        for (int cuti=0; cuti<2; cuti++)
        {
            leg_gates_dx[cuti]=singlet_tensor_from_basis_params(leg_gates_basis[cuti],leg_gate_params_dx);
            updated_site_tensors_dx[cuti]=kagome_normal_wf_distance_params_fu->evolved_site_tensors[cuti]*leg_gates_dx[cuti];
            updated_site_tensors_dx[cuti].noprime();

            //calculate time
            std::chrono::time_point<std::chrono::system_clock> start, middle, end;
            start=std::chrono::system_clock::now();
            auto updated_wf_norm_sq_dx_rel=kagome_normal_rdm.expect_val_from_replaced_tensors({updated_site_tensors_dx,updated_site_tensors})/wf_norm_sq;
            middle=std::chrono::system_clock::now();
            auto updated_evolved_wf_overlap_dx_rel=kagome_normal_wf_distance_params_fu->kagome_normal_rdm.expect_val_from_replaced_tensors({updated_site_tensors_dx,evolved_site_tensors},{bond_tensors,{bond_tensors[0],evolved_bond_tensor,bond_tensors[2]}})/wf_norm_sq;
            end=std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed0=middle-start, elapsed1=end-middle;
            //cout << "time to get derivative:" << endl;
            //Print(elapsed0.count());
            //Print(elapsed1.count());

            dfi+=(updated_wf_norm_sq_dx_rel-updated_wf_norm_sq_rel)/dxi-(updated_evolved_wf_overlap_dx_rel-updated_evolved_wf_overlap_rel)/dxi;

            //recover leg_gates_dx and updated_site_tensors_dx
            leg_gates_dx[cuti]=leg_gates[cuti];
            updated_site_tensors_dx[cuti]=updated_site_tensors[cuti];

            //Print((updated_wf_norm_sq_dx-updated_wf_norm_sq)/dxi);
            //Print((updated_evolved_wf_overlap_dx-updated_evolved_wf_overlap)/dxi);
        }

        //Print(i);
        //Print(dfi.real()*2.);
        //gsl_vector *x_plus_dxi;
        //x_plus_dxi=gsl_vector_alloc(x->size);
        //gsl_vector_memcpy(x_plus_dxi,x);
        //gsl_vector_set(x_plus_dxi,i,gsl_vector_get(x,i)+1E-10);
        //Print((kagome_normal_wf_distance_sq_fu_f(x_plus_dxi,params)-kagome_normal_wf_distance_sq_fu_f(x,params))/1E-10);
        //gsl_vector_set(x_plus_dxi,i,gsl_vector_get(x,i)+1E-12);
        //Print((kagome_normal_wf_distance_sq_fu_f(x_plus_dxi,params)-kagome_normal_wf_distance_sq_fu_f(x,params))/1E-12);
        //gsl_vector_free(x_plus_dxi);

        gsl_vector_set(df,i,dfi.real()*2.);

        df_sq+=pow(std::abs(dfi),2.);
    }
    Print(sqrt(df_sq));
}

void kagome_normal_wf_distance_sq_fu_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
    *f=kagome_normal_wf_distance_sq_fu_f(x,params);
    kagome_normal_wf_distance_sq_fu_df(x,params,df);
}


void obtain_kagome_normal_env_MPO(int env_option, std::vector<int> cutting_sites, std::vector<int> cutting_bonds, const std::array<std::vector<IQTensor>,2> &su_env_mats, const IQPEPS &kagome_rvb, std::vector<IQTensor> &env_tensors)
{
    env_tensors.clear();

    //simple update env (double layer)
    if (env_option==0)
    {
        for (int i=0; i<2; i++)
        {
            for (const auto &single_layer_env_mat : su_env_mats[i])
            {
                env_tensors.push_back(swapPrime(single_layer_env_mat,1,2)*prime(single_layer_env_mat).dag());
            }
        }
        return;
    }

    //obtain double layer env_mat, the primed leg has Out dir and the umprimed leg has In dir
    IQTensor env_mat=su_env_mats[0][0];
    if ((env_mat.indices()[0].dir()==In && env_mat.indices()[0].primeLevel()==1) || (env_mat.indices()[0].dir()==Out && env_mat.indices()[0].primeLevel()==0)) env_mat.dag();
    env_mat=swapPrime(env_mat,1,2)*prime(env_mat).dag();

    //triangle shape
    if (env_option==1)
    {
        //obtain env tensors connect to one leg
        for (int i=0; i<2; i++)
        {
            for (const auto &single_layer_env_mat : su_env_mats[i])
            {
                if (commonIndex(single_layer_env_mat,kagome_rvb.bond_tensors(cutting_bonds[0]))!=IQIndex::Null()) continue;
                if (commonIndex(single_layer_env_mat,kagome_rvb.bond_tensors(cutting_bonds[0]+2))!=IQIndex::Null()) continue;
                env_tensors.push_back(swapPrime(single_layer_env_mat,1,2)*prime(single_layer_env_mat).dag());
            }
        }

        //obtain env tensor connect two tensors
        std::vector<IQIndex> indices_contract;
        indices_contract.push_back(commonIndex(kagome_rvb.site_tensors(cutting_sites[0]+1),kagome_rvb.bond_tensors(cutting_bonds[0]+3)));
        Coordinate bond_coord_temp=kagome_rvb.lattice().bond_list_to_coord(cutting_bonds[0]);
        bond_coord_temp[0]-=1;
        bond_coord_temp[2]+=5;
        indices_contract.push_back(commonIndex(kagome_rvb.site_tensors(cutting_sites[0]+1),kagome_rvb.bond_tensors(bond_coord_temp)));

        env_tensors.push_back(env_site_combined_tensor(env_mat,kagome_rvb.site_tensors(cutting_sites[0]+1)*kagome_rvb.bond_tensors(cutting_bonds[0]+2),indices_contract));
        //Print(env_tensors);

        return;
    }

    //tree shape
    if (env_option==2)
    {
        //obtain the six "leaf" tensors
        Coordinate coord=kagome_rvb.lattice().site_list_to_coord(cutting_sites[0]);
        std::vector<IQTensor> leaf_tensors(6);
        std::vector<Coordinate> leaf_coords(6);

        leaf_coords[0]={coord[0]-1,coord[1],0};
        leaf_coords[1]={coord[0],coord[1]-1,0};
        leaf_coords[2]={coord[0]+1,coord[1]-1,0};
        leaf_coords[3]={coord[0]+1,coord[1],0};
        leaf_coords[4]={coord[0],coord[1]+1,0};
        leaf_coords[5]={coord[0]-1,coord[1]+1,0};

        for (int leafi=0; leafi<leaf_tensors.size(); leafi++)
        {
            std::vector<IQTensor> uc_site_tensors(3);
            uc_site_tensors[0]=kagome_rvb.site_tensors(leaf_coords[leafi])*kagome_rvb.bond_tensors(leaf_coords[leafi])*kagome_rvb.bond_tensors({leaf_coords[leafi][0],leaf_coords[leafi][1],1});
            uc_site_tensors[1]=kagome_rvb.site_tensors({leaf_coords[leafi][0],leaf_coords[leafi][1],1})*kagome_rvb.bond_tensors({leaf_coords[leafi][0],leaf_coords[leafi][1],2});
            uc_site_tensors[2]=kagome_rvb.site_tensors({leaf_coords[leafi][0],leaf_coords[leafi][1],2});

            if (leafi==0 || leafi==5) leaf_tensors[leafi]=env_kagome_normal_uc_combined_tensor(env_mat,uc_site_tensors,2); 
            if (leafi==1 || leafi==2) leaf_tensors[leafi]=env_kagome_normal_uc_combined_tensor(env_mat,uc_site_tensors,1); 
            if (leafi==3 || leafi==4) leaf_tensors[leafi]=env_kagome_normal_uc_combined_tensor(env_mat,uc_site_tensors,0); 
        }

        //obtain env for cutting_sites[0]
        IQTensor connect_bond_tensor;
        IQTensor left_env_tens;
        connect_bond_tensor=kagome_rvb.bond_tensors({coord[0]-1,coord[1]-1,5});
        left_env_tens=leaf_tensors[0]*(connect_bond_tensor*dag(connect_bond_tensor).prime())*leaf_tensors[1];
        connect_bond_tensor=kagome_rvb.bond_tensors({coord[0]-1,coord[1],4});
        left_env_tens*=connect_bond_tensor*dag(connect_bond_tensor).prime();
        connect_bond_tensor=kagome_rvb.bond_tensors({coord[0],coord[1]-1,3});
        left_env_tens*=connect_bond_tensor*dag(connect_bond_tensor).prime();
        //Print(left_env_tens.indices());
        env_tensors.push_back(left_env_tens);

        //obtain env for cutting_sites[1]
        IQTensor right_env_tens;
        connect_bond_tensor=kagome_rvb.bond_tensors({coord[0]+1,coord[1]-1,3});
        right_env_tens=leaf_tensors[2]*(connect_bond_tensor*dag(connect_bond_tensor).prime())*leaf_tensors[3];
        connect_bond_tensor=kagome_rvb.bond_tensors({coord[0],coord[1]-1,5});
        right_env_tens*=connect_bond_tensor*dag(connect_bond_tensor).prime();
        //Print(right_env_tens.indices());
        env_tensors.push_back(right_env_tens);

        //obtain env connect two tensors
        IQTensor lr_env_tens;
        connect_bond_tensor=kagome_rvb.bond_tensors({coord[0]-1,coord[1]+1,4});
        lr_env_tens=leaf_tensors[4]*(connect_bond_tensor*dag(connect_bond_tensor).prime())*leaf_tensors[5];
        connect_bond_tensor=kagome_rvb.bond_tensors({coord[0]-1,coord[1],5});
        lr_env_tens*=connect_bond_tensor*dag(connect_bond_tensor).prime();
        connect_bond_tensor=kagome_rvb.bond_tensors({coord[0],coord[1],3});
        lr_env_tens*=connect_bond_tensor*dag(connect_bond_tensor).prime();
        lr_env_tens=lr_env_tens*kagome_rvb.site_tensors(cutting_sites[0]+1)*dag(kagome_rvb.site_tensors(cutting_sites[0]+1)).prime(Link);
        lr_env_tens*=kagome_rvb.bond_tensors(cutting_bonds[0]+2)*prime(kagome_rvb.bond_tensors(cutting_bonds[0]+2)).dag();
        //Print(lr_env_tens.indices());
        env_tensors.push_back(lr_env_tens);
    }
}


IQTensor env_site_combined_tensor(const IQTensor &env_mat, const IQTensor &site_tensor, const std::vector<IQIndex> &indices_contract, bool contract_phys_ind)
{
    std::vector<IQTensor> env_mats_contract;
    for (const auto &indice : indices_contract)
    {
        IQTensor env_mat_contract(dag(indice),prime(indice));
        if (dag(indice).dir()==env_mat.indices()[0].dir())
        {
            tensor_assignment(env_mat_contract,env_mat);
        }
        else
        {
            tensor_assignment(env_mat_contract,dag(env_mat));
        }
        env_mats_contract.push_back(env_mat_contract);
    }

    IQTensor combined_tensor=site_tensor;
    for (const auto &env_mat_temp : env_mats_contract) combined_tensor*=env_mat_temp;
    if (contract_phys_ind)
    {
        return combined_tensor*dag(site_tensor).prime(Link);
    }
    else
    {
        return combined_tensor*dag(site_tensor).prime();
    }

    return IQTensor();
}


IQTensor env_kagome_normal_uc_combined_tensor(const IQTensor &env_mat, const std::vector<IQTensor> &site_tensors, int bulk_no)
{
    std::vector<std::vector<IQIndex>> indices_contract(site_tensors.size());
    std::vector<IQTensor> env_site_combined_tensors(site_tensors.size());
    IQTensor all_boundary_tensor;

    //get boundary tensors combined with env_tensor
    for (int i=0; i<site_tensors.size(); i++)
    {
        if (i==bulk_no) continue;
        for (const auto &indice : site_tensors[i].indices())
        {
            if (indice.type()==Site) continue;
            bool is_comm_indice=false;
            for (int j=0; j<site_tensors.size(); j++)
            {
                if (i==j) continue;
                if (hasindex(site_tensors[j],indice))
                {
                    is_comm_indice=true;
                    break;
                }
            }
            if (!is_comm_indice) indices_contract[i].push_back(indice);
        }

        env_site_combined_tensors[i]=env_site_combined_tensor(env_mat,site_tensors[i],indices_contract[i]);

        if (all_boundary_tensor.valid()==false) 
            all_boundary_tensor=env_site_combined_tensors[i];
        else
            all_boundary_tensor*=env_site_combined_tensors[i];
    }

    return (all_boundary_tensor*site_tensors[bulk_no]*dag(site_tensors[bulk_no]).prime(Link));
}
