
#include "tensor_update.h"


void kagome_normal_rvb_tensor_update(IQPEPS &kagome_rvb)
{
    //init params
    int Lx=kagome_rvb.lattice().n_uc()[0],
        Ly=kagome_rvb.lattice().n_uc()[1];
    std::vector<int> cutting_sites={3*(Lx+1),3*(Lx+1)+2}, 
                     cutting_bonds={6*(Lx+1),6*(Lx+1)+1,6*(Lx+1)+4};


    //symmetric basis is for cutting_sites[1] 
    Singlet_Tensor_Basis cutting_site_basis(kagome_rvb.site_tensors(cutting_sites[0]).indices());
    std::vector<IQTensor> symmetric_site_basis;
    obtain_kagome_rvb_normal_site_tensor_symmetric_basis(cutting_site_basis,symmetric_site_basis);
    //we should set basis such that projection basis get real numbers
    IQPEPS random_kagome_rvb(kagome_rvb.lattice(),kagome_rvb.indexset());
    random_init_kagome_rvb_normal_peps(random_kagome_rvb);
    for (auto &site_base : symmetric_site_basis)
    {
        Complex site_param=(dag(site_base)*random_kagome_rvb.site_tensors(cutting_sites[0])).toComplex();
        if (std::abs(std::imag(site_param))>std::abs(std::real(site_param))) site_base*=Complex_i;
    }

    std::vector<double> site_params;
    for(const auto &site_base : symmetric_site_basis)
    {
        Complex param=(dag(site_base)*kagome_rvb.site_tensors(cutting_sites[0])).toComplex();
        site_params.push_back(param.real());
    }


    //init symmetric basis for another cutting site
    std::vector<IQTensor> symmetric_site_basis_v2=symmetric_site_basis;
    for (auto &site_base : symmetric_site_basis_v2)
    {
        for (int indi=0; indi<site_base.r(); indi++)
        {
            IQIndex oind=site_base.indices()[indi],
                    nind=kagome_rvb.site_tensors(cutting_sites[1]).indices()[indi];
            site_base.replaceIndex(oind,nind);
        }
    }

    std::vector<std::vector<IQTensor>> cutting_sites_symmetric_basis{symmetric_site_basis,symmetric_site_basis_v2};


    //init environment
    //options: 0:mode. If mode=1: chi only; if mode=2: cutoff only; if mode=3, both chi and cutoff
    //1: if mode=1: this is chi. if mode=2: cut off is 1.E-(this number). eg. this number=14. if mode=3, this is chi
    //2: number of update steps, 
    //[3,4]: starting_corr (e.g.: could be [1,1])
    //[5,6]: unit cell size (e.g. could be [2,1] for pi_flux)
    //7: svd_method, 1 if using Itensor svd, 2 if using lapack svd.
    iVec env_options="3,40,5,1,1,1,1,2";
    cout << "chi=" << env_options[1] << endl;
    if (abs(kagome_psg::mu_12+1.)<EPSILON) env_options(5)=2;

    Tnetwork_Storage<IQTensor> kagome_rvb_tnetwork_storage=peps_to_tnetwork_storage(kagome_rvb);
    full_update_itebd<IQTensor> full_update_env_tensors(kagome_rvb_tnetwork_storage,env_options,"double layer X");
    auto env_tensors_temp=full_update_env_tensors.output_env_MPOs();
    //Print(env_tensors_temp);
    SzSz_measure(full_update_env_tensors);
    std::vector<IQTensor> env_tensors;
    for (int envi=0; envi<env_tensors_temp.size(); envi++) env_tensors.push_back(env_tensors_temp(envi));
    env_options(2)=2;

    //init kagome RDM use full update env
    PEPSt_RDM<IQTensor> kagome_normal_rdm("two sites shape",cutting_sites,cutting_bonds,{env_tensors[0]*env_tensors[1],env_tensors[5],env_tensors[2]*env_tensors[3],env_tensors[4]},kagome_rvb,{0,1,0,1});
    //Print(kagome_normal_rdm.wf_norm_sq());


    int iter=0, max_iter=1000;
    double tol=0.2,
           max_step_size=4e-3,
           step_size=max_step_size;
    while (iter<max_iter)
    {
        Print(iter);

        std::vector<double> energy_deriv=heisenberg_energy_derivative(kagome_normal_rdm,cutting_sites_symmetric_basis,site_params);
        double energy_deriv_norm=0.;
        for (double deriv : energy_deriv) energy_deriv_norm+=deriv*deriv;
        energy_deriv_norm=sqrt(energy_deriv_norm);
        Print(energy_deriv_norm);
        Print(energy_deriv);

        double energy=heisenberg_energy_from_RDM(kagome_normal_rdm),
               updated_energy=energy+1;
        IQPEPS updated_kagome_rvb=kagome_rvb;
        PEPSt_RDM<IQTensor> updated_kagome_normal_rdm=kagome_normal_rdm;

        //choose step size such that tensor does not changes a lot
        if (step_size>max_step_size) step_size=max_step_size;
        do
        {
            IQTensor updated_site_tensor=kagome_normal_rdm.cutting_site_tensors(0);
            for (int basei=0; basei<site_params.size(); basei++)
            {
                const auto &site_base=cutting_sites_symmetric_basis[0][basei];
                updated_site_tensor-=step_size*energy_deriv[basei]/energy_deriv_norm*site_base;
            }
            updated_kagome_rvb.generate_site_tensors({updated_site_tensor,updated_site_tensor,updated_site_tensor});
            updated_kagome_normal_rdm.update_RDM({env_tensors[0]*env_tensors[1],env_tensors[5],env_tensors[2]*env_tensors[3],env_tensors[4]},updated_kagome_rvb);
            updated_energy=heisenberg_energy_from_RDM(updated_kagome_normal_rdm);

            //Print(step_size);
            //Print(energy);
            //Print(updated_energy);

            if (updated_energy<energy) break;
            step_size*=tol;
        }
        while (updated_energy>energy);

        //update site_params, kagome_rvb, env_tensors and kagome_normal_rdm
        for (int parami=0; parami<site_params.size(); parami++) site_params[parami]-=step_size*energy_deriv[parami]/energy_deriv_norm;
        kagome_rvb=updated_kagome_rvb;

        kagome_rvb_tnetwork_storage=peps_to_tnetwork_storage(kagome_rvb);
        full_update_env_tensors.update_env(kagome_rvb_tnetwork_storage,env_options);
        env_tensors_temp=full_update_env_tensors.output_env_MPOs();
        env_tensors.clear();
        for (int envi=0; envi<env_tensors_temp.size(); envi++) 
        {
            env_tensors.push_back(env_tensors_temp(envi));
            Print(env_tensors.back().norm());
        }
        Print(SzSz_measure(full_update_env_tensors));

        kagome_normal_rdm.update_RDM({env_tensors[0]*env_tensors[1],env_tensors[5],env_tensors[2]*env_tensors[3],env_tensors[4]},kagome_rvb);

        updated_energy=heisenberg_energy_from_RDM(kagome_normal_rdm);

        Print(step_size);
        Print(energy);
        Print(updated_energy);
        Print(site_params);

        iter++;

        if (updated_energy<energy)
            step_size*=2.;
        else
            step_size*=tol;
        if (step_size<5e-8)
        {
            cout << "bad init state!" << endl;
            exit(0);
        }

        //update_environment accurately
        if (iter%10==0)
        {
            cout << "update environment accurately!" << endl;
            int env_step=env_options(2);
            env_options(2)=5;
            kagome_rvb_tnetwork_storage=peps_to_tnetwork_storage(kagome_rvb);
            full_update_env_tensors.update_env(kagome_rvb_tnetwork_storage,env_options);
            env_tensors_temp=full_update_env_tensors.output_env_MPOs();
            env_tensors.clear();
            for (int envi=0; envi<env_tensors_temp.size(); envi++) 
            {
                env_tensors.push_back(env_tensors_temp(envi));
                Print(env_tensors.back().norm());
            }
            Print(SzSz_measure(full_update_env_tensors));
            kagome_normal_rdm.update_RDM({env_tensors[0]*env_tensors[1],env_tensors[5],env_tensors[2]*env_tensors[3],env_tensors[4]},kagome_rvb);
            env_options(2)=env_step;
        }
    }
}



std::vector<double> heisenberg_energy_derivative(PEPSt_RDM<IQTensor> &kagome_normal_rdm, const std::vector<std::vector<IQTensor>> &cutting_sites_symmetric_basis, const std::vector<double> &site_params)
{
    NN_Heisenberg_Hamiltonian heisenberg_gate({kagome_normal_rdm.cutting_phys_legs(0),kagome_normal_rdm.cutting_phys_legs(1)});
    IQTensor hamiltonian=heisenberg_gate.site_tensors(0)*heisenberg_gate.bond_tensor()*heisenberg_gate.site_tensors(1);
    Complex normsq=kagome_normal_rdm.wf_norm_sq(),
            energy=(kagome_normal_rdm.RDM()*hamiltonian).toComplex();

    //Print(normsq);
    //Print(energy);
    //Print(energy/normsq);

    //evolve_legs_combiners are used to combine legs of peps site tensors and evolve gate site tensors
    std::vector<IQCombiner> evolve_legs_combiners;
    std::vector<IQTensor> site_tensors_evolved;
    for (int cuti=0; cuti<kagome_normal_rdm.cutting_sites_no(); cuti++)
    {
        IQIndex evolve_gate_virt_leg=commonIndex(heisenberg_gate.site_tensors(cuti),heisenberg_gate.bond_tensor()),
                cutting_virt_leg=commonIndex(kagome_normal_rdm.cutting_site_tensors(cuti),kagome_normal_rdm.cutting_bond_tensors(1));
        evolve_legs_combiners.push_back(IQCombiner(cutting_virt_leg,evolve_gate_virt_leg));
        site_tensors_evolved.push_back(kagome_normal_rdm.cutting_site_tensors(cuti)*heisenberg_gate.site_tensors(cuti)*evolve_legs_combiners[cuti]);
        site_tensors_evolved[cuti].noprime();
    }
    std::vector<IQTensor> bond_tensors_evolved=kagome_normal_rdm.cutting_bond_tensors();
    bond_tensors_evolved[1]=bond_tensors_evolved[1]*heisenberg_gate.bond_tensor()*dag(evolve_legs_combiners[0])*dag(evolve_legs_combiners[1]);
    //Print(site_tensors_evolved);
    //Print(bond_tensors_evolved);

    double dxi=1e-3;
    std::vector<double> df;
    for (int parami=0; parami<site_params.size(); parami++)
    {
        Complex normsq_dfi=0, energy_dfi=0;

        //contribution from cutting_sites[0]
        std::vector<IQTensor> site_tensors_dx={kagome_normal_rdm.cutting_site_tensors(0)+dxi*cutting_sites_symmetric_basis[0][parami],kagome_normal_rdm.cutting_site_tensors(1)};
        normsq_dfi+=(kagome_normal_rdm.expect_val_from_replaced_tensors({site_tensors_dx,kagome_normal_rdm.cutting_site_tensors()})-normsq)/dxi;
        energy_dfi+=(kagome_normal_rdm.expect_val_from_replaced_tensors({site_tensors_dx,site_tensors_evolved},{kagome_normal_rdm.cutting_bond_tensors(),bond_tensors_evolved})-energy)/dxi;

        //contribution from cutting_sites[1]
        std::vector<IQTensor> site_tensors_dx_v2={kagome_normal_rdm.cutting_site_tensors(0),kagome_normal_rdm.cutting_site_tensors(1)+dxi*cutting_sites_symmetric_basis[1][parami]};
        normsq_dfi+=(kagome_normal_rdm.expect_val_from_replaced_tensors({site_tensors_dx_v2,kagome_normal_rdm.cutting_site_tensors()})-normsq)/dxi;
        energy_dfi+=(kagome_normal_rdm.expect_val_from_replaced_tensors({site_tensors_dx_v2,site_tensors_evolved},{kagome_normal_rdm.cutting_bond_tensors(),bond_tensors_evolved})-energy)/dxi;

        //factor 2. counts for contribution from ket
        df.push_back((-energy/(normsq*normsq)*normsq_dfi+energy_dfi/normsq).real()*2.);
    }

    return df;
}
