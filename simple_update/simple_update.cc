
#include "simple_update.h"


void spin_square_peps_simple_update(IQPEPS &square_peps, const Evolution_Params &square_su_params)
{
    //Initialize trotter gate
    std::array<IQIndex,2> site01_legs{{square_peps.phys_legs(0),square_peps.phys_legs(1)}};
    NN_Heisenberg_Trotter_Gate evolve_gate(site01_legs);

    //Initialize env_tens
    //env_tens[i] stores environment of site i, which is direct product of matrices for the legs of site i 
    //there is no env matrix between A and B
    std::array<std::vector<IQTensor>,2> env_tens;

    int comm_bond=square_peps.lattice().comm_bond(0,1);
    auto comm_bond_tensor=square_peps.bond_tensors(comm_bond);
    for (int i=0; i<2; i++)
    {
        for (const auto &leg : square_peps.site_tensors(i).indices())
        {
            if (leg.type()==Site) continue;
            if (hasindex(comm_bond_tensor,dag(leg))) continue;

            //we set the initial value of env as diagonal matrix
            auto leg_dag=dag(leg);
            auto leg_prime=prime(leg);
            IQTensor leg_tensor(dag(leg),prime(leg));
            for (int val=1; val<=leg.m(); val++)
            {
                leg_tensor(leg_dag(val),leg_prime(val))=1.;
            }

            env_tens[i].push_back(leg_tensor);
        }
    }


    //leg gates is used to approx evolve_gate, which formed by three indices, two in one out (primed). 
    //two in legs are contracting to site tensor after applying trotter gate
    //one out leg is to contract comm bond tensor
    std::array<IndexSet<IQIndex>,2> leg_gates_indices;
    std::array<Singlet_Tensor_Basis,2> leg_gates_basis;
    for (int i=0; i<2; i++)
    {
        leg_gates_indices[i].addindex(commonIndex(dag(square_peps.site_tensors(i)),comm_bond_tensor));
        leg_gates_indices[i].addindex(commonIndex(dag(evolve_gate.site_tensors(i)),evolve_gate.bond_tensors(0)));
        leg_gates_indices[i].addindex(commonIndex(square_peps.site_tensors(i),dag(comm_bond_tensor)).prime());

        leg_gates_basis[i]=Singlet_Tensor_Basis(leg_gates_indices[i]);
    }

    //leg_gates_for_site0 is used to update site_tensors[0]
    std::vector<IQTensor> leg_gates_site0;
    auto indice_from_evolve_gate=commonIndex(dag(evolve_gate.site_tensors(0)),evolve_gate.bond_tensors(0));
    for (const auto &leg_indice : square_peps.site_tensors(0).indices())
    {
        if (leg_indice.type()==Site) continue;

        std::vector<IQIndex> leg_gates_site0_indices{{dag(leg_indice),indice_from_evolve_gate,prime(leg_indice)}};

        leg_gates_site0.push_back(IQTensor(leg_gates_site0_indices));
    }

    //leg_gate_params stores tunable parameters for the leg gate.
    std::vector<double> leg_gate_params(leg_gates_basis[0].dim(),1);
    
    for (int iter=0; iter<square_su_params.iter_nums; iter++)
    {
        evolve_gate.change_time(square_su_params.ts[iter]);

        for (int step=0; step<square_su_params.steps_nums[iter]; step++)
        {
            get_env_tensor(square_peps.site_tensors(0)*comm_bond_tensor,square_peps.site_tensors(1),env_tens);

            std::array<IQTensor,2> site_env_tens{{square_peps.site_tensors(0),square_peps.site_tensors(1)}};
            
            for (int sitei=0; sitei<2; sitei++)
            {
                for (const auto &env_leg_tensor : env_tens[sitei])
                {
                    site_env_tens[sitei]*=env_leg_tensor;
                }
                site_env_tens[sitei].noprime();
            }

            obtain_spin_sym_leg_gates_params(site_env_tens,comm_bond_tensor,evolve_gate,leg_gates_basis,leg_gate_params);


            //using leg_gate_params to generate all leg gates
            auto leg_gate_sample=singlet_tensor_from_basis_params(leg_gates_basis[0],leg_gate_params);

            for (auto &gate : leg_gates_site0)
            {
                auto gate_indices=gate.indices();
                gate=leg_gate_sample;

                int legi=0;
                for (const auto &sample_indice : leg_gate_sample.indices())
                {
                    gate.replaceIndex(sample_indice,gate_indices[legi]);
                    legi++;
                }
            }


            //updated site_tensors[0]
            auto &site_tens0=square_peps.site_tensors(0);

            int legi=0;
            for (const auto &leg_gate : leg_gates_site0)
            {
                site_tens0*=evolve_gate.site_tensors(0)*leg_gate;
                site_tens0.noprime();
                legi++;
            }

            //symmetrize site_tens0, and using site_tens0 to generate all site tensors of peps
            C4_symmetrized_tensor(site_tens0);
            square_peps.generate_site_tensors({site_tens0});

        }//trotter steps
    }//finish all simple update
}


void obtain_spin_sym_leg_gates_params(const std::array<IQTensor,2> &site_tensors, const IQTensor &bond_tensor, Trotter_Gate &evolve_gate, const std::array<Singlet_Tensor_Basis,2> &leg_gates_basis, std::vector<double> &leg_gate_params)
{
    //get time evloved tensors
    std::array<IQTensor,2> site_tensors_evolved(site_tensors);
    IQTensor bond_tensor_evolved(bond_tensor);

    for (int i=0; i<2; i++)
    {
        site_tensors_evolved[i]*=evolve_gate.site_tensors(i);
        site_tensors_evolved[i].noprime();

    }

    bond_tensor_evolved*=evolve_gate.bond_tensors(0);

    double evolved_tensor_norm=(site_tensors_evolved[0]*bond_tensor_evolved*site_tensors_evolved[1]).norm();


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
    }


    //fix another gate (gate v), and get updated site1
    int N_basis=leg_gates_basis[0].dim();
    arma::Col<double> vec_params(leg_gate_params);
    //matN_{ij}=\langle T_i|T_j\rangle
    //vecb_i=\langle T_i|\psi_0\rangle
    //where T_i=site0_updated_basis[i]*bond_tensor*site1_updated
    //\psi_0 is the wavefunction after time evolving
    //Then we get 
    //(abs(|\psi\rangle-|\psi_0\rangle))^2=\sum_{i,j}p_ip_j.matN_{ij}-\sum_ip_i.(vecb_i+vecb_i^*)+\langle\psi_0|\psi_0\rangle
    arma::Mat<double> matN(N_basis,N_basis);
    arma::Col<double> vecb(N_basis);

    int iter=0, max_iter=10;
    while (iter<max_iter)
    {
        IQTensor gate_v=singlet_tensor_from_basis_params(leg_gates_basis[1],vec_params);
        IQTensor site1_updated=site_tensors_evolved[1]*gate_v;
        IQTensor site1_updated_dag=dag(site1_updated);
        site1_updated.noprime();


        for (int i=0; i<N_basis; i++)
        {
            auto b_i=(site_tensors_evolved[0]*bond_tensor_evolved)*(site0_updated_basis_dag[i]*bond_tensor_dag)*(site_tensors_evolved[1]*site1_updated_dag);
            //we expect the result is a real number although the tensor is complex
            vecb(i)=b_i.toComplex().real();
            for (int j=0; j<N_basis; j++)
            {
                auto N_ij=(site0_updated_basis_dag[i]*bond_tensor_dag)*(site0_updated_basis[j]*bond_tensor)*(site1_updated_dag*site1_updated);
                //we expect the result is a real number although the tensor is complex
                matN(i,j)=N_ij.toComplex().real();
            }
        }

        double diff=std::sqrt(arma::dot(vec_params,matN*vec_params)-2*arma::dot(vec_params,vecb)+evolved_tensor_norm*evolved_tensor_norm);

        cout << "Step " << iter << ": distance of psi and psi_0 is " << diff << endl;


        vec_params=arma::inv(matN)*vecb;

        iter++;
    }

    for (int i=0; i<N_basis; i++)
    {
        leg_gate_params[i]=vec_params(i);
    }

}


void C4_symmetrized_tensor(IQTensor &site_tensor)
{
    IQIndex phys_indice;
    IndexSet<IQIndex> virt_indices;
    std::vector<int> virt_indices_dim;

    for (const auto &indice : site_tensor.indices())
    {
        if (indice.type()==Site) 
        {
            phys_indice=indice;
            continue;
        }

        virt_indices.addindex(indice);
        virt_indices_dim.push_back(indice.m());
    }

    int total_virt_dim=virt_indices.dim(),
        total_virt_legs_num=virt_indices.r(); //In this case, legs num = 4

    for (int phys_val_num=0; phys_val_num<phys_indice.m(); phys_val_num++)
    {
        auto phys_indice_val=phys_indice(phys_val_num+1);
            
        std::vector<bool> val_finished(total_virt_dim,false);

        for (int virt_val_num=0; virt_val_num<total_virt_dim; virt_val_num++)
        {
            std::vector<int> virt_val_list=list_from_num(virt_val_num,virt_indices_dim);

            if (val_finished[virt_val_num]) continue;
            val_finished[virt_val_num]=true;

            std::vector<IQIndexVal> virt_indices_val{{virt_indices[0](virt_val_list[0]+1),virt_indices[1](virt_val_list[1]+1),virt_indices[2](virt_val_list[2]+1),virt_indices[3](virt_val_list[3]+1)}};

            //check sz quantum number of indval
            if (div(site_tensor)!=(phys_indice_val.qn()+virt_indices_val[0].qn()+virt_indices_val[1].qn()+virt_indices_val[2].qn()+virt_indices_val[3].qn())) 
                continue;

            std::vector<Complex> rotation_related_elems(total_virt_legs_num);
            auto site_tensor_real_part=realPart(site_tensor);
            auto site_tensor_imag_part=imagPart(site_tensor);

            rotation_related_elems[0]=
                site_tensor_real_part(phys_indice_val,virt_indices_val[0],virt_indices_val[1],virt_indices_val[2],virt_indices_val[3])+
                Complex_i*site_tensor_imag_part(phys_indice_val,virt_indices_val[0],virt_indices_val[1],virt_indices_val[2],virt_indices_val[3]);
            for (int i=1; i<total_virt_legs_num; i++)
            {
                std::rotate(virt_val_list.begin(),virt_val_list.begin()+1,virt_val_list.end());
                val_finished[num_from_list(virt_val_list,virt_indices_dim)]=true;

                rotation_related_elems[i]=
                    site_tensor_real_part(phys_indice_val,virt_indices_val[0],virt_indices_val[1],virt_indices_val[2],virt_indices_val[3])+
                    Complex_i*site_tensor_imag_part(phys_indice_val,virt_indices_val[0],virt_indices_val[1],virt_indices_val[2],virt_indices_val[3]);
                rotation_related_elems[i]/=chi_c4*Theta_c4;
            }
            //rotate back to original val_list
            std::rotate(virt_val_list.begin(),virt_val_list.begin()+1,virt_val_list.end());

            Complex symmetrized_elem=0;
            for (const auto &elem : rotation_related_elems)
            {
                symmetrized_elem+=elem;
            }
            symmetrized_elem/=total_virt_legs_num*1.;


            for (int i=0; i<total_virt_legs_num; i++)
            {
                std::rotate(virt_val_list.begin(),virt_val_list.begin()+1,virt_val_list.end());

                auto ratio=(symmetrized_elem/rotation_related_elems[i]).real();

                site_tensor(phys_indice_val,virt_indices_val[0],virt_indices_val[1],virt_indices_val[2],virt_indices_val[3])*=ratio;
            }

        }
    }
}

