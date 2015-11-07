
#include "simple_update_env.h"

void get_env_tensor(const IQTensor &site_tensA, const IQTensor &site_tensB, std::array<std::vector<IQTensor>,2> &env_tens)
{
    //env_tens is initialized according to spin and flavor qn
    auto comm_leg=commonIndex(site_tensA,site_tensB);
    if (env_tens[0].empty())
    {
        std::vector<int> virt_leg_flavor_deg;
        std::vector<Spin_Basis> virt_leg_spin_basis;
        iqind_spin_rep(comm_leg,virt_leg_flavor_deg);
        iqind_to_spin_basis(comm_leg,virt_leg_flavor_deg,virt_leg_spin_basis);
        std::vector<double> env_tens_init_diag(comm_leg.m());
        std::vector<bool> diag_elems_inited(comm_leg.m(),false);
        for (int vali=0; vali<comm_leg.m(); vali++)
        {
            if (diag_elems_inited[vali]) continue;
            diag_elems_inited[vali]=true;
            env_tens_init_diag[vali]=rand_gen();
            for (int valj=vali+1; valj<comm_leg.m(); valj++)
            {
                if (diag_elems_inited[valj]) continue;
                if (virt_leg_spin_basis[vali].S()==virt_leg_spin_basis[valj].S() && virt_leg_spin_basis[vali].t()==virt_leg_spin_basis[valj].t())
                {
                    env_tens_init_diag[valj]=env_tens_init_diag[vali];
                    diag_elems_inited[valj]=true;
                }
            }
        }
        //Print(virt_leg_spin_basis);
        //Print(env_tens_init_diag);

        for (const auto &leg : site_tensA.indices())
        {
            if (leg.type()==Site) continue;
            if (leg==comm_leg) continue;
            IQTensor temp_tens(dag(leg),prime(leg));
            for (int val=1; val<=leg.m(); val++) temp_tens(dag(leg)(val),prime(leg)(val))=env_tens_init_diag[val-1];
            env_tens[0].push_back(temp_tens);
        }

        for (const auto &leg : site_tensB.indices())
        {
            if (leg.type()==Site) continue;
            if (leg==comm_leg) continue;
            IQTensor temp_tens(dag(leg),prime(leg));
            for (int val=1; val<=leg.m(); val++) temp_tens(dag(leg)(val),prime(leg)(val))=env_tens_init_diag[val-1];
            env_tens[1].push_back(temp_tens);
        }
    }

    //for (const auto &tensor : env_tens[0]) PrintDat(tensor);
    //for (const auto &tensor : env_tens[1]) PrintDat(tensor);

    int n_out_legs=env_tens[0].size();
    assert(n_out_legs==env_tens[1].size());

    //get env_mat iteratively
    int iter=0, max_iter=100;
    while (iter<max_iter)
    {
        IQTensor site_env_tensA=site_tensA,
                 site_env_tensB=site_tensB;
        for (int legi=0; legi<n_out_legs; legi++)
        {
            site_env_tensA*=env_tens[0][legi];
            site_env_tensB*=env_tens[1][legi];
        }
        site_env_tensA.noprime();
        site_env_tensB.noprime();
        site_env_tensA.clean();
        site_env_tensB.clean();

        auto env_tens_diag=nondeg_spin_sym_env_updated(site_env_tensA,site_env_tensB);

        //Print(iter);
        //PrintDat(site_env_tensA);
        //PrintDat(site_env_tensB);
        //Print(env_tens_diag);

        IndexSet<IQIndex> new_env_indset=env_tens[0][0].indices();
        IQTensor new_env_tensor(new_env_indset[0],new_env_indset[1]);

        for (int val=0; val<env_tens_diag.size(); val++)
        {
            new_env_tensor(new_env_indset[0](val+1),new_env_indset[1](val+1))=env_tens_diag[val];
        }

        //PrintDat(new_env_tensor);
        Print((new_env_tensor-env_tens[0][0]).norm());

        if ((new_env_tensor-env_tens[0][0]).norm()<1E-3) break;

        for (int legi=0; legi<n_out_legs; legi++)
        {
            tensor_assignment(env_tens[0][legi],new_env_tensor);
            tensor_assignment(env_tens[1][legi],new_env_tensor);
        }

        iter++;
    }

    if (iter==max_iter)
    {
        cout << "Unable to find good environment by iteration! Redo the iteration!" << endl;
        //exit(EXIT_FAILURE);
        for (auto &tensors : env_tens)
            tensors.clear();
        get_env_tensor(site_tensA,site_tensB,env_tens);
    }

}


std::vector<double> nondeg_spin_sym_env_updated(const IQTensor &tens_A, const IQTensor &tens_B)
{
    IQIndex comm_leg_A=commonIndex(tens_A,dag(tens_B)),
            comm_leg_B=dag(comm_leg_A);

    assert(comm_leg_A!=IQIndex::Null());

    IQTensor tens_A_dag=dag(tens_A).prime(dag(comm_leg_A)),
             tens_B_dag=dag(tens_B).prime(dag(comm_leg_B));

    IQTensor tens_A_norm_square=(tens_A_dag*tens_A).takeRealPart(),
             tens_B_norm_square=(tens_B_dag*tens_B).takeRealPart();

    //PrintDat(tens_A_norm_square);
    //PrintDat(tens_B_norm_square);

    std::vector<double> env_diag;

    for (int i=1; i<=comm_leg_A.m(); i++)
    {
        double normA=sqrt(tens_A_norm_square(comm_leg_A(i),dag(prime(comm_leg_A(i))))),
               normB=sqrt(tens_B_norm_square(comm_leg_B(i),dag(prime(comm_leg_B(i)))));
        env_diag.push_back(normA*normB);
    }

    //we fix the norm of env_tens_diag to be sqrt(dim)
    double env_norm=0;
    for (auto elem : env_diag) env_norm+=elem*elem;
    env_norm=sqrt(env_norm);
    double env_amp=sqrt(env_diag.size()*1.)/env_norm;
    for (auto &elem : env_diag) elem*=env_amp;
    //Print(env_diag);

    return env_diag;
}
