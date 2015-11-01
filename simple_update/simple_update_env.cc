
#include "simple_update_env.h"

void get_env_tensor(const IQTensor &site_tensA, const IQTensor &site_tensB, std::array<std::vector<IQTensor>,2> &env_tens)
{
    //init env_tens
    if (env_tens[0].empty())
    {
        auto comm_leg=commonIndex(site_tensA,site_tensB);
        for (const auto &leg : site_tensA.indices())
        {
            if (leg.type()==Site) continue;
            if (leg==comm_leg) continue;
            IQIndex temp_tens(dag(leg),prime(leg));
            for (int val=1; val<=leg.m(); val++) temp_tens(dag(leg)(val),prime(leg)(val))=1.;
            env_tens[0].push_back(temp_tens);
        }
    }

    if (env_tens[1].empty())
    {
        auto comm_leg=commonIndex(site_tensA,site_tensB);
        for (const auto &leg : site_tensB.indices())
        {
            if (leg.type()==Site) continue;
            if (leg==comm_leg) continue;
            IQIndex temp_tens(dag(leg),prime(leg));
            for (int val=1; val<=leg.m(); val++) temp_tens(dag(leg)(val),prime(leg)(val))=1.;
            env_tens[1].push_back(temp_tens);
        }
    }

    int n_out_legs=env_tens[0].size();
    assert(n_out_legs==env_tens[1].size());

    //get env_mat iteratively
    int iter=0, max_iter=10;
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

        auto env_tens_diag=nondeg_spin_sym_env_updated(site_env_tensA,site_env_tensB);
        IndexSet<IQIndex> new_env_indset=env_tens[0][0].indices();
        IQTensor new_env_tensor(new_env_indset[0],new_env_indset[1]);

        for (int val=0; val<env_tens_diag.size(); val++)
        {
            new_env_tensor(new_env_indset[0](val+1),new_env_indset[1](val+1))=env_tens_diag[val];
        }

        Print(iter);
        PrintDat(new_env_tensor);
        Print((new_env_tensor-env_tens[0][0]).norm());

        for (int legi=0; legi<n_out_legs; legi++)
        {
            tensor_assignment(env_tens[0][legi],new_env_tensor);
            tensor_assignment(env_tens[1][legi],new_env_tensor);
        }

        iter++;
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

    std::vector<double> env_diag;

    for (int i=1; i<=comm_leg_A.m(); i++)
    {
        double normA=sqrt(tens_A_norm_square(comm_leg_A(i),dag(prime(comm_leg_A(i))))),
               normB=sqrt(tens_B_norm_square(comm_leg_B(i),dag(prime(comm_leg_B(i)))));
        env_diag.push_back(normA*normB);
    }

    //we fix the norm of env_tens_diag to be sqrt(dim)
    double env_norm=0;
    std::accumulate(env_diag.begin(),env_diag.end(),env_norm,
            [](double norm, double elem){norm+=elem*elem});
    env_norm=sqrt(env_norm);
    double env_amp=sqrt(env_diag.size()*1.)/env_norm;
    for (auto &elem : env_diag) elem*=env_amp;

    return env_diag;
}
