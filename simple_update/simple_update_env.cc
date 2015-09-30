
#include "simple_update_env.h"

void get_env_tensor(const IQTensor &site_tensA, const IQTensor &site_tensB, std::array<std::vector<IQTensor>,2> &env_tens)
{
    int n_leg=env_tens[0].size();
    std::array<std::vector<IndexSet<IQIndex>>,2> env_tens_indices;

    for (int i=0; i<2; i++)
    {
        for (const auto &tensor:env_tens[i])
        {
            env_tens_indices[i].push_back(tensor.indices());
        }
    }

    //get env_mat iteratively
    int iter=0, max_iter=10;
    while (iter<max_iter)
    {
        IQTensor site_env_tensA=site_tensA,
                 site_env_tensB=site_tensB;
        for (int i=0; i<n_leg; i++)
        {
            site_env_tensA*=env_tens[0][i];
            site_env_tensB*=env_tens[1][i];
        }
        site_env_tensA.noprime();
        site_env_tensB.noprime();

        auto env_tens_diag=nondeg_spin_sym_env_updated(site_env_tensA,site_env_tensB);

        int leg_dim=env_tens_diag.size();
        for (int i=0; i<2; i++)
        {
            for (int j=0; j<n_leg; j++)
            {
                for (int k=0; k<leg_dim; k++)
                {
                    env_tens[i][j](env_tens_indices[i][j][0](k+1),env_tens_indices[i][j][1](k+1))=env_tens_diag[k];
                }
            }
        }
        
        iter++;
    }

}


std::vector<double> nondeg_spin_sym_env_updated(const IQTensor &tens_A, const IQTensor &tens_B)
{
    IQIndex comm_leg_A=commonIndex(tens_A,dag(tens_B)),
            comm_leg_B=dag(comm_leg_A);

       assert(comm_leg_A!=IQIndex::Null());

    IQCombiner iqcomb_A, iqcomb_B;


    IQTensor tens_A_dag=dag(tens_A).prime(dag(comm_leg_A),1),
             tens_B_dag=dag(tens_B).prime(dag(comm_leg_B),2);

    //env_tens_square equals to square root of norm_A*norm_B
    IQTensor env_tens_square=tens_A_dag*tens_A*tens_B*tens_B_dag;
    IQIndex env_leg1=prime(dag(comm_leg_A),1),
            env_leg2=prime(dag(comm_leg_B),2);

    std::vector<double> env_diag;

    for (int i=1; i<comm_leg_A.m(); i++)
    {

#ifndef NDEBUG
        //check for off-diagonal elements
        for(int j=1; j<comm_leg_A.m(); j++)
        {
            auto elem=env_tens_square(env_leg1(i),env_leg2(j));
            if (i!=j) 
            {
                assert(std::abs(elem)<EPSILON);
                continue;
            }
        }
#endif

        env_diag.push_back(std::sqrt(env_tens_square(env_leg1(i),env_leg2(i))));
    }

    return env_diag;
}
