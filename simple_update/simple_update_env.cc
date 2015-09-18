
#include "simple_update_env.h"

const IQTensor &nondeg_spin_sym_env(const IQTensor &tens_A, const IQTensor &tens_B)
{
    IQIndex comm_leg_A=commonIndex(tens_A,dag(tens_B)),
            comm_leg_B=dag(comm_leg_A);

       assert(comm_leg_A!=IQIndex::Null());

    IQCombiner iqcomb_A, iqcomb_B;

    //for (const auto &left_A : tens_A.indices())
    //{
    //    if (left_A==comm_leg_A) continue;
    //    iqcomb_A.addleft(left_A);
    //}
    ////The direction of combined leg is opposite to the remaining leg
    ////then the "cgtensor"is identity
    ////iqcomb_A.init("comb_leg_A",All,-comm_leg_A.dir());
    //IQTensor mat_A=iqcomb_A*tens_A;

    //for (const auto &left_B : tens_B.indices())
    //{
    //    if (left_B==comm_leg_B) continue;
    //    iqcomb_B.addleft(left_B);
    //}
    ////iqcomb_B.init("comb_leg_B",All,-comm_leg_B.dir());
    //IQTensor mat_B=iqcomb_B*tens_B;


    IQTensor tens_A_dag=dag(tens_A).prime(dag(comm_leg_A),1),
             tens_B_dag=dag(tens_B).prime(dag(comm_leg_B),2);

    //env_tens equals to square root of norm_A*norm_B
    IQTensor env_tens=tens_A_dag*tens_A*tens_B*tens_B_dag;
    IQIndex env_leg1=prime(dag(comm_leg_A),1),
            env_leg2=prime(dag(comm_leg_B),2);

    for (int i=1; i<comm_leg_A.m(); i++)
    {
        for(int j=1; j<comm_leg_A.m(); j++)
        {
            auto elem=env_tens(env_leg1(i),env_leg2(j));
            if (i!=j) 
            {
                assert(std::abs(elem)<EPSILON);
                continue;
            }
                
            env_tens(env_leg1(i),env_leg2(j))=std::sqrt(elem);
        }
    }

    return env_tens;
}
