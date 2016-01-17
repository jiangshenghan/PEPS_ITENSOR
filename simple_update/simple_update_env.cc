
#include "simple_update_env.h"

void get_env_tensor_minimization(const IQTensor &site_tensA, const IQTensor &site_tensB, std::array<std::vector<IQTensor>,2> &env_tens)
{
    if (env_tens[0].empty())
    {
        init_env_tensor(site_tensA,site_tensB,env_tens);
    }

    //using conjugate gradient minimization to minimize distance square between updated env_tens and original env_tens
    int find_min_status;
    int iter=0, max_iter=100;

    const gsl_multimin_fdfminimizer_type *minimize_T;
    gsl_multimin_fdfminimizer *s;

    //params to do minimization, including
    auto comm_leg=commonIndex(site_tensA,site_tensB);
    std::vector<int> flavor_deg;
    std::vector<Spin_Basis> spin_basis;
    std::vector<int> flavor_accumulate_deg;

    iqind_spin_rep(comm_leg,flavor_deg);
    iqind_to_spin_basis(comm_leg,flavor_deg,spin_basis);

    flavor_accumulate_deg.push_back(0);
    for (int spini=1; spini<=flavor_deg.size();spini++)
    {
        flavor_accumulate_deg.push_back(flavor_accumulate_deg[spini-1]+flavor_deg[spini-1]);
    }

    //Print(flavor_deg);
    //Print(flavor_accumulate_deg);
    //Print(spin_basis);

    Env_Tens_Params *updated_env_tens_diff_params=new Env_Tens_Params{flavor_deg,flavor_accumulate_deg,spin_basis,{site_tensA,site_tensB},env_tens}; 

    //x stores elems of env_tens according to spin and flavor
    //TODO:consider more about flavor_deg!
    gsl_vector *x;
    gsl_multimin_function_fdf updated_env_tens_diff_func;

    //init x
    int x_size=flavor_accumulate_deg.back();
    x=gsl_vector_alloc(x_size);
    auto env_elems=env_elems_from_env_tens(flavor_accumulate_deg,spin_basis,env_tens);
    for (int i=0; i<x_size; i++) gsl_vector_set(x,i,env_elems[i]);

    //init updated_env_tens_diff_func
    updated_env_tens_diff_func.n=x_size;
    updated_env_tens_diff_func.f=updated_env_tens_diff_f;
    updated_env_tens_diff_func.df=updated_env_tens_diff_df;
    updated_env_tens_diff_func.fdf=updated_env_tens_diff_fdf;
    updated_env_tens_diff_func.params=updated_env_tens_diff_params;

    minimize_T=gsl_multimin_fdfminimizer_conjugate_fr;
    s=gsl_multimin_fdfminimizer_alloc(minimize_T,x_size);

    gsl_multimin_fdfminimizer_set(s,&updated_env_tens_diff_func,x,0.1,0.1);

    do
    {
        iter++;
        find_min_status=gsl_multimin_fdfminimizer_iterate(s);
        if (find_min_status) break;
        find_min_status=gsl_multimin_test_gradient(s->gradient,1E-3);

        //Print(iter);
        //Print(s->f);
    }
    while (find_min_status=GSL_CONTINUE && iter<max_iter);
    
    //Print(iter);
    //Print(s->f);
    //Print(x);
    Print(updated_env_tens_diff_f(s->x,updated_env_tens_diff_params));

    //if the result is much larger than zero, we will redo the minimization
    if (s->f>1E-3)
    {
        gsl_multimin_fdfminimizer_free(s);
        gsl_vector_free(x);
        delete updated_env_tens_diff_params;

        env_tens[0].clear();
        env_tens[1].clear();
        get_env_tensor_minimization(site_tensA,site_tensB,env_tens);

        return;
    }

    std::vector<double> env_elems_final;
    for (int i=0; i<s->x->size; i++) env_elems_final.push_back(gsl_vector_get(s->x,i));
    Print(env_elems_final);
    obtain_env_tens_from_env_elems(flavor_accumulate_deg,spin_basis,env_elems_final,env_tens);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
    delete updated_env_tens_diff_params;
}


void init_env_tensor(const IQTensor &site_tensA, const IQTensor &site_tensB, std::array<std::vector<IQTensor>,2> &env_tens)
{
    //env_tens is initialized according to spin and flavor qn
    auto comm_leg=commonIndex(site_tensA,site_tensB);
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


std::vector<double> env_elems_from_env_tens(const std::vector<int> &flavor_accumulate_deg, const std::vector<Spin_Basis> &spin_basis, const std::array<std::vector<IQTensor>,2> &env_tens)
{
    std::vector<double> env_elems(flavor_accumulate_deg.back());
    const auto &temp_tensor=env_tens[0][0];
    auto indices=temp_tensor.indices();
    for (int vali=1; vali<=indices[0].m(); vali++)
    {
        env_elems[flavor_accumulate_deg[spin_basis[vali-1].S()]+spin_basis[vali-1].t()]=temp_tensor(indices[0](vali),indices[1](vali));
    }

    //cout << "check for function env_elems_from_env_tens:" << endl;
    //Print(env_elems);
    //PrintDat(temp_tensor);

    return env_elems;
}

void obtain_env_tens_from_env_elems(const std::vector<int> &flavor_accumulate_deg, const std::vector<Spin_Basis> &spin_basis, const std::vector<double> env_elems, std::array<std::vector<IQTensor>,2> &env_tens)
{
    //assign env_elems to temp_tensor
    IQTensor temp_tensor(env_tens[0][0]);
    temp_tensor-=temp_tensor; //set temp_tensor to zero
    auto temp_indices=temp_tensor.indices();
    for (int vali=1; vali<=temp_indices[0].m(); vali++)
    {
        double elem=env_elems[flavor_accumulate_deg[spin_basis[vali-1].S()]+spin_basis[vali-1].t()];
        temp_tensor(temp_indices[0](vali),temp_indices[1](vali))=elem;
    }

    //get all env_tens by replaceIndex
    //cout << "check for function obtain_env_tens_from_env_elems:" << endl;
    //Print(env_elems);
    for (auto &env_tens_one_side : env_tens)
    {
        for (auto &env_tens_one_leg : env_tens_one_side) 
        {
            if (temp_tensor.indices()[0].dir()!=env_tens_one_leg.indices()[0].dir()) temp_tensor.dag();
            //Print(env_tens_one_leg.indices());
            tensor_assignment(env_tens_one_leg,temp_tensor); 
            //PrintDat(env_tens_one_leg);
        }
    }
    //cout << "finish checking obtain_env_tens_from_env_elems!" << endl;
}


double updated_env_tens_diff_f(const gsl_vector *x, void *params)
{
    Env_Tens_Params *env_tens_params=(Env_Tens_Params *)params;
    std::vector<double> env_elems(env_tens_params->flavor_accumulate_deg.back());

    for (int i=0; i<env_elems.size(); i++)
        env_elems[i]=gsl_vector_get(x,i);

    auto env_tens=env_tens_params->env_tens;
    obtain_env_tens_from_env_elems(env_tens_params->flavor_accumulate_deg,env_tens_params->spin_basis,env_elems,env_tens);

    //get tensors combined site tensors and env tensors
    auto combined_tens=env_tens_params->site_tens;
    for (int sitei=0; sitei<2; sitei++)
    {
        for (const auto &env_tens_one_leg : env_tens[sitei])
            combined_tens[sitei]*=env_tens_one_leg;
    }

    const auto &env_tens_original=env_tens[0][0];
    auto env_tens_updated=env_tens_original;
    auto updated_env_diag=nondeg_spin_sym_env_updated(combined_tens[0],combined_tens[1]);
    auto &env_tens_indices=env_tens_updated.indices();
    for (int vali=1; vali<=env_tens_indices[0].m(); vali++)
    {
        env_tens_updated(env_tens_indices[0](vali),env_tens_indices[1](vali))=updated_env_diag[vali-1];
    }

    double diff_norm=(env_tens_updated-env_tens_original).norm();

    //cout << "check for updated_env_tens_diff_f:" << endl;
    //Print(env_elems);
    //Print(updated_env_diag);
    //Print(diff_norm);
    
    return diff_norm;
}

void updated_env_tens_diff_df(const gsl_vector *x, void *params, gsl_vector *df)
{
    int x_size=x->size;

    gsl_vector *x_plus_dx;
    x_plus_dx=gsl_vector_alloc(x_size);
    gsl_vector_memcpy(x_plus_dx,x);

    //x_plus_dx=gsl_vector_alloc(x_size);
    //for (int i=0; i<x_size; i++)
    //{
    //    gsl_vector_set(x_plus_dx,i,gsl_vector_get(x,i));
    //}

    //TODO:improve the numerical derivative?
    double f_x=updated_env_tens_diff_f(x,params);
    //cout << "check updated_env_tens_diff_df:" << endl;
    //Print(f_x);
    for (int i=0; i<x_size; i++)
    {
        double dxi=1E-10;
        gsl_vector_set(x_plus_dx,i,gsl_vector_get(x,i)+dxi);
        double f_x_plus_dxi=updated_env_tens_diff_f(x_plus_dx,params);
        gsl_vector_set(df,i,(f_x_plus_dxi-f_x)/dxi);
        gsl_vector_set(x_plus_dx,i,gsl_vector_get(x,i));

        //Print(i);
        //Print(f_x_plus_dxi);
        //Print(gsl_vector_get(df,i));
    }

    //for (int i=0; i<x_size; i++) Print(gsl_vector_get(df,i))

    gsl_vector_free(x_plus_dx);
}

void updated_env_tens_diff_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
    *f=updated_env_tens_diff_f(x,params);
    updated_env_tens_diff_df(x,params,df);
}


void get_env_tensor_iterative(const IQTensor &site_tensA, const IQTensor &site_tensB, std::array<std::vector<IQTensor>,2> &env_tens)
{
    if (env_tens[0].empty())
    {
        init_env_tensor(site_tensA,site_tensB,env_tens);
    }

    //for (const auto &tensor : env_tens[0]) PrintDat(tensor);
    //for (const auto &tensor : env_tens[1]) PrintDat(tensor);

    int n_out_legs=env_tens[0].size();
    //assert(n_out_legs==env_tens[1].size());

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

        for (int sitei=0; sitei<2; sitei++)
        {
            for (int legi=0; legi<n_out_legs; legi++)
            {
                if (new_env_indset[0].dir()==env_tens[sitei][legi].indices()[0].dir()) 
                    tensor_assignment(env_tens[sitei][legi],new_env_tensor);
                else
                    tensor_assignment(env_tens[sitei][legi],dag(new_env_tensor));
            }
        }

        iter++;
    }

    if (iter==max_iter)
    {
        cout << "Unable to find good environment by iteration! Redo the iteration!" << endl;
        //exit(EXIT_FAILURE);
        for (auto &tensors : env_tens)
            tensors.clear();
        get_env_tensor_iterative(site_tensA,site_tensB,env_tens);
    }

}

