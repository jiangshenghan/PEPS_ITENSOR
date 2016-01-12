
#include "utilities.h"

IQIndex Spin_leg(const std::vector<int> &flavor_deg, const std::string &iqind_name, Arrow dir, IndexType it, int qn_order)
{
    //spin_dim=2S+1, where S is the largest spin #
    int spin_dim=flavor_deg.size();
    //while (flavor_deg[spin_dim-1]==0) spin_dim--;

    //sz_deg[qn] stores the degeneracy of Sz=(spin_dim-qn-1)/2
    std::vector<int> sz_deg;
    spin_to_sz(flavor_deg,sz_deg);
    
    //for (const auto &deg : flavor_deg) cout << deg << " ";
    //cout << endl;

    //for (const auto &deg : sz_deg) cout << deg << " ";
    //cout << endl;

    //init legs support rep of spin rotation with flavor_deg
    std::vector<IndexQN> indqn;
    for (int qn=0; qn<2*spin_dim-1; qn++)
    {
        if (sz_deg[qn]==0) continue;

        //stores indqn from high sz to low sz
        if (qn_order==-1)
        {
            indqn.push_back(IndexQN(Index(nameint("2Sz=",spin_dim-qn-1),sz_deg[qn],it),QN(spin_dim-qn-1,0)));
        }
        //stores indqn from low sz to high sz
        if (qn_order==1)
        {
            indqn.push_back(IndexQN(Index(nameint("2Sz=",-spin_dim+qn+1),sz_deg[qn],it),QN(-spin_dim+qn+1,0)));
        }
        
        //cout << indqn[qn] << endl;
    }
    
    IQIndex iqind=IQIndex(iqind_name,indqn,dir);

    return iqind;
}


bool iqind_spin_rep(const IQIndex &sz_leg, std::vector<int> &flavor_deg)
{
    //spin_dim=2S_{max}+1
    int spin_dim=0;
    for (const auto &indqn : sz_leg.indices())
    {
        if (indqn.qn.sz()>spin_dim) 
            spin_dim=indqn.qn.sz();
    }
    spin_dim+=1;

    //sz_deg[qn] stores the degeneracy of 2Sz=spin_dim-qn-1
    std::vector<int> sz_deg(2*spin_dim-1,0);
    for (const auto &indqn : sz_leg.indices())
    {
        sz_deg[spin_dim-indqn.qn.sz()-1]=indqn.m();
    }

    return sz_to_spin(sz_deg,flavor_deg);
}


bool sz_to_spin(const std::vector<int> &sz_deg, std::vector<int> &flavor_deg)
{
    int spin_dim=(sz_deg.size()+1)/2;
    flavor_deg=std::vector<int>(spin_dim,0);

    //Notice, spin_S and spin_Sz are in the unit of 1/2
    //deg of -Sz must equal to deg of Sz
    for (int spin_Sz=0; spin_Sz<=spin_dim-1; spin_Sz++)
    {
        if (sz_deg[spin_dim-spin_Sz-1]!=sz_deg[spin_dim+spin_Sz-1])
            return false;
    }

    std::vector<int> sz_deg_prime(sz_deg);
    for (int spin_S=spin_dim-1; spin_S>=0; spin_S--)
    {
        flavor_deg[spin_S]=sz_deg_prime[spin_dim-spin_S-1];

        for (int spin_Sz=-spin_S; spin_Sz<=spin_S; spin_Sz+=2)
        {
            sz_deg_prime[spin_dim-spin_Sz-1]-=flavor_deg[spin_S];
            if (sz_deg_prime[spin_dim-spin_Sz-1]<0) return false;
        }
        
        //cout << spin_S/2.0 << ": " << flavor_deg[spin_S] << endl;
    }

    //all elements of sz_deg_prime must vanish to form a spin rep
    auto sz_iter=std::find_if_not(sz_deg_prime.begin(),sz_deg_prime.end(),[](int i){return i==0;});

    if (sz_iter!=sz_deg_prime.end())
    {
        flavor_deg=std::vector<int>();
        return false;
    }
    return true;
}


void spin_to_sz(const std::vector<int> &flavor_deg, std::vector<int> &sz_deg)
{
    int spin_dim=flavor_deg.size();
    sz_deg=std::vector<int>(2*spin_dim-1,0);

    for (int spin_qn=0; spin_qn<spin_dim; spin_qn++)
    {
        for (int qn=spin_dim-spin_qn-1; qn<=spin_dim+spin_qn-1; qn+=2)
        {
            sz_deg[qn]+=flavor_deg[spin_qn];
        }
    }
}


bool iqind_to_spin_basis(const IQIndex &sz_leg, std::vector<Spin_Basis> &spin_basis)
{
    std::vector<int> flavor_deg;

    if (!iqind_spin_rep(sz_leg,flavor_deg)) return false;

    return iqind_to_spin_basis(sz_leg,flavor_deg,spin_basis);
}

bool iqind_to_spin_basis(const IQIndex &sz_leg, const std::vector<int> &flavor_deg, std::vector<Spin_Basis> &spin_basis)
{
    for (const auto &indqn : sz_leg.indices())
    {
        int m=indqn.qn.sz();
        for (int S=std::abs(m); S<flavor_deg.size(); S+=2)
        {
            for (int t=0; t<flavor_deg[S]; t++)
            {
                spin_basis.push_back(Spin_Basis(S,m,t));
            }
        }
    }
    
    assert(spin_basis.size()==sz_leg.m());

    return true;
}


std::vector<int> list_from_num(int num, const std::vector<int> &max_nums)
{
    std::vector<int> list(max_nums.size(),0);

    for (int i=list.size()-1; i>=0; i--)
    {
        list[i]=num%max_nums[i];
        num/=max_nums[i];

        if (num==0) break;

    }

    if (num>0)
    {
        cout << "Invalid num!" << endl;
        exit (EXIT_FAILURE);
    }

    return list;
}

int num_from_list(const std::vector<int> &list, const std::vector<int> &max_nums)
{
    for (int i=0; i<list.size(); i++)
    {
        if (list[i]>=max_nums[i])
        {
            cout << "Invalid list!" << endl;
            exit (EXIT_FAILURE);
        }
    }

    int num=list[0];

    for (int i=1; i<list.size(); i++)
    {
        num*=max_nums[i];
        num+=list[i];
    }

    //Print(list);
    //Print(max_nums);
    //Print(num);

    return num;
}


IQTensor eta_from_mu(double mu, const std::vector<int> &flavor_deg)
{
    auto eta_leg=Spin_leg(flavor_deg,"eta_leg",In),
         eta_leg_prime=prime(dag(eta_leg));
    auto eta=IQTensor(eta_leg,eta_leg_prime);

    if (std::abs(mu-1)<EPSILON)
    {
        for (int val=1; val<=eta_leg.m(); val++)
        {
            eta(eta_leg(val),eta_leg_prime(val))=1;
        }
        return eta;
    }
    else
    {
        std::vector<Spin_Basis> eta_leg_spin_basis;
        iqind_to_spin_basis(eta_leg,flavor_deg,eta_leg_spin_basis);

        for (int val=1; val<=eta_leg.m(); val++)
        {
            eta(eta_leg(val),eta_leg_prime(val))=std::pow(mu,eta_leg_spin_basis[val-1].S());
        }
        return eta;
    }

}


IQTensor eta_from_mu(double mu, IQIndex eta_leg)
{
    //if (eta_leg.dir()==In) eta_leg.dag();
    auto eta_leg_prime=prime(dag(eta_leg));
    auto eta=IQTensor(eta_leg,eta_leg_prime);
    
    if (std::abs(mu-1)<EPSILON)
    {
        for (int val=1; val<=eta_leg.m(); val++)
        {
            eta(eta_leg(val),eta_leg_prime(val))=1;
        }
        return eta;
    }
    else
    {
        std::vector<Spin_Basis> eta_leg_spin_basis;
        iqind_to_spin_basis(eta_leg,eta_leg_spin_basis);

        for (int val=1; val<=eta_leg.m(); val++)
        {
            eta(eta_leg(val),eta_leg_prime(val))=std::pow(mu,eta_leg_spin_basis[val-1].S());
        }
        return eta;
    }

}

ITensor eta_from_mu(double mu, Index eta_leg)
{
    auto eta_leg_prime=prime(eta_leg);
    auto eta=ITensor(eta_leg,eta_leg_prime);
    for (int val=1; val<=eta_leg.m(); val++)
    {
        eta(eta_leg(val),eta_leg_prime(val))=1;
    }
    return eta;
}


Index isomorphic_legs(const Index &old_leg, const std::string &new_leg_name)
{
    return Index(new_leg_name,old_leg.m(),old_leg.type(),old_leg.primeLevel());
}

IQIndex isomorphic_legs(const IQIndex &old_leg, const std::string &new_leg_name)
{
    auto indices_qn=old_leg.indices();
    return IQIndex(new_leg_name,indices_qn,old_leg.dir(),old_leg.primeLevel());
}


void tensor_assignment(ITensor &TA, const ITensor &TB)
{
    ITensor tensor_tmp=TB;

    for (int leg_i=0; leg_i<TA.r(); leg_i++)
    {
        auto oind=tensor_tmp.indices()[leg_i],
             nind=TA.indices()[leg_i];
        tensor_tmp.replaceIndex(oind,nind);
    }
    TA=tensor_tmp;

    return;
}

void tensor_assignment(IQTensor &TA, const IQTensor &TB)
{
    //cout << Global::args() << endl;
    IQTensor tensor_tmp=TB;

    //cout << TB.indices() << endl;
    //cout << TA.indices();
    for (int leg_i=0; leg_i<TA.r(); leg_i++)
    {
        auto oind=tensor_tmp.indices()[leg_i],
             nind=TA.indices()[leg_i];
        tensor_tmp.replaceIndex(oind,nind);
        //cout << tensor_tmp.indices() << endl;
    }
    TA=tensor_tmp;


    //cout << "Check for tensor_assignment" << endl;
    //cout << "old tensor:" << endl; 
    //PrintDat(TB);
    //cout << "new tensor:" << endl;
    //PrintDat(TA);
    
    return;
}


template <class TensorT>
TensorT tensor_permutation(const std::vector<int> &permuted_indices, const TensorT &tensor_origin)
{
    if (tensor_origin.isComplex())
    {
        TensorT tensor_origin_real=realPart(tensor_origin),
                tensor_origin_imag=imagPart(tensor_origin);
        TensorT tensor_permutation_real=tensor_permutation(permuted_indices,tensor_origin_real),
                tensor_permutation_imag=tensor_permutation(permuted_indices,tensor_origin_imag);

        return (tensor_permutation_real+Complex_i*tensor_permutation_imag);
    }

    using IndexValT=typename TensorT::IndexValT;
    auto tensor_permutation(tensor_origin);
    tensor_permutation-=tensor_origin;

    std::vector<int> max_val_list;
    auto tensor_legs=tensor_origin.indices();
    for (const auto &leg : tensor_legs) max_val_list.push_back(leg.m());
    for (int val_num=0; val_num<tensor_legs.dim(); val_num++)
    {
        auto val_list=list_from_num(val_num,max_val_list);
        std::vector<IndexValT> leg_vals(NMAX,IndexValT::Null());
        for (int legi=0; legi<tensor_legs.r(); legi++) leg_vals[legi]=tensor_legs[legi](val_list[legi]+1);
        
        auto elem=tensor_origin(leg_vals[0],leg_vals[1],leg_vals[2],leg_vals[3],leg_vals[4],leg_vals[5],leg_vals[6],leg_vals[7]);
        if (std::abs(elem)<EPSILON) continue;

        std::vector<int> val_list_permuted(val_list.size());
        for (int posi=0; posi<val_list.size(); posi++)
            val_list_permuted[permuted_indices[posi]]=val_list[posi];
        std::vector<IndexValT> leg_vals_permuted(NMAX,IndexValT::Null());
        for (int legi=0; legi<tensor_legs.r(); legi++) leg_vals_permuted[legi]=tensor_legs[legi](val_list_permuted[legi]+1);

        //TODO:set complex number?
        tensor_permutation(leg_vals_permuted[0],leg_vals_permuted[1],leg_vals_permuted[2],leg_vals_permuted[3],leg_vals_permuted[4],leg_vals_permuted[5],leg_vals_permuted[6],leg_vals_permuted[7])=elem;
    }
    clean(tensor_permutation);

    return tensor_permutation;
}
template
ITensor tensor_permutation(const std::vector<int> &permuted_indices, const ITensor &tensor_origin);
template 
IQTensor tensor_permutation(const std::vector<int> &permuted_indices, const IQTensor &tensor_origin);


template <class TensorT>
void tensor_assignment_diff_order(TensorT &TA, const TensorT &TB)
{
    if (TB.isComplex())
    {
        TensorT TA_real=realPart(TA), TA_imag=imagPart(TA),
                TB_real=realPart(TB), TB_imag=imagPart(TB);
        tensor_assignment_diff_order(TA_real,TB_real);
        tensor_assignment_diff_order(TA_imag,TB_imag);
        TA=TA_real+Complex_i*TA_imag;
        return;
    }

    using IndexValT=typename TensorT::IndexValT;

    std::vector<int> max_val_list;
    auto tensor_legs=TA.indices();
    for (const auto &leg : tensor_legs) max_val_list.push_back(leg.m());

    for (int val_num=0; val_num<tensor_legs.dim(); val_num++)
    {
        auto val_list=list_from_num(val_num,max_val_list);
        std::vector<IndexValT> leg_vals(NMAX,IndexValT::Null());
        for (int legi=0; legi<tensor_legs.r(); legi++) leg_vals[legi]=tensor_legs[legi](val_list[legi]+1);
        
        auto elem=TB(leg_vals[0],leg_vals[1],leg_vals[2],leg_vals[3],leg_vals[4],leg_vals[5],leg_vals[6],leg_vals[7]);
        if (std::abs(elem)<EPSILON) continue;
        TA(leg_vals[0],leg_vals[1],leg_vals[2],leg_vals[3],leg_vals[4],leg_vals[5],leg_vals[6],leg_vals[7])=elem;
        
    }
    clean(TA);
}
template
void tensor_assignment_diff_order(ITensor &TA, const ITensor &TB);
template
void tensor_assignment_diff_order(IQTensor &TA, const IQTensor &TB);

