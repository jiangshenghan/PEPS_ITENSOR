
#include "utilities.h"

IQIndex Spin_leg(const std::vector<int> &spin_deg, const std::string &iqind_name, Arrow dir, IndexType it, int qn_order)
{
    //spin_dim=2S+1, where S is the largest spin #
    int spin_dim=spin_deg.size();
    //while (spin_deg[spin_dim-1]==0) spin_dim--;

    //sz_deg[qn] stores the degeneracy of Sz=(spin_dim-qn-1)/2
    std::vector<int> sz_deg;
    spin_to_sz(spin_deg,sz_deg);
    
    //for (const auto &deg : spin_deg) cout << deg << " ";
    //cout << endl;

    //for (const auto &deg : sz_deg) cout << deg << " ";
    //cout << endl;

    //init legs support rep of spin rotation with spin_deg
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


bool iqind_spin_rep(const IQIndex &sz_leg, std::vector<int> &spin_deg)
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

    return sz_to_spin(sz_deg,spin_deg);
}


bool sz_to_spin(const std::vector<int> &sz_deg, std::vector<int> &spin_deg)
{
    int spin_dim=(sz_deg.size()+1)/2;
    spin_deg=std::vector<int>(spin_dim,0);

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
        spin_deg[spin_S]=sz_deg_prime[spin_dim-spin_S-1];

        for (int spin_Sz=-spin_S; spin_Sz<=spin_S; spin_Sz+=2)
        {
            sz_deg_prime[spin_dim-spin_Sz-1]-=spin_deg[spin_S];
            if (sz_deg_prime[spin_dim-spin_Sz-1]<0) return false;
        }
        
        //cout << spin_S/2.0 << ": " << spin_deg[spin_S] << endl;
    }

    //all elements of sz_deg_prime must vanish to form a spin rep
    auto sz_iter=std::find_if_not(sz_deg_prime.begin(),sz_deg_prime.end(),[](int i){return i==0;});

    if (sz_iter!=sz_deg_prime.end())
    {
        spin_deg=std::vector<int>();
        return false;
    }
    return true;
}


void spin_to_sz(const std::vector<int> &spin_deg, std::vector<int> &sz_deg)
{
    int spin_dim=spin_deg.size();
    sz_deg=std::vector<int>(2*spin_dim-1,0);

    for (int spin_qn=0; spin_qn<spin_dim; spin_qn++)
    {
        for (int qn=spin_dim-spin_qn-1; qn<=spin_dim+spin_qn-1; qn+=2)
        {
            sz_deg[qn]+=spin_deg[spin_qn];
        }
    }
}


bool iqind_to_spin_basis(const IQIndex &sz_leg, std::vector<Spin_Basis> &spin_basis)
{
    std::vector<int> spin_deg;

    if (!iqind_spin_rep(sz_leg,spin_deg)) return false;

    return iqind_to_spin_basis(sz_leg,spin_deg,spin_basis);
}

bool iqind_to_spin_basis(const IQIndex &sz_leg, const std::vector<int> &spin_deg, std::vector<Spin_Basis> &spin_basis)
{
    for (const auto &indqn : sz_leg.indices())
    {
        int m=indqn.qn.sz();
        for (int S=std::abs(m); S<spin_deg.size(); S+=2)
        {
            for (int t=0; t<spin_deg[S]; t++)
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
        num*=max_nums[i-1];
        num+=list[i];
    }

    return num;
}


IQTensor eta_from_mu(double mu, const std::vector<int> &spin_deg)
{
    auto eta_leg=Spin_leg(spin_deg,"eta_leg",In),
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

    std::vector<Spin_Basis> eta_leg_spin_basis;
    iqind_to_spin_basis(eta_leg,spin_deg,eta_leg_spin_basis);

    for (int val=1; val<=eta_leg.m(); val++)
    {
        eta(eta_leg(val),eta_leg_prime(val))=std::pow(mu,eta_leg_spin_basis[val-1].S());
    }
    return eta;

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

    std::vector<Spin_Basis> eta_leg_spin_basis;
    iqind_to_spin_basis(eta_leg,eta_leg_spin_basis);

    for (int val=1; val<=eta_leg.m(); val++)
    {
        eta(eta_leg(val),eta_leg_prime(val))=std::pow(mu,eta_leg_spin_basis[val-1].S());
    }
    return eta;

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


