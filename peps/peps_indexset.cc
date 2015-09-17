
#include "peps_indexset.h"

template<class IndexT>
PEPSt_IndexSet_Base<IndexT>::PEPSt_IndexSet_Base(const int &d, const int &D, const int &n_sites_total, const int &n_bonds_to_one_site): 
    d_(d), 
    D_(D),
    phys_legs_(n_sites_total),
    virt_legs_(n_bonds_to_one_site*n_sites_total)
{ }
template 
PEPSt_IndexSet_Base<Index>::PEPSt_IndexSet_Base(const int &d, const int &D, const int &n_sites_total, const int &n_bonds_to_one_site);
template 
PEPSt_IndexSet_Base<IQIndex>::PEPSt_IndexSet_Base(const int &d, const int &D, const int &n_sites_total, const int &n_bonds_to_one_site);


//
//PEPS_IndexSet
//
PEPS_IndexSet::PEPS_IndexSet(const int &d, const int &D, const Lattice_Torus_Base &lattice): PEPSt_IndexSet_Base<Index>(d,D,lattice.n_sites_total(),lattice.n_bonds_to_one_site())
{
    init_phys_legs();
    init_virt_legs();
}

void PEPS_IndexSet::init_phys_legs()
{
    int site_i=0;
    for (auto &leg : phys_legs_)
    {
        leg=Index(nameint("phys_leg ",site_i),d_,Site);
        site_i++;
        
        //cout << leg << endl;
    }
}

void PEPS_IndexSet::init_virt_legs()
{
    int leg_i=0;
    for (auto &leg : virt_legs_)
    {
        leg=Index(nameint("virt_leg ",leg_i),D_,Link);
        leg_i++;

        //cout << leg << endl;
    }
}


//
//IQPEPS_IndexSet_SpinHalf
//
IQPEPS_IndexSet_SpinHalf::IQPEPS_IndexSet_SpinHalf(const int &D, const Lattice_Torus_Base &lattice): PEPSt_IndexSet_Base<IQIndex>(2,D,lattice.n_sites_total(),lattice.n_bonds_to_one_site())
{
    //construct physical legs
    init_phys_legs();

    //Decompose a virtual leg to spins, and check D
    //spin_dim=2S+1, which is used to store the largest spin #
    int Dprime=D, spin_dim=1;
    for (; Dprime>0; spin_dim++)
    {
        Dprime-=spin_dim;
    }
    spin_dim-=1;
    assert(Dprime==0);

    std::vector<int> virt_indqn_deg(2*spin_dim-1);

    for (int qn_i=0; qn_i<spin_dim; qn_i++)
    {
        virt_indqn_deg[qn_i]=qn_i/2+1;
        virt_indqn_deg[2*spin_dim-2-qn_i]=qn_i/2+1;
        
        //cout << "virt_indqn_deg[" << qn_i << "]=" << virt_indqn_deg[qn_i] << endl
        //     << "virt_indqn_deg[" << 2*spin_dim-1-qn_i << "]=" << virt_indqn_deg[2*spin_dim-1-qn_i] << endl;
    }

    //cout << "bond_iqindex_qn_deg: "; 
    //for (const auto &deg : virt_indqn_deg)
    //    cout << deg << " ";
    //cout << endl;

    //construct virtual legs
    init_virt_legs(spin_dim,virt_indqn_deg);
}

IQPEPS_IndexSet_SpinHalf::IQPEPS_IndexSet_SpinHalf(const int &D, const std::vector<int> &virt_leg_spin, const Lattice_Torus_Base &lattice): PEPSt_IndexSet_Base<IQIndex>(2,D,lattice.n_sites_total(),lattice.n_bonds_to_one_site())
{
    //Check the input of D
    int Dprime=0;
    //spin_dim=2S+1, where S is the largest spin #
    int spin_dim=virt_leg_spin.size();
    for (int spin_i=0; spin_i<spin_dim; spin_i++)
    {
        Dprime+=(spin_i+1)*virt_leg_spin[spin_i];
    }
    //cout << Dprime << endl;
    assert(Dprime==D);


    //Constructor physical legs
    init_phys_legs();

    //virt_indqn_deg[qn_i] stores the degeneracy of space Sz=(spin_dim-qn_i-1)/2
    std::vector<int> virt_indqn_deg(2*spin_dim-1,0);

    //construct virt_indqn_deg, which stores information about degeneracy associated with S_z quantum number
    for (int spin_i=0; spin_i<spin_dim; spin_i++)
    {
        for (int qn_i=spin_dim-spin_i-1; qn_i<=spin_dim+spin_i-1; qn_i+=2)
        {
            virt_indqn_deg[qn_i]+=virt_leg_spin[spin_i];
        }
    }

    //Construct virtual legs
    init_virt_legs(spin_dim,virt_indqn_deg);
}



void IQPEPS_IndexSet_SpinHalf::init_phys_legs()
{
    int site_i=0;
    for (auto &leg : phys_legs_)
    {
        leg=IQIndex(nameint("S=1/2 phys_leg ",site_i),
                Index(nameint("Up for phys_leg ",site_i),1,Site),QN(+1,0),
                Index(nameint("Down for phys_leg ",site_i),1,Site),QN(-1,0));
        site_i++;
    }
}

void IQPEPS_IndexSet_SpinHalf::init_virt_legs(const int &spin_dim, const std::vector<int> &virt_indqn_deg)
{
    int leg_i=0;
    for (auto &leg : virt_legs_)
    {
        //virt_indqn[qn_i] stores index with Sz=(spin_dim-qn_i-1)/2
        //with dimension virt_indqn_deg[qn_i]
        std::vector<IndexQN> virt_indqn;
        for (int qn_i=0; qn_i<2*spin_dim-1; qn_i++)
        {
            if (virt_indqn_deg[qn_i]==0) continue;

            std::stringstream ss;
            ss << "Sz=" << (spin_dim-qn_i-1)/2.0 << " for virt_leg " << leg_i;
            std::string ind_name=ss.str();
            virt_indqn.push_back(IndexQN(Index(ind_name,virt_indqn_deg[qn_i],Link),QN(spin_dim-qn_i-1,0)));
        }
        leg=IQIndex(nameint("virt_leg ",leg_i),virt_indqn);

        leg_i++;
    }
}
