
#include "peps_indexset.h"

template<class IndexT>
PEPSt_IndexSet_Base<IndexT>::PEPSt_IndexSet_Base(const int &d, const int &D, const int &n_sites_total, const int &n_bonds_total): 
    d_(d), 
    D_(D),
    phys_legs_(n_sites_total),
    virt_legs_(n_bonds_total)
{ }
template 
PEPSt_IndexSet_Base<Index>::PEPSt_IndexSet_Base(const int &d, const int &D, const int& n_sites_total, const int& n_bonds_total);
template 
PEPSt_IndexSet_Base<IQIndex>::PEPSt_IndexSet_Base(const int &d, const int &D, const int& n_sites_total, const int& n_bonds_total);


//
//PEPS_IndexSet
//
PEPS_IndexSet::PEPS_IndexSet(const int &d, const int &D, const Lattice_Torus_Base &lattice): PEPSt_IndexSet_Base<Index>(d,D,lattice.n_sites_total(),lattice.n_bonds_total())
{
    init_phys_legs(lattice.n_sites_total());
    init_virt_legs(lattice.n_bonds_total());
}

void PEPS_IndexSet::init_phys_legs(const int &n_sites_total)
{
    for (int site_i=0; site_i<n_sites_total; site_i++)
    {
        phys_legs_[site_i]=Index(nameint("site=",site_i),d_,Site);
    }
}

void PEPS_IndexSet::init_virt_legs(const int &n_bonds_total)
{
    for (int bond_i=0; bond_i<n_bonds_total; bond_i++)
    {
        virt_legs_[bond_i]=Index(nameint("bond=",bond_i),D_,Link);
    }
}


//
//IQPEPS_IndexSet_SpinHalf
//
IQPEPS_IndexSet_SpinHalf::IQPEPS_IndexSet_SpinHalf(const int &D, const Lattice_Torus_Base &lattice): PEPSt_IndexSet_Base<IQIndex>(2,D,lattice.n_sites_total(),lattice.n_bonds_total())
{
    //construct physical legs
    init_phys_legs(lattice.n_sites_total());

    //Decompose a virtual leg to spins, and check D
    //spin_dim=2S+1, which is used to store the largest spin #
    int Dprime=D, spin_dim=1;
    for (; Dprime>0; spin_dim++)
    {
        Dprime-=spin_dim;
    }
    spin_dim-=1;
    assert(Dprime==0);

    std::vector<int> bond_indqn_deg(2*spin_dim-1);

    for (int qn_i=0; qn_i<spin_dim; qn_i++)
    {
        bond_indqn_deg[qn_i]=qn_i/2+1;
        bond_indqn_deg[2*spin_dim-2-qn_i]=qn_i/2+1;
        
        //cout << "bond_indqn_deg[" << qn_i << "]=" << bond_indqn_deg[qn_i] << endl
        //     << "bond_indqn_deg[" << 2*spin_dim-1-qn_i << "]=" << bond_indqn_deg[2*spin_dim-1-qn_i] << endl;
    }

    //cout << "bond_iqindex_qn_deg: "; 
    //for (const auto &deg : bond_indqn_deg)
    //    cout << deg << " ";
    //cout << endl;

    //construct virtual legs
    init_virt_legs(lattice.n_bonds_total(),spin_dim,bond_indqn_deg);
}

IQPEPS_IndexSet_SpinHalf::IQPEPS_IndexSet_SpinHalf(const int &D, const std::vector<int> &virt_leg_spin, const Lattice_Torus_Base &lattice): PEPSt_IndexSet_Base<IQIndex>(2,D,lattice.n_sites_total(),lattice.n_bonds_total())
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
    for (int site_i=0; site_i<lattice.n_sites_total(); site_i++)
    {
        phys_legs_[site_i]=IQIndex(nameint("S=1/2 site=",site_i),
                Index(nameint("Up for site",site_i),1,Site),QN(+1,0),
                Index(nameint("Down for site",site_i),1,Site),QN(-1,0));
    }

    //bond_indqn_deg[qn_i] stores the degeneracy of space Sz=(spin_dim-qn_i-1)/2
    std::vector<int> bond_indqn_deg(2*spin_dim-1,0);

    //construct bond_indqn_deg
    for (int spin_i=0; spin_i<spin_dim; spin_i++)
    {
        for (int qn_i=spin_dim-spin_i-1; qn_i<=spin_dim+spin_i-1; qn_i+=2)
        {
            bond_indqn_deg[qn_i]+=virt_leg_spin[spin_i];
        }
    }

    //Construct virtual legs
    init_virt_legs(lattice.n_bonds_total(),spin_dim,bond_indqn_deg);
}



void IQPEPS_IndexSet_SpinHalf::init_phys_legs(const int &n_sites_total)
{
    for (int site_i=0; site_i<n_sites_total; site_i++)
    {
        phys_legs_[site_i]=IQIndex(nameint("S=1/2 site=",site_i),
                Index(nameint("Up for site",site_i),1,Site),QN(+1,0),
                Index(nameint("Down for site",site_i),1,Site),QN(-1,0));
    }
}

void IQPEPS_IndexSet_SpinHalf::init_virt_legs(const int &n_bonds_total, const int &spin_dim, const std::vector<int> &bond_indqn_deg)
{
    for (int bond_i=0; bond_i<n_bonds_total; bond_i++)
    {
        //bond_indqn[qn_i] stores index with Sz=(spin_dim-qn_i-1)/2
        //with dimension bond_indqn_deg[qn_i]
        std::vector<IndexQN> bond_indqn;
        for (int qn_i=0; qn_i<2*spin_dim-1; qn_i++)
        {
            std::stringstream ss;
            ss << "Sz=" << (spin_dim-qn_i-1)/2.0 << " for bond " << bond_i;
            std::string ind_name=ss.str();
            bond_indqn.push_back(IndexQN(Index(ind_name,bond_indqn_deg[qn_i],Link),QN(spin_dim-qn_i-1,0)));
        }
        virt_legs_[bond_i]=IQIndex(nameint("Bond",bond_i),bond_indqn);
    }
}
