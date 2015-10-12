
#include "peps_indexset.h"

template<class IndexT>
PEPSt_IndexSet_Base<IndexT>::PEPSt_IndexSet_Base(const int &d, const int &D, const int &n_sites_total, const int &n_bonds_to_one_site, const int &n_boundary_legs): 
    d_(d), 
    D_(D),
    phys_legs_(n_sites_total),
    //notice for we only add 1/2 boundary legs to virtual legs, which belongs to boundary bond virtual legs, since the other 1/2 boundary legs has been counted in boundary site virtual legs
    virt_legs_(n_bonds_to_one_site*n_sites_total+n_boundary_legs/2)
{ 
    //cout << n_boundary_legs << endl;
    //cout << virt_legs_.size() << endl;
}
template 
PEPSt_IndexSet_Base<Index>::PEPSt_IndexSet_Base(const int &d, const int &D, const int &n_sites_total, const int &n_bonds_to_one_site, const int &n_boundary_legs);
template 
PEPSt_IndexSet_Base<IQIndex>::PEPSt_IndexSet_Base(const int &d, const int &D, const int &n_sites_total, const int &n_bonds_to_one_site, const int &n_boundary_legs);


//
//PEPS_IndexSet
//
PEPS_IndexSet::PEPS_IndexSet(const int &d, const int &D, const Lattice_Base &lattice): PEPSt_IndexSet_Base<Index>(d,D,lattice.n_sites_total(),lattice.n_bonds_to_one_site(),lattice.n_boundary_legs())
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
IQPEPS_IndexSet_SpinHalf::IQPEPS_IndexSet_SpinHalf(const int &D, const Lattice_Base &lattice): 
    PEPSt_IndexSet_Base<IQIndex>(2,D,lattice.n_sites_total(),lattice.n_bonds_to_one_site(),lattice.n_boundary_legs())
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

    name_+="virt_legs: nondeg, ";
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

IQPEPS_IndexSet_SpinHalf::IQPEPS_IndexSet_SpinHalf(const int &D, const std::vector<int> &virt_leg_spin, const Lattice_Base &lattice): 
    PEPSt_IndexSet_Base<IQIndex>(2,D,lattice.n_sites_total(),lattice.n_bonds_to_one_site(),lattice.n_boundary_legs())
{
    //Check the input of D
    int Dprime=0;
    //spin_dim=2S+1, where S is the largest spin #
    int spin_dim=virt_leg_spin.size();
    bool extra_deg=false;
    for (int spin_i=0; spin_i<spin_dim; spin_i++)
    {
        Dprime+=(spin_i+1)*virt_leg_spin[spin_i];
        if (virt_leg_spin[spin_i]>1) extra_deg=true;
    }
    //cout << Dprime << endl;
    assert(Dprime==D);

    if (!extra_deg)
        name_+="virt_legs: nondeg, ";


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
    name_="phys_leg: half spin, ";

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


IQPEPS_IndexSet_Spin_Sym::IQPEPS_IndexSet_Spin_Sym(int d, int D, const std::vector<int> &phys_leg_spin_deg, const std::vector<int> &virt_leg_spin_deg, const Lattice_Base &lattice, int phys_legs_qn_order, int virt_legs_qn_order):
    PEPSt_IndexSet_Base<IQIndex>(d,D,lattice.n_sites_total(),lattice.n_bonds_to_one_site(),lattice.n_boundary_legs())
{
    name_="spin symmetric, ";

    init_phys_legs(phys_leg_spin_deg,phys_legs_qn_order);
    init_virt_legs(virt_leg_spin_deg,virt_legs_qn_order);
}


void IQPEPS_IndexSet_Spin_Sym::init_phys_legs(const std::vector<int> &phys_leg_spin_deg, int phys_legs_qn_order)
{
    int d_prime=0;
    int spini=0;
    for (int deg : phys_leg_spin_deg)
    {
        d_prime+=deg*(spini+1);
        spini++;
    }
    assert(d_prime==d_);

    int sitei=0;
    for (auto &leg : phys_legs_)
    {
        leg=Spin_leg(phys_leg_spin_deg,nameint("phys_leg ",sitei),Out,Site,phys_legs_qn_order);
        sitei++;
    }
}

void IQPEPS_IndexSet_Spin_Sym::init_virt_legs(const std::vector<int> &virt_leg_spin_deg, int virt_legs_qn_order)
{
    int D_prime=0;
    int spini=0;
    for (int deg : virt_leg_spin_deg)
    {
        D_prime+=deg*(spini+1);
        spini++;
    }
    assert(D_prime==D);

    int linki=0;
    for (auto &leg : virt_legs_)
    {
        leg=Spin_leg(virt_leg_spin_deg,nameint("virt_leg ",linki),Out,Link,virt_legs_qn_order);
        linki++;
    }
}

