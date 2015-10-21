
#include "featureless_honeycomb.h"

int main()
{
    double A1=0.9, A2=1-A1;

    std::vector<IQIndex> phys_legs(6), virt_legs(18);

    for (int sitei=0; sitei<phys_legs.size(); sitei++)
    {
        phys_legs[sitei]=Spin_leg(std::vector<int>{0,1},nameint("phys_leg ",sitei),Out,Site);
    }
    for (int linki=0; linki<virt_legs.size(); linki++)
    {
        virt_legs[linki]=Spin_leg(std::vector<int>{0,2},nameint("virt_leg ",linki),Out,Link);
    }

    std::vector<IQTensor> site_tensors(6), bond_tensors(5), boundary_tensors(8);

    //parameter for site tensors
    std::array<double,2> norms={2.,std::sqrt(12.)};
    std::vector<double> T_tilde={0,-0.5*A2*norms[0],-0.5*A2*norms[0],0,A2*norms[0],1.5*A1*norms[0],-1.5*A1*norms[0],0,0,0.5*A2*norms[1],-0.5*A2*norms[1],A1*norms[1],0,-0.5*A1*norms[1],-0.5*A1*norms[1],0};
    //create site tensors
    for (int sitei=0; sitei<site_tensors.size(); sitei++)
    {
        std::vector<IQIndex> sitei_indices={phys_legs[sitei],virt_legs[3*sitei],virt_legs[3*sitei+1],virt_legs[3*sitei+2]};
        Singlet_Tensor_Basis sitei_basis(sitei_indices);
        site_tensors[sitei]=singlet_tensor_from_basis_params(sitei_basis,T_tilde);
    }
    //create bond tensors
    std::vector<double> B_tilde={std::sqrt(2),0,0,std::sqrt(2)};
    std::vector<std::vector<IQIndex>> bond_tensors_indices={{dag(virt_legs[0]),dag(virt_legs[3])},{dag(virt_legs[4]),dag(virt_legs[7])},{dag(virt_legs[5]),dag(virt_legs[13])},{dag(virt_legs[10]),dag(virt_legs[14])},{dag(virt_legs[12]),dag(virt_legs[15])}};
    for (int bondi=0; bondi<bond_tensors.size(); bondi++)
    {
        Singlet_Tensor_Basis bondi_basis(bond_tensors_indices[bondi]);
        bond_tensors[bondi]=singlet_tensor_from_basis_params(bondi_basis,B_tilde);
    }
    //create boundary tensors
    boundary_tensors[0]=IQTensor(dag(virt_legs[1]));
    boundary_tensors[1]=IQTensor(dag(virt_legs[2]));
    boundary_tensors[2]=IQTensor(dag(virt_legs[6]));
    boundary_tensors[3]=IQTensor(dag(virt_legs[8]));
    boundary_tensors[4]=IQTensor(dag(virt_legs[9]));
    boundary_tensors[5]=IQTensor(dag(virt_legs[11]));
    boundary_tensors[6]=IQTensor(dag(virt_legs[16]));
    boundary_tensors[7]=IQTensor(dag(virt_legs[17]));

    std::vector<std::vector<double>> boundary_vector={{1,0,1,0},{0,1,0,0},{0,1,0,0},{1,0,1,0},{0,1,0,0},{1,0,1,0},{0,1,0,0},{1,0,1,0}};
    int boundaryi=0;
    for (auto &tensor : boundary_tensors)
    {
        IQIndex ind=tensor.indices()[0];
        for (int val=1; val<=ind.m(); val++)
        {
            tensor(ind(val))=boundary_vector[boundaryi][val-1];
        }
        boundaryi++;
    }

    cout << "Site tensors:" << endl;
    for (const auto &tensor : site_tensors) PrintDat(tensor);
    cout << "Bond tensors:" << endl;
    for (const auto &tensor : bond_tensors) PrintDat(tensor);
    cout << "Boundary tensors:" << endl;
    for (const auto &tensor : boundary_tensors) PrintDat(tensor);


    //obtain contracted tensors
    std::vector<IQTensor> combined_site_tensors={site_tensors[0]*boundary_tensors[0]*boundary_tensors[1]*bond_tensors[0],site_tensors[1],site_tensors[2]*bond_tensors[1]*boundary_tensors[2]*boundary_tensors[3],site_tensors[3]*boundary_tensors[4]*boundary_tensors[5]*bond_tensors[3],site_tensors[4]*bond_tensors[2],site_tensors[5]*boundary_tensors[6]*boundary_tensors[7]*bond_tensors[4]};
    IQTensor wf_tensor=combined_site_tensors[0]*combined_site_tensors[1]*combined_site_tensors[2]*combined_site_tensors[3]*combined_site_tensors[4]*combined_site_tensors[5];
    wf_tensor.clean();
    PrintDat(wf_tensor);

    //Spin operators
    std::vector<IQTensor> Sx(6), iSy(6), Sz(6);
    for (int sitei=0; sitei<phys_legs.size(); sitei++)
    {
        auto ind=phys_legs[sitei];
        Sx[sitei]=IQTensor(dag(ind)(1),prime(ind)(2));
        Sx[sitei](dag(ind)(2),prime(ind)(1))=1;
        iSy[sitei]=IQTensor(dag(ind)(1),prime(ind)(2));
        iSy[sitei](dag(ind)(2),prime(ind)(1))=-1;
        Sz[sitei]=IQTensor(dag(ind)(1),prime(ind)(1));
        Sz[sitei](dag(ind)(2),prime(ind)(2))=-1;

        //PrintDat(Sx[sitei]);
        //PrintDat(iSy[sitei]);
        //PrintDat(Sz[sitei]);
    }

    //calculate correlator
    double wf_norm=(wf_tensor*dag(wf_tensor)).toReal();
    cout << "wavefunction norm:" << endl << wf_norm << endl << endl;
    std::vector<int> bond0{0,1}, bond1{1,2}, bond2{1,4}, bond3{3,4},bond4{4,5};

    cout << "SxSx correlators:" << endl;
    cout << ((wf_tensor*Sx[bond0[0]]*Sx[bond0[1]]).noprime()*dag(wf_tensor)/wf_norm).toReal() << endl;
    cout << ((wf_tensor*Sx[bond1[0]]*Sx[bond1[1]]).noprime()*dag(wf_tensor)/wf_norm).toReal() << endl;
    cout << ((wf_tensor*Sx[bond2[0]]*Sx[bond2[1]]).noprime()*dag(wf_tensor)/wf_norm).toReal() << endl;
    cout << ((wf_tensor*Sx[bond3[0]]*Sx[bond3[1]]).noprime()*dag(wf_tensor)/wf_norm).toReal() << endl;
    cout << ((wf_tensor*Sx[bond4[0]]*Sx[bond4[1]]).noprime()*dag(wf_tensor)/wf_norm).toReal() << endl;
    return 0;
}
