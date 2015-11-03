
#include "square_rvb.h"
#include "simple_update_env.h"

int main()
{
    //Check for simple_update_env
    int Lx=4, Ly=4;

    Square_Lattice_Torus square_lattice({Lx,Ly});
    IQPEPS_IndexSet_SpinHalf index_set(10,square_lattice);
    IQPEPS square_peps(square_lattice,index_set);
    random_init_square_rvb_peps(square_peps);

    cout << "Check symmetry operation:" << endl;

    //IQPEPS square_peps=square_srvb_peps(Lx,Ly);
    
    //for (const auto &tensor : square_peps.site_tensors())
    //{
    //    PrintDat(tensor);
    //}

    //for (const auto &tensor : square_peps.bond_tensors())
    //{
    //    PrintDat(tensor);
    //}

    //Check environment update
    std::array<std::vector<IQTensor>,2> env_tensors;
    auto site_tensorA=square_peps.site_tensors(0)*square_peps.bond_tensors(1),
         site_tensorB=square_peps.site_tensors(1);
    get_env_tensor(site_tensorA,site_tensorB,env_tensors);
    
    //Check for Hamiltonian
    //std::array<IQIndex,2> sites{{Spin_leg({0,1},"site1",Out),Spin_leg({0,1},"site2",Out)}};
    //NN_Heisenberg_Hamiltonian hamiltonian_gate(sites,0.2);

    //cout << "Site tensors:" << endl;
    //PrintDat(hamiltonian_gate.site_tensors(0));
    //PrintDat(hamiltonian_gate.site_tensors(1));

    //cout << "Bond tensor:" << endl; 
    //PrintDat(hamiltonian_gate.bond_tensor());
    //cout << "J: " << hamiltonian_gate.J() << endl;


    //Check for trotter_gate
    //std::array<IQIndex,2> sites{{Spin_leg({0,1},"site1",Out),Spin_leg({0,1},"site2",Out)}};
    //NN_Heisenberg_Trotter_Gate evolve_gate(sites,0.1);

    //cout << "Site tensors:" << endl;
    //PrintDat(evolve_gate.site_tensors(0));
    //PrintDat(evolve_gate.site_tensors(1));

    //cout << "Bond tensor:" << endl; 
    //PrintDat(evolve_gate.bond_tensors(0));
    //cout << "Imag time: " << evolve_gate.t() << endl;


    return 0;
}
