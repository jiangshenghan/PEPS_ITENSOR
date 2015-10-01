//TODO: 
//2. check simple_update_env
//3. check simple_update algorithm


#include "simple_update.h"

int main()
{
    //Check for simple_update_env
    Square_Lattice_Torus square_lattice{std::array<int,2>{2,2}};
   

    IQPEPS_IndexSet_SpinHalf index_set(3,square_lattice);
    IQPEPS square_peps(square_lattice,index_set);
    
    randomize_spin_sym_square_peps(square_peps);

    for (const auto &tensor : square_peps.site_tensors())
    {
        PrintDat(tensor);
    }

    for (const auto &tensor : square_peps.bond_tensors())
    {
        PrintDat(tensor);
    }
    
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
