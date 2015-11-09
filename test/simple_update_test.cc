
//#include "square_rvb.h"
#include "simple_update.h"

int main()
{
    int Lx=4, Ly=4;

    Square_Lattice_Torus square_lattice({Lx,Ly});
    IQPEPS_IndexSet_SpinHalf index_set(3,square_lattice);
    IQPEPS square_peps(square_lattice,index_set);
    random_init_square_rvb_peps(square_peps);

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
    //std::array<std::vector<IQTensor>,2> env_tensors;
    //auto site_tensorA=square_peps.site_tensors(0)*square_peps.bond_tensors(1),
    //     site_tensorB=square_peps.site_tensors(1);
    //get_env_tensor_minimization(site_tensorA,site_tensorB,env_tensors);
    
    //Check for trotter_gate
    //std::array<IQIndex,2> site_legs{Spin_leg({0,1},"site1",Out),Spin_leg({0,1},"site2",Out)};
    //NN_Heisenberg_Trotter_Gate evolve_gate(site_legs,0.1);

    //cout << "Site tensors:" << endl;
    //PrintDat(evolve_gate.site_tensors(0));
    //PrintDat(evolve_gate.site_tensors(1));

    //cout << "Bond tensor:" << endl; 
    //PrintDat(evolve_gate.bond_tensors(0));
    //cout << "Imag time: " << evolve_gate.t() << endl;

    //Check for optimazation
    //Evolution_Params square_su_params(1,{10},{0.1});
    //spin_square_peps_simple_update(square_peps,square_su_params);

    return 0;
}
