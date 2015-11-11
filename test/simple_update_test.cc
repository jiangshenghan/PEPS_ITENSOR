
//#include "square_rvb.h"
#include "simple_update.h"

int main()
{
    int Lx=8, Ly=8;

    //Square_Lattice_Torus square_lattice({Lx,Ly});
    //IQPEPS_IndexSet_SpinHalf index_set(6,square_lattice);
    //IQPEPS square_peps(square_lattice,index_set);
    //random_init_square_rvb_peps(square_peps);

    IQPEPS square_peps=square_srvb_peps(Lx,Ly);
    
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
    //auto site_tensorA=square_peps.site_tensors(0),
    //     site_tensorB=square_peps.site_tensors(1);
    //for (int neighi=0; neighi<square_peps.n_bonds_to_one_site(); neighi++)
    //{
    //    int bondi=square_peps.lattice().site_neighbour_bonds(0,neighi);
    //    site_tensorA*=square_peps.bond_tensors(bondi);
    //}
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
    Evolution_Params square_su_params(1,{1},{0.1});
    spin_square_peps_simple_update(square_peps,square_su_params);

    //Check the output result
    //std::vector<IQTensor> tens_step;
    //for (int step=0; step<20; step++)
    //{
    //    Tnetwork_Storage<IQTensor> square_rvb_from_file;
    //    std::stringstream ss;
    //    ss << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_rvb_D=3" << "_Lx=" << Lx << "_Ly=" << Ly << "_iter=0_step=" << step << ".txt";
    //    std::string file_name=ss.str();
    //    readFromFile(file_name,square_rvb_from_file);
    //    //Print(square_rvb_from_file._Lx);
    //    //Print(square_rvb_from_file._Ly);
    //    //PrintDat(square_rvb_from_file._tensor_list);
    //    tens_step.push_back(square_rvb_from_file._tensor_list(0));
    //    Print(step);
    //    //PrintDat(tens_step.back());
    //    Print(tens_step.back().norm());
    //    if (step>0) 
    //    {
    //        Print(tens_step[step-1].norm());
    //        Print((tens_step[step]-tens_step[step-1]).norm());
    //        //Print((tens_step[step]+tens_step[step-1]).norm());
    //    }
    //}

    return 0;
}
