
//#include "square_rvb.h"
#include "simple_update.h"

using namespace square_psg;


int main()
{
    //control PSG parameters
    mu_12=1; 


    //init lattice
    int Lx=8, Ly=8;
    Square_Lattice_Torus square_lattice({Lx,Ly});


    //construct random peps with a good initial state
    int D=6;
    IQPEPS_IndexSet_SpinHalf index_set(D,square_lattice);
    IQPEPS square_peps(square_lattice,index_set);

    double init_energy=0;
    do
    {
        //get "energy" of init state by env_tens
        random_init_square_rvb_peps(square_peps);
        std::array<IQIndex,2> site01_legs{square_peps.phys_legs(0),square_peps.phys_legs(1)};
        NN_Heisenberg_Hamiltonian hamiltonian_gate(site01_legs);
        std::array<std::vector<IQTensor>,2> env_tens;

        int comm_bond=square_peps.lattice().comm_bond(0,1);
        auto comm_bond_tensor=square_peps.bond_tensors(comm_bond);

        auto combined_site_tens0=square_peps.site_tensors(0);
        for (int neighi=0; neighi<square_peps.n_bonds_to_one_site(); neighi++)
        {
            int bondi=square_peps.lattice().site_neighbour_bonds(0,neighi);
            if (bondi==comm_bond) continue;
            combined_site_tens0*=square_peps.bond_tensors(bondi);
        }
        get_env_tensor_minimization(combined_site_tens0*comm_bond_tensor,square_peps.site_tensors(1),env_tens);

        std::array<IQTensor,2> site_env_tens{{combined_site_tens0,square_peps.site_tensors(1)}}; 
        for (int sitei=0; sitei<2; sitei++)
        {
            for (const auto &env_leg_tensor : env_tens[sitei])
            {
                site_env_tens[sitei]*=env_leg_tensor;
            }
            site_env_tens[sitei].noprime();
        }

        init_energy=heisenberg_energy_from_site_env_tensors(site_env_tens,comm_bond_tensor,hamiltonian_gate);
        Print(init_energy);
    }
    while (init_energy>-0.1);

    //construct short-range rvb
    //IQPEPS square_peps=square_srvb_peps(Lx,Ly);
    

    //construct PEPS from file
    //std::stringstream ss;

    ////zero-flux state
    //ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_rvb_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_optimized";
    ////pi flux state
    ////ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_pi_rvb_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_optimized";

    //std::string file_name=ss.str();
    //IQPEPS square_peps(square_lattice);
    //readFromFile(file_name,square_peps);


    //Check symmetry
    //Print(diff_tensor_and_singlet_projection(square_peps.site_tensors(0)));
    //Print(diff_tensor_and_singlet_projection(square_peps.bond_tensors(0)));
    //Print(diff_tensor_and_singlet_projection(square_peps.bond_tensors(1)));


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
    

    //Check for trotter_gate and hamiltonian gate
    //std::array<IQIndex,2> site_legs{Spin_leg({0,1},"site1",Out),Spin_leg({0,1},"site2",Out)};
    //NN_Heisenberg_Hamiltonian hamiltonian_gate(site_legs);
    //NN_Heisenberg_Trotter_Gate evolve_gate(site_legs,1);

    //PrintDat(hamiltonian_gate.site_tensors(0)*hamiltonian_gate.bond_tensor()*hamiltonian_gate.site_tensors(1));
    //PrintDat(evolve_gate.site_tensors(0)*evolve_gate.bond_tensors(0)*evolve_gate.site_tensors(1));

    //cout << "Trotter gate site tensors:" << endl;
    //PrintDat(evolve_gate.site_tensors(0));
    //PrintDat(evolve_gate.site_tensors(1));

    //cout << "Trotter gate bond tensor:" << endl; 
    //PrintDat(evolve_gate.bond_tensors(0));
    //cout << "Imag time: " << evolve_gate.t() << endl;


    //optimazation
    Evolution_Params square_su_params(6,{20,50,300,1000,3000,100000},{1,1e-1,1e-2,1e-3,1e-4,1e-5});
    //Evolution_Params square_su_params(1,{25},{1e-0});
    spin_square_peps_simple_update(square_peps,square_su_params);


    //Check the output result
    //for (int step=0; step<20; step++)
    //{
    //    Print(step);
    //    Tnetwork_Storage<IQTensor> square_rvb_from_file;
    //    std::stringstream ss;
    //    ss << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_rvb_D=6" << "_Lx=" << Lx << "_Ly=" << Ly << "_iter=0_step=" << step;
    //    std::string file_name=ss.str();
    //    readFromFile(file_name,square_rvb_from_file);
    //    //Print(square_rvb_from_file._Lx);
    //    //Print(square_rvb_from_file._Ly);
    //    //PrintDat(square_rvb_from_file._tensor_list);

    //    IQTensor site_tensA=square_rvb_from_file._tensor_list(0),
    //             site_tensB=square_rvb_from_file._tensor_list(1);
    //    std::array<IQIndex,2> site_legs{findtype(site_tensA,Site),findtype(site_tensB,Site)};
    //    std::array<std::vector<IQTensor>,2> env_tens;
    //    get_env_tensor_minimization(site_tensA,site_tensB,env_tens);
    //    auto site_env_tens=std::array<IQTensor,2>{site_tensA,site_tensB};
    //    NN_Heisenberg_Hamiltonian hamiltonian_gate(site_legs);

    //    for (int sitei=0; sitei<2; sitei++) 
    //    {
    //        for (auto &env_tensor_one_leg : env_tens[sitei])
    //            site_env_tens[sitei]*=env_tensor_one_leg;
    //    }

    //    Print(heisenberg_energy_from_site_env_tensors(site_env_tens,hamiltonian_gate));
    //}

    return 0;
}
