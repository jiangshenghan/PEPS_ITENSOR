
#include "transfer_to_square.h"
#include "square_double_layer_peps.h"

int main()
{
    //system size
    int Ly=3,Lx=24*Ly;

    //tunable parameter for wavefunctions
    std::default_random_engine generator(std::time(0));
    std::uniform_real_distribution<double> distribution(-1,1);
    auto rand_gen = std::bind(distribution,generator);
    //double A1=rand_gen(), A2=rand_gen();
    double A1=0.9, A2=1-A1;

    //Output file name for double_layer_peps
    std::stringstream ss;
    //file name for cylinder
    ss << "/home/jiangsb/code/peps_itensor/result/featureless_honeycomb_peps_Ly=" << Ly << "/Lx=" << Lx << "_A1=" << A1 << "_A2=" << A2 << "_double_layer_peps.txt";
    //file name for ribbon
    //ss << "/home/jiangsb/code/peps_itensor/result/featureless_honeycomb_peps_ribbon" << "/Ly=" << Ly << "_Lx=" << Lx << "_A1=" << A1 << "_A2=" << A2 << "_double_layer_peps.txt";
    std::string file_name=ss.str();

    //Input honeycomb tensor in one uc and transfer to square peps
    std::vector<IQTensor> honeycomb_site_tensors_uc, honeycomb_bond_tensors_uc;
    generate_featureless_honeycomb_ansatz_uc(A1,A2,honeycomb_site_tensors_uc,honeycomb_bond_tensors_uc);

    //peps on square cylinder
    Square_Lattice_Cylinder square_cylinder(std::array<int,2>{Lx,Ly});
    auto square_peps=spin_sym_square_peps_from_honeycomb_tensor_uc(honeycomb_site_tensors_uc,honeycomb_bond_tensors_uc,square_cylinder);

    //peps on square ribbon
    //Square_Lattice_Open square_ribbon(std::array<int,2>{Lx,Ly});
    //auto square_peps=spin_sym_square_peps_from_honeycomb_tensor_uc(honeycomb_site_tensors_uc,honeycomb_bond_tensors_uc,square_ribbon);

    cout << "\n========================================\n" << endl;
    cout << square_peps.name() << endl;
    cout << "System Size: " << Lx << "x" << Ly << endl;
    cout << "Wavefunction Params: " << "A1=" << A1 << ", A2=" << A2 << endl;
    cout << "Output file: " << endl << file_name << endl;
    cout << "\n========================================\n" << endl;


    //set ferromagnet boundary condition
    IQIndex boundary_leg=Spin_leg(std::vector<int>{0,2},"boundary leg",In,Link);
    IQTensor boundary_tensor(boundary_leg);
    boundary_tensor(boundary_leg(1))=1;
    boundary_tensor(boundary_leg(2))=1-boundary_tensor(boundary_leg(1));
    //boundary_tensor(boundary_leg(3))=rand_gen();
    //boundary_tensor(boundary_leg(4))=rand_gen();
    cout << "Boundary condition:" << endl;
    PrintDat(boundary_tensor);
    for (auto &tensor : square_peps.boundary_tensors())
    {
        //cout << tensor;
        auto oind=boundary_tensor.indices()[0],
             nind=tensor.indices()[0];
        if (oind.dir()==-nind.dir())
        {
            oind.dag();
            boundary_tensor.dag();
        }
        boundary_tensor.replaceIndex(oind,nind);
        tensor=boundary_tensor;
        //PrintDat(tensor);
    }

    //set anti-ferromagnet boundary condition
    //int boundaryi=0;
    //for (auto &tensor : square_peps.boundary_tensors())
    //{
    //    IQIndex ind=tensor.indices()[0];
    //    if (boundaryi%2==0)
    //    {
    //        tensor(ind(1))=1;
    //        tensor(ind(2))=1;
    //    }
    //    else
    //    {
    //        tensor(ind(3))=1;
    //        tensor(ind(4))=1;
    //    }
    //    //PrintDat(tensor);
    //    boundaryi++;
    //}

    //set random boundary condition for short boundary for ribbon geometry
    //for (int boundaryi=0; boundaryi<Ly; boundaryi++)
    //{
    //    auto &left_tensor=square_peps.boundary_tensors(boundaryi);
    //    auto &right_tensor=square_peps.boundary_tensors(Lx+Ly+boundaryi);
    //    auto left_ind=left_tensor.indices()[0],
    //         right_ind=right_tensor.indices()[0];
    //    for (int val=1; val<=left_ind.m(); val++)
    //    {
    //        left_tensor(left_ind(val))=rand_gen();
    //        right_tensor(right_ind(val))=rand_gen();
    //    }
    //}

    
    //set random boundary condition
    //for (auto &tensor: square_peps.boundary_tensors())
    //{
    //    IQIndex ind=tensor.indices()[0];
    //    for (int val=1; val<=ind.m(); val++)
    //    {
    //        tensor(ind(val))=rand_gen();
    //    }
    //    //PrintDat(tensor);
    //}


    Square_Double_Layer_PEPSt<IQTensor> square_double_layer_peps(square_peps);

    //check transfer matrix
    //square_double_layer_peps.obtain_transfer_matrix(2);
    //auto transfer_mat_eigvals=square_double_layer_peps.transfer_matrix_eigvals();
    //for (auto &eigval : transfer_mat_eigvals)
    //    cout << eigval << " ";
    //cout << endl;
    //return 0;

    //calculate sigma_lr_
    square_double_layer_peps.obtain_sigma_lr_iterative(Lx/2-1,Lx/2);

    //cout << "\n========================================\n" << endl;
    //cout << "Data for original double layer PEPS:" << endl;
    //cout << "col_left=" << square_double_layer_peps.col_lr(0) << ", col_right=" << square_double_layer_peps.col_lr(1) << endl;
    //PrintDat(square_double_layer_peps.sigma_lr(0));
    //PrintDat(square_double_layer_peps.sigma_lr(1));
    //cout << "\n========================================\n" << endl;

    //we should always decombine sigma_lr_ before writing to file
    square_double_layer_peps.decombine_sigma_lr();
    writeToFile(file_name,square_double_layer_peps);
    square_double_layer_peps.recombine_sigma_lr();

    //get entanglement property
    square_double_layer_peps.move_sigma_lr({Lx/2-1,Lx/2});
    square_double_layer_peps.from_sigma_lr_to_sigma_b();
    square_double_layer_peps.obtain_density_matrix_spectrum();
    //cout << "Density Matrix spectrum: " << endl; 
    //for (double eigval : square_double_layer_peps.density_mat_spectrum()) cout << eigval << " ";
    //cout << endl;
    cout << "Entanglement entropy: " << square_double_layer_peps.entanglement_entropy_vN() << endl;

    //Square_Double_Layer_PEPSt<IQTensor> double_layer_peps_from_file(square_cylinder);
    //readFromFile(file_name,double_layer_peps_from_file);

    //cout << "\n========================================\n" << endl;
    //cout << "Data for double layer PEPS reading from file:" << endl;
    //cout << "col_left=" << double_layer_peps_from_file.col_lr(0) << ", col_right=" << double_layer_peps_from_file.col_lr(1) << endl;
    //PrintDat(double_layer_peps_from_file.sigma_lr(0));
    //PrintDat(double_layer_peps_from_file.sigma_lr(1));
    //cout << "\n========================================\n" << endl;


    //square_double_layer_peps.obtain_boundary_theory_iterative();
    //cout << "\n========================================\n" << endl;
    //cout << "Square lattice with Lx=" << square_cylinder.n_uc()[0] << " and Ly=" << square_cylinder.n_uc()[1] << " cylinder " << endl;
    //cout << "A1=" << A1 << ", A2=" << A2 << endl;
    //cout << "Density Matrix spectrum: " << endl; 
    //for (double eigval : square_double_layer_peps.density_mat_spectrum()) cout << eigval << " ";
    //cout << endl;
    //cout << "Entanglement entropy: " << square_double_layer_peps.entanglement_entropy_vN() << endl;
    //cout << "\n========================================\n" << endl;

    cout << "Finish Successfully!" << endl;
    return 0;
}


