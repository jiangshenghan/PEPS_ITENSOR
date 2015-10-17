
#include "transfer_to_square.h"
#include "double_layer_peps.h"

int main()
{
    //system size
    int Ly=3, Lx=12*Ly;

    double A1=1, A2=1;

    //Output file name for double_layer_peps
    std::stringstream ss;
    ss << "/home/jiangsb/code/peps_itensor/result/featureless_honeycomb_peps_Ly=" << Ly << "/A1=" << A1 << "_A2=" << A2 << "_double_layer_peps.txt";
    std::string file_name=ss.str();

    cout << "\n========================================\n" << endl;
    cout << "System Size: " << Lx << "x" << Ly << endl;
    cout << "Wavefunction Params: " << "A1=" << A1 << ", A2=" << A2 << endl;
    cout << "Input file: " << endl << file_name << endl;

    //Input honeycomb tensor in one uc and transfer to square peps
    std::vector<IQTensor> honeycomb_site_tensors_uc, honeycomb_bond_tensors_uc;
    generate_featureless_honeycomb_ansatz_uc(A1,A2,honeycomb_site_tensors_uc,honeycomb_bond_tensors_uc);

    Square_Lattice_Cylinder square_cylinder(std::array<int,2>{Lx,Ly});
    auto square_peps=spin_sym_square_peps_from_honeycomb_tensor_uc(honeycomb_site_tensors_uc,honeycomb_bond_tensors_uc,square_cylinder);

    //set boundary tensors
    IQIndex boundary_leg=Spin_leg(std::vector<int>{0,2},"boundary leg",In,Link);
    IQTensor boundary_tensor(boundary_leg(2));
    boundary_tensor(boundary_leg(1))=1;

    int boundary_i=0;
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
    }

    cout << "Boundary condition:" << endl;
    PrintDat(boundary_tensor);
    cout << "========================================\n" << endl;

    //reading square_cylinder_double_peps
    Cylinder_Square_Double_Layer_PEPSt<IQTensor> square_cylinder_double_peps(square_cylinder);
    readFromFile(file_name,square_cylinder_double_peps);

    //calculate correlators
    std::array<Coordinate,2> acting_sites_coord;
    acting_sites_coord[0]=Coordinate{1,0,0};
    acting_sites_coord[1]=Coordinate{2,0,0};

    cout << honeycomb_cylinder_SzSz_correlator(acting_sites_coord,square_cylinder_double_peps) << endl;

    return 0;
}
