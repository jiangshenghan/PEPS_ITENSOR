
#include "transfer_to_square.h"
#include "double_layer_peps.h"
//#include "boundary_theory.h"

int main()
{
    //std::default_random_engine generator(std::time(0));
    //std::uniform_real_distribution<double> distribution(-1,1);
    //auto rand_gen = std::bind(distribution,generator);
    //double A1=rand_gen(), A2=rand_gen();
    double A1=1, A2=1;

    std::vector<IQTensor> honeycomb_site_tensors_uc, honeycomb_bond_tensors_uc;
    generate_featureless_honeycomb_ansatz_uc(A1,A2,honeycomb_site_tensors_uc,honeycomb_bond_tensors_uc);

    Square_Lattice_Cylinder square_cylinder(std::array<int,2>{38,2});
    auto square_peps=spin_sym_square_peps_from_honeycomb_tensor_uc(honeycomb_site_tensors_uc,honeycomb_bond_tensors_uc,square_cylinder);

    //set boundary tensors
    IQIndex boundary_leg=Spin_leg(std::vector<int>{0,2},"boundary leg",In,Link);
    IQTensor boundary_tensor(boundary_leg(2));
    boundary_tensor(boundary_leg(1))=1;
    //cout << "Boundary condition:" << endl;
    //PrintDat(boundary_tensor);
    //TODO: div(boundary_tensor) may cause inconsistency for generic boundary condition?
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

    //calculate entanglement entropy and related quantities
    Cylinder_Square_Double_Layer_PEPSt<IQTensor> square_double_layer_peps(square_peps);

    square_double_layer_peps.obtain_sigma_lr_iterative(17,24);
    //we should always decombine sigma_lr_ before writing to file
    cout << "\n========================================\n" << endl;
    cout << "Data for original double layer PEPS:" << endl;
    cout << "col_left=" << square_double_layer_peps.col_lr(0) << ", col_right=" << square_double_layer_peps.col_lr(1) << endl;
    PrintDat(square_double_layer_peps.sigma_lr(0));
    PrintDat(square_double_layer_peps.sigma_lr(1));
    cout << "\n========================================\n" << endl;

    square_double_layer_peps.decombine_sigma_lr();
    writeToFile("/home/jiangsb/code/peps_itensor/result/test/iotest.txt",square_double_layer_peps);
    Cylinder_Square_Double_Layer_PEPSt<IQTensor> double_layer_peps_from_file(square_cylinder);
    readFromFile("/home/jiangsb/code/peps_itensor/result/test/iotest.txt",double_layer_peps_from_file);

    cout << "\n========================================\n" << endl;
    cout << "Data for double layer PEPS reading from file:" << endl;
    cout << "col_left=" << double_layer_peps_from_file.col_lr(0) << ", col_right=" << double_layer_peps_from_file.col_lr(1) << endl;
    PrintDat(double_layer_peps_from_file.sigma_lr(0));
    PrintDat(double_layer_peps_from_file.sigma_lr(1));
    cout << "\n========================================\n" << endl;


    //square_double_layer_peps.obtain_boundary_theory_iterative();
    //cout << "\n========================================\n" << endl;
    //cout << "Square lattice with Lx=" << square_cylinder.n_uc()[0] << " and Ly=" << square_cylinder.n_uc()[1] << " cylinder " << endl;
    //cout << "Density Matrix spectrum: " << endl; 
    //for (double eigval : square_double_layer_peps.density_mat_spectrum()) cout << eigval << " ";
    //cout << endl;
    //cout << "Entanglement entropy: " << square_double_layer_peps.entanglement_entropy_vN() << endl;
    //cout << "\n========================================\n" << endl;


    return 0;
}

//int main()
//{
//    
//    //construct phys legs and virt legs for site tensors and bond tensors of honeycomb lattice in one uc
//    std::vector<IQIndex> phys_legs={Spin_leg(std::vector<int>{0,1},"phys_leg 1",Out,Site), Spin_leg(std::vector<int>{0,1},"phys_leg 2",Out,Site)};
//
//    std::vector<IQIndex> virt_legs;
//    for (int i=0; i<8; i++)
//    {
//        virt_legs.push_back(Spin_leg(std::vector<int>{0,2},nameint("virt_leg ",i)));
//    }
//
//    //obtain singlet basis for two site tensors, then construct symmetric site tensors using these basis
//    //in our gauge, these two site tensors are the same
//    std::vector<std::vector<IQIndex>> site_tensors_indices={{phys_legs[0],virt_legs[0],virt_legs[1],virt_legs[2]},{phys_legs[1],virt_legs[3],virt_legs[4],virt_legs[5]}};
//
//    std::vector<Singlet_Tensor_Basis> site_tensors_basis={Singlet_Tensor_Basis(site_tensors_indices[0]),Singlet_Tensor_Basis(site_tensors_indices[1])};
//
//    //PrintDat(site_tensors_basis[0]);
//
//    //In current convention \widetilde{T}^i_{\alpha\beta\gamma} transform to [8*(i-1)+BinarytoDecimal(\alpha\beta\gamma)]'s elem of T_tilde
//    //Here, i=1,2 labels the fusion channel, \alpha=0,1 labels deg in flavor space
//    //std::default_random_engine generator(std::time(0));
//    //std::uniform_real_distribution<double> distribution(-1,1);
//    //auto rand_gen = std::bind(distribution,generator);
//    //double A1=rand_gen(), A2=rand_gen();
//    //double A1=1, A2=1;
//    double A1=1, A2=1;
//    cout << "params: A1=" << A1 << ", A2=" << A2 << endl;
//
//    std::array<double,2> norms={2.,std::sqrt(12.)};
//    //peps suppose to have Z2 IGG and independ on A1/A2
//    //std::vector<double> T_tilde={0,-0.5*A2*norms[0],-0.5*A2*norms[0],A1*norms[0],A2*norms[0],-0.5*A1*norms[0],-0.5*A1*norms[0],0,0,0.5*A2*norms[1],-0.5*A2*norms[1],0,0,-0.5*A1*norms[1],0.5*A1*norms[1],0};
//    //trivial IGG peps
//    std::vector<double> T_tilde={0,-0.5*A2*norms[0],-0.5*A2*norms[0],0,A2*norms[0],1.5*A1*norms[0],-1.5*A1*norms[0],0,0,0.5*A2*norms[1],-0.5*A2*norms[1],A1*norms[1],0,-0.5*A1*norms[1],-0.5*A1*norms[1],0};
//    
//    std::vector<IQTensor> honeycomb_site_tensors_uc={singlet_tensor_from_basis_params(site_tensors_basis[0],T_tilde),singlet_tensor_from_basis_params(site_tensors_basis[1],T_tilde)};
//
//    cout << "\n---------------------------------\n" << endl;
//    cout << "honeycomb site tensor:" << endl;
//    for (const auto &tensor : honeycomb_site_tensors_uc)
//    {
//        PrintDat(tensor);
//    }
//    cout << "\n---------------------------------\n" << endl;
//
//
//    //obtain singlet basis for three bond tensors, and then construct bond tensors
//    std::vector<std::vector<IQIndex>> bond_tensors_indices={{dag(virt_legs[0]),dag(virt_legs[3])},{dag(virt_legs[5]),dag(virt_legs[6])},{dag(virt_legs[4]),dag(virt_legs[7])}};
//
//    std::vector<Singlet_Tensor_Basis> bond_tensors_basis={Singlet_Tensor_Basis(bond_tensors_indices[0]),Singlet_Tensor_Basis(bond_tensors_indices[1]),Singlet_Tensor_Basis(bond_tensors_indices[2])};
//
//    std::vector<double> B_tilde={std::sqrt(2),0,0,std::sqrt(2)};
//
//    std::vector<IQTensor> honeycomb_bond_tensors_uc;
//
//    for (const auto &basis : bond_tensors_basis)
//    {
//        honeycomb_bond_tensors_uc.push_back(singlet_tensor_from_basis_params(basis,B_tilde));
//    }
//
//    //for (const auto &tensor : honeycomb_bond_tensors_uc)
//    //{
//    //    PrintDat(tensor);
//    //}
//
//    //Combine two site tensors and bond tensors connect them to make a single site tensor for square lattice
//    //the remaining two bond tensors are square bond tensors
//
//    IQTensor square_site_tensor_no_order;
//    square_site_tensor_no_order=honeycomb_site_tensors_uc[0]*honeycomb_bond_tensors_uc[0]*honeycomb_site_tensors_uc[1];
//
//    IQCombiner phys_legs_combiner(phys_legs[0],phys_legs[1]);
//    phys_legs_combiner.init("square_phys_leg",Site,Out);
//    square_site_tensor_no_order=square_site_tensor_no_order*phys_legs_combiner;
//
//    int phys_legs_qn_order=0;
//    std::vector<int> sz_qns;
//    for (const auto &indqn : phys_legs_combiner.right().indices())
//    {
//        sz_qns.push_back(indqn.qn.sz());
//    }
//    for (int i=0; i<sz_qns.size()-1; i++)
//    {
//        if (phys_legs_qn_order==0)
//        {
//            if (sz_qns[i]<sz_qns[i+1]) phys_legs_qn_order=1;
//            if (sz_qns[i]>sz_qns[i+1]) phys_legs_qn_order=-1;
//        }
//
//        if (phys_legs_qn_order==1) assert(sz_qns[i]<sz_qns[i+1]);
//        if (phys_legs_qn_order==-1) assert(sz_qns[i]>sz_qns[i+1]);
//    }
//    
//
//    //get the right order of square_site_tensor, which should set to be phys,1,5,4,2! Now, the order for square_site_tensor_no_order is 4,5,1,2,phys
//    //the right order of site tensor is used to initialize square lattice peps, where we use replaceIndex
//    IQTensor square_site_tensor_uc(phys_legs_combiner.right(),virt_legs[1],virt_legs[5],virt_legs[4],virt_legs[2]);
//    std::vector<int> leg_dims={phys_legs_combiner.right().m(),virt_legs[1].m(),virt_legs[5].m(),virt_legs[4].m(),virt_legs[2].m()};
//    int max_vals=phys_legs_combiner.right().m()*virt_legs[1].m()*virt_legs[5].m()*virt_legs[4].m()*virt_legs[2].m();
//
//    for (int vals=0; vals<max_vals; vals++)
//    {
//        auto val_list=list_from_num(vals,leg_dims);
//        square_site_tensor_uc(phys_legs_combiner.right()(val_list[0]+1),virt_legs[1](val_list[1]+1),virt_legs[5](val_list[2]+1),virt_legs[4](val_list[3]+1),virt_legs[2](val_list[4]+1))=square_site_tensor_no_order(phys_legs_combiner.right()(val_list[0]+1),virt_legs[1](val_list[1]+1),virt_legs[5](val_list[2]+1),virt_legs[4](val_list[3]+1),virt_legs[2](val_list[4]+1));
//    }
//    square_site_tensor_uc.clean();
//
//    //we do not care about order of bond tensor since it at most contributes overall -1
//    std::vector<IQTensor> square_bond_tensors_uc={honeycomb_bond_tensors_uc[1],honeycomb_bond_tensors_uc[2]};
//
//    //for (const auto &bond_tensor : square_bond_tensors_uc)
//    //{
//    //    PrintDat(bond_tensor);
//    //}
//
//    //PrintDat(square_site_tensor_uc);
//    //cout << square_site_tensor_uc.indices();
//    //cout << "norm: " << (square_site_tensor_uc-square_site_tensor_no_order).norm() << endl;
//
//    //construct square peps on cylinder
//    Square_Lattice_Cylinder square_cylinder(std::array<int,2>{40,4});
//    IQPEPS_IndexSet_Spin_Sym square_indexset(4,4,std::vector<int>{1,0,1},std::vector<int>{0,2},square_cylinder,phys_legs_qn_order);
//    IQPEPS square_peps(square_cylinder,square_indexset,std::vector<IQTensor>{square_site_tensor_uc},square_bond_tensors_uc);
//
//    //square_cylinder.print_lattice_inf();
//
//    //for (const auto &leg : square_indexset.phys_legs()) cout << leg;
//    //for (const auto &leg : square_indexset.virt_legs()) cout << leg;
//
//    //set boundary tensors
//    IQIndex boundary_leg=Spin_leg(std::vector<int>{0,2},"boundary leg",In,Link);
//    IQTensor boundary_tensor(boundary_leg(1));
//    boundary_tensor(boundary_leg(2))=1;
//    cout << "Boundary condition:" << endl;
//    PrintDat(boundary_tensor);
//    //TODO: div(boundary_tensor) may cause inconsistency for generic boundary condition?
//    int boundary_i=0;
//    for (auto &tensor : square_peps.boundary_tensors())
//    {
//        //cout << tensor;
//        auto oind=boundary_tensor.indices()[0],
//             nind=tensor.indices()[0];
//        if (oind.dir()==-nind.dir())
//        {
//            oind.dag();
//            boundary_tensor.dag();
//        }
//        boundary_tensor.replaceIndex(oind,nind);
//        tensor=boundary_tensor;
//    }
//
//    //cout << "square lattice site tensors: " << endl;
//    //for (const auto &tensor : square_peps.site_tensors()) PrintDat(tensor); 
//    //cout << "square lattice bond tensors: " << endl;
//    //for (const auto &tensor : square_peps.bond_tensors()) PrintDat(tensor);
//    //cout << "square lattice boundary tensors: " << endl;
//    //for (const auto &tensor : square_peps.boundary_tensors()) PrintDat(tensor);
//
//
//    //Double_Layer_PEPSt<IQTensor> square_double_layer_peps(square_peps);
//    //cout << "\n---------------------------------\n" << endl;
//    //cout << "double layer tensors of suqare lattice: " << endl;
//    //for (const auto &tensor: square_double_layer_peps.layered_site_tensors()) PrintDat(tensor);
//    //cout << div(square_double_layer_peps.layered_site_tensors(0));
//    //cout << "\n---------------------------------\n" << endl;
//    //return 0;
//
//    //calculate entanglement entropy
//    Boundary_Theory<IQTensor> square_boundary_theory(square_peps);
//
//    cout << "Square lattice with Lx=" << square_cylinder.n_uc()[0] << " and Ly=" << square_cylinder.n_uc()[1] << " cylinder " << endl;
//    cout << "Density Matrix spectrum: " << endl; 
//    for (double eigval : square_boundary_theory.density_mat_spectrum()) cout << eigval << " ";
//    cout << endl;
//    cout << "Entanglement entropy: " << square_boundary_theory.entanglement_entropy_vN() << endl;
//
//    return 0;
//}

