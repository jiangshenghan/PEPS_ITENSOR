
#include "square_double_layer_peps.h"
#include "simple_update.h"

using namespace square_psg;

int main()
{
    mu_12=1;
    int D=6;
    int Lx=8, Ly=8;
    Square_Lattice_Torus square_lattice({Lx,Ly});

    //reading PEPS from file
    std::stringstream ss;
    if (std::abs(mu_12-1)<EPSILON)
    {
        ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_rvb_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_iter=2.txt";
    }
    if (std::abs(mu_12+1)<EPSILON)
    {
        ss << "/home/jiangsb/code/peps_itensor/result/peps_storage/square_pi_rvb_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_iter=2.txt";
    }
    std::string file_name=ss.str();
    IQPEPS square_peps_old(square_lattice);
    readFromFile(file_name,square_peps_old);

    //PrintDat(square_peps_old.bond_tensors(0));
    //PrintDat(square_peps_old.bond_tensors(2));

    //translate PEPS to another size
    Ly=2;
    Lx=32*Ly;
    Square_Lattice_Cylinder square_cylinder({Lx,Ly});
    IQPEPS_IndexSet_SpinHalf index_set(D,square_cylinder);

    std::vector<IQTensor> site_tensors_uc, bond_tensors_uc;
    for (int sitei=0; sitei<square_lattice.n_sites_uc(); sitei++) site_tensors_uc.push_back(square_peps_old.site_tensors(sitei).takeRealPart());
    for (int bondi=0; bondi<square_lattice.n_bonds_uc(); bondi++) bond_tensors_uc.push_back(square_peps_old.bond_tensors(bondi).takeRealPart());

    IQPEPS square_peps_cylinder(square_cylinder,index_set,site_tensors_uc,bond_tensors_uc,mu_12);

    //set boundary condition
    IQIndex boundary_leg=Spin_leg({1,1,1},"boundary leg",In,Link);
    IQTensor boundary_tensor(boundary_leg);
    boundary_tensor(boundary_leg(1))=1;
    for (auto &tensor : square_peps_cylinder.boundary_tensors())
    {
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
    cout << "Successfully construct square_peps_cylinder!" << endl;

    //obtain double layer peps and calculate sigma_lr_
    Square_Double_Layer_PEPSt<IQTensor> square_double_layer_peps(square_peps_cylinder);
    square_double_layer_peps.obtain_sigma_lr_iterative(Lx/2-1,Lx/2);
    cout << "Obtain boundary theory!" << endl;

    //construct SiSi operator
    IQIndex ope_leg=Spin_leg({0,1},"ope_leg0",Out,Site);
    IQTensor Sx(dag(ope_leg(1)),prime(ope_leg(2))),
            iSy(dag(ope_leg(2)),prime(ope_leg(1))),
            Sz(dag(ope_leg(1)),prime(ope_leg(1)));

    Sx(dag(ope_leg(2)),prime(ope_leg(1)))=1;
    iSy(dag(ope_leg(1)),prime(ope_leg(2)))=-1;
    Sz(dag(ope_leg(2)),prime(ope_leg(2)))=-1;

    PrintDat(Sx);
    PrintDat(iSy);
    PrintDat(Sz);

    //single bond correlation
    int origin_site=4;
    Print(square_double_layer_peps.obtain_correlators({origin_site,origin_site+Lx},{Sx,Sx}));
    Print(square_double_layer_peps.obtain_correlators({origin_site,origin_site+Lx},{iSy,-iSy}));
    Print(square_double_layer_peps.obtain_correlators({origin_site,origin_site+Lx},{Sz,Sz}));

    //for (const auto &tensor : square_peps_new.bond_tensors()) PrintDat(tensor);

    //stores square_peps_new as tnetwork_storage
    //Tnetwork_Storage<IQTensor> square_peps_storage=peps_to_tnetwork_storage(square_peps_new);
    //ss.str(std::string());
    //ss.clear();
    //if (std::abs(mu_12-1)<EPSILON)
    //{
    //    ss << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_rvb_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_iter=2_step=999.txt";
    //}
    //if (std::abs(mu_12+1)<EPSILON)
    //{
    //    ss << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_pi_rvb_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_iter=2_step=999.txt";
    //}
    //file_name=ss.str();
    //writeToFile(file_name,square_peps_storage);

    return 0;
}
