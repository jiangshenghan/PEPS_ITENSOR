
#include "transfer_to_square.h"
#include "square_double_layer_peps.h"

int main()
{
    //system size
    int Ly=2, Lx=40*Ly;

    double A1=0.9, A2=1-A1;

    //Output file name for double_layer_peps
    std::stringstream ss;
    //file name for cylinder
    //ss << "/home/jiangsb/code/peps_itensor/result/featureless_honeycomb_peps_Ly=" << Ly << "/Lx=" << Lx << "_A1=" << A1 << "_A2=" << A2 << "_double_layer_peps.txt";
    //file name for ribbon
    ss << "/home/jiangsb/code/peps_itensor/result/featureless_honeycomb_peps_ribbon" << "/Ly=" << Ly << "_Lx=" << Lx << "_A1=" << A1 << "_A2=" << A2 << "_double_layer_peps.txt";

    std::string file_name=ss.str();

    //cylinder geometry peps
    //Square_Lattice_Cylinder square_cylinder(std::array<int,2>{Lx,Ly});
    //Square_Double_Layer_PEPSt<IQTensor> square_double_layer_peps(square_cylinder);

    //ribbon geometry peps
    Square_Lattice_Open square_ribbon(std::array<int,2>{Lx,Ly});
    Square_Double_Layer_PEPSt<IQTensor> square_double_layer_peps(square_ribbon);

    readFromFile(file_name,square_double_layer_peps);
    
    cout << "\n========================================\n" << endl;
    cout << square_double_layer_peps.lattice().name() << endl;
    cout << "System Size: " << Lx << "x" << Ly << endl;
    cout << "Wavefunction Params: " << "A1=" << A1 << ", A2=" << A2 << endl;
    cout << "Input file: " << endl << file_name << endl;

    cout << "Successfully reading from file!" << endl;
    cout << "========================================\n" << endl;

    //cout << "\n========================================\n" << endl;
    //cout << "Data for double layer PEPS reading from file:" << endl;
    //cout << "col_left=" << square_double_layer_peps.col_lr(0) << ", col_right=" << square_double_layer_peps.col_lr(1) << endl;
    //PrintDat(square_double_layer_peps.sigma_lr(0));
    //PrintDat(square_double_layer_peps.sigma_lr(1));
    //cout << "\n========================================\n" << endl;


    //construct Si (i=x,y,z) operators of honeycomb lattice and then transfer to square
    //there are three cases: SiId, IdSi, SiSi
    std::array<IQIndex,2> honeycomb_legs={Spin_leg(std::vector<int>{0,1}, "honeycomb_leg 0", Out, Site), Spin_leg(std::vector<int>{0,1},"honeycomb_leg 1",Out,Site)};
    std::array<IQTensor,2> Id_honeycomb_uc={IQTensor(dag(honeycomb_legs[0]),prime(honeycomb_legs[0])),IQTensor(dag(honeycomb_legs[1]),prime(honeycomb_legs[1]))}; 
    auto Sx_honeycomb_uc=Id_honeycomb_uc, 
         iSy_honeycomb_uc=Id_honeycomb_uc,
         Sz_honeycomb_uc=Id_honeycomb_uc;
    for (int sitei=0; sitei<2; sitei++)
    {
        Id_honeycomb_uc[sitei](dag(honeycomb_legs[sitei])(1),prime(honeycomb_legs[sitei])(1))=1;
        Id_honeycomb_uc[sitei](dag(honeycomb_legs[sitei])(2),prime(honeycomb_legs[sitei])(2))=1;

        Sx_honeycomb_uc[sitei](dag(honeycomb_legs[sitei])(1),prime(honeycomb_legs[sitei])(2))=1;
        Sx_honeycomb_uc[sitei](dag(honeycomb_legs[sitei])(2),prime(honeycomb_legs[sitei])(1))=1;

        iSy_honeycomb_uc[sitei](dag(honeycomb_legs[sitei])(1),prime(honeycomb_legs[sitei])(2))=1;
        iSy_honeycomb_uc[sitei](dag(honeycomb_legs[sitei])(2),prime(honeycomb_legs[sitei])(1))=-1;

        Sz_honeycomb_uc[sitei](dag(honeycomb_legs[sitei])(1),prime(honeycomb_legs[sitei])(1))=1;
        Sz_honeycomb_uc[sitei](dag(honeycomb_legs[sitei])(2),prime(honeycomb_legs[sitei])(2))=-1;

        Sx_honeycomb_uc[sitei]/=2;
        iSy_honeycomb_uc[sitei]/=2;
        Sz_honeycomb_uc[sitei]/=2;

        //PrintDat(Id_honeycomb_uc[sitei]);
        //PrintDat(Sx_honeycomb_uc[sitei]);
        //PrintDat(iSy_honeycomb_uc[sitei]);
        //PrintDat(Sz_honeycomb_uc[sitei]);
    }

    IQCombiner honeycomb_legs_combiner(honeycomb_legs[0],honeycomb_legs[1]);
    honeycomb_legs_combiner.init("honeycomb_uc_legs_to_square",Site,Out);

    IQTensor IdId=Id_honeycomb_uc[0]*Id_honeycomb_uc[1]*dag(honeycomb_legs_combiner)*prime(honeycomb_legs_combiner),
             SxId=Sx_honeycomb_uc[0]*Id_honeycomb_uc[1]*dag(honeycomb_legs_combiner)*prime(honeycomb_legs_combiner),
             IdSx=Id_honeycomb_uc[0]*Sx_honeycomb_uc[1]*dag(honeycomb_legs_combiner)*prime(honeycomb_legs_combiner),
             SxSx=Sx_honeycomb_uc[0]*Sx_honeycomb_uc[1]*dag(honeycomb_legs_combiner)*prime(honeycomb_legs_combiner),
             iSyId=iSy_honeycomb_uc[0]*Id_honeycomb_uc[1]*dag(honeycomb_legs_combiner)*prime(honeycomb_legs_combiner),
             IdiSy=Id_honeycomb_uc[0]*iSy_honeycomb_uc[1]*dag(honeycomb_legs_combiner)*prime(honeycomb_legs_combiner),
             SySy=-iSy_honeycomb_uc[0]*iSy_honeycomb_uc[1]*dag(honeycomb_legs_combiner)*prime(honeycomb_legs_combiner),
             SzId=Sz_honeycomb_uc[0]*Id_honeycomb_uc[1]*dag(honeycomb_legs_combiner)*prime(honeycomb_legs_combiner),
             IdSz=Id_honeycomb_uc[0]*Sz_honeycomb_uc[1]*dag(honeycomb_legs_combiner)*prime(honeycomb_legs_combiner),
             SzSz=Sz_honeycomb_uc[0]*Sz_honeycomb_uc[1]*dag(honeycomb_legs_combiner)*prime(honeycomb_legs_combiner);

    //PrintDat(IdId);
    //PrintDat(SxId);
    //PrintDat(IdSx);
    //PrintDat(SxSx);
    //PrintDat(iSyId);
    //PrintDat(IdiSy);
    //PrintDat(SySy);
    //PrintDat(SzId);
    //PrintDat(IdSz);
    //PrintDat(SzSz);

    //get correlation function
    //std::vector<int> acting_sites_list={3};
    //cout << "SxSx correlator: " << square_double_layer_peps.obtain_correlators(acting_sites_list,std::vector<IQTensor>{SxSx}) << endl;
    //cout << "SySy correlator: " << square_double_layer_peps.obtain_correlators(acting_sites_list,std::vector<IQTensor>{SySy}) << endl;
    //cout << "SzSz correlator: " << square_double_layer_peps.obtain_correlators(acting_sites_list,std::vector<IQTensor>{SzSz}) << endl;

    //std::vector<int> acting_sites_list={3,3+1};
    //cout << "SxSx correlator: " << square_double_layer_peps.obtain_correlators(acting_sites_list,std::vector<IQTensor>{IdSx,SxId}) << endl;
    //cout << "SySy correlator: " << square_double_layer_peps.obtain_correlators(acting_sites_list,std::vector<IQTensor>{IdiSy,-iSyId}) << endl;
    //cout << "SzSz correlator: " << square_double_layer_peps.obtain_correlators(acting_sites_list,std::vector<IQTensor>{IdSz,SzId}) << endl;

    int sitei=2+(Ly-1)/2*Lx;
    //int sitei=2;

    //single bond correlation
    //Sx correlators
    cout << "SxSx correlator for bond a: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei},std::vector<IQTensor>{SxSx}) << endl;
    cout << "SxSx correlator for bond b: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei,sitei+1},std::vector<IQTensor>{IdSx,SxId}) << endl;
    cout << "SxSx correlator for bond c: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei,sitei+Lx},std::vector<IQTensor>{IdSx,SxId}) << endl;
    //Sy correlators
    cout << "SySy correlator for bond a: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei},std::vector<IQTensor>{SySy}) << endl;
    cout << "SySy correlator for bond b: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei,sitei+1},std::vector<IQTensor>{IdiSy,-iSyId}) << endl;
    cout << "SySy correlator for bond c: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei,sitei+Lx},std::vector<IQTensor>{IdiSy,-iSyId}) << endl;
    //Sz correlators
    cout << "SzSz correlator for bond a: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei},std::vector<IQTensor>{SzSz}) << endl;
    cout << "SzSz correlator for bond b: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei,sitei+1},std::vector<IQTensor>{IdSz,SzId}) << endl;
    cout << "SzSz correlator for bond c: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei,sitei+Lx},std::vector<IQTensor>{IdSz,SzId}) << endl;
    

    //get entanglement property
    square_double_layer_peps.move_sigma_lr({Lx/2-1,Lx/2});
    square_double_layer_peps.from_sigma_lr_to_sigma_b();
    square_double_layer_peps.obtain_density_matrix_spectrum();
    cout << "Density Matrix spectrum: " << endl; 
    for (double eigval : square_double_layer_peps.density_mat_spectrum()) cout << eigval << " ";
    cout << endl;
    cout << "Entanglement entropy: " << square_double_layer_peps.entanglement_entropy_vN() << endl;

    return 0;
}
