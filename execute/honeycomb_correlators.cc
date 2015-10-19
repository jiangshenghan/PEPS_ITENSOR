
#include "transfer_to_square.h"
#include "square_cylinder.h"

int main()
{
    //system size
    int Ly=6, Lx=16*Ly;

    double A1=0.9, A2=1-A1;

    //Output file name for double_layer_peps
    std::stringstream ss;
    ss << "/home/jiangsb/code/peps_itensor/result/featureless_honeycomb_peps_Ly=" << Ly << "/Lx=" << Lx << "_A1=" << A1 << "_A2=" << A2 << "_double_layer_peps.txt";
    std::string file_name=ss.str();

    cout << "\n========================================\n" << endl;
    cout << "System Size: " << Lx << "x" << Ly << endl;
    cout << "Wavefunction Params: " << "A1=" << A1 << ", A2=" << A2 << endl;
    cout << "Input file: " << endl << file_name << endl;
    cout << "========================================\n" << endl;

    Square_Lattice_Cylinder square_cylinder(std::array<int,2>{Lx,Ly});

    //reading square_cylinder_double_peps
    Cylinder_Square_Double_Layer_PEPSt<IQTensor> square_cylinder_double_peps(square_cylinder);
    readFromFile(file_name,square_cylinder_double_peps);
    
    cout << "Successfully reading from file!" << endl;

    //cout << "\n========================================\n" << endl;
    //cout << "Data for double layer PEPS reading from file:" << endl;
    //cout << "col_left=" << square_cylinder_double_peps.col_lr(0) << ", col_right=" << square_cylinder_double_peps.col_lr(1) << endl;
    //PrintDat(square_cylinder_double_peps.sigma_lr(0));
    //PrintDat(square_cylinder_double_peps.sigma_lr(1));
    //cout << "\n========================================\n" << endl;

    //calculate correlators
    //std::array<Coordinate,2> acting_sites_coord;
    //acting_sites_coord[0]=Coordinate{1,0,0};
    //acting_sites_coord[1]=Coordinate{2,0,0};

    //cout << honeycomb_cylinder_SzSz_correlator(acting_sites_coord,square_cylinder_double_peps) << endl;


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
    //cout << "SxSx correlator: " << square_cylinder_double_peps.obtain_correlators(acting_sites_list,std::vector<IQTensor>{SxSx}) << endl;
    //cout << "SySy correlator: " << square_cylinder_double_peps.obtain_correlators(acting_sites_list,std::vector<IQTensor>{SySy}) << endl;
    //cout << "SzSz correlator: " << square_cylinder_double_peps.obtain_correlators(acting_sites_list,std::vector<IQTensor>{SzSz}) << endl;

    //std::vector<int> acting_sites_list={3,3+1};
    //cout << "SxSx correlator: " << square_cylinder_double_peps.obtain_correlators(acting_sites_list,std::vector<IQTensor>{IdSx,SxId}) << endl;
    //cout << "SySy correlator: " << square_cylinder_double_peps.obtain_correlators(acting_sites_list,std::vector<IQTensor>{IdiSy,-iSyId}) << endl;
    //cout << "SzSz correlator: " << square_cylinder_double_peps.obtain_correlators(acting_sites_list,std::vector<IQTensor>{IdSz,SzId}) << endl;

    cout << "SzSz correlator for bond a (weak bond): " << square_cylinder_double_peps.obtain_correlators(std::vector<int>{3},std::vector<IQTensor>{SzSz}) << endl;
    cout << "SzSz correlator for bond b (strong bond): " << square_cylinder_double_peps.obtain_correlators(std::vector<int>{3,4},std::vector<IQTensor>{IdSz,SzId}) << endl;
    cout << "SzSz correlator for bond c (weak bond): " << square_cylinder_double_peps.obtain_correlators(std::vector<int>{3,3+Lx},std::vector<IQTensor>{IdSz,SzId}) << endl;

    //get entanglement property
    square_cylinder_double_peps.move_sigma_lr({Lx/2-1,Lx/2});
    square_cylinder_double_peps.from_sigma_lr_to_sigma_b();
    square_cylinder_double_peps.obtain_density_matrix_spectrum();
    cout << "Density Matrix spectrum: " << endl; 
    for (double eigval : square_cylinder_double_peps.density_mat_spectrum()) cout << eigval << " ";
    cout << endl;
    cout << "Entanglement entropy: " << square_cylinder_double_peps.entanglement_entropy_vN() << endl;

    return 0;
}
