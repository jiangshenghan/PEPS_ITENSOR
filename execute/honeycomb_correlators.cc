
#include "transfer_to_square.h"
#include "square_double_layer_peps.h"

int main()
{
    //system size
    int Ly=6, Lx=30*Ly;

    double A1=0.9, A2=1-A1;

    //Output file name for double_layer_peps
    std::stringstream ss;
    //file name for cylinder
    ss << "/home/jiangsb/code/peps_itensor/result/featureless_honeycomb_peps_Ly=" << Ly << "/Lx=" << Lx << "_A1=" << A1 << "_A2=" << A2 << "_double_layer_peps.txt";
    //file name for ribbon
    //ss << "/home/jiangsb/code/peps_itensor/result/featureless_honeycomb_peps_ribbon" << "/Ly=" << Ly << "_Lx=" << Lx << "_A1=" << A1 << "_A2=" << A2 << "_double_layer_peps.txt";

    std::string file_name=ss.str();

    //cylinder geometry peps
    Square_Lattice_Cylinder square_cylinder(std::array<int,2>{Lx,Ly});
    Square_Double_Layer_PEPSt<IQTensor> square_double_layer_peps(square_cylinder);

    //ribbon geometry peps
    //Square_Lattice_Open square_ribbon(std::array<int,2>{Lx,Ly});
    //Square_Double_Layer_PEPSt<IQTensor> square_double_layer_peps(square_ribbon);

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

    //test transfer matrix
    //int col_no=2, transfer_site_no=col_no;
    //IQTensor transfer_mat;
    //for (int coli=0; coli<Ly; coli++) 
    //{
    //    if (coli==0) transfer_mat=square_double_layer_peps.double_layer_tensors(transfer_site_no);
    //    if (coli>0) transfer_mat*=square_double_layer_peps.double_layer_tensors(transfer_site_no);

    //    //IQIndex left_indice=commonIndex(square_double_layer_peps.double_layer_tensors(transfer_site_no),dag(square_double_layer_peps.double_layer_tensors(transfer_site_no-1))),
    //    //        right_indice=commonIndex(square_double_layer_peps.double_layer_tensors(transfer_site_no),dag(square_double_layer_peps.double_layer_tensors(transfer_site_no+1)));
    //    //transfer_mat.replaceIndex(left_indice,dag(right_indice).prime());
    //    transfer_site_no+=Lx;
    //}
    ////PrintDat(transfer_mat);

    //square_double_layer_peps.move_sigma_lr(std::array<int,2>{col_no-1,col_no+1});
    //square_double_layer_peps.decombine_sigma_lr();
    //PrintDat(square_double_layer_peps.sigma_lr(0).norm());
    //PrintDat((square_double_layer_peps.sigma_lr(0)*transfer_mat).norm());

    //std::vector<int> sorted_indices(eigvals_transfer.size());
    //for (int i=0; i<sorted_indices.size(); i++) sorted_indices[i]=i;
    //std::sort(sorted_indices.begin(),sorted_indices.end(),
    //        [&eigvals_transfer](int ind1, int ind2){ return eigvals[ind1]<eigvals[ind2]; });

    //return 0;

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

    //single bond correlation
    //int sitei=2+(Ly-1)/2*Lx;
    //Sx correlators
    //cout << "SxSx correlator for bond a: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei},std::vector<IQTensor>{SxSx}) << endl;
    //cout << "SxSx correlator for bond b: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei,sitei+1},std::vector<IQTensor>{IdSx,SxId}) << endl;
    //cout << "SxSx correlator for bond c: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei,sitei+Lx},std::vector<IQTensor>{IdSx,SxId}) << endl;
    ////Sy correlators
    //cout << "SySy correlator for bond a: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei},std::vector<IQTensor>{SySy}) << endl;
    //cout << "SySy correlator for bond b: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei,sitei+1},std::vector<IQTensor>{IdiSy,-iSyId}) << endl;
    //cout << "SySy correlator for bond c: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei,sitei+Lx},std::vector<IQTensor>{IdiSy,-iSyId}) << endl;
    ////Sz correlators
    //cout << "SzSz correlator for bond a: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei},std::vector<IQTensor>{SzSz}) << endl;
    //cout << "SzSz correlator for bond b: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei,sitei+1},std::vector<IQTensor>{IdSz,SzId}) << endl;
    //cout << "SzSz correlator for bond c: " << square_double_layer_peps.obtain_correlators(std::vector<int>{sitei,sitei+Lx},std::vector<IQTensor>{IdSz,SzId}) << endl;
    
    int origin_site=2, distance=15;

    //obtain long-range spin-spin correlator
    cout << "benchmark with previous result about bond SzSz: " << square_double_layer_peps.obtain_correlators(std::vector<int>{origin_site},std::vector<IQTensor>{SzSz}) << endl << endl;
    for (int sitei=origin_site+1; sitei<origin_site+distance; sitei++)
    {
        cout << "correlator between col " << origin_site%Lx << " and col " << sitei%Lx << endl;
        //cout << "SxSx00: " << square_double_layer_peps.obtain_correlators(std::vector<int>{origin_site,sitei},std::vector<IQTensor>{SxId,SxId}) << endl;
        //cout << "SxSx01: " << square_double_layer_peps.obtain_correlators(std::vector<int>{origin_site,sitei},std::vector<IQTensor>{SxId,IdSx}) << endl;
        //cout << "SySy00: " << square_double_layer_peps.obtain_correlators(std::vector<int>{origin_site,sitei},std::vector<IQTensor>{iSyId,-iSyId}) << endl;
        //cout << "SySy01: " << square_double_layer_peps.obtain_correlators(std::vector<int>{origin_site,sitei},std::vector<IQTensor>{iSyId,-IdiSy}) << endl;
        cout << "SzSz00: " << square_double_layer_peps.obtain_correlators({origin_site,sitei},{SzId,SzId}) << endl;
        cout << "SzSz01: " << square_double_layer_peps.obtain_correlators({origin_site,sitei},{SzId,IdSz}) << endl;
        cout << endl;
    }

    //obtain long-range bond-bond correlator
    //double bond_a_amplitude=square_double_layer_peps.obtain_correlators({origin_site},{SzSz}),
    //       bond_b_amplitude=square_double_layer_peps.obtain_correlators({origin_site,origin_site+1},{IdSz,SzId}),
    //       bond_c_amplitude=square_double_layer_peps.obtain_correlators({origin_site,origin_site+Lx},{IdSz,SzId});
    //cout << "bond a amplitude: " << bond_a_amplitude << endl;
    //cout << "bond b amplitude: " << bond_b_amplitude << endl;
    //cout << "bond c amplitude: " << bond_c_amplitude << endl;

    //for (int sitei=origin_site; sitei<origin_site+distance; sitei++)
    //{
    //    cout << "bond-bond correlator between col " << origin_site%Lx << " and col " << sitei%Lx << endl;

    //    double bond_ba_correlator=square_double_layer_peps.obtain_correlators({origin_site-1,origin_site,sitei},{IdSz,SzId,SzSz}),
    //           bond_bb_correlator=square_double_layer_peps.obtain_correlators({origin_site-1,origin_site,sitei,sitei+1},{IdSz,SzId,IdSz,SzId}),
    //           bond_bc_correlator=square_double_layer_peps.obtain_correlators({origin_site-1,origin_site,sitei,sitei+Lx},{IdSz,SzId,IdSz,SzId});

    //    cout << "bond_ba correlator: " << bond_ba_correlator-bond_b_amplitude*bond_a_amplitude << endl;
    //    cout << "bond_bb correlator: " << bond_bb_correlator-bond_b_amplitude*bond_b_amplitude << endl;
    //    cout << "bond_bc correlator: " << bond_bc_correlator-bond_b_amplitude*bond_c_amplitude << endl;
    //}

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
