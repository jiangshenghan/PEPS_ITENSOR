
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

    int origin_site=2, distance=15;
    //correlators for S.S bonds
    //we choose bonds inside one uc (bond a) for simplicity
    //IQTPO SdotS_honeycomb=SpinSpin();
    //IQCombiner phys_legs_combiner(SdotS_honeycomb.phys_legs(0),SdotS_honeycomb.phys_legs(1));
    //phys_legs_combiner.init("honeycomb_uc_legs_to_square", Site,Out);
    //IQTensor SdotS_square_site=SdotS_honeycomb.site_tensors(0)*SdotS_honeycomb.bond_tensors(0)*SdotS_honeycomb.site_tensors(1)*dag(phys_legs_combiner)*prime(phys_legs_combiner);
    //double SdotS_bonda_expect_val=square_double_layer_peps.obtain_correlators({origin_site},{SdotS_square_site});
    //cout << "correlator for single bond a: " << SdotS_bonda_expect_val << endl;
    //for (int sitei=origin_site+1; sitei<origin_site+distance; sitei++)
    //{
    //    double bond_aa_correlator=square_double_layer_peps.obtain_correlators({origin_site,sitei},{SdotS_square_site,SdotS_square_site});
    //    cout << "correlator between bond a in col " << origin_site%Lx << " and col " << sitei%Lx << " is: " 
    //         << bond_aa_correlator-SdotS_bonda_expect_val*SdotS_bonda_expect_val << endl;
    //}


    //correlators for SxS bonds
    //we choose bonds inside on uc (bond a) for simpicity
    //we only measure (SxS)_y component, others are the same due to spin rotation symmetry
    auto SxS_honeycomb=vectorSpinChirality();
    IQCombiner phys_legs_combiner(SxS_honeycomb[0].phys_legs(0),SxS_honeycomb[0].phys_legs(1));
    phys_legs_combiner.init("honeycomb_uc_legs_to_square",Site,Out);
    std::array<IQTensor,3> SxS_square_site;
    for (int i=0; i<3; i++)
    {
        SxS_square_site[i]=SxS_honeycomb[i].site_tensors(0)*SxS_honeycomb[i].site_tensors(1)*dag(phys_legs_combiner)*prime(phys_legs_combiner);
        SxS_square_site[i].clean();
    }
    //TODO: SxS_square_site[2] do not vanish??
    //cout << "SxS correlator for single bond a: " << endl 
    //    << square_double_layer_peps.obtain_correlators({origin_site},{SxS_square_site[0]}) << endl
    //    << square_double_layer_peps.obtain_correlators({origin_site},{SxS_square_site[1]}) << endl
    //    << square_double_layer_peps.obtain_correlators({origin_site},{SxS_square_site[2]}) << endl;

    for (int sitei=origin_site+1; sitei<origin_site+distance; sitei++)
    {
        double bond_aa_correlator=square_double_layer_peps.obtain_correlators({origin_site,sitei},{SxS_square_site[1],SxS_square_site[1]});
        cout << "correlator between bonda in col " << origin_site%Lx << " and col " << sitei%Lx << " is: " << bond_aa_correlator << endl;
    }

    return 0;
}
