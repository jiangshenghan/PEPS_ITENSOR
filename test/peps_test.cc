
//#include "double_layer_peps.h"
//#include "square_rvb.h"
#include "kagome_rvb.h"

int main()
{
    kagome_psg::mu_12=1;
    kagome_psg::mu_c6=1;

    auto kagome_srvb=kagome_cirac_srvb_peps(4,4);

    Kagome_Cirac_Lattice_Torus kagome_cirac_lattice({4,4});
    IQPEPS_IndexSet_SpinHalf index_set(3,kagome_cirac_lattice);
    //IQPEPS kagome_rvb(kagome_cirac_lattice,index_set);
    IQPEPS kagome_rvb(kagome_cirac_lattice,kagome_srvb.indexset());
    random_init_kagome_rvb_cirac_peps(kagome_rvb,{sqrt(2),sqrt(6)});

    for (int sitei=0; sitei<kagome_srvb.n_sites_total(); sitei++)
    {
        PrintDat(kagome_rvb.site_tensors(sitei));
        PrintDat(kagome_srvb.site_tensors(sitei));
        PrintDat((kagome_srvb.site_tensors(sitei)-kagome_rvb.site_tensors(sitei)).norm());
    }
    for (int bondi=0; bondi<kagome_srvb.n_bonds_total(); bondi++)
    {
        //PrintDat(kagome_rvb.bond_tensors(bondi));
        //PrintDat(kagome_srvb.bond_tensors(bondi));
        Print((kagome_srvb.bond_tensors(bondi)-kagome_rvb.bond_tensors(bondi)).norm());
    }

    Tnetwork_Storage<IQTensor> kagome_rvb_storage=peps_to_tnetwork_storage(kagome_rvb);

    //Print(kagome_rvb_storage._tnetwork_type);
    //Print(kagome_rvb_storage._Lx);
    //Print(kagome_rvb_storage._Ly);
    //Print(kagome_rvb_storage._boundary_condition);
    //Print(kagome_rvb_storage._n_subl);
    //Print(kagome_rvb_storage._coor_to_siteind);
    //Print(kagome_rvb_storage._tensor_list);

    std::stringstream ss;
    ss << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/kagome_srvb_Lx=4_Ly=4_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6;
    std::string file_name=ss.str();
    writeToFile(file_name,kagome_rvb_storage);

    //PrintDat(kagome_rvb.site_tensors(0));
    //PrintDat(kagome_rvb.site_tensors(1));
    //PrintDat(kagome_rvb.site_tensors(2));

    //PrintDat(kagome_rvb.bond_tensors(0));
    //PrintDat(kagome_rvb.bond_tensors(1));

    //Square_Lattice_Torus square_lattice{std::array<int,2>{2,2}};
    //Square_Lattice_Cylinder square_lattice{{3,3}};
    //Square_Lattice_Open square_lattice{{3,3}};
    
    //square_lattice.print_lattice_inf();

    //PEPS_IndexSet index_set(2,2,square_lattice);
    //PEPS peps_test(square_lattice,index_set);

    //IQPEPS_IndexSet_SpinHalf index_set(6,square_lattice);
    //IQPEPS_IndexSet_SpinHalf index_set(6,kagome_cirac);
    //IQPEPS_IndexSet_SpinHalf index_set(17,std::vector<int>{1,3,0,0,2},square_lattice);

    //cout << "\n========================================\n" << endl;
    //cout << "Original PEPS index set: " << endl;
    //cout << index_set.name() << endl << endl;
    //cout << "Physical Legs d=" << index_set.d() <<  endl;
    //for (const auto &leg : index_set.phys_legs())
    //{
    //    cout << leg << endl;
    //}
    //cout << "Virtual Legs D=" << index_set.D() << endl;
    //for (const auto &leg : index_set.virt_legs())
    //{
    //    cout << leg << endl;
    //}
    //cout << "\n========================================\n" << endl;

    //writeToFile("/home/jiangsb/code/peps_itensor/result/test/iotest.txt",index_set);
    //IQPEPS_IndexSet_SpinHalf index_set_from_file;
    //readFromFile("/home/jiangsb/code/peps_itensor/result/test/iotest.txt",index_set_from_file);

    //cout << "\n========================================\n" << endl;
    //cout << "PEPS index set from iotest.txt: " << endl;
    //cout << index_set_from_file.name() << endl << endl;
    //cout << "Physical Legs d=" << index_set.d() <<  endl;
    //for (const auto &leg : index_set_from_file.phys_legs())
    //{
    //    cout << leg << endl;
    //}
    //cout << "Virtual Legs D=" << index_set.D() << endl;
    //for (const auto &leg : index_set_from_file.virt_legs())
    //{
    //    cout << leg << endl;
    //}
    //cout << "\n========================================\n" << endl;

    //IQPEPS peps_test(square_lattice,index_set);
    //IQPEPS peps_test(kagome_cirac,index_set);

    //cout << "\n========================================\n" << endl;
    //cout << "Original PEPS:" << endl;
    //cout << peps_test.name() << endl << endl;

    //cout << "Site Tensors: " << endl;
    //for (const auto &tensor : peps_test.site_tensors())
    //{
    //    cout << tensor << endl;
    //}

    //cout << "Bond Tensors: " << endl;
    //for (const auto &tensor : peps_test.bond_tensors())
    //{
    //    cout << tensor << endl;
    //}

    //cout << "Boundary Tensors: " << endl;
    //for (const auto &tensor : peps_test.boundary_tensors())
    //{
    //    cout << tensor << endl;
    //}
    //cout << "\n========================================\n" << endl;


    //IQPEPS peps_from_file(square_lattice);
    //writeToFile("/home/jiangsb/code/peps_itensor/result/test/iotest.txt",peps_test);
    //readFromFile("/home/jiangsb/code/peps_itensor/result/test/iotest.txt",peps_from_file);

    //cout << "\n========================================\n" << endl;
    //cout << "PEPS reading from file:" << endl;
    //cout << peps_from_file.name() << endl << endl;

    //cout << "Site Tensors: " << endl;
    //for (const auto &tensor : peps_from_file.site_tensors())
    //{
    //    cout << tensor << endl;
    //}

    //cout << "Bond Tensors: " << endl;
    //for (const auto &tensor : peps_from_file.bond_tensors())
    //{
    //    cout << tensor << endl;
    //}

    //cout << "Boundary Tensors: " << endl;
    //for (const auto &tensor : peps_from_file.boundary_tensors())
    //{
    //    cout << tensor << endl;
    //}
    //cout << "\n========================================\n" << endl;

    //Double_Layer_IQPEPS layered_peps_test(peps_test);

    //cout << "\n========================================\n" << endl;
    //cout << "Original double layer PEPS:" << endl;
    //cout << "Virtual Leg Combiners:" << endl;
    //int sitei=0;
    //for (const auto &site_combiners : layered_peps_test.virt_leg_combiners())
    //{
    //    cout << "Cominers for virtual legs of site " << sitei << endl;
    //    for (const auto &combiner : site_combiners)
    //    {
    //        cout << combiner;
    //    }
    //    cout << endl;
    //    sitei++;
    //}

    //cout << "Layered Tensors:" << endl;
    //for (const auto &tensor : layered_peps_test.double_layer_tensors())
    //{
    //    cout << tensor << endl;
    //}
    //cout << "\n========================================\n" << endl;

    //writeToFile("/home/jiangsb/code/peps_itensor/result/test/iotest.txt",layered_peps_test);
    //Double_Layer_IQPEPS layered_peps_from_file(square_lattice);
    //readFromFile("/home/jiangsb/code/peps_itensor/result/test/iotest.txt",layered_peps_from_file);

    //cout << "\n========================================\n" << endl;
    //cout << "Double layer PEPS reading from file:" << endl;
    //cout << "Virtual Leg Combiners:" << endl;
    //sitei=0;
    //for (const auto &site_combiners : layered_peps_from_file.virt_leg_combiners())
    //{
    //    cout << "Cominers for virtual legs of site " << sitei << endl;
    //    for (const auto &combiner : site_combiners)
    //    {
    //        cout << combiner;
    //    }
    //    cout << endl;
    //    sitei++;
    //}

    //cout << "Layered Tensors:" << endl;
    //for (const auto &tensor : layered_peps_from_file.double_layer_tensors())
    //{
    //    cout << tensor << endl;
    //}
    //cout << "\n========================================\n" << endl;

    return 0;
}
