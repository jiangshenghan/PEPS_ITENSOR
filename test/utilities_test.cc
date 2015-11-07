
#include "peps.h"

int main()
{
    //Construct square lattice sRVB
    int Lx=16, Ly=16;
    Square_Lattice_Torus square_torus({Lx,Ly});
    IQPEPS_IndexSet_SpinHalf indexset(3,square_torus);
    IQPEPS square_srvb(square_torus,indexset);

    //input sqaure srvb ansatz
    Singlet_Tensor_Basis site_tensor_basis(square_srvb.site_tensors(0).indices()),
                         bond_tensor_basis(square_srvb.bond_tensors(0).indices());
    //check tensor basis
    //cout << site_tensor_basis.dim() << " site tensor basis:" << endl;
    //PrintDat(site_tensor_basis);
    //cout << "bond tensor basis:" << endl;
    //PrintDat(bond_tensor_basis);
    //return 0;

    std::vector<double> site_tensor_params(site_tensor_basis.dim()),
                        bond_tensor_params={1,sqrt(2)};
    site_tensor_params[0]=sqrt(2);
    site_tensor_params[1]=sqrt(2);
    site_tensor_params[2]=sqrt(2);
    site_tensor_params[5]=sqrt(2);
    IQTensor site_tensor=singlet_tensor_from_basis_params(site_tensor_basis,site_tensor_params),
             bond_tensor=singlet_tensor_from_basis_params(bond_tensor_basis,bond_tensor_params);
    square_srvb.generate_site_tensors({site_tensor});
    square_srvb.generate_bond_tensors({bond_tensor,bond_tensor});
    //for (const auto &site_tensor : square_srvb.site_tensors()) PrintDat(site_tensor);
    //for (const auto &bond_tensor : square_srvb.bond_tensors()) PrintDat(bond_tensor);

    
    //store the peps ansatz for measurement
    Tnetwork_Storage<IQTensor> square_srvb_storage;

    square_srvb_storage._tnetwork_type=1;
    square_srvb_storage._boundary_condition=1;
    square_srvb_storage._n_subl=1;
    square_srvb_storage._Lx=Lx;
    square_srvb_storage._Ly=Ly;

    //combine bond tensors to site tensors and then store in _tensor_list
    square_srvb_storage._tensor_list.set_size(square_torus.n_sites_total());
    auto tensor_list=square_srvb.combined_site_tensors();
    //for (const auto &tensor : tensor_list) PrintDat(tensor);
    for (int sitei=0; sitei<square_torus.n_sites_total(); sitei++) square_srvb_storage._tensor_list(sitei)=tensor_list[sitei];

    square_srvb_storage._coor_to_siteind.set_size(Lx,Ly);
    for (int x=0; x<Lx; x++)
    {
        for (int y=0; y<Ly; y++)
        {
            square_srvb_storage._coor_to_siteind(x,y).set_size(square_srvb_storage._n_subl);
            for (int subli=0; subli<square_srvb_storage._n_subl; subli++) square_srvb_storage._coor_to_siteind(x,y)(subli)=square_torus.site_coord_to_list(x,y,subli);
        }
    }
    //cout << square_srvb_storage._coor_to_siteind << endl;

    //Output file name for square_srvb_storage
    std::stringstream ss;
    ss << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_srvb_Lx=" << Lx << "_Ly=" << Ly << ".txt";
    std::string file_name=ss.str();
    writeToFile(file_name,square_srvb_storage);

    //Tnetwork_Storage<IQTensor> square_srvb_from_file;
    //readFromFile(file_name,square_srvb_from_file);
    //cout << square_srvb_from_file._Lx << " " << square_srvb_from_file._Ly << endl << square_srvb_from_file._tensor_list << endl;

    return 0;
}
