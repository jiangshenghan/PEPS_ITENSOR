
#include "tnetwork_storage.h"
#include "corner_transfer_matrix.h"

int main()
{
    //reading tnetwork from file
    std::stringstream ss;
    ss << "/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_srvb_Lx=4_Ly=4.txt";
    std::string file_name=ss.str();
    Tnetwork_Storage<IQTensor> square_srvb_tnetwork;
    readFromFile(file_name,square_srvb_tnetwork);
    //check tnetwork from file
    Print(square_srvb_tnetwork._Lx);
    Print(square_srvb_tnetwork._Ly);
    PrintDat(square_srvb_tnetwork._tensor_list(0));

    //stores ordered virtual indices 
    const auto &tensor_list=square_srvb_tnetwork._tensor_list;
    int Lx=square_srvb_tnetwork._Lx,
        Ly=square_srvb_tnetwork._Ly;
    std::array<IQIndex,4> tens0_ordered_virt_indices;
    tens0_ordered_virt_indices[0]=commonIndex(tensor_list(0),tensor_list(Lx-1));
    tens0_ordered_virt_indices[1]=commonIndex(tensor_list(0),tensor_list(Lx));
    tens0_ordered_virt_indices[1]=commonIndex(tensor_list(0),tensor_list(1));
    tens0_ordered_virt_indices[1]=commonIndex(tensor_list(0),tensor_list(Lx*(Ly-1)));

    Corner_Transfer_Matrix<IQTensor> square_srvb_ctm({tensor_list(0)},{tens0_ordered_virt_indices});

    IQTPO heisenberg_hamiltonian=SpinSpin();
    PrintDat(heisenberg_hamiltonian);

    square_srvb_ctm.obtain_all_bonds_energy(heisenberg_hamiltonian);
    for (int stepi=0; stepi<500; stepi++)
    {
        square_srvb_ctm.update_env_tensor_one_step();
        square_srvb_ctm.obtain_all_bonds_energy(heisenberg_hamiltonian);
    }

    return 0;
}
