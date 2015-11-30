
//#include "peps.h"
//#include "square_rvb.h"
#include "kagome_rvb.h"

//using namespace square_psg;

int main()
{
    //random_init_square_rvb_peps(square_peps);
    //Construct square lattice sRVB
    int Lx=4, Ly=4;
    IQPEPS kagome_peps=kagome_srvb_cirac_peps(Lx,Ly);

    //Square_Lattice_Torus square_lattice({Lx,Ly});
    //IQPEPS_IndexSet_SpinHalf index_set(6,square_lattice);
    //IQPEPS square_peps(square_lattice,index_set);
    //random_init_square_rvb_peps(square_peps);

    //mu_12=1;
    //IQPEPS square_srvb_zero=square_srvb_peps(Lx,Ly);
    //Tnetwork_Storage<IQTensor> square_srvb_zero_storage=peps_to_tnetwork_storage(square_srvb_zero);
    //writeToFile("/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_srvb_zero_flux.txt",square_srvb_zero_storage);

    //mu_12=-1;
    //IQPEPS square_srvb_pi=square_srvb_peps(Lx,Ly);
    //Tnetwork_Storage<IQTensor> square_srvb_pi_storage=peps_to_tnetwork_storage(square_srvb_pi);
    //writeToFile("/home/jiangsb/code/peps_itensor/result/tnetwork_storage/square_srvb_pi_flux.txt",square_srvb_pi_storage);

    //cout << "Finish writing as tnetwork_storage!" << endl;

    return 0;
}
