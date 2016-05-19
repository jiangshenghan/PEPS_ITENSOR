
#include "kagome_rvb.h"
#include "tensor_vmc.h"

int main(int argc, char *argv[])
{
    //initialize mpi
    MPI_Init(&argc,&argv);
    int mpi_id, mpi_procnum;
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_id);
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_procnum);

    kagome_psg::mu_12=-1;
    kagome_psg::mu_c6=1;
    
    if (mpi_id==0)
    {
        Print(kagome_psg::mu_12);
        Print(kagome_psg::mu_c6);
        Print(mpi_procnum);
    }

    int Lx=8, Ly=8, D=6;
    Kagome_Normal_Lattice_Torus kagome_lattice({Lx,Ly});
    IQPEPS_IndexSet_SpinHalf index_set(D,kagome_lattice);
    //IQPEPS kagome_rvb(kagome_lattice,index_set);
    //random_init_kagome_rvb_normal_peps(kagome_rvb);
    //IQPEPS kagome_rvb=kagome_normal_srvb_peps(Lx,Ly);

    //construct PEPS from file
    std::stringstream ss;
    ss << "/home/jiangsb/tn_ying/tensor_network/result/peps_storage/kagome_simple_update/kagome_rvb_normal_D=" << D << "_Lx=" << Lx << "_Ly=" << Ly << "_mu12=" << kagome_psg::mu_12 << "_muc6=" << kagome_psg::mu_c6 << "_step=1e-2";
    std::string file_name=ss.str();
    IQPEPS kagome_rvb(kagome_lattice);
    readFromFile(file_name,kagome_rvb);
    cout << "Reading files successfully!" << endl;

    tensor_vmc_parallel<IQTensor>(kagome_rvb,{"Operator","SzSz","Maxm",15,"InitSpins","antiferro","ThermalSteps",5,"MeasureSteps",100});

    MPI_Finalize();

    return 0;
}
