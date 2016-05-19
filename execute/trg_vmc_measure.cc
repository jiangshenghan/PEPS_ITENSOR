
#include "kagome_rvb.h"
#include "tensor_vmc.h"

int main(int argc, char *argv[])
{
    //initialize mpi
    MPI_Init(&argc,&argv);
    int mpi_id, mpi_procnum;
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_id);
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_procnum);

    //to construct peps from file, we need to specify lattice first
    int Lx=8, Ly=8;
    Kagome_Normal_Lattice_Torus kagome_lattice({Lx,Ly});
    IQPEPS kagome_rvb(kagome_lattice);

    if (mpi_id==0) 
    { 
        Print(mpi_procnum); 
        Print(kagome_lattice.name());
        Print(Lx);
        Print(Ly);
        Print(argv[1]);
    }

    readFromFile(argv[1],kagome_rvb);
    cout << "Reading files successfully!" << endl;

    tensor_vmc_parallel<IQTensor>(kagome_rvb,{"Operator","SzSz","Maxm",8,"InitSpins","antiferro","ThermalSteps",10,"MeasureSteps",3000});

    MPI_Finalize();

    return 0;
}

