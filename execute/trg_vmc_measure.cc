
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

    Args measure_args={"Operator","SzSz","Maxm",20,"InitSpins","antiferro","ThermalSteps",10,"MeasureSteps",300,"SpinFlipPrint",true};

    tensor_vmc_parallel<IQTensor>(kagome_rvb,measure_args);

    //using ITensor instead of IQTensor
    //kagome_rvb.obtain_combined_site_tensors();
    //std::vector<ITensor> combined_tensors;
    //for (const auto &tensor : kagome_rvb.combined_tensors()) combined_tensors.push_back(tensor.toITensor());
    //tensor_vmc_parallel<ITensor>(combined_tensors,kagome_lattice,measure_args);

    MPI_Finalize();

    return 0;
}

