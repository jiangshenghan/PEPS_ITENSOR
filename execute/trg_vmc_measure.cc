
#include "kagome_rvb.h"
#include "tensor_vmc.h"

int main(int argc, char *argv[])
{
    //to run this file, type e.g.:
    //mpirun -np 10 ./trg_vmc_measure.cc file_path/file_name operator(Heisenberg/SzSz) maxm(an int)

    //initialize mpi
    MPI_Init(&argc,&argv);
    int mpi_id, mpi_procnum;
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_id);
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_procnum);

    // construct tnetwork_storage from file
    Tnetwork_Storage<IQTensor> tnetwork_storage;
    readFromFile(argv[1],tnetwork_storage);
    cout << "Reading files successfully!" << endl;

    //kagome lattice
    int Lx=tnetwork_storage._Lx,
        Ly=tnetwork_storage._Ly;
    if (tnetwork_storage._tnetwork_type!=8 || tnetwork_storage._boundary_condition!=1) 
    {
        Print(tnetwork_storage._tnetwork_type);
        Print(tnetwork_storage._boundary_condition);
        cout << "Incorrect tnetwork_storage" << endl;
        exit(1);
    }
    Kagome_Normal_Lattice_Torus kagome_lattice({Lx,Ly});

    if (mpi_id==0) 
    { 
        Print(mpi_procnum); 
        Print(kagome_lattice.name());
        Print(Lx);
        Print(Ly);
        Print(argv[1]);
    }

    std::string ope=argv[2];
    int maxm=std::atoi(argv[3]);
    Args measure_args={"Operator",ope,"Maxm",maxm,"InitSpins","antiferro","ThermalSteps",20,"MeasureSteps",500,"SpinFlipPrint",true};

    std::vector<IQTensor> combined_tensors;
    for (int y=0; y<Ly; y++)
        for (int x=0; x<Lx; x++)
            for (int subi=0; subi<tnetwork_storage._n_subl; subi++) 
                combined_tensors.push_back(tnetwork_storage._tensor_list(tnetwork_storage._coor_to_siteind(x,y)(subi)));
    tensor_vmc_parallel<IQTensor>(combined_tensors,kagome_lattice,measure_args);

    MPI_Finalize();

    return 0;
}

