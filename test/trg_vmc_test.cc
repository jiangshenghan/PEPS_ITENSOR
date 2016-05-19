
#include "kagome_rvb.h"
#include "tensor_vmc.h"

int main()
{
    kagome_psg::mu_12=-1;
    kagome_psg::mu_c6=-1;
    
    Print(kagome_psg::mu_12);
    Print(kagome_psg::mu_c6);

    int Lx=8, Ly=8, D=3;
    //Kagome_Normal_Lattice_Torus kagome_lattice({Lx,Ly});
    //IQPEPS_IndexSet_SpinHalf index_set(3,kagome_lattice);
    //IQPEPS kagome_rvb(kagome_lattice,index_set);
    //random_init_kagome_rvb_normal_peps(kagome_rvb);
    IQPEPS kagome_rvb=kagome_normal_srvb_peps(Lx,Ly);

    tensor_vmc<IQTensor>(kagome_rvb,{"Operator","Heisenberg","Maxm",20,"InitSpins","antiferro","ThermalSteps",20,"MeasureSteps",200});

    return 0;
}
