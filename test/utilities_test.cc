
//#include "peps.h"
//#include "square_rvb.h"
#include "kagome_rvb.h"
//#include "simple_update.h"

int main()
{
    auto kagome_rvb=kagome_cirac_srvb_peps(4,4);
    //IQTPO SdotS_kagome_cirac=SpinSpin_kagome_cirac({kagome_rvb.phys_legs(0),kagome_rvb.phys_legs(2),kagome_rvb.phys_legs(1)});
    //PrintDat(SdotS_kagome_cirac);
    //PrintDat(SdotS_kagome_cirac.site_tensors(0)*SdotS_kagome_cirac.bond_tensors(0)*SdotS_kagome_cirac.site_tensors(1)*SdotS_kagome_cirac.site_tensors(2));

    IQTPO trotter_gate=trotter_gate_kagome_cirac({kagome_rvb.phys_legs(0),kagome_rvb.phys_legs(2),kagome_rvb.phys_legs(1)},1);
    PrintDat(trotter_gate);
    PrintDat(trotter_gate.site_tensors(0)*trotter_gate.bond_tensors(0)*trotter_gate.site_tensors(1)*trotter_gate.site_tensors(2));
    return 0;
}

