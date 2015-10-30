
#include "TPO.h"


int main()
{
    IQTPO T=SpinSpin();
    PrintDat(T);
    PrintDat(T.site_tensors(0)*T.bond_tensors(0)*T.site_tensors(1));

    return 0;
}
