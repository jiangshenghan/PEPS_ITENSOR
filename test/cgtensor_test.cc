
#include "cgtensor.h"

int main()
{
    //std::vector<Spin> spins{Spin(1,Out), Spin(1,Out)};
    //std::vector<Spin> spins{Spin(1,Out),Spin(0,Out),Spin(1,Out),Spin(1,Out),Spin(1,Out)};
    std::vector<Spin> spins{Spin(1,Out),Spin(1,Out),Spin(1,In),Spin(1,In)};
    CGTensors spin_singlets(spins);

    if (!spin_singlets.valid()) cout << "Inconsistent spins!" << endl;
    for (const auto &tensor : spin_singlets.K())
    {
        PrintDat(tensor);
    }

    return 0;
}
