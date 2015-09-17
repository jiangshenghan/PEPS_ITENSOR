
#include "singlet_tensor_basis.h"


int main()
{
    IQIndex i1=Spin_leg({0,1},"leg1",Out),
            i2=Spin_leg({1,1},"leg2",Out),
            i3=Spin_leg({1,1},"leg3",Out),
            i4=Spin_leg({1,1},"leg4",Out),
            i5=Spin_leg({1,1},"leg5",Out);

    Singlet_Tensor_Basis tensor_basis(std::vector<IQIndex>{i1,i2,i3,i4,i5});

    int i=0;
    for (const auto &tensor : tensor_basis.tensors())
    {
        cout << tensor_basis.spin_configs(i) << endl;
        PrintDat(tensor);
        i++;
    }
    cout << "Number of singlet basis: " << tensor_basis.tensors().size() << endl << endl;


    return 0;
}
