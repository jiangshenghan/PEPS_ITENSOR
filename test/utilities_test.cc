
#include "singlet_tensor_basis.h"


int main()
{
    auto sigma_left=readFromFile<IQTensor>("/home/jiangsb/code/peps_itensor/result/test/sigma_left.txt");
    auto sigma_right=readFromFile<IQTensor>("/home/jiangsb/code/peps_itensor/result/test/sigma_right.txt");
    PrintDat(sigma_left);
    PrintDat(sigma_right);

    return 0;
}
