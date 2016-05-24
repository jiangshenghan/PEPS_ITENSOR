
#include "spin_sym_decomp.h"

int main()
{
    IQIndex leg_a=Spin_leg({1,2,2},"leg_a",In),
            leg_b=Spin_leg({2,2,2,1},"leg_b",Out);
    Singlet_Tensor_Basis tensor_basis(std::vector<IQIndex>{leg_a,leg_b});
    std::vector<double> params;
    for (int basei=0; basei<tensor_basis.dim(); basei++) params.push_back(5.*rand_gen());
    IQTensor T=singlet_tensor_from_basis_params(tensor_basis,params);
    IQTensor U(leg_a),D,V;
    spin_sym_svdRank2(T,leg_a,leg_b,U,D,V);
    PrintDat(T);
    Print(norm(T-U*D*V));
    return 0;
}

