
#include "kagome_rvb.h"
#include "tensor_update.h"

const dP imps_normalization_factor=1.E2;

int main()
{
    IQIndex leg_a=Spin_leg({1,1,1},"leg_a",In),
            leg_b=Spin_leg({2,0,1},"leg_b",In),
            leg_c=Spin_leg({1,2},"leg_c",Out),
            leg_d=Spin_leg({0,1,1},"leg_d",Out);
    IQTensor T(leg_a(2),leg_b(4),leg_c(2),leg_d(3)),
             A(leg_a,leg_b), B;
    T.randomize();
    Print(T);
    factor(T,A,B,{"ShowEigs",true,"Maxm",30});
    Print(norm(T-A*B));
    return 0;
}

