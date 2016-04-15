
#include "peps.h"
//#include "kagome_rvb.h"


int main()
{
    IQIndex spin_one_phys_leg=Spin_leg({0,0,1},"phys_leg");
    std::vector<IQIndex> spin_half_virt_legs(4);
    int legi=0;
    for (auto &leg : spin_half_virt_legs)
    {
        leg=IQIndex(nameint("S=1/2 virt_leg ",legi),
                Index(nameint("Up",legi),1,Link),QN(+1,0),
                Index(nameint("Down",legi),1,Link),QN(-1,0));
        legi++;
    }

    std::vector<IQIndex> spin_zero_one_virt_legs(2);
    legi=0;
    for (auto &leg : spin_zero_one_virt_legs)
    {
        leg=Spin_leg({1,0,1},nameint("virt_leg",legi));
        legi++;
    }
    
    Singlet_Tensor_Basis tripletA(std::vector<IQIndex>{spin_half_virt_legs[0],spin_half_virt_legs[1],spin_one_phys_leg}), 
                         tripletB(std::vector<IQIndex>{spin_half_virt_legs[2],spin_half_virt_legs[3],spin_one_phys_leg}),
                         singlet_A(std::vector<IQIndex>{spin_half_virt_legs[2],spin_half_virt_legs[3]}),
                         singlet_B(std::vector<IQIndex>{spin_half_virt_legs[0],spin_half_virt_legs[1]}),
                         fuse_left(std::vector<IQIndex>{dag(spin_half_virt_legs[0]),dag(spin_half_virt_legs[2]),spin_zero_one_virt_legs[0]}),
                         fuse_right(std::vector<IQIndex>{dag(spin_half_virt_legs[1]),dag(spin_half_virt_legs[3]),spin_zero_one_virt_legs[1]}),
                         site_basis(std::vector<IQIndex>{spin_zero_one_virt_legs[0],spin_zero_one_virt_legs[1],spin_one_phys_leg});

    //PrintDat(fuse_left[0]);
    //PrintDat(fuse_left[1]);
    PrintDat(site_basis);

    IQTensor siteA=tripletA[0]*singlet_A[0]*(fuse_left[0]+sqrt(3.)*fuse_left[1])*(fuse_right[0]+sqrt(3.)*fuse_right[1]),
             siteB=tripletB[0]*singlet_B[0]*(fuse_left[0]+sqrt(3.)*fuse_left[1])*(fuse_right[0]+sqrt(3.)*fuse_right[1]),
             bond=singlet_A[0]*singlet_B[0]*(fuse_left[0]+sqrt(3.)*fuse_left[1])*(fuse_right[0]+sqrt(3.)*fuse_right[1]);

    PrintDat(siteA);
    PrintDat(siteB);
    PrintDat(bond);

    PrintDat(siteA*dag(site_basis[0]));
    PrintDat(siteA*dag(site_basis[1]));
    PrintDat(siteA*dag(site_basis[2]));

    PrintDat(siteB*dag(site_basis[0]));
    PrintDat(siteB*dag(site_basis[1]));
    PrintDat(siteB*dag(site_basis[2]));


    return 0;
}

