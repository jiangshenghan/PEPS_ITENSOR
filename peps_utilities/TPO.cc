
#include "TPO.h"

IQTPO SpinSpin()
{
    std::array<IQIndex,2> phys_legs={Spin_leg({0,1},"phys_leg 0",Out,Site),Spin_leg({0,1},"phys_leg 1",Out,Site)};
    return SpinSpin(phys_legs);
}

IQTPO SpinSpin(const std::array<IQIndex,2> &phys_legs)
{
    std::array<IQIndex,2> virt_legs{Spin_leg({0,0,1},"virt_leg 0"),Spin_leg({0,0,1},"virt_leg 1")};
    IQTPO SdotS(2,1);

    //Init site tensors
    for (int sitei=0; sitei<2; sitei++)
    {
        auto &tensor=SdotS.site_tensors(sitei);
        Singlet_Tensor_Basis tensor_basis(std::vector<IQIndex>{dag(phys_legs[sitei]),prime(phys_legs[sitei]),virt_legs[sitei]});
        tensor=tensor_basis[0]*std::sqrt(3.)/2.;
    }
    //Init bond tensor
    Singlet_Tensor_Basis bond_tensor_basis(std::vector<IQIndex>{dag(virt_legs[0]),dag(virt_legs[1])});
    SdotS.bond_tensors(0)=-std::sqrt(3.)*bond_tensor_basis[0];

    return SdotS;
}


std::array<IQTPO,3> vectorSpinChirality()
{
    std::array<IQIndex,2> phys_legs={Spin_leg({0,1},"phys_leg 0",Out,Site),Spin_leg({0,1},"phys_leg 1",Out,Site)};
    return vectorSpinChirality(phys_legs);
}

std::array<IQTPO,3> vectorSpinChirality(const std::array<IQIndex,2> &phys_legs)
{
    //In fact, quantum number of virt_leg is not well defined
    IQIndex virt_leg("virt_leg",Index("virt_leg_ind",2),QN(0,0),Out);
    std::array<IQTPO,3> SxS{IQTPO(2,0),IQTPO(2,0),IQTPO(2,0)};

    //Init site tensors
    for (auto &tensor_operator : SxS)
    {
        tensor_operator.site_tensors(0)=IQTensor(dag(phys_legs[0]),prime(phys_legs[0]),virt_leg);
        tensor_operator.site_tensors(1)=IQTensor(dag(phys_legs[1]),prime(phys_legs[1]),dag(virt_leg));
    }
    //(SxS)_x=S_yS_z-S_zS_y
    SxS[0].site_tensors(0).set(dag(phys_legs[0])(1),prime(phys_legs[0])(2),virt_leg(1),0.5*Complex_i);
    SxS[0].site_tensors(0).set(dag(phys_legs[0])(2),prime(phys_legs[0])(1),virt_leg(1),-0.5*Complex_i);
    SxS[0].site_tensors(0).set(dag(phys_legs[0])(1),prime(phys_legs[0])(1),virt_leg(2),-0.5);
    SxS[0].site_tensors(0).set(dag(phys_legs[0])(2),prime(phys_legs[0])(2),virt_leg(2),0.5);
    SxS[0].site_tensors(1).set(dag(phys_legs[1])(1),prime(phys_legs[1])(1),dag(virt_leg(1)),0.5);
    SxS[0].site_tensors(1).set(dag(phys_legs[1])(2),prime(phys_legs[1])(2),dag(virt_leg(1)),-0.5);
    SxS[0].site_tensors(1).set(dag(phys_legs[1])(1),prime(phys_legs[1])(2),dag(virt_leg(2)),0.5*Complex_i);
    SxS[0].site_tensors(1).set(dag(phys_legs[1])(2),prime(phys_legs[1])(1),dag(virt_leg(2)),-0.5*Complex_i);
    //PrintDat(SxS[0]);
    //(SxS)_y=S_zS_x-S_xS_z
    SxS[1].site_tensors(0).set(dag(phys_legs[0])(1),prime(phys_legs[0])(1),virt_leg(1),0.5);
    SxS[1].site_tensors(0).set(dag(phys_legs[0])(2),prime(phys_legs[0])(2),virt_leg(1),-0.5);
    SxS[1].site_tensors(0).set(dag(phys_legs[0])(1),prime(phys_legs[0])(2),virt_leg(2),-0.5);
    SxS[1].site_tensors(0).set(dag(phys_legs[0])(2),prime(phys_legs[0])(1),virt_leg(2),-0.5);
    SxS[1].site_tensors(1).set(dag(phys_legs[1])(1),prime(phys_legs[1])(2),dag(virt_leg(1)),0.5);
    SxS[1].site_tensors(1).set(dag(phys_legs[1])(2),prime(phys_legs[1])(1),dag(virt_leg(1)),0.5);
    SxS[1].site_tensors(1).set(dag(phys_legs[1])(1),prime(phys_legs[1])(1),dag(virt_leg(2)),0.5);
    SxS[1].site_tensors(1).set(dag(phys_legs[1])(2),prime(phys_legs[1])(2),dag(virt_leg(2)),-0.5);
    //PrintDat(SxS[1]);
    //(SxS)_z=S_xS_y-S_yS_x
    SxS[2].site_tensors(0).set(dag(phys_legs[0])(1),prime(phys_legs[0])(2),virt_leg(1),0.5);
    SxS[2].site_tensors(0).set(dag(phys_legs[0])(2),prime(phys_legs[0])(1),virt_leg(1),0.5);
    SxS[2].site_tensors(0).set(dag(phys_legs[0])(1),prime(phys_legs[0])(2),virt_leg(2),-0.5*Complex_i);
    SxS[2].site_tensors(0).set(dag(phys_legs[0])(2),prime(phys_legs[0])(1),virt_leg(2),0.5*Complex_i);
    SxS[2].site_tensors(1).set(dag(phys_legs[1])(1),prime(phys_legs[1])(2),dag(virt_leg(1)),0.5*Complex_i);
    SxS[2].site_tensors(1).set(dag(phys_legs[1])(2),prime(phys_legs[1])(1),dag(virt_leg(1)),-0.5*Complex_i);
    SxS[2].site_tensors(1).set(dag(phys_legs[1])(1),prime(phys_legs[1])(2),dag(virt_leg(2)),0.5);
    SxS[2].site_tensors(1).set(dag(phys_legs[1])(2),prime(phys_legs[1])(1),dag(virt_leg(2)),0.5);

    for (const auto &tensor_operator : SxS)
    {
        auto tensor=(tensor_operator.site_tensors(0)*tensor_operator.site_tensors(1));
        tensor.clean();
        PrintDat(tensor);
    }

    return SxS;
}
