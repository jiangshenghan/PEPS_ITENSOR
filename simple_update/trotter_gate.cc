
#include "trotter_gate.h"

NN_Heisenberg_Trotter::NN_Heisenberg_Trotter(std::array<IQIndex,2> nn_sites, int t): 
    phy_legs_{{dag(nn_sites[0]),prime(nn_sites[0])},{dag(nn_sites[1]),prime(nn_sites[1])}},
    virt_legs_{Spin_leg(std::vector<int>{1,0,1},"virt_leg1",Out),
               Spin_leg(std::vector<int>{1,0,1},"virt_leg2",Out)},
    t_(t)
{
    for (int sitei=0; sitei<2; sitei++)
    {
        site_tensors[sitei]=IQTensor(phys_leg_[sitei][0],phys_leg_[sitei][1],virt_legs_[sitei]);
    }
    bond_tensor_=IQTensor(dag(virt_legs_[0]),dag(virt_legs[1]));

    init_site_tensors();
    init_bond_tensor();
}
