
#include "transfer_to_square.h"

IQPEPS spin_sym_square_peps_from_honeycomb_tensor_uc(const std::vector<IQTensor> &honeycomb_site_tensors_uc, const std::vector<IQTensor> &honeycomb_bond_tensors_uc, const Lattice_Base &square_lattice)
{
    IQTensor square_site_tensor_no_order;
    square_site_tensor_no_order=honeycomb_site_tensors_uc[0]*honeycomb_bond_tensors_uc[0]*honeycomb_site_tensors_uc[1];

    std::vector<IQIndex> phys_legs, virt_legs;

    for (int sitei=0; sitei<2; sitei++)
    {
        for (const auto &ind : honeycomb_site_tensors_uc[sitei].indices())
        {
            if (ind.type()==Site) phys_legs.push_back(ind);
            if (ind.type()==Link) virt_legs.push_back(ind);
        }
    }

    IQCombiner phys_legs_combiner(phys_legs[0],phys_legs[1]);
    phys_legs_combiner.init("square_phys_leg",Site,Out);
    auto square_phys_leg=phys_legs_combiner.right();
    square_site_tensor_no_order=square_site_tensor_no_order*phys_legs_combiner;

    //get the right order of square_site_tensor, which should set to be phys,1,5,4,2! Now, the order for square_site_tensor_no_order is 4,5,1,2,phys
    //the right order of site tensor is used to initialize square lattice peps, where we use replaceIndex
    IQTensor square_site_tensor_uc(square_phys_leg,virt_legs[1],virt_legs[5],virt_legs[4],virt_legs[2]);
    std::vector<int> leg_dims={square_phys_leg.m(),virt_legs[1].m(),virt_legs[5].m(),virt_legs[4].m(),virt_legs[2].m()};
    int max_vals=square_phys_leg.m()*virt_legs[1].m()*virt_legs[5].m()*virt_legs[4].m()*virt_legs[2].m();

    for (int vals=0; vals<max_vals; vals++)
    {
        auto val_list=list_from_num(vals,leg_dims);
        square_site_tensor_uc(square_phys_leg(val_list[0]+1),virt_legs[1](val_list[1]+1),virt_legs[5](val_list[2]+1),virt_legs[4](val_list[3]+1),virt_legs[2](val_list[4]+1))=square_site_tensor_no_order(square_phys_leg(val_list[0]+1),virt_legs[1](val_list[1]+1),virt_legs[5](val_list[2]+1),virt_legs[4](val_list[3]+1),virt_legs[2](val_list[4]+1));
    }
    square_site_tensor_uc.clean();

    //TODO:fix the order of bond tensor, which will be important for IGG nontrivial case
    std::vector<IQTensor> square_bond_tensors_uc={honeycomb_bond_tensors_uc[1],honeycomb_bond_tensors_uc[2]};

    //get the phys_leg_qn_order
    int phys_leg_qn_order=0;
    std::vector<int> sz_qns;
    for (const auto &indqn : square_phys_leg.indices())
    {
        sz_qns.push_back(indqn.qn.sz());
    }
    for (int i=0; i<sz_qns.size()-1; i++)
    {
        if (phys_leg_qn_order==0)
        {
            if (sz_qns[i]<sz_qns[i+1]) phys_leg_qn_order=1;
            if (sz_qns[i]>sz_qns[i+1]) phys_leg_qn_order=-1;
        }

        if (phys_leg_qn_order==1) assert(sz_qns[i]<sz_qns[i+1]);
        if (phys_leg_qn_order==-1) assert(sz_qns[i]>sz_qns[i+1]);
    }

    //get spin rep for phys leg and virt leg
    std::vector<int> phys_leg_spin_deg, virt_leg_spin_deg;
    iqind_spin_rep(square_phys_leg,phys_leg_spin_deg);
    iqind_spin_rep(virt_legs[0],virt_leg_spin_deg);

    //create square peps with spin sym
    IQPEPS_IndexSet_Spin_Sym square_indexset(square_phys_leg.m(),virt_legs[0].m(),phys_leg_spin_deg,virt_leg_spin_deg,square_lattice,phys_leg_qn_order);

    return IQPEPS(square_lattice,square_indexset,std::vector<IQTensor>{square_site_tensor_uc},square_bond_tensors_uc);
}
