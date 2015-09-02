
#ifndef _SIMPLE_UPDATE_H_
#define _SIMPLE_UPDATE_H_

#include "ipeps.h"

simple_update_one_step(iPEPS& psi, vector<ITensor>& trotter_gate, int& link_ind)
{
    //site_ind_A/B are indices of two sites to be acted by the trotter gate
    int site_ind_A=psi.ilattice().link_to_sites(link_ind,0)[2], 
        site_ind_B=psi.ilattice().link_to_sites(link_ind,1)[2];
    //AB_link_tensor connected two sites under trotter gate action
    ITensor &AB_link_tensor=psi.link_tensors(link_ind);
    //T_A and T_B are two site tensors acted by current trotter gate
    ITensor &T_A=psi.site_tensors(site_ind_A),
            &T_B=psi.site_tensors(site_ind_B);
    //A/B_link_tensors store link tensors connected to A/B except the one connects both A and B (AB_link_tensor), those tensors are invariant in this function
    vector<ITensor> A_link_tensors, B_link_tensors;
    //We also stores inverse of A/B_link_tensors in order to update the site tensor
    vector<ITensor> A_link_tensors_inverse, B_link_tensors_inverse;

    //Define environment link tensors as well as their inverse, we use the fact that the link tensors are diag.
    for (int i=0; i<psi.ilattice().n_site_to_links(); i++)
    {
        if (i!=link_ind)
        {
            A_link_tensors.push_back(psi.link_tensors(i));
            B_link_tensors.push_back(psi.link_tensors(i));

            IndexSet<ITensor> link_indices(psi.link_tensors(i).indices);
            ITensor link_tensor_inverse(link_indices);
            
            for (int indval=1; indval<=link_indices[0].m(); inval++)
            {
                link_tensor_inverse(link_indices[0](indval),link_indices[1](indval))=1.0/psi.link_tensors(i)(link_indices[0](indval),link_indices[1](indval));
            }

            A_link_tensors_inverse.push_back(link_tensor_inverse);
            B_link_tensors_inverse.push_back(link_tensor_inverse);
        }
    }

    //set prime level of tensors of A_link_tensors to be 1 and B_link_tensors to be 2
    for (auto& tensor_i : A_link_tensors)
    {
        tensor_i.prime(1);
    }

    for (auto& tensor_i : B_link_tensors)
    {
        tensor_i.prime(2);
    }

    //set the prime level of of T_A to be 1 and T_B to be 2 and leaving legs connect to AB_link_tensor invariant. Also, we will store the primed leg as well as physical leg for future use
    IndexSet<Index> UA_legs, VB_legs;

    for (const auto& leg_i : T_A.indices())
    {
        if (leg_i.type()==Site)
        {
            UA_legs.addindex(leg_i);
            continue;
        }

        if ((leg_i!=AB_link_tensor.indices()[0]) && (leg_i!=AB_link_tensor.indices()[1]))
        {
            T_A.prime(leg_i,1);
            UA_legs.addindex(leg_i);
        }
    }

    for (const auto& leg_i : T_B.indices())
    {
        if (leg_i.type()==Site)
        {
            VB_legs.addindex(leg_i);
            continue;
        }

        if ((leg_i!=AB_link_tensor.indices()[0]) && (leg_i!=AB_link_tensor.indices[1]))
        {
            T_B.prime(leg_i,2);
            VB_legs.addindex(leg_i);
        }
    }

    //The combined site tensor are obtained by multiplying the Trotter gate as well as the environment link tensors
    ITensor TA_Combined=T_A*trotter_gate[0], TB_Combined=T_B*trotter_gate[1];
    for (const auto& tensor_i : A_link_tensors)
    {
        TA_Combined*=tensor_i;
    }

    for (const auto& tensor_i : B_link_tensors)
    {
        TB_Combined*=tensor_i;
    }

    //Do density matrix decomposition on both TA_Combined & TB_Combined, TA_Combined=U_A.R_A, TB=L_B.V_B, where U and V are isometries. 
    OptSet decomp_opts;
    ITensor U_A(UA_legs), VB(VB_legs);
    ITensor R_A, L_B;

    denmatDecomp(TA_Combined,U_A,R_A,Fromleft);
    denmatDecomp(TB_Combined,L_B,V_B,Fromright);

    //Then do svd on R.AB_link_tensor.L, truncate the singular value to update AB_link_tensor. Finally, update TA & TB.
    decomp_opts.add("Maxm",psi.D());

}

#endif
