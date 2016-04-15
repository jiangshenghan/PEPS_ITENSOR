
#include "single_layer_env.h"

SL_Env_Tensors::SL_Env_Tensors(const std::vector<IQTensor> &site_tensors, double cut_off): site_tensors_(site_tensors), cut_off_(cut_off)
{
    comm_leg_=commonIndex(site_tensors_[0],site_tensors_[1]);

    //init env_tensors_ 
    for (const auto &tensor : site_tensors_)
    {
        for (const auto &leg : tensor.indices())
        {
            if (leg.type()==Site || leg==comm_leg_) continue;

            IQTensor env_tensor(dag(leg),prime(leg));
            for (int val=1; val<=leg.m(); val++)
            {
                env_tensor(dag(leg)(val),prime(leg)(val))=1.;
            }

            env_tensors_.push_back(env_tensor);
        }
    }

    //init spectrum_sqrt_
    int dim=comm_leg_.m();
    spectrum_sqrt_=std::vector<double>(dim,1.);
}

void SL_Env_Tensors::update_env_tensor_one_step()
{
    //get env & site combined_tensors
    std::vector<IQTensor> combined_tensors=site_tensors_;
    for (auto &env_tensor : env_tensors_)
    {
        if (commonIndex(env_tensor,combined_tensors[0])!=IQIndex::Null())
        {
            combined_tensors[0]*=env_tensor;
        }
        else
        {
            combined_tensors[1]*=env_tensor;
        }
    }

    //Perform density decomposition on two combined tensor separately
    //combined_tensors[i]=U[i]*DV_tensors[i], where U[i] is unitary and DV_tensors are square matrix
    //We should also compute inverse of DV_tensors for further use
    std::vector<IQTensor> DV_tensors, DV_inv_tensors;
    for (const auto &tensor : combined_tensors)
    {
        std::vector<IQIndex> U_legs;
        for (const auto leg : tensor.indices())
        {
            if (leg==comm_leg_) continue;
            U_legs.push_back(leg);
        }

        //IQTensor U(U_legs), DV;
        //denmatDecomp(tensor,U,DV,Fromleft);
        //PrintDat(U);

        IQTensor U(U_legs),D,V,DV,DV_inv;
        svd(tensor,U,D,V);
        DV=D*V;
        DV_inv=inverse_rank2_tensor_by_arma_mat(DV);

        //get inverse of DV
        
        DV_tensors.push_back(DV);
        DV_inv_tensors.push_back(DV_inv);
    }
    PrintDat(DV_inv_tensors);
    PrintDat(DV_tensors[0]*DV_tensors[1]);

    //obtain new spectrum_sqrt and new env_tensors
    IQTensor spectrum_sqrt_tensor,left_tensor,right_tensor;
    for (const auto &leg : DV_tensors[0].indices())
    {
        if (leg==comm_leg_) continue;
        left_tensor=IQTensor(leg);
    }
    svd(DV_tensors[0]*DV_tensors[1],left_tensor,spectrum_sqrt_tensor,right_tensor);
    //TODO: take care of qn direction on spectrum_sqrt_tensor
    IQTensor update_env_tensor=spectrum_sqrt_tensor*left_tensor*DV_inv_tensors[0];

    PrintDat(spectrum_sqrt_tensor);
    PrintDat(left_tensor);
    PrintDat(right_tensor);
    PrintDat(update_env_tensor);
}
