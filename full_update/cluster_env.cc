
#include "cluster_env.h"

Cluster_Env::Cluster_Env(const std::string &name, const std::vector<int> &cutting_sites, const std::vector<IQTensor> &site_tensors, double cut_off): 
    name_(name),
    cutting_sites_(cutting_sites),
    site_tensors_(site_tensors),
    cut_off_(cut_off)
{
    cutting_leg_=commonIndex(site_tensors_[cutting_sites_[0]],site_tensors_[cutting_sites_[1]]);

    init_boundary_legs();
    init_sl_env();
}


void Cluster_Env::init_boundary_legs()
{
    for (int i=0; i<site_tensors_.size(); i++)
    {
        const auto &tensA=site_tensors_[i];
        std::vector<IQIndex> tensA_bulk_legs;
        for (int j=0; j<site_tensors_.size(); j++)
        {
            if (i==j) continue;
            const auto &tensB=site_tensors_[j];
            //we assume tensA and tensB at most share one common leg
            IQIndex temp_leg=commonIndex(tensA,tensB);
            if (temp_leg!=IQIndex::Null())
            {
                tensA_bulk_legs.push_back(temp_leg);
            }
        }
        //Print(tensA_bulk_legs);
        std::vector<IQIndex> tensA_boundary_legs;
        for (const auto &tensA_leg : tensA.indices())
        {
            auto leg_iter=std::find(tensA_bulk_legs.begin(),tensA_bulk_legs.end(),tensA_leg);
            if (leg_iter==tensA_bulk_legs.end() && tensA_leg.type()==Link)
            {
                tensA_boundary_legs.push_back(tensA_leg);
            }
        }
        boundary_legs_.push_back(tensA_boundary_legs);
    }
    //Print(boundary_legs_);
}

void Cluster_Env::init_sl_env()
{
    sl_env_tensor_=IQTensor(dag(cutting_leg_),prime(cutting_leg_));
    for (int val=1; val<=cutting_leg_.m(); val++)
    {
        sl_env_tensor_(dag(cutting_leg_)(val),prime(cutting_leg_)(val))=std::abs(rand_gen());
    }
    //normalize to sqrt(m)
    double env_norm=sl_env_tensor_.norm();
    sl_env_tensor_*=sqrt(cutting_leg_.m()*1.)/env_norm;
}

void Cluster_Env::obtain_sl_env_iterative_nodeg()
{
    std::vector<IQTensor> dl_boundary_tensors;
    IQTensor updated_sl_env_tensor=sl_env_tensor_;
    int iter=0;
    do
    {
        sl_env_tensor_=updated_sl_env_tensor;
        obtain_dl_boundary_tensors(dl_boundary_tensors);
        updated_sl_env_tensor=update_sl_env_one_step_nodeg(dl_boundary_tensors);
        //Print(iter);
        //PrintDat(sl_env_tensor_);
        iter++;
        if (iter>100) init_sl_env();
    }
    while ((updated_sl_env_tensor-sl_env_tensor_).norm()>cut_off_);
    Print(iter);
    PrintDat(sl_env_tensor_);
}

void Cluster_Env::obtain_dl_boundary_tensors(std::vector<IQTensor> &dl_boundary_tensors)
{
    dl_boundary_tensors.clear();
    for (int sitei=0; sitei<site_tensors_.size(); sitei++)
    {
        if (boundary_legs_[sitei].empty()) continue;
        //multiply double layer env tensor of each boundary legs of sitei tensor to sitei tensor
        IQTensor temp_dl_tensor=site_tensors_[sitei];
        for (const auto &boundary_leg : boundary_legs_[sitei])
        {
            IQTensor temp_env_tensor=sl_env_tensor_;
            if (cutting_leg_.dir()==boundary_leg.dir())
            {
                temp_env_tensor.replaceIndex(dag(cutting_leg_),dag(boundary_leg));
                temp_env_tensor.replaceIndex(prime(cutting_leg_),prime(boundary_leg));
            }
            else
            {
                temp_env_tensor.replaceIndex(dag(cutting_leg_),boundary_leg);
                temp_env_tensor.replaceIndex(prime(cutting_leg_),dag(prime(boundary_leg)));
                temp_env_tensor.dag();
            }
            temp_env_tensor=temp_env_tensor*dag(swapPrime(temp_env_tensor,0,2));
            temp_env_tensor.mapprime(2,1);
            temp_dl_tensor*=temp_env_tensor;
        }
        //get double layer boundary tensor of sitei
        temp_dl_tensor*=dag(site_tensors_[sitei]).prime(Link);
        //PrintDat(temp_dl_tensor);
        dl_boundary_tensors.push_back(temp_dl_tensor);
    }
}

IQTensor Cluster_Env::update_sl_env_one_step_nodeg(const std::vector<IQTensor> &dl_boundary_tensors)
{
    IQTensor updated_sl_env_tensor=sl_env_tensor_;
    if (name_.find("kagome normal")!=std::string::npos && name_.find("triangle shape")!=std::string::npos)
    {
        int env_site;
        for (int sitei=0; sitei<site_tensors_.size(); sitei++)
        {
            if (sitei==cutting_sites_[0] || sitei==cutting_sites_[1]) continue;
            env_site=sitei;
            break;
        }
        //obtain tensA_normsq
        IQTensor tensA_normsq=dl_boundary_tensors[env_site]*dl_boundary_tensors[cutting_sites_[0]]*(trace(dl_boundary_tensors[cutting_sites_[1]],dag(cutting_leg_),prime(cutting_leg_))),
                 tensB_normsq=dl_boundary_tensors[env_site]*dl_boundary_tensors[cutting_sites_[1]]*(trace(dl_boundary_tensors[cutting_sites_[0]],cutting_leg_,dag(prime(cutting_leg_))));
        //PrintDat(tensA_normsq);
        //PrintDat(tensB_normsq);
        
        for (int diagi=1; diagi<=cutting_leg_.m(); diagi++)
        {
            double normA=sqrt(tensA_normsq(cutting_leg_(diagi),dag(prime(cutting_leg_))(diagi))),
                   normB=sqrt(tensB_normsq(dag(cutting_leg_)(diagi),prime(cutting_leg_)(diagi)));
            //updated_sl_env_tensor(dag(cutting_leg_)(diagi),prime(cutting_leg_)(diagi))=sqrt(normA*normB);
            updated_sl_env_tensor(dag(cutting_leg_)(diagi),prime(cutting_leg_)(diagi))=normA*normB;
        }
    }

    double env_norm=updated_sl_env_tensor.norm();
    updated_sl_env_tensor*=sqrt(cutting_leg_.m()*1.)/env_norm;
    //PrintDat(updated_sl_env_tensor);

    return updated_sl_env_tensor;
}
