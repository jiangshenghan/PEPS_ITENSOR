
#ifndef _CLUSTER_ENV_H_
#define _CLUSTER_ENV_H_

#include "peps.h"

//algorithm to obtain direct product environment for a small cluster of tensor wavefunction. The cluster contains loop in general
//0. we want to get entanglement information between leg 'a' and leg 'b'
//1. Contract cluster tensor T, with legs beyond 'a' and 'b' grounped as leg 's', get a T_{s,ab}
//2. Do QR decompostion for T_{sb,a}=Q_{sb,a'}R_{a'a}
//3. Do similar thing for leg 'b', we get T_{b,sa}=L_{b'b}Q_{b',sa}
//4. Insert R^{-1}RLL^{-1} between leg 'a' and 'b': T_{s,ab}=T^{s,ab}R^{-1}_{aa'}R_{a'c}L_{cb'}L^{-1}_{b'b} 
//5. Do SVD on RL=U.D.V. 
//6. FIXME: put sqrt(D) on boundary leg
//7. Repeat the above procedure until D (normalized) converges
//
//For legs with no spin deg, R,L is diag, and the diag element is norm of col (row) of T_{sb,a} (T_{b,sa}). Then, sqrt(D)=sqrt(normA*normB)
//
//TODO: consider flavor deg case!

class Cluster_Env
{
    public:
        //
        //Constructor
        //
        Cluster_Env(const std::string &name, const std::vector<int> &cutting_sites, const std::vector<IQTensor> &site_tensors, double cut_off=1e-8);

        //Acess Methods
        const IQTensor &sl_env_tensor() const
        {
            return sl_env_tensor_;
        }
        
        //
        //Methods
        //
        void init_boundary_legs();
        void init_sl_env();
        void obtain_sl_env_iterative_nodeg();
        void obtain_dl_boundary_tensors(std::vector<IQTensor> &dl_boundary_tensors);
        IQTensor update_sl_env_one_step_nodeg(const std::vector<IQTensor> &dl_boundary_tensors);
        //we require the updated site tensors share the same indices as original ones
        void update_site_tensors(const std::vector<IQTensor> &site_tensors) { site_tensors_=site_tensors; }

    private:
        //name_ characterize the shape of cluster
        //1. kagome normal lattice
        //  a. triangle shape
        //          \ /
        //           v
        //          / \
        //         /   \
        //      --u-  --w--
        //       /       \
        //
        std::string name_;
        std::vector<int> cutting_sites_;
        IQIndex cutting_leg_;
        std::vector<std::vector<IQIndex>> boundary_legs_;
        //we already absorb bond tensors to site_tensors_
        std::vector<IQTensor> site_tensors_;
        //sl_env_tensor_ stores single layer env tensor 
        IQTensor sl_env_tensor_;
        double cut_off_;

};

#endif
