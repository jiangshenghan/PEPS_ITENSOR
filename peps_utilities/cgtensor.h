
//This file defines tensors store information about CG coefficient

#ifndef _CGTENSOR_H_
#define _CGTENSOR_H_

#include "utilities.h"

class CGTensors
{
    public:
        //
        //Constructor
        //
        //spin_qn stores total spin quantum number (times 2) for each legs
        CGTensors() {}
        CGTensors(const std::vector<Spin> &spins);

        //
        //Access method
        //
        bool valid() { return valid_; }
        const std::vector<IQTensor> &K() const { return K_; }
        
        //
        //Construct helper
        //
        //V_{a_1a_2...}^{b_1b_2...}=\oplus_{c}(V_{a_1a_2...}^u\otimes V_c^{b_1b_2...})
        //Further, we have V_{a_1a_2...}^c=\oplus_{d_a,1d_2,...}(V_{a_1a_2}^{d_1}\otimes V_{d_1a_3}^{d_2}...V_{d_{n-2}a_n}^c)
        //If there are only out legs or only in legs, we will make a spin singlet
        //Since cg coefficient is real, we have K_{a_1a_2...}^c=K_c^{a_1a_2...}
        void init();
        //to obtain all sets of different mediate spins
        //# of sets = # of fusion channel
        //mediate spins are used to decompose K's to CG coefficient
        //bool obtain_mediate_spins(const std::vector<int> &out_spins, const std::vector<int> &in_spins, std::vector<std::vector<int> > &mediate_spins_sets);
        
        //obtain all possible sets of K_'s by performing CG decomposition
        //K_{a_1a_2...}^{b_1b_2...}=K_{a_1a_2...}^c.K_c^{b_1b_2...}
        //K_{a_1a_2...}^c=K_{a_1a_2}^{d_1}...K_{d_{n-2}a_n}^c
        //different choice of c and d_i corresponds to different fusion channel
        //Here we use recursive method to obtain K
        bool obtain_K(const std::vector<IndexSpin> &out_spin_legs, const std::vector<IndexSpin> &in_spin_legs, std::vector<IQTensor> &K);

        //functions to get tensor K_{S_1S_2}^{S_3}
        IQTensor obtain_CG(const IndexSpin &S1, const IndexSpin &S2, const IndexSpin &S3);

    private:
        std::vector<IndexSpin> spin_legs_;
        //there may be more than one fusion channel
        std::vector<IQTensor> K_;
        //valid_=false means K_ is empty. Namely, the spin_qns are inconsistent
        bool valid_;
};

#endif
