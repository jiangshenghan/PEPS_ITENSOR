
//This file defines tensors store information about CG coefficient

#ifndef _CGTENSOR_H_
#define _CGTENSOR_H_

#include "utilities.h"

//class IndexSpin is an IQIndex which stores indices with fixed spin quantum number
class IndexSpin 
{
    public:
        //
        //Constructor
        //
        IndexSpin() {}

        explicit IndexSpin(const IQIndex &iqind) : spin_(iqind.nindex()-1), sz_legs_(iqind) { assert(valid()); }

        IndexSpin(const IQIndex &iqind, int spin) : spin_(spin), sz_legs_(iqind) { assert(valid()); }

        //construct IndexSpin from spin quantum number, no extra deg
        explicit IndexSpin(int spin, Arrow dir=Out): spin_(spin)
        {
            std::vector<IndexQN> sz_indqns;
            for (int sz=spin_; sz>=-spin_; sz-=2)
            {
                std::stringstream ss;
                ss << "S=" << spin_/2.0 << ", Sz=" << sz/2.0 << " leg";
                std::string sz_leg_name=ss.str();
                sz_indqns.push_back(IndexQN(Index(sz_leg_name,1),QN(sz)));
            }
            sz_legs_=IQIndex(nameint("leg with 2S=", spin_),sz_indqns,dir);
            
            assert(valid());
        }

        //
        //Access method
        //
        Arrow dir() const { return sz_legs_.dir(); }

        int spin_qn() const { return spin_; }
        
        const IQIndex &leg() const { return sz_legs_; }

        //return the first found IQIndexVal with qn=sz
        IQIndexVal indval (int sz) const
        {
            int n=0;
            for (const auto &indqn : sz_legs_.indices())
            {
                if (indqn.qn.sz()==sz) return sz_legs_.operator()(n+1);
                n+=indqn.m();
            }
            return IQIndexVal();
        }


        IQIndexVal operator()(int n) const
        {
            return sz_legs_.operator()(n);
        }
        
        int m() const { return sz_legs_.m(); }

        //
        //Method
        //
        IndexSpin &dag() { sz_legs_.dag(); return *this; }

        //
        //Constructor helpers
        //
        //To determine if the input IQIndex and spin are valid for a SU(2) representation \mathbb{D}\otime\mathbb{V}_S
        bool valid()
        {
            
            int dim=sz_legs_[0].m();
            for (const auto &indqn : sz_legs_.indices())
            {
                int sz=indqn.qn.sz();
                if (dim!=indqn.m() || //consistent extra deg
                    std::abs(sz)>spin_ ||
                    std::abs(sz)%2!=spin_%2) //check sz qn
                {
                    return false;
                }
            }
            return true;   //every sz has same extra deg, all int/half-int, and go from -spin,...,spin
        }

    private:
        int spin_;
        IQIndex sz_legs_;
};

class CGTensors
{
    public:
        //
        //Constructor
        //
        //spin_qn stores total spin quantum number (times 2) for each legs
        CGTensors() {}
        CGTensors(const std::vector<int> &spins, const std::vector<Arrow> &dirs);
        CGTensors(const std::vector<IQIndex> &sz_legs);
        CGTensors(const std::vector<IndexSpin> &spin_legs);

        //
        //Access method
        //
        bool valid() { return valid_; }
        const std::vector<IQTensor> &K() const { return K_; }
        const std::vector<IndexSpin> &spin_indices() const { return spin_legs_; }
        const IQIndex &iqindice(int i) const { return spin_legs_[i].leg(); }
        
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


inline std::ostream &operator<<(std::ostream &s, const IndexSpin &indspin)
{
    return s << indspin.spin_qn() << indspin.leg();
}

inline std::ostream &operator<<(std::ostream &s, const CGTensors &cg_tensors)
{
    for (const auto &tensor : cg_tensors.K())
    {
        s << tensor << endl;
    }
    return s;
}

#endif
