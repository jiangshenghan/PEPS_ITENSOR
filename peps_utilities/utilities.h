
#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include "predef.h"

class Spin_Basis
{
    public:
        Spin_Basis() {}

        Spin_Basis(const std::array<int,3> &basis_label) : basis_label_(basis_label) {}

        Spin_Basis(int S, int m, int t=0) : basis_label_{{S,m,t}} {}

        int S() const { return basis_label_[0]; }

        int m() const { return basis_label_[1]; }

        int t() const { return basis_label_[2]; }

        const std::array<int,3> &basis_label() const
        {
            return basis_label_;
        }

        bool operator==(const Spin_Basis &other)
        {
            return (basis_label_==other.basis_label());
        }

    private:
        std::array<int,3> basis_label_;
};

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


inline std::ostream &operator<<(std::ostream &s, const Coordinate &coord)
{
    return s << "(" << coord[0] << "," << coord[1] << "," << coord[2] << ")";
}

inline std::ostream &operator<<(std::ostream &s, const IndexSpin &indspin)
{
    return s << indspin.spin_qn() << indspin.leg();
}

inline std::ostream &operator<<(std::ostream &s, const Spin_Basis &spin_basis)
{
    return s << "(" << spin_basis.S()/2.0 << "," << spin_basis.m()/2.0 << "," << spin_basis.t() << ")";
}

inline std::ostream &operator<<(std::ostream &s, const std::vector<int> &ivec)
{
    for (const auto &i : ivec) s<< i << " ";
    return s;
}

//This function transfers a number to a list
std::vector<int> list_from_num(int num, const std::vector<int> &max_nums);

//Given sz quantum numbers, get corresponding spin representation
//sz_deg[qn] stores the degeneracy of 2Sz=spin_dim-qn-1
bool sz_to_spin(const std::vector<int> &sz_deg, std::vector<int> &spin_deg);

//Given spin quantum numbers, get corresponding sz quantum numbers
void spin_to_sz(const std::vector<int> &spin_deg, std::vector<int> &sz_deg);


//Given spin_deg this function creates an IQIndex leg,
//where the leg accommodates rep of spin_qn with deg equals to spin_deg[spin_qn] 
IQIndex Spin_leg(const std::vector<int> &spin_deg, const std::string &iqind_name, Arrow dir=Out);

//Given an IQIndex (with sz quantum number), this function determines the corresponding SU(2) rep
//spin_deg[n] stores extra deg for spin n/2
bool iqind_spin_rep(const IQIndex &sz_leg, std::vector<int> &spin_deg);

//Given an IQIndex sz_leg (with sz quantum number) which accommodates SU(2) rep, this function relates integer i (0<=i<sz_leg.m()) and the spin basis (|S,m,t\rangle)
bool iqind_to_spin_basis(const IQIndex &sz_leg, std::vector<Spin_Basis> &spin_basis);

bool iqind_to_spin_basis(const IQIndex &sz_leg, const std::vector<int> &spin_deg, std::vector<Spin_Basis> &spin_basis);


#endif
