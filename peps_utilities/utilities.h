
#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include "peps_class_params.h"

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


inline std::ostream &operator<<(std::ostream &s, const Coordinate &coord)
{
    return s << "(" << coord[0] << "," << coord[1] << "," << coord[2] << ")";
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
//This function transfers a list to a num
int num_from_list(const std::vector<int> &list, const std::vector<int> &max_nums);


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


//This function construct eta from mu
//the Out leg has prime while In leg has no prime
IQTensor eta_from_mu(double mu, const std::vector<int> &spin_deg);
//construct eta by given leg
IQTensor eta_from_mu(double mu, IQIndex eta_leg);
//define the function for Index for convience
ITensor eta_from_mu(double mu, Index eta_leg);

#endif
