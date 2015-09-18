
#ifndef _SINGLET_TENSOR_BASIS_
#define _SINGLET_TENSOR_BASIS_

#include "cgtensor.h"

//Singlet_Tensor_Basis generates all spin singlet basis for a given set of legs which accommodate rep for SU(2) symmetry
class Singlet_Tensor_Basis
{
    public:
        //
        //Type Aliens
        //

        //
        //Constructor
        //
        Singlet_Tensor_Basis() {}

        explicit Singlet_Tensor_Basis(const std::vector<IQIndex> &iqinds);

        explicit Singlet_Tensor_Basis(const IndexSet<IQIndex> &iqinds_set);

        //
        //Access Methods
        //
        const IndexSet<IQIndex> &indices() const { return is_; }
        const std::vector<IQTensor> &tensors() const {return singlet_tensors_; }
        const std::vector<int> &spin_configs(int i) const { return spin_configs_[i]; }
        const IQTensor &operator()(int i) const { return singlet_tensors_[i]; }

        //
        //Constructor Helper
        //
        void init_spin_deg_and_basis();
        void init_singlet_tensors();

    private:
        IndexSet<IQIndex> is_;
        //is_spin_degs_[i] stores spin deg for IQIndex is_[i]
        std::vector<std::vector<int>> is_spin_degs_;
        //is_spin_basis_[i][j] stores spin basis |S,m,t\rangle for is_[i][j]
        std::vector<std::vector<Spin_Basis>> is_spin_basis_;
        std::vector<IQTensor> singlet_tensors_;
        //spin_configs_ stores the spin_list for singlet_tensors_
        std::vector<std::vector<int>> spin_configs_;


        //
        //friend functions
        //
        friend bool iqind_spin_rep(const IQIndex &sz_leg, std::vector<int> &spin_deg);

        friend bool iqind_to_spin_basis(const IQIndex &sz_leg, std::vector<Spin_Basis> &spin_basis);

        friend bool iqind_to_spin_basis(const IQIndex &sz_leg, const std::vector<int> &spin_deg, std::vector<Spin_Basis> &spin_basis);

};

inline std::ostream &operator<<(std::ostream &s, const Singlet_Tensor_Basis &tensor_basis)
{
    for (const auto &tensor : tensor_basis.tensors())
    {
        s << tensor << endl;
    }

    return s;
}


#endif
