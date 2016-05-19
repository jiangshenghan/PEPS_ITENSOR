
#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include "predef.h"
#include "tnetwork_storage.h"

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

template <class T>
inline std::ostream &operator<<(std::ostream &s, const std::vector<T> &vec)
{
    for (const auto &elem : vec) s << elem << " ";
    return s;
}

//This function transfers a number to a list
std::vector<int> list_from_num(int num, const std::vector<int> &max_nums);
//This function transfers a list to a num
int num_from_list(const std::vector<int> &list, const std::vector<int> &max_nums);


//Given sz quantum numbers, get corresponding spin representation
//sz_deg[qn] stores the degeneracy of 2Sz=spin_dim-qn-1
bool sz_to_spin(const std::vector<int> &sz_deg, std::vector<int> &flavor_deg);

//Given spin quantum numbers, get corresponding sz quantum numbers
void spin_to_sz(const std::vector<int> &flavor_deg, std::vector<int> &sz_deg);


//Given flavor_deg this function creates an IQIndex leg,
//where the leg accommodates rep of spin_qn with deg equals to flavor_deg[spin_qn] 
//qn_order denotes the way to store IndexQN in IQIndex
//qn_order=-1 means from high sz to low sz
//while qn_order=1 means from low sz to high sz
IQIndex Spin_leg(const std::vector<int> &flavor_deg, const std::string &iqind_name, Arrow dir=Out, IndexType it=Link, int qn_order=-1);

//Given an IQIndex (with sz quantum number), this function determines the corresponding SU(2) rep
//flavor_deg[n] stores extra deg for spin n/2
bool iqind_spin_rep(const IQIndex &sz_leg, std::vector<int> &flavor_deg);

//Given an IQIndex sz_leg (with sz quantum number) which accommodates SU(2) rep, this function relates integer i (0<=i<sz_leg.m()) and the spin basis (|S,m,t\rangle)
bool iqind_to_spin_basis(const IQIndex &sz_leg, std::vector<Spin_Basis> &spin_basis);

bool iqind_to_spin_basis(const IQIndex &sz_leg, const std::vector<int> &flavor_deg, std::vector<Spin_Basis> &spin_basis);


//This function construct eta from mu
//the Out leg has prime while In leg has no prime
IQTensor eta_from_mu(double mu, const std::vector<int> &flavor_deg);
//construct eta by given leg
IQTensor eta_from_mu(double mu, IQIndex eta_leg);
//define the function for Index for convience
ITensor eta_from_mu(double mu, Index eta_leg);


//This function creates isomorphic leg, which share the same prime level
Index isomorphic_legs(const Index &old_leg, const std::string &new_leg_name);
IQIndex isomorphic_legs(const IQIndex &old_leg, const std::string &new_leg_name);


//Assign value of tensor TB to tensor TA without changing the index.
//Hilbert space of TB should be morphism to that of TA 
void tensor_assignment(ITensor &TA, const ITensor &TB);
//For IQTensor, arrows should also match
//IMPORTANT: to replace iqindex, we should also make sure the order of indexqn is the same.
void tensor_assignment(IQTensor &TA, const IQTensor &TB);

//Assign index IB to index IA, with different index id to avoid contraction
inline void index_assignment(Index &IA, const Index &IB, const std::string &IA_name="IndexA")
{
    IA=Index(IA_name,IB.m(),IB.type(),IB.primeLevel());
}
inline void index_assignment(IQIndex &IA, const IQIndex &IB, const std::string &IA_name="IndexA")
{
    std::vector<IndexQN> ind_qn(IB.indices());
    IA=IQIndex(IA_name,ind_qn,IB.dir(),IB.primeLevel());
}


//these two functions can solve the inconsistent interface to get left_[i] of Combiner and IQCombiner
inline Index left_leg_of_combiners(const Combiner &combiner, int i) { return combiner.left(i); }
inline IQIndex left_leg_of_combiners(const IQCombiner &combiner, int i) { return combiner.left()[i]; }


//these two functions solve inconsistent interface to get qn preserve block of iqtensor
inline void clean(ITensor &tensor) {}
inline void clean(IQTensor &tensor) {return tensor.clean(); }


inline ITensor toITensor(const ITensor &tensor) { return tensor; }
inline ITensor toITensor(const IQTensor &tensor) { return tensor.toITensor(); }


//given a permuted indices, we obtain a permuted tensor such that
//(tensor_permute)_{\alpha\beta...}=(tensor)_{P(\alpha)P(\beta)...}
template <class TensorT>
TensorT tensor_permutation(const std::vector<int> &permuted_indices, const TensorT &tensor_origin);

//given two tensors with the same indices but different ordering, assign TB to TA without changing the order of TA
template <class TensorT>
void tensor_assignment_diff_order(TensorT &TA, const TensorT &TB);


//Multiply the corresponding eta of mu on one leg, and leave the indices order invariant
template <class TensorT>
inline void obtain_tensor_after_eta_action(double mu, TensorT &tens, const typename TensorT::IndexT &act_leg)
{
    auto eta=eta_from_mu(mu,dag(act_leg));
    TensorT eta_tens_unordered=eta*tens;
    eta_tens_unordered.replaceIndex(prime(act_leg),act_leg);
    tensor_assignment_diff_order(tens,eta_tens_unordered);
}

template <class TensorT>
inline TensorT tensor_after_eta_action(double mu, TensorT tens, const typename TensorT::IndexT &act_leg)
{
    obtain_tensor_after_eta_action(mu,tens,act_leg);
    return tens;
}


//contract tensor by identifying given indices
template <class TensorT, class IndexT>
inline TensorT tensor_contraction(TensorT TA, TensorT TB, const std::vector<IndexT> &TA_inds, const std::vector<IndexT> &TB_inds)
{
    for (int indi=0; indi<TA_inds.size(); indi++)
    {
        TB.replaceIndex(TB_inds[indi],dag(TA_inds[indi]));
    }
    return TA*TB;
}

//contract multiple tensors in order
//template <class TensorT, class IndexT>
//inline TensorT tensor_contracion(std::list<TensorT> Tens, std::list<std::vector<IndexT>> Tens_inds)
//{
//    if (Tens.size()<2) 
//    {
//        cout << "Tensor contraction requires at least two tensors!" << endl;
//        exit(EXIT_FAILURE);
//    }
//    TensorT temp_tensor=tensor_contracion(Tens[0],Tens[1],Tens_inds[0],Tens_inds[1]);
//    if (Tens.size()==2) 
//    {
//        return temp_tensor;
//    }
//    else
//    {
//        Tens.pop_front();
//        Tens.pop_front();
//        Tens.push_front(temp_tensor);
//        Tens_inds.pop_front();
//        Tens_inds.pop_front();
//        return tensor_contraction(Tens,Tens_inds);
//    }
//}

//Translate Combiner/IQCombiner to ITensor/IQTensor 
inline void combiner_to_tensor(Combiner &combiner, ITensor &tensor)
{
    combiner.init();
    tensor=combiner.toITensor();
}
inline void combiner_to_tensor(IQCombiner &iqcombiner, IQTensor &iqtensor)
{
    iqcombiner.init();
    iqtensor=iqcombiner.toIQTensor();
}

//count number of legs of tensors after contraction
//we exclude the case that more than two tensors share the same leg
template <class TensorT>
inline int legs_num_after_contraction(const std::vector<TensorT> &curr_tensors)
{
    int legs_num=0;
    std::vector<typename TensorT::IndexT> uncontract_indices;

    for (const auto &tens : curr_tensors)
    {
        for (const auto &ind : tens.indices())
        {
            auto ind_iter=std::find(uncontract_indices.begin(),uncontract_indices.end(),ind);
            if (ind_iter==uncontract_indices.end())
            {
                legs_num++;
                uncontract_indices.push_back(ind);
            }
            else
            {
                legs_num--;
            }
        }
    }

    return legs_num;
}


//transfer an rank 2 ITensor to armadillo matrix. T can be either real or cplx
arma::Mat<Complex> arma_mat_from_rank2_itensor(const ITensor &tensor);
arma::Mat<Complex> arma_mat_from_rank2_itensor(const ITensor &tensor, Index leg1, Index leg2);

//transfer an armadillo matrix to a rank 2 ITensor
ITensor rank2_itensor_from_arma_mat(arma::Mat<double> matrix, const Index &leg1, const Index &leg2);
ITensor rank2_itensor_from_arma_mat(arma::Mat<Complex> matrix, const Index &leg1, const Index &leg2);

//get inverse tensor of rank 2 tensor by transform to arma_mat
ITensor inverse_rank2_tensor_by_arma_mat(const ITensor &tensor);
//for IQTensor, we get the inverse of each block separately, and then combine to the total inverse tensor
IQTensor inverse_rank2_tensor_by_arma_mat(const IQTensor &tensor);


//
// The "factor" decomposition is based on the SVD, but factorizes a tensor T into only two tensors T=A*B where A and B share a single common index.
//
// If the SVD of T is T=U*S*V where S is a diagonal matrix of singular values, then A and B are schematically A=U*sqrt(S) and B=sqrt(S)*V.
//
// In addition to the named Args recognized by the svd routine, factor accepts an Arg "IndexName" which wil be the name of the common index connecting A and B.
template <typename Tensor>
void factor(Tensor const& T, Tensor &A, Tensor &B, Args const &args = Args::global());

#endif
