
//
//define tensor_matrix class with size and product, used in arnoldi method etc
//
//we get a set of tensors
//---output_indices-contract_tensors-input_indices--- . ---phi
//the multiplications sequence follows the contract_seq_ 
//where contract_seq_[i]=-1 denotes multiply phi 
//To prevent exceeds of leg numbers, we use combiners to combine/decombine legs
//TODO: the contract_seq should be a tree form instead of a vector
//

#ifndef _TENSOR_MATRIX_H_
#define _TENSOR_MATRIX_H_

#include "utilities.h"

template <class TensorT>
class TensorT_Matrix
{
    public:
        //
        //type alias
        //
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;

        //
        //constructor
        //
        TensorT_Matrix(const std::vector<IndexT> &input_indices, const std::vector<IndexT> &output_indices, const std::vector<TensorT> &contract_tensors, const std::vector<int> &contract_seq=std::vector<int>(), const std::vector<CombinerT> &input_combiners=std::vector<CombinerT>(), const std::vector<CombinerT> &output_combiners=std::vector<CombinerT>(), const std::vector<int> &mult_input_combiners=std::vector<int>(), const std::vector<int> &mult_output_combiners=std::vector<int>()):
            input_indices_(input_indices),
            output_indices_(output_indices),
            input_combiners_(input_combiners),
            contract_tensors_(contract_tensors),
            contract_seq_(contract_seq),
            mult_input_combiners_(mult_input_combiners),
            mult_output_combiners_(mult_output_combiners_)
        {
            size_=1;
            for (const auto &ind: input_indices_) size_*=ind.m();

            //the default contract seq is to obtain sparse matrix first and then multiply phi at last
            if (contract_seq_.empty())
            {
                for (int i=0; i<contract_tensors.size(); i++) contract_seq_.push_back(i);
                contract_seq_.push_back(-1);
            }

            if (input_combiners_.empty() || mult_input_combiners_.empty())
                mult_input_combiners_=std::vector<int>(contract_tensors.size()+1,-1);
            if (output_combiners_.empty() || mult_output_combiners_.empty())
                mult_output_combiners_=std::vector<int>(contract_tensors.size()+1,-1);
        }

        //
        //Sparse matrix methods
        //
        int size() const { return size_; }

        void product(const TensorT &phi, TensorT &phip) const
        {
            for (int i=0; i<contract_seq_.size(); i++)
            {
                TensorT temp_tensor;
                if (contract_seq_[i]==-1) 
                {
                    temp_tensor=phi;
                    if (mult_input_combiners_[i]>=0) temp_tensor=temp_tensor*input_combiners_[mult_input_combiners_[i]];
                }
                else
                {
                    temp_tensor=contract_tensors_[contract_seq_[i]];
                    if (mult_input_combiners_[i]>=0) phip=phip*input_combiners_[mult_input_combiners_[i]]; //phi should contain in phip
                }

                if (i==0) phip=temp_tensor;
                else phip*=temp_tensor;

                if (mult_output_combiners_[i]>=0) phip=phip*output_combiners_[mult_output_combiners_[i]];
            }

            //TODO: consider the case where input_dir==output_dir
            for (int indi=0; indi<input_indices_.size(); indi++)
            {
                if (input_indices_[indi].dir()==output_indices_[indi].dir()) 
                {
                    cout << "Input indices and output indices of sparse matrix have conflict directions!" << endl;
                    exit(1);
                }
                phip.replaceIndex(output_indices_[indi],dag(input_indices_[indi]));
            }
        }

        //TODO: implement diag() method

    private:
        int size_;
        //output_indices and input_indices should have one-to-one correpondance
        std::vector<IndexT> input_indices_, output_indices_;
        std::vector<CombinerT> input_combiners_, output_combiners_;
        std::vector<TensorT> contract_tensors_;
        std::vector<int> contract_seq_, mult_input_combiners_, mult_output_combiners_;
};

#endif
