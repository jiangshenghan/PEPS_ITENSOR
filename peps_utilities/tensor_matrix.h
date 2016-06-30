

#ifndef _TENSOR_MATRIX_H_
#define _TENSOR_MATRIX_H_

#include "utilities.h"

//TODO: debug this file
//
//define tensor_matrix_arnoldi class with size and product, used in arnoldi method 
//
//we get a set of tensors
//---output_indices-contract_tensors-input_indices--- . ---phi
//the multiplications sequence follows the contract_seq_ 
//contract_seq_[i]==-1 denotes multiply phi, and contract_seq_[i]==n denotes multiply contract_tensors_[n]
//To prevent exceeds of leg numbers, we use combiners to combine/decombine legs
//leg_combiners_ stores the combiners. Notice, we always put empty combiners in leg_combiners_[0]
//the multiplication sequence to combiners are stored in combiner_seq_
//combiner_seq_[i] is empty means doing nothing, while combiner_seq_[i]=={n1,-n2} means multiply leg_combiners_[n1] and dag(leg_combiners[n2]) to phip. 
//TODO: In general, the contract_seq should be a tree instead of a vector
//
template <class TensorT>
class TensorT_Matrix_Arnoldi
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
        TensorT_Matrix_Arnoldi(const std::vector<IndexT> &input_indices, const std::vector<IndexT> &output_indices, const std::vector<TensorT> &contract_tensors, const std::vector<int> &contract_seq=std::vector<int>(), const std::vector<CombinerT> &leg_combiners=std::vector<CombinerT>(), std::vector<std::vector<int>> combiner_seq=std::vector<std::vector<int>>()):
            input_indices_(input_indices),
            output_indices_(output_indices),
            leg_combiners_(leg_combiners),
            contract_tensors_(contract_tensors),
            contract_seq_(contract_seq),
            combiner_seq_(combiner_seq)
        {
            size_=1;
            for (const auto &ind: input_indices_) size_*=ind.m();

            //the default contract seq is to multiply phi first and then obtain sparse matrix 
            if (contract_seq_.empty())
            {
                contract_seq_.push_back(-1);
                for (int i=0; i<contract_tensors_.size(); i++) contract_seq_.push_back(i);
            }

            if (leg_combiners_.empty() || combiner_seq_.empty()) combiner_seq_=std::vector<std::vector<int>>(contract_tensors_.size()+1,std::vector<int>());
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
                    temp_tensor=phi;
                else
                    temp_tensor=contract_tensors_[contract_seq_[i]];

                if (i==0) 
                    phip=temp_tensor;
                else
                    phip*=temp_tensor;

                //Print(i);
                //Print(phip.indices());
                for (int combineri: combiner_seq_[i])
                {
                    //Print(leg_combiners_[std::abs(combineri)]);
                    if (combineri>0) 
                        phip=phip*leg_combiners_[combineri];
                    else if (combineri<0)
                        phip=phip*dag(leg_combiners_[-combineri]);
                }

                //Print(i);
                //Print(temp_tensor.indices());
                //Print(phip.indices());
                //Print(temp_tensor.norm());
                //Print(phip.norm());
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
        //matrix size
        int size_;
        //output_indices and input_indices should have one-to-one correpondance
        std::vector<IndexT> input_indices_, output_indices_;
        std::vector<CombinerT> leg_combiners_;
        std::vector<TensorT> contract_tensors_;
        std::vector<int> contract_seq_;
        std::vector<std::vector<int>>combiner_seq_;
};

#endif
