
#ifndef _DOUBLE_LAYER_PEPS_
#define _DOUBLE_LAYER_PEPS_

#include "peps.h"

//
//Double_Layer_PEPSt is obtained by contracting physical leg of PEPS
//We only stores site tensors while bond tensors and boundary tensors are contracted to sites
//We also combine virtual legs of up layer and down layer, which can be decomposed by combiners
//
template <class TensorT>
class Double_Layer_PEPSt
{
    public:
        //
        //type alias
        //
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;

        //
        //Constructors
        //
        //Double_Layer_PEPSt() {}

        Double_Layer_PEPS () {}
        Double_Layer_PEPSt(const Lattice_Base &lattice);  

        Double_Layer_PEPSt(const PEPSt<TensorT> &peps);

        //
        //Access Method
        //
        const Lattice_Base &lattice() const
        {
            return single_layer_peps_.lattice();
        }

        const TensorT &double_layer_tensors(int sitei) const
        {
            return double_layer_tensors_[sitei];
        }
        const std::vector<TensorT> &double_layer_tensors() const
        {
            return double_layer_tensors_;
        }

        const CombinerT &virt_leg_combiners(int sitei, int j) const
        {
            return virt_leg_combiners_[sitei][j];
        }
        //return combiner at sitei with right index equal to double_virt_ind
        const CombinerT &virt_leg_combiners(int sitei, const IndexT &double_virt_ind) const
        {
            for (const auto &combiner : virt_leg_combiners_[sitei])
            {
                if (combiner.right()==double_virt_ind) return combiner;
            }
            cout << "Fail to find combiner!" << endl;
            exit(EXIT_FAILURE);
        }
        //return combiner shared by double layer tensors at sitei and sitej
        //if they share more than one combiners, return arbitrary one
        const CombinerT &virt_leg_combiners(const std::array<int,2> &sites_no) const
        {
            for (const auto &combineri : virt_leg_combiners_[sites_no[0]])
            {
                for (const auto &combinerj : virt_leg_combiners_[sites_no[1]])
                {
                    if (combineri.right()==combinerj.right()) return combineri;
                }
            }
            cout << "Fail to find combiner!" << endl;
            exit(EXIT_FAILURE);
        }
        const std::vector<CombinerT> &virt_leg_combiners(int sitei) const
        { 
            return virt_leg_combiners_[sitei]; 
        }
        const std::vector<std::vector<CombinerT>> &virt_leg_combiners() const
        {
            return virt_leg_combiners_;
        }

        const IndexT &phys_leg(int sitei) const
        {
            return single_layer_peps_.phys_legs(sitei);
        }

        //
        //Other Methods
        //
        //get the original lower_indice of virt_ind of sitei. upper_indic=dag(prime(lower_indice))
        CombinerT decombine_double_virt_indice(int sitei, const IndexT &double_virt_ind, IndexT &lower_ind);
        CombinerT decombine_double_virt_indice(std::array<int,2> sites_no, IndexT &lower_ind);
        //Sandwich operators between double layer PEPS
        //Insert direct product operators: the operator is formed by two indices, one with prime and one without prime
        void obtain_peps_sandwich_single_site_operators(std::vector<TensorT> direct_prod_operators, const std::vector<int> &acting_sites_list, std::vector<TensorT> &sandwiched_tensors);
        //Insert tensor product operators
        void obtain_peps_sandwich_tensor_prod_operators(std::vector<TPOt<TensorT>> tensor_prod_operators, const std::vector<std::vector<int>> &acting_sites_list, std::vector<TensorT> &sandwiched_tensors);

        //
        //Constructor Helpers
        //
        //Absorb bond tensors and boundary tensors to site tensors
        void obtain_single_layer_tensors();
        //from lower_tensors to layered_tensors with all pairs of virtual legs combined
        void obtain_layered_tensors_with_combined_legs();
        //when combine to double layered tensors, the number of legs may exceed NMAX, so we need to first combine virtual legs in one layer and then decombine and recombine in pairs
        //we do not require phys leg to contract
        TensorT double_layer_tensor_from_lower_upper_tensors(TensorT lower_tensor, TensorT upper_tensor, const std::vector<CombinerT> &pair_combiners);

        //read/write from/to file
        //Notice we need to reconstruct double_layer_tensors_ and virt_leg_combiners_ after reading from file
        void read(std::istream &s);
        void write(std::ostream &s) const;

    protected:
        PEPSt<TensorT> single_layer_peps_;
        //single layer tensors absorb bond tensors and boundary tensors to nearby site tensors
        std::vector<TensorT> single_layer_tensors_;
        //double_layer_tensors obtained from contracting physical legs of single_layer_tensors_ and combine pairs of virtual legs
        std::vector<TensorT> double_layer_tensors_;
        //Notice there is no order in virt_leg_combiners_[sitei]. To combine/decombine a special leg, we should use hasindex() to select the particular combiner we want
        std::vector<std::vector<CombinerT>> virt_leg_combiners_;

        //bool double_layer_initted_;
        
};
using Double_Layer_PEPS=Double_Layer_PEPSt<ITensor>;
using Double_Layer_IQPEPS=Double_Layer_PEPSt<IQTensor>;


#endif
