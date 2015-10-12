
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
        Double_Layer_PEPSt() {}
        Double_Layer_PEPSt(const PEPSt<TensorT> &peps);
        //Construct double layer peps on square cylinder using a single site tensor where the boundary tensor is direct product of single indice tensors
        //the virtual leg of site tensor we input is in the order as left, up, right, down
        //Double_Layer_PEPSt(const Lattice_Base &square_lattice_cylinder, const TensorT &square_site_tensor, const TensorT &boundary_tensor);

        //
        //Access Method
        //
        int DD() const
        {
            return DD_;
        }

        const Lattice_Base &lattice() const
        {
            return lattice_;
        }

        const TensorT &layered_site_tensors(int sitei) const
        {
            return layered_site_tensors_[sitei];
        }
        const std::vector<TensorT> &layered_site_tensors() const
        {
            return layered_site_tensors_;
        }

        const CombinerT &virt_leg_combiners(int sitei, int j) const
        {
            return virt_leg_combiners_[sitei][j];
        }
        const std::vector<CombinerT> &virt_leg_combiners(int sitei) const
        { 
            return virt_leg_combiners_[sitei]; 
        }
        const std::vector<std::vector<CombinerT>> &virt_leg_combiners() const
        {
            return virt_leg_combiners_;
        }

        //
        //Constructor Helpers
        //
        //Absorb bond tensors and boundary tensors to site tensors
        void obtain_combined_site_tensors(const PEPSt<TensorT> &peps, std::vector<TensorT> &combined_site_tensors);
        //from lower_tensors to layered_tensors with all pairs of virtual legs combined
        void obtain_layered_tensors_with_combined_legs(const std::vector<TensorT> &lower_tensors);


    private:
        int DD_;

        const Lattice_Base &lattice_;

        std::vector<TensorT> layered_site_tensors_;

        //Notice there is no order in virt_leg_combiners_[sitei]. To combine/decombine a special leg, we should use hasindex() to select the particular combiner we want
        std::vector<std::vector<CombinerT>> virt_leg_combiners_;
        
};
using Double_Layer_PEPS=Double_Layer_PEPSt<ITensor>;
using Double_Layer_IQPEPS=Double_Layer_PEPSt<IQTensor>;

#endif
