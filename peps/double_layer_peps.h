
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
        void obtain_single_layer_tensors();
        //from lower_tensors to layered_tensors with all pairs of virtual legs combined
        void obtain_layered_tensors_with_combined_legs();

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


//
//class for double layered square peps on cylinder geometry
//can be used to obtain entanglement property and measure correlation functions
//
template<class TensorT>
class Cylinder_Square_Double_Layer_PEPSt : public Double_Layer_PEPSt<TensorT>
{
    public:
        //
        //type alias
        //
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;

        //
        //Constructor
        //
        Cylinder_Square_Double_Layer_PEPSt() {}
        Cylinder_Square_Double_Layer_PEPSt(const Lattice_Base &square_cylinder);
        Cylinder_Square_Double_Layer_PEPSt(const PEPSt<TensorT> &square_peps, int cutting_col=-1);

        //
        //Acess Method
        //
        int col_lr(int i) const { return col_lr_[i]; }
        const std::array<TensorT,2> &sigma_lr() const { return sigma_lr_; }
        const TensorT &sigma_lr(int i) const { return sigma_lr_[i]; }

        const std::vector<double> &density_mat_spectrum() const { return density_mat_spectrum_; }

        //
        //Method to obtain sigma
        //
        //using iterative method to obtain boundary theories
        //when do_decombine==true, we will expand indices of sigma_lr_
        void obtain_boundary_theory_iterative();
        void obtain_sigma_lr_iterative(int left_end_col, int right_end_col);

        //decombine/recombine sigma_lr_, where the decombined sigma_lr_ shares the same indices as double_layer_tensors_
        void decombine_sigma_lr();
        void recombine_sigma_lr();
        //match indices of current decombined sigma_lr_ and double_layer_tensors_
        void match_indices_sigma_lr();

        //snake walking for one bulk or boundary col 
        //vertical_dir==1(-1) denotes walking from left(right) to right(left),
        //horizontal_dir==1(-1) denotes walking from down(up) to up(down)
        void snake_walking_boundary_col();
        void snake_walking_bulk_col(int coli, int horizontal_dir, int vertical_dir);
        //turn a one indice tensor (vector) to a two indices (upper leg and lower leg) matrix
        //we need to decombine the index, and recombine the upper ones and lower ones separetely
        void from_sigma_vec_to_mat();
        //sigma_b=\sqrt(\sigma_l).\sigma_r.sqrt(sigma_l)
        void from_sigma_lr_to_sigma_b();
        //get spectrum of density matrix, which is identical to that of sigma_b (normalized)
        void obtain_density_matrix_spectrum();
        
        //Other Methods
        std::vector<double> entanglement_spectrum();
        double entanglement_entropy_vN();
        double entanglement_entropy_Renyi(double renyi_n);

        //void obtain_transfer_matrix(int coli=1);

        //print first several nonzero elems of vector
        void print_vector_nonzero_elem(TensorT vector, int nonzero_elem_num)
        {
            assert(vector.r()==1);

            Print(vector);
            auto vec_ind=vector.indices()[0];
            
            int val=1;
            int nonzeroi=0;
            while (nonzeroi<nonzero_elem_num && val<=vec_ind.m())
            {
                auto elem=vector(vec_ind(val));
                if (elem>EPSILON)
                {
                    cout << "(" << val << "): " << elem << endl;
                    nonzeroi++;
                }
                val++;
            }
            cout << endl << endl;
        }

        //Method to read/write from/to file
        //Before writing simga_lr_ to file, we should decombine their index. (vector -> tensor)
        //After reading sigma_lr_, we should to replace their index since virt_leg_combiners_ have been reconstructed.
        //We should also reconstruct iterative_combiners_
        void read(std::istream &s);
        void write(std::ostream &s) const;

    private:
        //cutting_col_ is used to obtain boundary theory
        int cutting_col_; 
        //col_lr_stores col no. for current sigma_lr
        std::array<int,2> col_lr_;
        //sigma_lr_[0(1)] stores sigma_left(right)
        std::array<TensorT,2> sigma_lr_;
        //iterative_combiners_ is used to contract cols to a big tensor
        std::array<std::vector<CombinerT>,2> iterative_combiners_;

        //eigval of sigma_b_ is identical to those of RDM
        TensorT sigma_b_;
        std::vector<double> density_mat_spectrum_;

        //TensorT transfer_mat_;
};


#endif
