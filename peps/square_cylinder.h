
#ifndef _SQUARE_CYLINDER_H_
#define _SQUARE_CYLINDER_H_

#include "double_layer_peps.h"
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
        using Parent=Double_Layer_PEPSt<TensorT>; 
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;

        //
        //Constructor
        //
        //Cylinder_Square_Double_Layer_PEPSt() {}
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

        //snake walking for one bulk or boundary col 
        //vertical_dir==1(-1) denotes walking from left(right) to right(left),
        //horizontal_dir==1(-1) denotes walking from down(up) to up(down)
        void snake_walking_boundary_col();
        void snake_walking_bulk_col(int coli, int horizontal_dir, int vertical_dir, bool do_normalize=false);
        //make indices of (meeting) sigma_lr_ match each other, so their multiplication gives wavefunction norm
        //keep sigma_lr_[1-lr_no] invariant, while change sigma_lr_[lr_no]
        void match_sigma_left_right(int lr_no=1);
        //sigma_b=\sqrt(\sigma_l).\sigma_r.sqrt(sigma_l)
        void from_sigma_lr_to_sigma_b();
        //get spectrum of density matrix, which is identical to that of sigma_b (normalized)
        void obtain_density_matrix_spectrum();


        //Methods relate to entanglement properties
        std::vector<double> entanglement_spectrum();
        double entanglement_entropy_vN();
        double entanglement_entropy_Renyi(double renyi_n);
        //void obtain_transfer_matrix(int coli=1);

        //
        //calculate correlators on square cylinder notice square_cylinder_double_peps should provide converged sigma_lr
        //this function only applies to those operators can be written as direct product form. Each operator is a Tensor formed by two Site indices (noprimed and primed)
        double obtain_correlators(const std::vector<int> &acting_sites_list, std::vector<TensorT> direct_prod_operators);
        //get value of <O>=\langle\psi|\hat{O}|\psi\rangle iteratively //where <O> is obtained by contracting sandwiched_tensors
        //we will start from col_lr, where data outside col_lr are stored in sigma_lr(decombined)
        //horizontal_dir is the contraction direction
        double sandwiched_peps_norm(const std::vector<TensorT> &sandwiched_tensors, int horizontal_dir=1);

        
        //
        //Other Methods
        //
        //decombine/recombine sigma_lr_ from/to matrix, where the decombined sigma_lr_ shares the same indices as double_layer_tensors_
        void decombine_sigma_lr(int lr_no);
        void decombine_sigma_lr();
        void recombine_sigma_lr();
        //match indices of current decombined sigma_lr_ and double_layer_tensors_
        void match_sigma_lr_double_peps();
        //move the converged sigma_lr_ to specific cols by replacing indices
        void move_sigma_lr(const std::array<int,2> &new_col_lr);


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
        //upper_combiners_ and lower_combiners_ are used to combine sigma_lr_ from tensor to matrix
        //upper_combiners_=prime(dag(lower_combiners_))
        std::array<std::vector<CombinerT>,2> lower_combiners_, upper_combiners_; ;

        //eigval of sigma_b_ is identical to those of RDM
        //we always choose sigma_b_ as ITensor instead
        ITensor sigma_b_;
        std::vector<double> density_mat_spectrum_;

        //TensorT transfer_mat_;
};



#endif
