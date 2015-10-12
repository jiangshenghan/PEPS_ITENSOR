
#ifndef _BOUNDARY_THEORY_
#define _BOUNDARY_THEORY_

#include "double_layer_peps.h"


//we obtain the boundary theory \sigma_b=\sqrt{\sigma_L}\cdot\sigma_R\cdot\sqrt{\sigma_L}
//where \sigma_{L/R} is obtained via contraction the left/right part of double layer PEPS on cylinder
//and \sqrt{\sigma_L}=U_L\cdot\sqrt(\lambda_L)\cdot V_L, where \lambda_L is the singular value of \sigma_L
//Here, we focus on square lattice on cylinder, other lattice should be transfer to square first
template<class TensorT>
class Boundary_Theory
{
    public:
        //
        //Type Alias
        //
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;

        //
        //Constructor
        //
        Boundary_Theory() {}
        Boundary_Theory(const PEPSt<TensorT> &square_peps, int cutting_col=-1, const std::string &method="iterative");

        //
        //Acess Method
        //
        const std::vector<double> &density_mat_spectrum() const { return density_mat_spectrum_; }

        //
        //Method to obtain \sigma_b
        //
        //using iterative method to obtain boundary theories
        void obtain_boundary_theory_iterative();
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

        void obtain_transfer_matrix(int coli=1);
        //void obtain_sigma_lr_from_transfer_matrix();
        //void apply_transfer_matrix_to_sigma_lr_vec(int lr_no);

    private:
        int n_rows_, n_cols_, cutting_col_;

        Double_Layer_PEPSt<TensorT> square_layered_peps_;

        //sigma_lr_[0(1)] stores sigma_left(right)
        std::array<TensorT,2> sigma_lr_;
        //iterative_combiners_ is used to contract cols to a big tensor
        std::array<std::vector<CombinerT>,2> iterative_combiners_;

        TensorT sigma_b_;
        std::vector<double> density_mat_spectrum_;

        TensorT transfer_mat_;
};


#endif
