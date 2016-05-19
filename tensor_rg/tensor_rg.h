
#ifndef _TENSOR_RG_H_
#define _TENSOR_RG_H_

#include "utilities.h"
#include "lattice.h"

//
//class TensorT_RG
//we do not assume translation invariant here
//we first transform the tensor network to a square shape
//and then do trg for this square
//
template <class TensorT>
class TensorT_RG
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
        TensorT_RG(const Lattice_Base &lattice): lattice_(lattice) {}
        TensorT_RG(const Lattice_Base &lattice, const std::vector<TensorT> &input_tensors, int maxm=100);

        //
        //Acess Methods
        //
        const Lattice_Base &lattice() const { return lattice_; }
        Complex trg_result() const { return trg_result_; }
        bool is_zero() const { return iszero_; }

        //
        //Methods to obtain/update layered trg_tensors and factor_tensors
        //
        //factor the input tensor, used when transform to square geometry for some lattices (e.g. kagome normal)
        void obtain_factor_input_tensor(int input_tensor_no);
        void obtain_zeroth_layer_trg_tensor(int tensor_no);
        //obtain trg tensor from factor tensors of last layer, should consider even and odd layer separately
        //the order of tensors to be cut:
        //   |    |            \ /
        // --B----C--           B
        //   |    |          \ / \ / 
        //   |    |           A   C
        // --A----D--        / \ / \
        //   |    |             D
        //                     / \ 
        //   
        void obtain_trg_tensor(int layer_no, int tensor_no);
        //
        //                |                    |
        //    |         --B                    B--
        //  --T-- ==>      \   (odd site),    /    (even site)
        //    |             A--            --A
        //                  |                |
        //
        //               \ /
        //    \ /         B                   \    /
        //     T   =>     |  (even rowy),      A--B  (odd rowy)
        //    / \         A                   /    \
        //               / \
        //
        void obtain_factor_tensors(int layer_no, int tensor_no);
        void obtain_trg_result();

        //update all layers of the tensor rg with some updated_tensors
        void update_trg_network(const std::vector<int> &input_inds, const std::vector<TensorT> &update_input_tensors);
        
        //get update inds of layer_no=0
        void obtain_zeroth_layer_update_inds(std::vector<int> &update_inds, const std::vector<int> &input_inds);
        //get update inds of layer_no>=1 from that of last layer
        void obtain_update_inds(std::vector<int> &update_inds_curr, const std::vector<int> &update_inds_last, int layer_no);

        //
        //Other Methods
        //
        //translate between odd layer tensor ind and plaq ind
        inline int plaq_ind_from_odd_layer_tensor_ind(int tensor_ind, int layer_no) const { return tensor_ind*2+(tensor_ind*2/layered_lattice_dim_[layer_no][0])%2; }

        inline std::vector<int> coord_from_ind(int ind, std::vector<int> lattice_dim) const { return {ind%lattice_dim[0],ind/lattice_dim[0]}; }
        inline int ind_from_coord(std::vector<int> coord, std::vector<int> lattice_dim) const 
        {
            coord[0]=(coord[0]+lattice_dim[0])%lattice_dim[0];
            coord[1]=(coord[1]+lattice_dim[1])%lattice_dim[1];
            return (coord[0]+coord[1]*lattice_dim[0]);
        }
        //obtain neighbouring tensor ind
        //neigh     =   0    1    2    3
        //even layer: left   up  right down
        //odd layer:    ld   lu  ru   rd
        inline int neigh_tensor_ind(int tensor_ind, int neigh, int layer_no) const
        {
            if (layer_no%2==0)
            {
                std::vector<int> neigh_coord=coord_from_ind(tensor_ind,layered_lattice_dim_[layer_no]);
                neigh_coord[0]+=(neigh-1)%2;
                neigh_coord[1]+=-(neigh-2)%2;
                return ind_from_coord(neigh_coord,layered_lattice_dim_[layer_no]);
            }
            else
            {
                std::vector<int> neigh_plaq_coord=coord_from_ind(plaq_ind_from_odd_layer_tensor_ind(tensor_ind,layer_no),layered_lattice_dim_[layer_no]);
                neigh_plaq_coord[0]+=(neigh<=1) ? -1:1;
                neigh_plaq_coord[1]+=(neigh==0||neigh==3) ? -1:1;
                return ind_from_coord(neigh_plaq_coord,layered_lattice_dim_[layer_no])/2;
            }
        }
        //obtain site inds around a plaquette
        //order: ld, lu, ru, rd
        //plad_ind is the same as its ld site ind
        inline int site_ind_surr_plaq(int plaq_ind, int surri, int layer_no) const
        {
            if (surri==0) return plaq_ind;
            const auto &lattice_dim=layered_lattice_dim_[layer_no];
            std::vector<int> plaq_coord=coord_from_ind(plaq_ind,lattice_dim);
            if (surri==1) return ind_from_coord({plaq_coord[0],plaq_coord[1]+1},lattice_dim);
            if (surri==2) return ind_from_coord({plaq_coord[0]+1,plaq_coord[1]+1},lattice_dim);
            if (surri==3) return ind_from_coord({plaq_coord[0]+1,plaq_coord[1]},lattice_dim);
        }
        //obtain neigh plaq inds around a center plaq
        //order: lurd
        inline int plaq_ind_surr_plaq(int plaq_ind, int surri, int layer_no) const
        {
            const auto &lattice_dim=layered_lattice_dim_[layer_no];
            std::vector<int> surr_plaq_coord=coord_from_ind(plaq_ind,lattice_dim);
            surr_plaq_coord[0]+=(surri-1)%2;
            surr_plaq_coord[1]+=-(surri-2)%2;
            return ind_from_coord(surr_plaq_coord,lattice_dim);
        }
        //obtain neigh plaq inds around a site
        //order: ld,lu,ru,rd
        inline int plaq_ind_surr_site(int site_ind, int surri, int layer_no) const
        {
            const auto &lattice_dim=layered_lattice_dim_[layer_no];
            std::vector<int> site_coord=coord_from_ind(site_ind,lattice_dim);
            if (surri==0) return ind_from_coord({site_coord[0]-1,site_coord[1]-1},lattice_dim);
            if (surri==1) return ind_from_coord({site_coord[0]-1,site_coord[1]},lattice_dim);
            if (surri==2) return site_ind;
            if (surri==3) return ind_from_coord({site_coord[0],site_coord[1]-1},lattice_dim);
        }

        //
        //Construct Helpers
        //
        //transform to square shape and init layered_lattice_dim_ as well as 0th layer of layered_trg_tensors_
        void init_to_square_network();


    private:
        //1. square lattice: just transform directly
        //2. normal kagome lattice: we first transform to a honeycomb lattice by SVD, and then to square lattice
        //        1             |
        //       /u\            a
        //  2---0---2   ==>  \ / \  ==> square
        //   \d/              b
        //    1               |
        //                               \ /
        //    /            /     \ /      B      \         \
        // --0-- ==> --B--A-- ,   1  ==>  |  ,  --2-- ==> --A--B-- 
        //  /         /          / \      A        \            \
        //                               / \ 
        //
        Lattice_Base lattice_;
        //steps of tensor rg
        int N_layer_;
        std::vector<TensorT> input_tensors_;
        //factor of input tensors, used for kagome normal case
        std::vector<std::vector<TensorT>> factor_input_tensors_;
        //we set the convention as follows: 
        //layered_lattice_dim_[2k+1]=layered_lattice_dim_[2k]
        //layered_lattice_dim_[2k+2]=layered_lattice_dim_[2k+1]/2;
        //for even layer, tensors live at site;
        //for odd layer, tensor i live at plaquette center i*2+i/Lx%2
        std::vector<std::vector<int>> layered_lattice_dim_;
        //layered tensors stores tensors obtained from intermediate steps of tensor rg
        //layered_trg_tensors_[n].size()=layered_trg_tensors_[n-1].size()/2;
        //layered_trg_tensors_[0] stores tensors of original lattice
        //for square lattice, the last layer stores four tensors
        std::vector<std::vector<TensorT>> layered_trg_tensors_;
        //layered_factor_tensors_ stores tensors from SVD decomposition 
        //the range is from 0 to N_layer_-1
        //order: down,left,right,up
        std::vector<std::vector<std::vector<TensorT>>> layered_factor_tensors_;
        Complex trg_result_;
        bool iszero_;

        //parameters for factor tensor
        Args factor_args_;
};

#endif
