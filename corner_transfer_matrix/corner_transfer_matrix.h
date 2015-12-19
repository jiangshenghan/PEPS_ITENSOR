
#ifndef _CORNER_TRANSFER_MATRIX_H_
#define _CORNER_TRANSFER_MATRIX_H_

//#include "utilities.h"
#include "TPO.h"

//
//class to obtain and store effective environment for iPEPS using CTM
//the algorithm follows arXiv:1402.2859
//
template <class TensorT>
class Corner_Transfer_Matrix
{
    public:

        //
        //type alias
        //
        using IndexT=typename TensorT::IndexT;
        using IndexValT=typename TensorT::IndexValT;
        using CombinerT=typename TensorT::CombinerT;
        enum CTM_Direction {Left=0, Up=1, Right=2, Down=3};

        //
        //Constructor
        //
        Corner_Transfer_Matrix() {}
        //construct class using building blocks of ipeps
        //ordered_virt_indices is virtual indices for each tensors, in Left/Up/Right/Down order
        Corner_Transfer_Matrix(const std::vector<TensorT> &single_layer_tensors, const std::vector<std::array<IndexT,4>> &ordered_virt_indices, int Lx=1, int Ly=1);

        //
        //Access Method
        //
        int Lx() const { return Lx_; }
        int Ly() const { return Ly_; }
        int N() const { return N_; }

        const TensorT &single_layer_tensors(int x, int y) const
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return single_layer_tensors_[x+y*Lx_];
        }
        const TensorT &single_layer_tensors(const std::array<int,2> &coord) const
        {
            return single_layer_tensors(coord[0],coord[1]);
        }
        
        const std::array<IndexT,4> &ordered_virt_indices(int x, int y) const
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return ordered_virt_indices_[x+y*Lx_];
        }
        const std::array<IndexT,4> &ordered_virt_indices(const std::array<int,2> &coord) const
        {
            return ordered_virt_indices(coord[0],coord[1]);
        }

        const TensorT &double_layer_uncontracted_tensors(int x, int y) const
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return double_layer_uncontracted_tensors_[x+y*Lx_];
        }
        const TensorT &double_layer_uncontracted_tensors(const std::array<int,2> &coord) const
        {
            return double_layer_uncontracted_tensors(coord[0],coord[1]);
        }

        const std::array<IndexT,4> &ordered_combined_virt_indices(int x, int y) const
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return ordered_combined_virt_indices_[x+y*Lx_];
        }
        const std::array<IndexT,4> &ordered_combined_virt_indices(const std::array<int,2> &coord) const
        {
            return ordered_combined_virt_indices(coord[0],coord[1]);
        }

        TensorT &bulk_tensors(int x, int y) 
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return bulk_tensors_[x+y*Lx_];
        }
        TensorT &bulk_tensors(const std::array<int,2> &coord)
        {
            return bulk_tensors(coord[0],coord[1]);
        }

        const std::array<TensorT,4> &edge_tensors(int x, int y) const
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return edge_tensors_[x+y*Lx_];
        }
        TensorT &edge_tensors(int x, int y, int edgei)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return edge_tensors_[x+y*Lx_][edgei];
        }
        TensorT &edge_tensors(const std::array<int,3> &coord)
        {
            return edge_tensors(coord[0],coord[1],coord[2]);
        }

        const std::array<TensorT,4> &corner_tensors(int x, int y) const
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return corner_tensors_[x+y*Lx_];
        }
        TensorT &corner_tensors(int x, int y, int corneri)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return corner_tensors_[x+y*Lx_][corneri];
        }
        TensorT &corner_tensors(const std::array<int,3> &coord)
        {
            return corner_tensors(coord[0],coord[1],coord[2]);
        }

        std::array<TensorT,2> &proj_tensors(int x, int y, int proj_dir)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return proj_tensors_[x+y*Lx_][proj_dir];
        }
        std::array<TensorT,2> &proj_tensors(const std::array<int,3> &coord)
        {
            return proj_tensors(coord[0],coord[1],coord[2]);
        }
        TensorT &proj_tensors(int x, int y, int proj_dir, int j)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return proj_tensors_[x+y*Lx_][proj_dir][j];
        }
        TensorT &proj_tensors(const std::array<int,4> &coord)
        {
            return proj_tensors(coord[0],coord[1],coord[2],coord[3]);
        }

        std::vector<double> &singular_vals(int x, int y, int spec_dir)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return singular_vals_[x+y*Lx_][spec_dir];
        }
        std::vector<double> &singular_vals(const std::array<int,3> &coord)
        {
            return singular_vals(coord[0],coord[1],coord[2]);
        }

        TensorT &networks_tensors(int network_no, std::array<int,2> coord)
        {
            coord[0]=(coord[0]+4)%4;
            coord[1]=(coord[1]+4)%4;
            return networks_tensors_[network_no][coord[0]+coord[1]*4];
        }


        std::array<IndexT,4> &bulk_tensors_indices(int x, int y)
        { 
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return bulk_tensors_indices_[x+y*Lx_];
        }
        std::array<IndexT,4> &bulk_tensors_indices(const std::array<int,2> &coord)
        {
            return bulk_tensors_indices(coord[0],coord[1]);
        }

        std::array<IndexT,4> &edge_tensors_indices(int x, int y, int edgei)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return edge_tensors_indices_[x+y*Lx_][edgei];
        }
        std::array<IndexT,4> &edge_tensors_indices(const std::array<int,3> &coord)
        {
            return edge_tensors_indices(coord[0],coord[1],coord[2]);
        }

        std::array<IndexT,4> &corner_tensors_indices(int x, int y, int corneri)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return corner_tensors_indices_[x+y*Lx_][corneri];
        }
        std::array<IndexT,4> &corner_tensors_indices(const std::array<int,3> &coord)
        {
            return corner_tensors_indices(coord[0],coord[1],coord[2]);
        }

        std::array<std::array<IndexT,3>,2> &proj_tensors_indices(int x, int y, int proj_dir)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return proj_tensors_indices_[x+y*Lx_][proj_dir];
        }
        std::array<std::array<IndexT,3>,2> &proj_tensors_indices(const std::array<int,3> &coord)
        {
            return proj_tensors_indices(coord[0],coord[1],coord[2]);
        }
        std::array<IndexT,3> &proj_tensors_indices(int x, int y, int proj_dir, int j)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return proj_tensors_indices_[x+y*Lx_][proj_dir][j];
        }
        std::array<IndexT,3> &proj_tensors_indices(const std::array<int,4> &coord)
        {
            return proj_tensors_indices(coord[0],coord[1],coord[2],coord[3]);
        }

        IndexT &networks_indices(int network_no, int linki) 
        {
            return networks_indices_[network_no][linki]; 
        }
        //return indice connect network site coord_a and coord_b
        IndexT networks_indices(int network_no, std::array<int,2> coord_a, std::array<int,2> coord_b)
        {
            for (int i=0; i<2; i++)
            {
                coord_a[i]=(coord_a[i]+4)%4;
                coord_b[i]=(coord_b[i]+4)%4;
            }

            int site0=coord_a[0]+coord_a[1]*4,
                site1=coord_b[0]+coord_b[1]*4;

            if (site1<site0) 
                return dag(networks_indices(network_no,coord_b,coord_a));

            if (coord_a[0]==(coord_b[0]-1) && coord_a[1]==coord_b[1])
            {
                return networks_indices_[network_no][coord_a[0]+coord_a[1]*3];
            }
            if (coord_a[0]==coord_b[0] && coord_a[1]==(coord_b[1]-1))
            {
                return networks_indices_[network_no][12+coord_a[0]+coord_a[1]*4];
            }

            return IndexT::Null();
        }


        //
        //Methods to get environment tensors
        //
        void obtain_env_tensors();

        //Update env tensor by left/right moves for one col or up/down moves for one row
        void left_move(int left_x0);
        void up_move(int up_y0);
        void right_move(int right_x0);
        void down_move(int down_y0);

        //update (single) corner tensor
        //      C---T---
        //      |   |
        //C'=    \ /
        //        P
        //        |
        void update_corner_tensor(const std::array<int,3> &corner_coord, const std::array<int,3> &edge_coord, const std::array<int,4> &proj_coord);
        //update (single) edge tensor
        //      |
        //      P0
        //     / \
        //    |   |
        //T'= T---a---
        //    |   |
        //     \ /
        //      P1
        //      |
        void update_edge_tensor(const std::array<int,3> &edge_coord, const std::array<int,2> &bulk_coord, const std::array<int,4> &proj_coord0, const std::array<int,4> &proj_coord1);

        //update left/up/right/down boundary of a given network
        void update_network_boundary(int network_no, int boundary_dir);
        void update_network_boundary(const std::array<int,2> &network_coord, int boundary_dir)
        {
            int x=(network_coord[0]%Lx_+Lx_)%Lx_,
                y=(network_coord[1]%Ly_+Ly_)%Ly_;
            int network_no=x+y*Lx_;
            return update_network_boundary(network_no,boundary_dir);
        }

        //obtain proj_tensor P and their indices for partiuclar network and particular move direction
        void obtain_proj_tensor(int network_no, const std::array<int,3> &proj_coord);
        void obtain_proj_tensor(const std::array<int,2> &network_coord , const std::array<int,3> &proj_coord)
        {
            int x=(network_coord[0]%Lx_+Lx_)%Lx_,
                y=(network_coord[1]%Ly_+Ly_)%Ly_;
            int network_no=x+y*Lx_;
            return obtain_proj_tensor(network_no,proj_coord);
        }

        void obtain_R_from_network(int network_no, int dir, std::array<TensorT,2> &R, std::array<IndexT,2> &outside_indices, std::array<IndexT,2> &inside_indices);


        //
        //Methods for measurement
        //
        Complex network_norm(int network_no);
        //measure expectation value for particular operator.
        //notice, we assume this only works for operator act within one network
        //bulk_position are position for bulk tensor in one network 
        //0:(1,1), 1:(2,1), 2:(1,2), 3:(2,2)
        void obtain_ctm_sandwich_operator(int network_no, const std::vector<int> &bulk_position, TPOt<TensorT> tensor_operator);

        //
        //Constructor helpers
        //
        //init tensors and their indices
        void init_double_layer_uncontracted_tensors();
        void init_bulk();
        void init_edge();
        void init_corner();
        void init_network(int network_no);

        //assign network tensors
        void assign_network_tensors(int network_no, int network_site, const TensorT &tens, const std::array<IndexT,4> &ordered_inds)
        {
            auto &network_tensor=networks_tensors_[network_no][network_site];
            network_tensor=tens;
            std::array<int,2> coord={network_site%4,network_site/4};
            if (ordered_inds[Left].valid())
                network_tensor.replaceIndex(ordered_inds[Left],networks_indices(network_no,coord,{coord[0]-1,coord[1]}));
            if (ordered_inds[Up].valid())
                network_tensor.replaceIndex(ordered_inds[Up],networks_indices(network_no,coord,{coord[0],coord[1]+1}));
            if (ordered_inds[Right].valid())
                network_tensor.replaceIndex(ordered_inds[Right],networks_indices(network_no,coord,{coord[0]+1,coord[1]}));
            if (ordered_inds[Down].valid())
                network_tensor.replaceIndex(ordered_inds[Down],networks_indices(network_no,coord,{coord[0],coord[1]-1}));
        }

    private:
        //
        //Data Member
        //
        int Lx_, Ly_, N_;
        //single layer
        std::vector<TensorT> single_layer_tensors_;
        std::vector<std::array<IndexT,4>> ordered_virt_indices_;
        //double_layer_tensors_ stores tensors are tensors formed by bottom and top single_layer_tensors_
        //top tensor is obtained by dag(prime(bottom_tensor))
        //We keep all indices uncontracted and combine virtual legs using CombinerT
        std::vector<TensorT> double_layer_uncontracted_tensors_;
        std::vector<std::array<IndexT,4>> ordered_combined_virt_indices_;
        //we store phys_leg contracted double layer tensors of iPEPS in one u.c. as bulk tensors. The number of sites in one uc equals Lx*Ly
        //we also store their indices, in the order of lurd
        std::vector<TensorT> bulk_tensors_;
        std::vector<std::array<IndexT,4>> bulk_tensors_indices_;
        //corner_tensors_ and edge_tensors_ can be viewed as effective environment of bulk tensors. 
        //At any position, we get four different kind of corner_tensors_(edge_tensors_), where we follow the convention in PRB 84, 041108.
        //we also store the indices of these tensors
        std::vector<std::array<TensorT,4>> edge_tensors_, corner_tensors_;
        std::vector<std::array<std::array<IndexT,4>,4>> edge_tensors_indices_, corner_tensors_indices_;
        //proj_tensors_ are used to update env tens
        //we label as proj_tensors_[network_no][direction][0/1]
        //the proj_tensors_indices_ is ordered as outside_ind, inside_ind and forward_indice 
        std::vector<std::array<std::array<TensorT,2>,4>> proj_tensors_;
        std::vector<std::array<std::array<std::array<IndexT,3>,2>,4>> proj_tensors_indices_;
        //we store svd eigenvalues to check the convergence 
        std::vector<std::array<std::vector<double>,4>> singular_vals_;
        double cutoff_, spec_diff_, converge_diff_;

        //we store tensors in a 4 by 4 network, which is used to update single step CTM of particular site (x,y)
        //
        // C1---T1---T1---C2
        // |    |    |    |
        // |    |    |    |
        // T0---a----a----T2
        // |    |    |    |
        // |    |    |    |
        // T0---a----a----T2
        // |    |    |    |
        // |    |    |    |
        // C0---T3---T3---C3
        //
        //where coordinate of C0 is [x,y]
        //
        //networks_indices are ordered as horizontal and vertical
        //we always choose indices from smaller site to larger site (left to right and down to up)
        std::vector<std::array<TensorT,16>> networks_tensors_;
        std::vector<std::array<IndexT,24>> networks_indices_;
};



#endif
