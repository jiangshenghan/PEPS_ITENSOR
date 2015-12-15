
#ifndef _CORNER_TRANSFER_MAT_H_
#define _CORNER_TRANSFER_MAT_H_

#include ""

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
        enum Direction {Left=0, Up=1, Right=2, Down=3};

        //
        //Constructor
        //
        Corner_Transfer_Matrix() {}
        //construct class using building blocks of ipeps
        //ordered_virt_indices is virtual indices for each tensors, in lurd order
        Corner_Transfer_Matrix(const std::vector<TensorT> &single_layer_tensors, const std::vector<std::array<IndexT,4>> &ordered_virt_indices, int Lx=1, int Ly=1);

        //
        //Access Method
        //
        int Lx() const { return Lx_; }
        int Ly() const { return Ly_; }
        int N() const { return N_; }

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
        TensorT &edge_tensors(const std::array<int,3> &coord);
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
            return corner_tensors_[x+y*Lx_][corneri]
        }
        TensorT &corner_tensors(const std::array<int,3> &coord);
        {
            return corner_tensors(coord[0],coord[1],coord[2]);
        }

        std::array<TensorT,2> &proj_tensors(int x, int y, int proj_dir)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return proj_tensors_[x+y*Lx_][proj_dir];
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

        TensorT &networks_tensors(int network_no, const std::array<int,2> &coord)
        {
            return networks_tensors_[network_no][coord[0]+coord[1]*4];
        }


        std::array<IndexT,4> &bulk_tensors_indices(int x, int y)
        { 
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return bulk_tensors_indices_[x+y*Lx_];
        }
        IndexT &bulk_tensors_indices(int x, int y, int bondi) 
        { 
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return bulk_tensors_indices_[x+y*Lx_][bondi]; 
        }
        IndexT &bulk_tensors_indices(const std::array<int,2> &coord, int bondi)
        {
            return bulk_tensors_indices(coord[0],coord[1],bondi);
        }

        std::array<IndexT,4> &edge_tensors_indices(int x, int y, int edgei)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return edge_tensors_indices_[x+y*Lx_][edgei];
        }
        IndexT &edge_tensors_indices(int x, int y, int edgei, int bondi)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return edge_tensors_indices_[x+y*Lx_][edgei][bondi];
        }
        IndexT &edge_tensors_indices(const std::array<int,3> &coord, int bondi)
        {
            return edge_tensors_indices(coord[0],coord[1],coord[2],bondi);
        }

        std::array<IndexT,4> &corner_tensors_indices(int x, int y, int corneri)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return corner_tensors_indices_[x+y*Lx_][corneri];
        }
        IndexT &corner_tensors_indices(int x, int y, int corneri, int bondi)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return corner_tensors_indices_[x+y*Lx_][corneri][bondi];
        }
        IndexT &corner_tensors_indices(const std::array<int,3> &coord, int bondi)
        {
            return corner_tensors_indices(coord[0],coord[1],coord[2],bondi);
        }

        std::array<std::array<IndexT,3>,2> &proj_tensors_indices(int x, int y, int proj_dir)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return proj_tensors_indices_[x+y*Lx_][proj_dir];
        }
        std::array<IndexT,3> &proj_tensors_indices(int x, int y, int proj_dir, int j)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return proj_tensors_indices_[x+y*Lx_][proj_dir][j];
        }
        IndexT &proj_tensors_indices(int x, int y, int proj_dir, int j int bondi)
        {
            x=(x%Lx_+Lx_)%Lx_;
            y=(y%Ly_+Ly_)%Ly_;
            return proj_tensors_indices_[x+y*Lx_][proj_dir][j][bondi];
        }
        IndexT &proj_tensors_indices(const std::array<int,4> &coord, int bondi)
        {
            return proj_tensors_indices(coord[0],coord[1],coord[2],coord[3],bondi);
        }

        IndexT &networks_indices(int network_no, int linki) 
        {
            return networks_indices_[network_no][linki]; 
        }
        IndexT networks_indices(int network_no, const std::array<int,2> &coord0, const std::array<int,2> &coord1)
        {
            int site0=coord0[0]+coord0[1]*4,
                site1=coord1[0]+coord1[1]*4;

            if (site1<site0) 
                return dag(networks_indices(network_no,coord1,coord0));

            if (coord0[0]==(coord1[0]-1) && coord0[1]==coord1[1])
            {
                return networks_indices_[network_no][coord0[0]+coord0[1]*4];
            }
            if (coord0[0]==coord1[0] && coord0[1]==(coord1[1]-1))
            {
                return networks_indices_[network_no][12+coord0[0]+coord0[1]*3];
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
        void update_corner_tensor(const std:array<int,3> &corner_coord, const std:array<int,3> &edge_coord, const std::array<int,3> &proj_coord);
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
        void update_edge_tensor(const std::array<int,3> &edge_coord, const std::array<int,3> &bulk_coord, const std::array<int,4> &proj_coord0, const std::array<int,4> &proj_coord1);

        //update left/up/right/down boundary of a given network
        void update_network_boundary(int network_no, Direction boundary_dir);
        void update_network_boundary(const std::array<int,2> &network_coord, Direction boundary_dir)
        {
            int x=(network_coord[0]%Lx_+Lx_)%Lx_,
                y=(network_coord[1]%Ly_+Ly_)%Ly_;
            network_no=x+y*Lx_;
            return update_network_boundary(network_no,boundary_dir);
        }

        //obtain proj_tensor P and their indices for partiuclar network and particular move direction
        void obtain_proj_tensor(int network_no, Direction dir, std::array<TensorT,2> &P, std::array<std::array<IndexT,3>,2> &P_indices, double cut_off=1E-4);
        void obtain_proj_tensor(const std::array<int,2> &network_coord , Direction dir, std::array<TensorT,2> &P, std::array<std::array<IndexT,3>,2> &P_indices, double cut_off=1E-4)
        {
            int x=(network_coord[0]%Lx_+Lx_)%Lx_,
                y=(network_coord[1]%Ly_+Ly_)%Ly_;
            network_no=x+y*Lx_;
            return obtain_proj_tensor(network_no,dir,P,P_indices,cut_off);
        }
        void obtain_R_from_network(int network_no, Direction dir, std::array<TensorT,2> &R, std::array<IndexT,2> &outside_indices, std::array<IndexT,2> &inside_indices);


        //
        //Constructor helpers
        //
        //init tensors and their indices
        void init_bulk(const std::vector<TensorT> &single_layer_tensors, const std::vector<std::array<IndexT,4>> &ordered_virt_indices);
        void init_edge(const std::vector<TensorT> &single_layer_tensors, const std::vector<std::array<IndexT,4>> &ordered_virt_indices);
        void init_corner(const std::vector<TensorT> &single_layer_tensors, const std::vector<std::array<IndexT,4>> &ordered_virt_indices);
        void init_network(int x, int y);

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
        //we store double layer tensors of iPEPS in one u.c. as bulk tensors. The number of sites in one uc equals Lx*Ly
        //we also store their indices, in the order of lurd
        std::vector<TensorT> bulk_tensors_;
        std::vector<std::array<IndexT,4>> bulk_tensors_indices_;
        //corner_tensors_ and edge_tensors_ can be viewed as effective environment of bulk tensors. 
        //At any position, we get four different kind of corner_tensors_(edge_tensors_), where we follow the convention in PRB 84, 041108.
        //we also store the indices of these tensors
        std::vector<std::array<TensorT,4>> edge_tensors_, corner_tensors_;
        std::vector<std::array<std::array<IndexT,4>,4>> edge_tensors_indices_, corner_tensors_indices_;
        //proj_tensors_ are used to update env tens
        //the proj_tensors_indices_ is ordered as outside_ind, inside_ind and indice after projection
        std::vector<std::array<std::array<TensorT,2>,4>> proj_tensors_;
        std::vector<std::array<std::array<std::array<IndexT,3>,2>4>> proj_tensors_indices_;
        //
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
        //we always choose indices from smaller site to larger site
        std::vector<std::array<Tensor,16>> networks_tensors_;
        std::vector<std::array<IndexT,24>> networks_indices_;
};



#endif
