
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
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return bulk_tensors_[x+y*Lx];
        }
        TensorT &bulk_tensors(const std::array<int,2> &coord)
        {
            return bulk_tensors(coord[0],coord[1]);
        }

        const std::array<TensorT,4> &edge_tensors(int x, int y) const
        {
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return edge_tensors_[x+y*Lx];
        }
        TensorT &edge_tensors(int x, int y, int edgei)
        {
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return edge_tensors_[x+y*Lx][edgei];
        }
        TensorT &edge_tensors(const std::array<int,3> &coord);
        {
            return edge_tensors(coord[0],coord[1],coord[2]);
        }

        const std::array<TensorT,4> &corner_tensors(int x, int y) const
        {
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return corner_tensors_[x+y*Lx];
        }
        TensorT &corner_tensors(int x, int y, int corneri)
        {
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return corner_tensors_[x+y*Lx][corneri]
        }
        TensorT &corner_tensors(const std::array<int,3> &coord);
        {
            return corner_tensors(coord[0],coord[1],coord[2]);
        }

        TensorT &proj_tensors(int x, int y, int i, int j)
        {
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return proj_tensors_[x+y*Lx][i][j];
        }
        TensorT &proj_tensors(const std::array<int,4> &coord)
        {
            return proj_tensors(coord[0],coord[1],coord[2],coord[4]);
        }

        std::array<IndexT,4> &bulk_tensors_indices(int x, int y)
        { 
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return bulk_tensors_indices_[x+y*Lx];
        }
        IndexT &bulk_tensors_indices(int x, int y, int bondi) 
        { 
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return bulk_tensors_indices_[x+y*Lx][bondi]; 
        }
        IndexT &bulk_tensors_indices(const std::array<int,2> &coord, int bondi)
        {
            return bulk_tensors_indices(coord[0],coord[1],bondi);
        }

        std::array<IndexT,4> &edge_tensors_indices(int x, int y, int edgei)
        {
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return edge_tensors_indices_[x+y*Lx][edgei];
        }
        IndexT &edge_tensors_indices(int x, int y, int edgei, int bondi)
        {
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return edge_tensors_indices_[x+y*Lx][edgei][bondi];
        }
        IndexT &edge_tensors_indices(const std::array<int,3> &coord, int bondi)
        {
            return edge_tensors_indices(coord[0],coord[1],coord[2],bondi);
        }

        std::array<IndexT,4> &corner_tensors_indices(int x, int y, int corneri)
        {
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return corner_tensors_indices_[x+y*Lx][corneri];
        }
        IndexT &corner_tensors_indices(int x, int y, int corneri, int bondi)
        {
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return corner_tensors_indices_[x+y*Lx][corneri][bondi];
        }
        IndexT &corner_tensors_indices(const std::array<int,3> &coord, int bondi)
        {
            return corner_tensors_indices(coord[0],coord[1],coord[2],bondi);
        }

        std::array<IndexT,3> &proj_tensors_indices(int x, int y, int i, int j)
        {
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return proj_tensors_indices_[x+y*Lx][i][j];
        }
        IndexT &proj_tensors_indices(int x, int y, int i, int j int bondi)
        {
            x=(x+Lx)%Lx;
            y=(y+Ly)%Ly;
            return proj_tensors_indices_[x+y*Lx][i][j][bondi];
        }
        IndexT &proj_tensors_indices(const std::array<int,4> &coord, int bondi)
        {
            return proj_tensors_indices(coord[0],coord[1],coord[2],coord[3],bondi);
        }

        //
        //Methods to get environment tensors
        //
        void obtain_env_tensors();

        //update env tensors
        void update_edge_tensors();

        //update (single) corner tensor
        //      C---T---
        //      |   |
        //C'=    \ /
        //        P
        //        |
        void update_corner_tensors(const std:array<int,3> &corner_coord, const std:array<int,3> &edge_coord, const std::array<int,3> &proj_coord);
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
        void update_edge_tensors(const std::array<int,3> &edge_coord, const std::array<int,3> &bulk_coord, const std::array<int,4> &proj_coord0, const std::array<int,4> &proj_coord1);

        //
        //Constructor helpers
        //
        //void init_bulk_tensors(const std::vector<TensorT> &single_layer_tensors, const std::vector<std::array<IndexT,4>> &ordered_virt_indices);
        //void init_edge_tensors();
        //void init_corner_tensors():

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
};

//
//class CTM_Network stores 16 tensors as a 4 by 4 network, which is used a single step CTM update of particular (x,y)
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
template <class TensorT>
class CTM_Network
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
        CTM_Network() {};
        //(x,y) are coordinate of C0 
        CTM_Network(const Corner_Transfer_Matrix<TensorT> &current_ctm, int x, int y);

        //
        //Access Method
        //
        IndexT &network_indices(int linki) { return network_indices_[linki]; }
        IndexT network_indices(const std::array<int,2> &coord0, const std::array<int,2> &coord1)
        {
            int site0=coord0[0]+coord0[1]*4,
                site1=coord1[0]+coord1[1]*4;

            if (site1<site0) 
                return dag(network_indices(coord1,coord0));

            if (coord0[0]==(coord1[0]-1) && coord0[1]==coord1[1])
            {
                return network_indices_[coord0[0]+coord0[1]*3];
            }
            if (coord0[0]==coord1[0] && coord0[1]==(coord1[1]-1))
            {
                return network_indices_[12+coord0[0]+coord0[1]*4];
            }

            return IndexT::Null();
        }

        TensorT network_tensors(int x, int y)
        {
            return network_tensors_[x+y*4];
        }

        //assign network tensors
        void assign_network_tensors(const std::array<int,2> &coord, const TensorT &tens, const std::array<IndexT,4> &ordered_inds)
        {
            int site=coord[0]+coord[1]*4;
            network_tensors_[site]=tens;

            if (ordered_inds[Left].valid())
                network_tensors_[site].replaceIndex(ordered_inds[Left],network_indices(coord,{coord[0]-1,coord[1]}));
            if (ordered_inds[Up].valid())
                network_tensors_[site].replaceIndex(ordered_inds[Up],network_indices(coord,{coord[0],coord[1]+1}));
            if (ordered_inds[Right].valid())
                network_tensors_[site].replaceIndex(ordered_inds[Right],network_indices(coord,{coord[0]+1,coord[1]}));
            if (ordered_inds[Down].valid())
                network_tensors_[site].replaceIndex(ordered_inds[Down],network_indices(coord,{coord[0],coord[1]-1}));
        }

        //obtain tensor R's from CTM_Network
        void obtain_R_from_network(Direction dir, std::array<TensorT,2> &R);


    private:
        //network_indices are ordered as horizontal and vertical
        //we always choose indices from smaller site to larger site
        std::array<IndexT,24> network_indices_;
        std::array<Tensor,16> network_tensors_;
};


//obtain projectors P's from R's
//---R[0]===R[1]--- = ---U--s--V--- 
//where s is (truncated) singluar value. So, we get the resolution
//===R[0]^{-1}---R===R[1]---R[1]^{-1}===
//= ===R[1]---V^{\dagger}--s^{-1/2}--s^{-1/2}--U^{\dagger}---R[0]===
//we define projector
//--P[0]=== = --s^{-1/2}--U^{\dagger}---R[0]===
//===P[1]-- = ===R[1]---V^{\dagger}--s^{-1/2}--
template <class TensorT>
void obtain_projectors(const std::array<TensorT,2> &R, std::array<TensorT,2> &P, double cutoff);

#endif
