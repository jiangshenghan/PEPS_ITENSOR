
#include "corner_transfer_mat.h"

//Class Corner_Transfer_Matrix
template <class TensorT>
Corner_Transfer_Matrix<TensorT>::Corner_Transfer_Matrix(const std::vector<TensorT> &single_layer_tensors, const std::vector<std::array<IndexT,4>> &ordered_virt_indices, int Lx, int Ly): 
    Lx_(Lx), Ly_(Ly), N_(Lx*Ly),
    bulk_tensors_(N_),
    bulk_tensors_indices_(N_),
    edge_tensors_(N_),
    corner_tensors_(N_),
    edge_tensors_indices_(N_),
    corner_tensors_indices_(N_),
    proj_tensors_(N_),
    proj_tensors_indices_(N_)
{
    //init bulk tensors
    for (int sitei=0; sitei<N_; sitei++)
    {
        auto temp_bulk_tensor=single_layer_tensors[sitei]*dag(single_layer_tensors[sitei]).prime(Link);
        std::array<IndexT,4> temp_bulk_tensor_indices;
        for (int bondi=0; bondi<4; bondi++)
        {
            CombinerT temp_combiner(ordered_virt_indices[sitei][bondi],dag(ordered_virt_indices[sitei][bondi]).prime());
            temp_bulk_tensor=temp_bulk_tensor*temp_combiner;
            temp_bulk_tensor_indices[bondi]=temp_combiner.right();
        }
        
        bulk_tensors_[sitei]=temp_bulk_tensor;
        bulk_tensors_indices_[sitei]=temp_bulk_tensor_indices;
    }

    //init edge tensors (three indices) 
    for (int sitei=0; sitei<N_; sitei++)
    {
        for (int edgei=0; edgei<4; edgei++)
        {
            auto temp_edge_tensor=single_layer_tensors[sitei]*dag(single_layer_tensors[sitei]).prime(Link).noprime(ordered_virt_indices[sitei][edgei]);
            std::array<IndexT,4> temp_edge_tensor_indices;
            for (int bondi=0; bondi<4; bondi++)
            {
                if (bondi==edgei) 
                    temp_edge_tensor_indices[bondi]=IndexT::Null();

                CombinerT temp_combiner(ordered_virt_indices[sitei][bondi],dag(ordered_virt_indices[sitei][bondi]).prime());
                temp_edge_tensor=temp_edge_tensor*temp_combiner;
                temp_edge_tensor_indices[bondi]=temp_combiner.right();
            }

            edge_tensors_[sitei][edgei]=temp_edge_tensor;
            edge_tensors_indices_[sitei][edgei]=temp_edge_tensor_indices;
        }
    }

    //init corner tensors (two indices)
    for (int sitei=0; sitei<N_; sitei++)
    {
        std::array<CombinerT,2> temp_combiners;

        //init C0
        corner_tensors_[sitei][0]=single_layer_tensors[sitei]*dag(single_layer_tensors[sitei]).prime(ordered_virt_indices[sitei][1]).prime(ordered_virt_indices[sitei][2]);
        temp_combiners[0]=CombinerT(ordered_virt_indices[sitei][1],dag(ordered_virt_indices[sitei][1]).prime());
        temp_combiners[1]=CombinerT(ordered_virt_indices[sitei][2],dag(ordered_virt_indices[sitei][2]).prime());
        corner_tensors_[sitei][0]=corner_tensors_[sitei][0]*temp_combiners[0]*temp_combiners[1];
        corner_tensors_indices_[sitei][0]=std::array<IndexT,4>{IndexT::Null(),temp_combiners[0].right(),temp_combiners[1].right(),IndexT::Null()};

        //init C1
        corner_tensors_[sitei][1]=single_layer_tensors[sitei]*dag(single_layer_tensors).prime(ordered_virt_indices[sitei][2]).prime(ordered_virt_indices[sitei][3]);
        temp_combiners[0]=CombinerT(ordered_virt_indices[sitei][2],dag(ordered_virt_indices[sitei][2]).prime());
        temp_combiners[1]=CombinerT(ordered_virt_indices[sitei][3],dag(ordered_virt_indices[sitei][3]).prime());
        corner_tensors_[sitei][1]=corner_tensors_[sitei][1]*temp_combiners[0]*temp_combiners[1];
        corner_tensors_indices_[sitei][1]=std::array<IndexT,4>{IndexT::Null(),IndexT::Null(),temp_combiners[0].right(),temp_combiners[1].right()};

        //init C2
        corner_tensors_[sitei][2]=single_layer_tensors[sitei]*dag(single_layer_tensors[sitei]).prime(ordered_virt_indices[sitei][0]).prime(ordered_virt_indices[sitei][3]);
        temp_combiners[0]=CombinerT(ordered_virt_indices[sitei][0],dag(ordered_virt_indices[sitei][0]).prime());
        temp_combiners[1]=CombinerT(ordered_virt_indices[sitei][3],dag(ordered_virt_indices[sitei][3]).prime());
        corner_tensors_[sitei][2]=corner_tensors_[sitei][2]*temp_combiners[0]*temp_combiners[1];
        corner_tensors_indices_[sitei][2]=std::array<IndexT,4>{temp_combiners[0].right(),IndexT::Null(),IndexT::Null(),temp_combiners[1].right()};

        //init C3
        corner_tensors_[sitei][3]=single_layer_tensors[sitei]*dag(single_layer_tensors[sitei]).prime(ordered_virt_indices[sitei][0]).prime(ordered_virt_indices[sitei][1]);
        temp_combiners[0]=CombinerT(ordered_virt_indices[sitei][0],dag(ordered_virt_indices[sitei][0]).prime());
        temp_combiners[1]=CombinerT(ordered_virt_indices[sitei][1],dag(ordered_virt_indices[sitei][1]).prime());
        corner_tensors_[sitei][3]=corner_tensors_[sitei][3]*temp_combiners[0]*temp_combiners[1];
        corner_tensors_indices_[sitei][3]=std::array<IndexT,4>{temp_combiners[0].right(),temp_combiners[1].right(),IndexT::Null(),IndexT::Null()};
    }

}
template 
Corner_Transfer_Matrix<ITensor>::Corner_Transfer_Matrix(const std::vector<ITensor> &single_layer_tensors, const std::vector<std::array<Index,4>> &ordered_virt_indices, int Lx, int Ly);
template 
Corner_Transfer_Matrix<IQTensor>::Corner_Transfer_Matrix(const std::vector<IQTensor> &single_layer_tensors, const std::vector<std::array<IQIndex,4>> &ordered_virt_indices, int Lx, int Ly);


template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::update_corner_tensors(const std:array<int,3> &corner_coord, const std:array<int,3> &edge_coord, const std::array<int,4> &proj_coord)
{
    //dir denotes the moving direction, proj_dir denotes the relative position of proj_tensor
    int dir=proj_coord[2], 
        inv_dir=(dir+2)%4, 
        proj_dir=(edge_coord[2]+2)%4;

    const IndexT &corner_contract_ind=corner_tensors_indices(corner_coord,inv_dir),
                 &edge_contract_ind=edge_tensors_indices(edge_coord,dir),
                 &edge_untouched_ind=edge_tensors_indices(edge_coord,inv_dir);
    //contract original corner tensor and edge tensor, get corner_edge_tensor
    //C---T--
    //|   |
    TensorT corner_edge_tensor=tensor_contraction(corner_tensors(corner_coord),edge_tensors(edge_coord),{corner_contract_ind},{edge_contract_ind});

    //act projector on corner_edge_tensor, get the corner_tensor_prime
    //C---T--
    //|   |
    // \ /
    //  P
    //  |
    //we get corner_tensor_prime with edge_untouched_ind and proj_tensors_indice(2)
    IndexT outside_ind=corner_tensors_indices(corner_coord,proj_dir),
           inside_ind=edge_tensors_indices(edge_coord,proj_dir);
    TensorT corner_tensor_prime=tensor_contraction(corner_edge_tensor,proj_tensors(proj_coord),{outside_ind,inside_ind},{proj_tensors_indices(proj_coord,0),proj_tensors_indices(proj_coord,1)});
    
    //update the corner_tensor as well as its indices
    auto updated_corner_coord=edge_coord;
    updated_corner_coord[2]=corner_coord[2];
    auto &updated_indices=corner_tensors_indices(updated_corner_coord);
    //index_assignment(updated_indices(inv_dir),edge_untouched_ind);
    index_assignment(updated_indices(proj_dir),proj_tensors_indices(proj_coord,2));

    corner_tensors(updated_corner_coord)=corner_tensor_prime;
    corner_tensors(updated_corner_coord).replaceIndex(edge_untouched_ind,updated_indices(inv_dir));
    corner_tensors(updated_corner_coord).replaceIndex(proj_tensors_indices(proj_coord,2),updated_indices(proj_dir));
}
template
void Corner_Transfer_Matrix<ITensor>::update_corner_tensors(const std:array<int,3> &corner_coord, const std:array<int,3> &edge_coord, const std::array<int,3> &proj_coord);
template
void Corner_Transfer_Matrix<IQTensor>::update_corner_tensors(const std:array<int,3> &corner_coord, const std:array<int,3> &edge_coord, const std::array<int,3> &proj_coord);

template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::update_edge_tensors(const std::array<int,3> &edge_coord, const std::array<int,3> &bulk_coord, const std::array<int,4> &proj_coord0, const std::array<int,4> &proj_coord1)
{
    //dir denotes the moving direction, proj_diri denotes the relateve
    //position of proj_tensori
    int dir=proj_coord0[2],
        inv_dir=(dir+2)%4,
        proj_dir1=proj_coord1[0]-edge_coord[0]+edge_coord[1]-proj_coord1[1]+3,
        proj_dir0=(proj_dir1+2)%4;

    const IndexT &edge_contract_ind=edge_tensors_indices(edge_coord,inv_dir),
                 &bulk_contract_ind=bulk_tensors_indices(bulk_coord,dir),
                 &bulk_untouched_ind=bulk_tensors_indices(bulk_coord,inv_dir),
                 &outside_ind0=edge_tensors_indices(edge_coord,proj_dir0),
                 &inside_ind0=bulk_tensors_indices(bulk_coord,proj_dir0),
                 &outside_ind1=edge_tensors_indices(edge_coord,proj_dir1),
                 &inside_ind1=bulk_tensors_indices(bulk_coord,proj_dir1);

    //      |
    //      P0
    //     / \
    //    |   |
    //T'= T---a---
    //    |   |
    //     \ /
    //      P1
    //      |
    const TensorT &proj_tens0=proj_tensors(proj_coord0),
                  &proj_tens1=proj_tensors(proj_coord1),
                  &curr_edge_tens=edge_tensors(edge_coord),
                  &bulk_tens=bulk_tensors(bulk_coord);
    TensorT edge_tensor_prime=tensor_contraction({curr_edge_tens,proj_tens0,bulk_tens,proj_tens1},{{outside_ind0},{proj_tensors_indices(proj_coord0,0)},{edge_contract_ind,proj_tensors_indices(proj_coord0,1)},{bulk_contract_ind,inside_ind0},{outside_ind1,inside_ind1},{proj_tensors_indices(proj_coord1,0),proj_tensors_indices(proj_coord1,1)}});

    //update edge tensors as well as its indices
    auto updated_edge_coord=bulk_coord;
    updated_edge_coord[2]=edge_coord[2];

    auto &updated_indices=edge_tensors_indices(updated_edge_coord);
    //index_assignment(updated_indices(inv_dir),bulk_untouched_ind);
    index_assignment(updated_indices(proj_dir0),proj_tensors_indices(proj_coord0,2));
    index_assignment(updated_indices(proj_dir1),proj_tensors_indices(proj_coord1,2));

    TensorT &updated_edge_tensor=edge_tensors(updated_edge_coord);
    updated_edge_tensor=edge_tensor_prime;
    updated_edge_tensor.replaceIndex(bulk_untouched_ind,updated_indices(inv_dir));
    updated_edge_tensor.replaceIndex(proj_tensors_indices(proj_coord0,2),updated_indices(proj_dir0));
    updated_edge_tensor.replaceIndex(proj_tensors_indices(proj_coord1,2),updated_indices(proj_dir1));
}
template
void Corner_Transfer_Matrix<ITensor>::update_edge_tensors(const std::array<int,3> &edge_coord, const std::array<int,3> &bulk_coord, const std::array<int,4> &proj_coord0, const std::array<int,4> &proj_coord1)
template
void Corner_Transfer_Matrix<IQTensor>::update_edge_tensors(const std::array<int,3> &edge_coord, const std::array<int,3> &bulk_coord, const std::array<int,4> &proj_coord0, const std::array<int,4> &proj_coord1)

//Class CTM_Network
template <class TensorT>
CTM_Network<TensorT>::CTM_Network(const Corner_Transfer_Matrix<TensorT> &current_ctm, int x, int y)
{

    //init network_indices_
    index_assignment(network_indices_[0],current_ctm.corner_tensors_indices(x,y,0,Right),"network_ind0");
    index_assignment(network_indices_[1],current_ctm.edge_tensors_indices(x+1,y,3,Right),"network_ind1");
    index_assignment(network_indices_[2],current_ctm.edge_tensors_indices(x+2,y,3,Right),"network_ind2");
    index_assignment(network_indices_[3],current_ctm.edge_tensors_indices(x,y+1,0,Right),"network_ind3");
    index_assignment(network_indices_[4],current_ctm.bulk_tensors_indices(x+1,y+1,Right),"network_ind4");
    index_assignment(network_indices_[5],current_ctm.bulk_tensors_indices(x+2,y+1,Right),"network_ind5");
    index_assignment(network_indices_[6],current_ctm.edge_tensors_indices(x,y+2,0,Right),"network_ind6");
    index_assignment(network_indices_[7],current_ctm.bulk_tensors_indices(x+1,y+2,Right),"network_ind7");
    index_assignment(network_indices_[8],current_ctm.bulk_tensors_indices(x+2,y+2,Right),"network_ind8");
    index_assignment(network_indices_[9],current_ctm.corner_tensors_indices(x,y+3,1,Right),"network_ind9");
    index_assignment(network_indices_[10],current_ctm.edge_tensors_indices(x+1,y+3,1,Right),"network_ind10");
    index_assignment(network_indices_[11],current_ctm.edge_tensors_indices(x+2,y+3,1,Right),"network_ind11");
    index_assignment(network_indices_[12],current_ctm.corner_tensors_indices(x,y,0,Up),"network_ind12");
    index_assignment(network_indices_[13],current_ctm.edge_tensors_indices(x+1,y,3,Up),"network_ind13");
    index_assignment(network_indices_[14],current_ctm.edge_tensors_indices(x+2,y,3,Up),"network_ind14");
    index_assignment(network_indices_[15],current_ctm.corner_tensors_indices(x+3,y,3,Up),"network_ind15");
    index_assignment(network_indices_[16],current_ctm.edge_tensors_indices(x,y+1,0,Up),"network_ind16");
    index_assignment(network_indices_[17],current_ctm.bulk_tensors_indices(x+1,y+1,Up),"network_ind17");
    index_assignment(network_indices_[18],current_ctm.bulk_tensors_indices(x+2,y+1,Up),"network_ind18");
    index_assignment(network_indices_[19],current_ctm.edge_tensors_indices(x+3,y+1,2,Up),"network_ind19");
    index_assignment(network_indices_[20],current_ctm.edge_tensors_indices(x,y+2,0,Up),"network_ind20");
    index_assignment(network_indices_[21],current_ctm.bulk_tensors_indices(x+1,y+2,Up),"network_ind21");
    index_assignment(network_indices_[22],current_ctm.bulk_tensors_indices(x+2,y+2,Up),"network_ind22");
    index_assignment(network_indices_[23],current_ctm.edge_tensors_indices(x+3,y+2,2,Up),"network_ind23");

    //init network_tensors_
    assign_network_tensors({0,0},current_ctm.corner_tensors(x,y,0),current_ctm.corner_tensors_indices(x,y,0));
    assign_network_tensors({1,0},current_ctm.edge_tensors(x+1,y,3),current_ctm.edge_tensors_indices(x+1,y,3));
    assign_network_tensors({2,0},current_ctm.edge_tensors(x+2,y,3),current_ctm.edge_tensors_indices(x+2,y,3));
    assign_network_tensors({3,0},current_ctm.corner_tensors(x+3,y,3),current_ctm.corner_tensors_indices(x+3,y,3));
    assign_network_tensors({0,1},current_ctm.edge_tensors(x,y+1,0),current_ctm.edge_tensors_indices(x,y+1,0));
    assign_network_tensors({1,1},current_ctm.bulk_tensors(x+1,y+1),current_ctm.bulk_tensors_indices(x+1,y+1));
    assign_network_tensors({2,1},current_ctm.bulk_tensors(x+2,y+1),current_ctm.bulk_tensors_indices(x+2,y+1));
    assign_network_tensors({3,1},current_ctm.edge_tensors(x+3,y+1,2),current_ctm.edge_tensors_indices(x+3,y+1,2));
    assign_network_tensors({0,2},current_ctm.edge_tensors(x,y+2,0),current_ctm.edge_tensors_indices(x,y+2,0));
    assign_network_tensors({1,2},current_ctm.bulk_tensors(x+1,y+2),current_ctm.bulk_tensors_indices(x+1,y+2));
    assign_network_tensors({2,2},current_ctm.bulk_tensors(x+2,y+2),current_ctm.bulk_tensors_indices(x+2,y+2));
    assign_network_tensors({3,2},current_ctm.edge_tensors(x+3,y+2,2),current_ctm.edge_tensors_indices(x+3,y+2,2));
    assign_network_tensors({0,3},current_ctm.corner_tensors(x,y+3,1),current_ctm.corner_tensors_indices(x,y+3,1));
    assign_network_tensors({1,3},current_ctm.edge_tensors(x+1,y+3,1),current_ctm.edge_tensors_indices(x+1,y+3,1));
    assign_network_tensors({2,3},current_ctm.edge_tensors(x+2,y+3,1),current_ctm.edge_tensors_indices(x+2,y+3,1));
    assign_network_tensors({3,3},current_ctm.corner_tensors(x+3,y+3,2),current_ctm.corner_tensors_indices(x+3,y+3,2));

}
template
CTM_Network<ITensor>::CTM_Network(const Corner_Transfer_Matrix<ITensor> &current_ctm, int x, int y);
template
CTM_Network<IQTensor>::CTM_Network(const Corner_Transfer_Matrix<IQTensor> &current_ctm, int x, int y);

template <class TensorT>
void CTM_Network<TensorT>::obtain_R_from_network(Direction dir, std::array<TensorT,2> &R)
{
    std::array<TensorT,2> Q;

    TensorT upper_left_tensor=network_tensors(0,2)*network_tensors(0,3)*network_tensors(1,3)*network_tensors(1,2),
            upper_right_tensor=network_tensors(3,2)*network_tensors(3,3)*network_tensors(2,3)*network_tensors(2,2),
            lower_left_tensor=network_tensors(0,1)*network_tensors(0,0)*network_tensors(1,0)*network_tensors(1,1),
            lower_right_tensor=network_tensors(3,1)*network_tensors(3,0)*network_tensors(2,0)*network_tensors(2,1);
    //TODO: modify following
    if (dir==Left || dir==Right) 
    {
        TensorT upper_half_tensor=upper_left_tensor*upper_right_tensor,
                lower_half_tensor=lower_left_tensor*lower_right_tensor;

        if (dir==Left) 
        {
            R[0]=TensorT(network_indices({0,2},{0,1}),network_indices({1,2},{1,1}));
            R[1]=dag(R[0]);
            denmatDecomp(upper_half_tensor,R[0],Q[0],FromRight);
            denmatDecomp(lower_half_tensor,R[1],Q[1],FromRight);
        }
        if (dir==Right)
        {
            R[0]=TensorT(network_indices({3,1},{3,2}),network_indices({2,1},{2,2}));
            R[1]=dag(R[0]);
            denmatDecomp(lower_half_tensor,R[0],Q[0],FromRight);
            denmatDecomp(upper_half_tensor,R[1],Q[1],FromRight);
        }
    }

    if (dir==Up || dir==Down)
    {
        TensorT left_half_tensor=upper_left_tensor*lower_left_tensor,
                right_half_tensor=upper_right_tensor*lower_right_tensor;

        if (dir==Up)
        {
            R[0]=TensorT(network_indices({2,3},{1,3}),network_indices({2,2},{1,2}));
            R[1]=dag(R[0]);
            denmatDecomp(right_half_tensor,R[0],Q[0],FromRight);
            denmatDecomp(left_half_tensor,R[1],Q[1],FromRight);
        }
        if (dir==Down)
        {
            R[0]=TensorT(network_indices({1,0},{2,0}),network_indices({1,1},{2,1}));
            R[1]=dag(R[0]);
            denmatDecomp(left_half_tensor,R[0],Q[0],FromRight);
            denmatDecomp(right_half_tensor,R[1],Q[1],FromRight);
        }
    }

}
template
void CTM_Network<ITensor>::obtain_R_from_network(Direction dir, std::array<ITensor,2> &R)
template
void CTM_Network<IQTensor>::obtain_R_from_network(Direction dir, std::array<IQTensor,2> &R)



//
//Other methods
//
template <class TensorT>
void obtain_projectors(const std::array<TensorT,2> &R, std::array<TensorT,2> &P, double cutoff)
{
    TensorT::IndexT U_ind=uniqueIndex(R[0],R[1]);
    TensorT U(U_ind),s,V,s_inv_sqrt;
    OptSet opts;
    opts.add("Cutoff",cutoff);
    //TODO: be careful about small singular value, use SVDThreshold in opts?
    svd(R[0]*R[1],U,s,V,opts);

    //get s^{-1/2}
    s_inv_sqrt=dag(s);
    s_inv_sqrt.pseudoInvert();
    std::function<double(double)> Sqrt=sqrt;
    s_inv_sqrt.mapElems(Sqrt);

    P[0]=s_inv_sqrt*dag(U)*R[0];
    P[1]=R[1]*dag(V)*s_inv_sqrt;

    //TODO: change index of P?
}
template
void obtain_projectors(const std::array<ITensor,2> &R, std::array<ITensor,2> &P, double cutoff);
template
void obtain_projectors(const std::array<IQTensor,2> &R, std::array<IQTensor,2> &P, double cutoff);
