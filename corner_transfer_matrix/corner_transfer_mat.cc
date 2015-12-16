
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
    proj_tensors_indices_(N_),
    singular_vals_(N_),
    networks_tensors_(N_),
    networks_indices_(N_)
{
    init_bulk(single_layer_tensors,ordered_virt_indices);
    init_edge(single_layer_tensors,ordered_virt_indices);
    init_corner(single_layer_tensors,ordered_virt_indices);

    for (int sitei=0; sitei<N_; sitei++) init_network(sitei);

    for (int sitei=0; sitei<N_; sitei++)
    {
        for (int diri=0; diri<4; dir++)
        {
            singular_vals_[sitei][diri]=std::vector<double>(200,0);
        }
    }

    obtain_env_tensors();
}
template 
Corner_Transfer_Matrix<ITensor>::Corner_Transfer_Matrix(const std::vector<ITensor> &single_layer_tensors, const std::vector<std::array<Index,4>> &ordered_virt_indices, int Lx, int Ly);
template 
Corner_Transfer_Matrix<IQTensor>::Corner_Transfer_Matrix(const std::vector<IQTensor> &single_layer_tensors, const std::vector<std::array<IQIndex,4>> &ordered_virt_indices, int Lx, int Ly);


template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::obtain_env_tensors()
{
    do
    {
        spec_diff_=0;

        for (int x=0; x<Lx_; x++) left_move(x);
        for (int y=Ly_-1; y>=0; y++) up_move(y);
        for (int x=Lx_-1; x>=0; x--) right_move(x);
        for (int y=0; y<Ly_; y++) down_move(y);
    }
    while (spec_diff_<1E-3);
}
template
void Corner_Transfer_Matrix<ITensor>::obtain_env_tensors();
template
void Corner_Transfer_Matrix<IQTensor>::obtain_env_tensors();


template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::left_move(int left_x0)
{
    //obtain proj_tensors at (left_x0,y)
    for (int y=0; y<Ly_; y++)
    {
        obtain_proj_tensor({left_x0,y-1},{Left_x0,y,Left});
    }

    //update env_tensors using proj_tensors obtained above
    for (int y=0; y<Ly_; y++)
    {
        update_corner_tensor({left_x0,y,0},{left_x0+1,y,3},{left_x0,y,Left,0});
        update_corner_tensor({left_x0,y,1},{left_x0+1,y,1},{left_x0,y-1,Left,0});
        update_edge_tensor({left_x0,y,0},{left_x0+1,y},{left_x0,y,Left,0},{left_x0,y-1,1});
    }

    
    //replace networks_tensors with updated env_tensors
    for (int y=0; y<Ly_; y++)
    {
        update_network_boundary({left_x0+1,y},Left);
    }

}
template
void Corner_Transfer_Matrix<ITensor>::left_move(int left_x0);
template
void Corner_Transfer_Matrix<IQTensor>::left_move(int left_x0);


template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::up_move(int up_y0)
{
    //obtain proj_tensors 
    for (int x=0; x<Lx_; x++)
    {
        obtain_proj_tensor({x-1,up_y0-3},{x,up_y0,Up});
    }

    //update env_tensors using proj_tensors obtained above
    for (int x=0; x<Lx_; x++)
    {
        update_corner_tensor({x,up_y0,1},{x,up_y0-1,0},{x,up_y0,Up,0});
        update_corner_tensor({x,up_y0,2},{x,up_y0-1,2},{x-1,up_y0,Up,1});
        update_edge_tensor({x,up_y0,1},{x,up_y0-1},{x,up_y0,Up,0},{x-1,up_y0,Up,1});
    }
    
    //replace networks_tensors with updated env_tensors
    for (int x=0; x<Lx_; x++)
    {
        update_network_boundary({x,up_y0-1},Up);
    }

}
template
void Corner_Transfer_Matrix<ITensor>::up_move(int up_y0);
template
void Corner_Transfer_Matrix<IQTensor>::up_move(int up_y0);


template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::right_move(int right_x0)
{
    //obtain proj_tensors
    for (int y=0; y<Ly_; y++)
    {
        obtain_proj_tensor({right_x0-3,y-2},{right_x0,y,Right});
    }

    //update env_tensors using proj_tensors obtained above
    for (int y=0; y<Ly_; y++)
    {
        update_corner_tensor({right_x0,y,2},{right_x0-1,y,1},{right_x0,y,Right,0});
        update_corner_tensor({right_x0,y,3},{right_x0-1,y,3},{right_x0,y+1,Right,1});
        update_edge_tensor({right_x0,y,2},{right_x0-1,y},{right_x0,y,Right,0},{right_x0,y+1,Right,1});
    }

    //replace networks_tensors with updated env_tensors
    for (int y=0; y<Ly_; y++)
    {
        update_network_boundary({right_x0-1,y},Right);
    }

}
template
void Corner_Transfer_Matrix<ITensor>::right_move(int right_x0);
template
void Corner_Transfer_Matrix<IQTensor>::right_move(int right_x0);


template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::down_move(int down_y0)
{
    //obtain proj_tensors
    for (int x=0; x<Lx_; x++)
    {
        obtain_proj_tensor({x-2,down_y0},{x,down_y0,Down});
    }

    //update env_tensors_using proj_tensors obtained above
    for (int x=0; x<Lx_; x++)
    {
        update_corner_tensor({x,down_y0,0},{x,down_y0+1,0},{x+1,down_y0,Down,1});
        update_corner_tensor({x,down_y0,3},{x,down_y0+1,2},{x,down_y0,Down,0});
        update_edge_tensor({x,down_y0,3},{x,down_y0+1},{x,down_y0,Down,0},{x+1,down_y0,Down,1});
    }
    
    //replace networks_tensors with updated env_tensors
    for (int x=0; x<Lx_; x++)
    {
        update_network_boundary({x,down_y0+1},Down);
    }

}
template
void Corner_Transfer_Matrix<ITensor>::down_move(down_y0);
template
void Corner_Transfer_Matrix<IQTensor>::down_move(down_y0);


template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::update_corner_tensor(const std:array<int,3> &corner_coord, const std:array<int,3> &edge_coord, const std::array<int,4> &proj_coord)
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
    TensorT corner_tensor_prime=tensor_contraction(corner_edge_tensor,proj_tensors(proj_coord),{outside_ind,inside_ind},{proj_tensors_indices(proj_coord)[0],proj_tensors_indices(proj_coord)[1]});
    
    //update the corner_tensor as well as its indices
    auto updated_corner_coord=edge_coord;
    updated_corner_coord[2]=corner_coord[2];
    auto &updated_indices=corner_tensors_indices(updated_corner_coord);
    //index_assignment(updated_indices(inv_dir),edge_untouched_ind);
    index_assignment(updated_indices[proj_dir],proj_tensors_indices(proj_coord)[2]);

    corner_tensors(updated_corner_coord)=corner_tensor_prime;
    corner_tensors(updated_corner_coord).replaceIndex(edge_untouched_ind,updated_indices[inv_dir]);
    corner_tensors(updated_corner_coord).replaceIndex(proj_tensors_indices(proj_coord)[2],updated_indices[proj_dir]);
}
template
void Corner_Transfer_Matrix<ITensor>::update_corner_tensor(const std:array<int,3> &corner_coord, const std:array<int,3> &edge_coord, const std::array<int,3> &proj_coord);
template
void Corner_Transfer_Matrix<IQTensor>::update_corner_tensor(const std:array<int,3> &corner_coord, const std:array<int,3> &edge_coord, const std::array<int,3> &proj_coord);

template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::update_edge_tensor(const std::array<int,3> &edge_coord, const std::array<int,3> &bulk_coord, const std::array<int,4> &proj_coord_a, const std::array<int,4> &proj_coord_b)
{
    //dir denotes the moving direction, proj_diri denotes the relateve
    //position of proj_tensori
    int dir=proj_coord_a[2],
        inv_dir=(dir+2)%4,
        proj_dir_a=(dir+1)%4,
        proj_dir_b=(dir+3)%4;

    const IndexT &edge_contract_ind=edge_tensors_indices(edge_coord)[inv_dir],
                 &bulk_contract_ind=bulk_tensors_indices(bulk_coord)[dir],
                 &bulk_untouched_ind=bulk_tensors_indices(bulk_coord)[inv_dir],
                 &outside_ind_a=edge_tensors_indices(edge_coord)[proj_dir_a],
                 &inside_ind_a=bulk_tensors_indices(bulk_coord)[proj_dir_a],
                 &outside_ind_b=edge_tensors_indices(edge_coord)[proj_dir_b],
                 &inside_ind_b=bulk_tensors_indices(bulk_coord)[proj_dir_b];

    //      |
    //      P0
    //     / \
    //    |   |
    //T'= T---a---
    //    |   |
    //     \ /
    //      P1
    //      |
    const TensorT &proj_tens_a=proj_tensors(proj_coord_a),
                  &proj_tens_b=proj_tensors(proj_coord_b),
                  &curr_edge_tens=edge_tensors(edge_coord),
                  &bulk_tens=bulk_tensors(bulk_coord);
    TensorT edge_tensor_prime=tensor_contraction({curr_edge_tens,proj_tens_a,bulk_tens,proj_tens_b},{{outside_ind_a},{proj_tensors_indices(proj_coord_a)[0]},{edge_contract_ind,proj_tensors_indices(proj_coord_a)[1]},{bulk_contract_ind,inside_ind_a},{outside_ind_b,inside_ind_b},{proj_tensors_indices(proj_coord_b)[0],proj_tensors_indices(proj_coord_b)[1]}});

    //update edge tensors as well as its indices
    auto updated_edge_coord=edge_coord;
    updated_edge_coord[0]=bulk_coord[0];
    updated_edge_coord[1]=bulk_coord[1];

    auto &updated_indices=edge_tensors_indices(updated_edge_coord);
    //index_assignment(updated_indices(inv_dir),bulk_untouched_ind);
    index_assignment(updated_indices(proj_dir_a),proj_tensors_indices(proj_coord_a)[2]);
    index_assignment(updated_indices(proj_dir_b),proj_tensors_indices(proj_coord_b)[2]);

    TensorT &updated_edge_tensor=edge_tensors(updated_edge_coord);
    updated_edge_tensor=edge_tensor_prime;
    updated_edge_tensor.replaceIndex(bulk_untouched_ind,updated_indices[inv_dir]);
    updated_edge_tensor.replaceIndex(proj_tensors_indices(proj_coord_a)[2],updated_indices[proj_dir_a]);
    updated_edge_tensor.replaceIndex(proj_tensors_indices(proj_coord_b)[2],updated_indices[proj_dir_b]);
}
template
void Corner_Transfer_Matrix<ITensor>::update_edge_tensor(const std::array<int,3> &edge_coord, const std::array<int,3> &bulk_coord, const std::array<int,4> &proj_coord_a, const std::array<int,4> &proj_coord_b);
template
void Corner_Transfer_Matrix<IQTensor>::update_edge_tensor(const std::array<int,3> &edge_coord, const std::array<int,3> &bulk_coord, const std::array<int,4> &proj_coord_a, const std::array<int,4> &proj_coord_b);


template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::update_network_boundary(int network_net, Direction boundary_dir)
{
    int x=network_no%Lx_,
        y=network_no/Lx_;
    auto &network_indices=networks_indices_[network_no];

    if (boundary_dir==Left)
    {
        index_assignment(network_indices[12],corner_tensors_indices(x,y,0)[Up],"network_ind12");
        index_assignment(network_indices[16],edge_tensors_indices(x,y+1,0)[Up],"network_ind16");
        index_assignment(network_indices[20],edge_tensors_indices(x,y+2,0)[Up],"network_ind20");

        assign_network_tensors(network_no,0,corner_tensors(x,y,0),corner_tensors_indices(x,y,0));
        assign_network_tensors(network_no,4,edge_tensors(x,y+1,0),edge_tensors_indices(x,y+1,0));
        assign_network_tensors(network_no,8,edge_tensors(x,y+2,0),edge_tensors_indices(x,y+2,0));
        assign_network_tensors(network_no,12,corner_tensors(x,y+3,1),corner_tensors_indices(x,y+3,1));

        return;
    }

    if (boundary_dir==Up)
    {
        index_assignment(network_indices[9],corner_tensors_indices(x,y+3,1)[Right],"network_ind9");
        index_assignment(network_indices[10],edge_tensors_indices(x+1,y+3,1)[Right],"network_ind10");
        index_assignment(network_indices[11],edge_tensors_indices(x+2,y+3,1)[Right],"network_ind11");

        assign_network_tensors(network_no,12,corner_tensors(x,y+3,1),corner_tensors_indices(x,y+3,1));
        assign_network_tensors(network_no,13,edge_tensors(x+1,y+3,1),edge_tensors_indices(x+1,y+3,1));
        assign_network_tensors(network_no,14,edge_tensors(x+2,y+3,1),edge_tensors_indices(x+2,y+3,1));
        assign_network_tensors(network_no,15,corner_tensors(x+3,y+3,2),corner_tensors_indices(x+3,y+3,2));
        
        return;
    }

    if (boundary_dir==Right)
    {
        index_assignment(network_indices[15],corner_tensors_indices(x+3,y,3)[Up],"network_ind15");
        index_assignment(network_indices[19],edge_tensors_indices(x+3,y+1,2)[Up],"network_ind19");
        index_assignment(network_indices[23],edge_tensors_indices(x+3,y+2,2)[Up],"network_ind23");

        assign_network_tensors(network_no,3,corner_tensors(x+3,y,3),corner_tensors_indices(x+3,y,3));
        assign_network_tensors(network_no,7,edge_tensors(x+3,y+1,2),edge_tensors_indices(x+3,y+1,2));
        assign_network_tensors(network_no,11,edge_tensors(x+3,y+2,2),edge_tensors_indices(x+3,y+2,2));
        assign_network_tensors(network_no,15,corner_tensors(x+3,y+3,2),corner_tensors_indices(x+3,y+3,2));

        return;
    }

    if (boundary_dir==Down)
    {
        index_assignment(network_indices[0],corner_tensors_indices(x,y,0)[Right],"network_ind0");
        index_assignment(network_indices[1],edge_tensors_indices(x+1,y,3)[Right],"network_ind1");
        index_assignment(network_indices[2],edge_tensors_indices(x+2,y,3)[Right],"network_ind2");

        assign_network_tensors(network_no,0,corner_tensors(x,y,0),corner_tensors_indices(x,y,0));
        assign_network_tensors(network_no,1,edge_tensors(x+1,y,3),edge_tensors_indices(x+1,y,3));
        assign_network_tensors(network_no,2,edge_tensors(x+2,y,3),edge_tensors_indices(x+2,y,3));
        assign_network_tensors(network_no,3,corner_tensors(x+3,y,3),corner_tensors_indices(x+3,y,3));

        return;
    }

}
template
void Corner_Transfer_Matrix<ITensor>::update_network_boundary(int network_net, Direction boundary_dir);
template
void Corner_Transfer_Matrix<IQTensor>::update_network_boundary(int network_net, Direction boundary_dir);


template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::obtain_proj_tensor(int network_no, const std::array<int,3> &proj_coord, double cut_off)
{
    //obtain R by contracting network and QR decomp
    std::array<TensorT,2> R;
    std::array<IndexT,2> outside_indices, inside_indices;
    obtain_R_from_network(network_no,proj_coord[2],R,outside_indices,inside_indices);

    //obtain projectors P's from R's
    //---R[0]===R[1]--- = ---U--s--V--- 
    //where s is (truncated) singluar value. So, we get the resolution
    //===R[0]^{-1}---R===R[1]---R[1]^{-1}===
    //= ===R[1]---V^{\dagger}--s^{-1/2}--s^{-1/2}--U^{\dagger}---R[0]===
    //we define projector
    //--P[0]=== = --s^{-1/2}--U^{\dagger}---R[0]===
    //===P[1]-- = ===R[1]---V^{\dagger}--s^{-1/2}--
    //We always make s normalized in order to check convergence
    IndexT U_ind=uniqueIndex(R[0],R[1]);
    TensorT U(U_ind),s,V,s_inv_sqrt;
    //TODO: be careful about small singular value, use SVDThreshold in opts?
    OptSet opts;
    opts.add("Cutoff",cutoff);
    svd(R[0]*R[1],U,s,V,opts);
    
    //obtain the distance between (normalized) singular value of current step with that of the last step. Then, add the distance to spec_diff_ to test convergence
    auto temp_singular_val=std::vector<double>(singular_vals(proj_coord).size(),0);
    auto s_diag=s.diag();
    int length=std::min(temp_singular_val.size(),s_diag.Length());
    double temp_singular_val_norm=0;
    for (int diagi=0; diagi<length; diagi++)
    {
        auto diag_elem=s_diag(diagi+1);
        temp_sigular_val[diagi]=diag_elem;
        temp_singular_val_norm+=diag_elem*diag_elem;
    }
    temp_singular_val_norm=sqrt(temp_singular_val_norm);

    double update_val_distance_sq=0;
    for (int diagi=0; diagi<temp_singular_val.size(); diagi++)
    {
        temp_singular_val[diagi]/=temp_singular_val_norm;
        update_val_distance_sq+=pow(singular_vals(proj_coord)[diagi]-temp_singular_val[diagi],2)
    }
    spec_diff_+=sqrt(update_val_distance_sq);
    singular_vals(proj_coord)=temp_singular_val;

    //get s^{-1/2}
    auto s_inv_sqrt=dag(s);
    s_inv_sqrt.pseudoInvert();
    std::function<double(double)> Sqrt=sqrt;
    s_inv_sqrt.mapElems(Sqrt);

    auto &P=proj_tensors(proj_coord);
    auto &P_indices=proj_tensors_indices(proj_coord);
    P[0]=s_inv_sqrt*dag(U)*R[0];
    P[1]=R[1]*dag(V)*s_inv_sqrt;

    //change index of P to avoid wrong contraction
    std::array<IndexT,2> forward_indices;
    for (int i=0; i<2; i++)
    {
        forward_indices[i]=uniqueIndex(P[i],TensorT(outside_indices[i][0],inside_indices[i][0]));

        index_assignment(P_indices[i][0],outside_indices[i]);
        index_assignment(P_indices[i][1],inside_indices[i]);
        index_assignment(P_indices[i][2],forward_indices[i]);
        
        P[i].replaceIndex(outside_indices[i],P_indices[i][0]);
        P[i].replaceIndex(inside_indices[i],P_indices[i][1]);
        P[i].replaceIndex(forward_indices[i],P_indices[i][2]);
    }

}
template
void Corner_Transfer_Matrix<ITensor>::obtain_proj_tensor(int network_no, const std::array<int,3> &proj_coord, double cut_off);
template
void Corner_Transfer_Matrix<IQTensor>::obtain_proj_tensor(int network_no, const std::array<int,3> &proj_coord, double cut_off);


//TODO: we should write QR decomposition rather than use denmatDecomp
template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::obtain_R_from_network(int network_no, Direction dir, std::array<TensorT,2> &R, std::array<IndexT,2> &outside_indices, std::array<IndexT,2> &inside_indices)
{
    std::array<TensorT,2> Q;

    TensorT upper_left_tensor=networks_tensors(network_no,{0,2})*networks_tensors(network_no,{0,3})*networks_tensors(network_no,{1,3})*networks_tensors(network_no,{1,2}),
            upper_right_tensor=networks_tensors(network_no,{3,2})*networks_tensors(network_no,{3,3})*networks_tensors(network_no,{2,3})*networks_tensors(network_no,{2,2}),
            lower_left_tensor=networks_tensors(network_no,{0,1})*networks_tensors(network_no,{0,0})*networks_tensors(network_no,{1,0})*networks_tensors(network_no,{1,1}),
            lower_right_tensor=networks_tensors(network_no,{3,1})*networks_tensors(network_no,{3,0})*networks_tensors(network_no,{2,0})*networks_tensors(network_no,{2,1});


    if (dir==Left) 
    {
        TensorT upper_half_tensor=upper_left_tensor*upper_right_tensor,
                lower_half_tensor=lower_left_tensor*lower_right_tensor;

        outside_indices[0]=networks_indices(network_no,{0,2},{0,1});
        inside_indices[0]=networks_indices(network_no,{1,2},{1,1});
        R[0]=TensorT(outside_indices[0],inside_indices[0]);

        outside_indices[1]=dag(outside_indices[0]);
        inside_indices[1]=dag(inside_indices[0]);
        R[1]=dag(R[0]);

        denmatDecomp(upper_half_tensor,R[0],Q[0],FromRight);
        denmatDecomp(lower_half_tensor,R[1],Q[1],FromRight);
    }
    if (dir==Right)
    {
        TensorT upper_half_tensor=upper_left_tensor*upper_right_tensor,
                lower_half_tensor=lower_left_tensor*lower_right_tensor;

        outside_indices[0]=networks_indices(network_no,{3,1},{3,2});
        inside_indices[0]=networks_indices(network_no,{2,1},{2,2});
        R[0]=TensorT(outside_indices[0],inside_indices[0]);

        outside_indices[1]=dag(outside_indices[0]);
        inside_indices[1]=dag(inside_indices[0]);
        R[1]=dag(R[0]);

        denmatDecomp(lower_half_tensor,R[0],Q[0],FromRight);
        denmatDecomp(upper_half_tensor,R[1],Q[1],FromRight);
    }

    if (dir==Up)
    {
        TensorT left_half_tensor=upper_left_tensor*lower_left_tensor,
                right_half_tensor=upper_right_tensor*lower_right_tensor;

        outside_indices[0]=networks_indices(network_no,{2,3},{1,3});
        inside_indices[0]=networks_indices(network_no,{2,2},{1,2});
        R[0]=TensorT(outside_indices[0],inside_indices[0]);

        outside_indices[1]=dag(outside_indices[0]);
        inside_indices[1]=dag(inside_indices[0]);
        R[1]=dag(R[0]);

        denmatDecomp(right_half_tensor,R[0],Q[0],FromRight);
        denmatDecomp(left_half_tensor,R[1],Q[1],FromRight);
    }
    if (dir==Down)
    {
        TensorT left_half_tensor=upper_left_tensor*lower_left_tensor,
                right_half_tensor=upper_right_tensor*lower_right_tensor;

        outside_indices[0]=networks_indices(network_no,{1,0},{2,0});
        inside_indices[0]=networks_indices(network_no,{1,1},{2,1});
        R[0]=TensorT(outside_indices[0],inside_indices[0]);

        outside_indices[1]=dag(outside_indices[0]);
        inside_indices[1]=dag(inside_indices[0]);
        R[1]=dag(R[0]);

        denmatDecomp(left_half_tensor,R[0],Q[0],FromRight);
        denmatDecomp(right_half_tensor,R[1],Q[1],FromRight);
    }

}
template
void Corner_Transfer_Matrix<ITensor>::obtain_R_from_network(int network_no, Direction dir, std::array<ITensor,2> &R, std::array<Index,2> &outside_indices, std::array<Index,2> &inside_indices);
template
void Corner_Transfer_Matrix<IQTensor>::obtain_R_from_network(int network_no, Direction dir, std::array<IQTensor,2> &R, std::array<IQIndex,2> &outside_indices, std::array<IQIndex,2> &inside_indices);


template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::init_bulk(const std::vector<TensorT> &single_layer_tensors, const std::vector<std::array<IndexT,4>> &ordered_virt_indices)
{
    for (int sitei=0; sitei<N_; sitei++)
    {
        auto temp_bulk_tensor=single_layer_tensors[sitei]*(dag(single_layer_tensors[sitei]).prime(Link));
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
}
template
void Corner_Transfer_Matrix<ITensor>::init_bulk(const std::vector<ITensor> &single_layer_tensors, const std::vector<std::array<Index,4>> &ordered_virt_indices);
template
void Corner_Transfer_Matrix<IQTensor>::init_bulk(const std::vector<IQTensor> &single_layer_tensors, const std::vector<std::array<IQIndex,4>> &ordered_virt_indices);

template <class TensorT>
void Corner_Transfer_Matrix<TensorT> init_edge(const std::vector<TensorT> &single_layer_tensors, const std::vector<std::array<IndexT,4>> &ordered_virt_indices)
{
    for (int sitei=0; sitei<N_; sitei++)
    {
        for (int edgei=0; edgei<4; edgei++)
        {
            //contract ordered_virt_indices[sitei][edgei]
            auto temp_edge_tensor=single_layer_tensors[sitei]*(dag(single_layer_tensors[sitei]).prime(Link).noprime(ordered_virt_indices[sitei][edgei]));
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
}
template
void Corner_Transfer_Matrix<ITensor> init_edge(const std::vector<ITensor> &single_layer_tensors, const std::vector<std::array<Index,4>> &ordered_virt_indices);
template
void Corner_Transfer_Matrix<IQIndex> init_edge(const std::vector<IQIndex> &single_layer_tensors, const std::vector<std::array<IQIndex,4>> &ordered_virt_indices);

template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::init_corner(const std::vector<TensorT> &single_layer_tensors, const std::vector<std::array<IndexT,4>> &ordered_virt_indices)
{
    for (int sitei=0; sitei<N_; sitei++)
    {
        std::array<CombinerT,2> temp_combiners;

        //init C0, leave up and right ind uncontracted
        corner_tensors_[sitei][0]=single_layer_tensors[sitei]*(dag(single_layer_tensors[sitei]).prime(ordered_virt_indices[sitei][1]).prime(ordered_virt_indices[sitei][2]));
        temp_combiners[0]=CombinerT(ordered_virt_indices[sitei][1],dag(ordered_virt_indices[sitei][1]).prime());
        temp_combiners[1]=CombinerT(ordered_virt_indices[sitei][2],dag(ordered_virt_indices[sitei][2]).prime());
        corner_tensors_[sitei][0]=corner_tensors_[sitei][0]*temp_combiners[0]*temp_combiners[1];
        corner_tensors_indices_[sitei][0]=std::array<IndexT,4>{IndexT::Null(),temp_combiners[0].right(),temp_combiners[1].right(),IndexT::Null()};

        //init C1, leave right and down ind uncontracted
        corner_tensors_[sitei][1]=single_layer_tensors[sitei]*(dag(single_layer_tensors).prime(ordered_virt_indices[sitei][2]).prime(ordered_virt_indices[sitei][3]));
        temp_combiners[0]=CombinerT(ordered_virt_indices[sitei][2],dag(ordered_virt_indices[sitei][2]).prime());
        temp_combiners[1]=CombinerT(ordered_virt_indices[sitei][3],dag(ordered_virt_indices[sitei][3]).prime());
        corner_tensors_[sitei][1]=corner_tensors_[sitei][1]*temp_combiners[0]*temp_combiners[1];
        corner_tensors_indices_[sitei][1]=std::array<IndexT,4>{IndexT::Null(),IndexT::Null(),temp_combiners[0].right(),temp_combiners[1].right()};

        //init C2, leave left and down ind uncontracted
        corner_tensors_[sitei][2]=single_layer_tensors[sitei]*(dag(single_layer_tensors[sitei]).prime(ordered_virt_indices[sitei][0]).prime(ordered_virt_indices[sitei][3]));
        temp_combiners[0]=CombinerT(ordered_virt_indices[sitei][0],dag(ordered_virt_indices[sitei][0]).prime());
        temp_combiners[1]=CombinerT(ordered_virt_indices[sitei][3],dag(ordered_virt_indices[sitei][3]).prime());
        corner_tensors_[sitei][2]=corner_tensors_[sitei][2]*temp_combiners[0]*temp_combiners[1];
        corner_tensors_indices_[sitei][2]=std::array<IndexT,4>{temp_combiners[0].right(),IndexT::Null(),IndexT::Null(),temp_combiners[1].right()};

        //init C3, leave left and up ind uncontracted
        corner_tensors_[sitei][3]=single_layer_tensors[sitei]*(dag(single_layer_tensors[sitei]).prime(ordered_virt_indices[sitei][0]).prime(ordered_virt_indices[sitei][1]));
        temp_combiners[0]=CombinerT(ordered_virt_indices[sitei][0],dag(ordered_virt_indices[sitei][0]).prime());
        temp_combiners[1]=CombinerT(ordered_virt_indices[sitei][1],dag(ordered_virt_indices[sitei][1]).prime());
        corner_tensors_[sitei][3]=corner_tensors_[sitei][3]*temp_combiners[0]*temp_combiners[1];
        corner_tensors_indices_[sitei][3]=std::array<IndexT,4>{temp_combiners[0].right(),temp_combiners[1].right(),IndexT::Null(),IndexT::Null()};
    }
}
template
void Corner_Transfer_Matrix<ITensor>::init_corner(const std::vector<ITensor> &single_layer_tensors, const std::vector<std::array<Index,4>> &ordered_virt_indices);
template
void Corner_Transfer_Matrix<IQTensor>::init_corner(const std::vector<IQTensor> &single_layer_tensors, const std::vector<std::array<IQIndex,4>> &ordered_virt_indices);


template <class TensorT>
void Corner_Transfer_Matrix<TensorT>::init_network(int network_no)
{ 
    //init networks_indices
    int x=network_no%Lx_,
        y=network_no/Lx_;
    auto &network_indices=networks_indices_[network_no];
    index_assignment(network_indices[0],corner_tensors_indices(x,y,0)[Right],"network_ind0");
    index_assignment(network_indices[1],edge_tensors_indices(x+1,y,3)[Right],"network_ind1");
    index_assignment(network_indices[2],edge_tensors_indices(x+2,y,3)[Right],"network_ind2");
    index_assignment(network_indices[3],edge_tensors_indices(x,y+1,0)[Right],"network_ind3");
    index_assignment(network_indices[4],bulk_tensors_indices(x+1,y+1)[Right],"network_ind4");
    index_assignment(network_indices[5],bulk_tensors_indices(x+2,y+1)[Right],"network_ind5");
    index_assignment(network_indices[6],edge_tensors_indices(x,y+2,0)[Right],"network_ind6");
    index_assignment(network_indices[7],bulk_tensors_indices(x+1,y+2)[Right],"network_ind7");
    index_assignment(network_indices[8],bulk_tensors_indices(x+2,y+2)[Right],"network_ind8");
    index_assignment(network_indices[9],corner_tensors_indices(x,y+3,1)[Right],"network_ind9");
    index_assignment(network_indices[10],edge_tensors_indices(x+1,y+3,1)[Right],"network_ind10");
    index_assignment(network_indices[11],edge_tensors_indices(x+2,y+3,1)[Right],"network_ind11");
    index_assignment(network_indices[12],corner_tensors_indices(x,y,0)[Up],"network_ind12");
    index_assignment(network_indices[13],edge_tensors_indices(x+1,y,3)[Up],"network_ind13");
    index_assignment(network_indices[14],edge_tensors_indices(x+2,y,3)[Up],"network_ind14");
    index_assignment(network_indices[15],corner_tensors_indices(x+3,y,3)[Up],"network_ind15");
    index_assignment(network_indices[16],edge_tensors_indices(x,y+1,0)[Up],"network_ind16");
    index_assignment(network_indices[17],bulk_tensors_indices(x+1,y+1)[Up],"network_ind17");
    index_assignment(network_indices[18],bulk_tensors_indices(x+2,y+1)[Up],"network_ind18");
    index_assignment(network_indices[19],edge_tensors_indices(x+3,y+1,2)[Up],"network_ind19");
    index_assignment(network_indices[20],edge_tensors_indices(x,y+2,0)[Up],"network_ind20");
    index_assignment(network_indices[21],bulk_tensors_indices(x+1,y+2)[Up],"network_ind21");
    index_assignment(network_indices[22],bulk_tensors_indices(x+2,y+2)[Up],"network_ind22");
    index_assignment(network_indices[23],edge_tensors_indices(x+3,y+2,2)[Up],"network_ind23");

    //init networks_tensors_
    assign_network_tensors(network_no,0,corner_tensors(x,y,0),corner_tensors_indices(x,y,0));
    assign_network_tensors(network_no,1,edge_tensors(x+1,y,3),edge_tensors_indices(x+1,y,3));
    assign_network_tensors(network_no,2,edge_tensors(x+2,y,3),edge_tensors_indices(x+2,y,3));
    assign_network_tensors(network_no,3,corner_tensors(x+3,y,3),corner_tensors_indices(x+3,y,3));
    assign_network_tensors(network_no,4,edge_tensors(x,y+1,0),edge_tensors_indices(x,y+1,0));
    assign_network_tensors(network_no,5,bulk_tensors(x+1,y+1),bulk_tensors_indices(x+1,y+1));
    assign_network_tensors(network_no,6,bulk_tensors(x+2,y+1),bulk_tensors_indices(x+2,y+1));
    assign_network_tensors(network_no,7,edge_tensors(x+3,y+1,2),edge_tensors_indices(x+3,y+1,2));
    assign_network_tensors(network_no,8,edge_tensors(x,y+2,0),edge_tensors_indices(x,y+2,0));
    assign_network_tensors(network_no,9,bulk_tensors(x+1,y+2),bulk_tensors_indices(x+1,y+2));
    assign_network_tensors(network_no,10,bulk_tensors(x+2,y+2),bulk_tensors_indices(x+2,y+2));
    assign_network_tensors(network_no,11,edge_tensors(x+3,y+2,2),edge_tensors_indices(x+3,y+2,2));
    assign_network_tensors(network_no,12,corner_tensors(x,y+3,1),corner_tensors_indices(x,y+3,1));
    assign_network_tensors(network_no,13,edge_tensors(x+1,y+3,1),edge_tensors_indices(x+1,y+3,1));
    assign_network_tensors(network_no,14,edge_tensors(x+2,y+3,1),edge_tensors_indices(x+2,y+3,1));
    assign_network_tensors(network_no,15,corner_tensors(x+3,y+3,2),corner_tensors_indices(x+3,y+3,2));

}
template
void Corner_Transfer_Matrix<ITensor>::init_network(int x, int y);
template
void Corner_Transfer_Matrix<IQTensor>::init_network(int x, int y);

