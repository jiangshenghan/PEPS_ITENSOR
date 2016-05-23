
#include "tensor_svd.h"
#include "tensor_rg.h"
#include <limits>

template <class TensorT>
TensorT_RG<TensorT>::TensorT_RG(const Lattice_Base &lattice, const std::vector<TensorT> &input_tensors, int maxm):
    lattice_(lattice),
    input_tensors_(input_tensors),
    factor_input_tensors_(input_tensors.size(),std::vector<TensorT>(2)),
    //factor_args_("Maxm",maxm,"ShowEigs"),
    factor_args_("Maxm",maxm),
    iszero_(false),
    almost_zero_trg_(false)
{
    N_layer_=std::round(std::log2(lattice_.n_uc()[0]*lattice_.n_uc()[1]*1./4.))+1;
    //TODO: consider case where we get zero tensor 
    init_to_square_network();

    //obtain layered trg_tensors and factor_tensors
    for (int layeri=1; layeri<N_layer_; layeri++)
    {
        //obtain factor tensors of layeri-1
        if (iszero_==true) 
        {
            break;
            Print(layeri);
        }
        for (int tensori=0; tensori<layered_trg_tensors_[layeri-1].size(); tensori++) obtain_factor_tensors(layeri-1,tensori);
        for (int tensori=0; tensori<layered_trg_tensors_[layeri].size(); tensori++) 
        {
            obtain_trg_tensor(layeri,tensori);
        }
    }

    obtain_trg_result();
}
template
TensorT_RG<ITensor>::TensorT_RG(const Lattice_Base &lattice, const std::vector<ITensor> &input_tensors, int maxm);
template
TensorT_RG<IQTensor>::TensorT_RG(const Lattice_Base &lattice, const std::vector<IQTensor> &input_tensors, int maxm);


template <class TensorT>
void TensorT_RG<TensorT>::init_to_square_network()
{
    std::vector<int> lattice_dim=std::vector<int>{lattice_.n_uc()[0],lattice_.n_uc()[1]};
    int N_tensors=lattice_dim[0]*lattice_dim[1];
    for (int layeri=0; layeri<N_layer_; layeri++)
    {
        layered_lattice_dim_.push_back(lattice_dim);
        layered_trg_tensors_.push_back(std::vector<TensorT>(N_tensors));
        layered_trg_tensors_norm_.push_back(std::vector<double>(N_tensors));
        layered_factor_tensors_.push_back(std::vector<std::vector<TensorT>>(N_tensors,std::vector<TensorT>(2))); 
        if ((layeri+1)%2==0) { lattice_dim[0]/=2; lattice_dim[1]/=2; }
        N_tensors/=2;
    }

    //Print(N_layer_);
    //Print(layered_lattice_dim_);

    //factor input tensors
    for (int tensori=0; tensori<input_tensors_.size(); tensori++) obtain_factor_input_tensor(tensori);
    //obtain trg tensor of the 0th layer
    for (int tensori=0; tensori<layered_trg_tensors_[0].size(); tensori++) obtain_zeroth_layer_trg_tensor(tensori); 
}
template
void TensorT_RG<ITensor>::init_to_square_network();
template
void TensorT_RG<IQTensor>::init_to_square_network();


template <class TensorT>
void TensorT_RG<TensorT>::obtain_factor_input_tensor(int input_tensor_no)
{
    //we do not need to factor input tensor of square network
    if (lattice_.name().find("square")!=std::string::npos) return;
    if (lattice_.name().find("kagome normal")!=std::string::npos)
    {
        auto tensor_coord=lattice_.site_list_to_coord(input_tensor_no);
        std::vector<int> neigh_tensors_no;
        neigh_tensors_no.push_back(lattice_.site_coord_to_list(tensor_coord[0],tensor_coord[1],(tensor_coord[2]+1)%3));
        neigh_tensors_no.push_back(lattice_.site_coord_to_list(tensor_coord[0],tensor_coord[1],(tensor_coord[2]+2)%3));

        std::vector<IndexT> tensorA_indices;
        for (auto neigh_no: neigh_tensors_no) tensorA_indices.push_back(commonIndex(input_tensors_[input_tensor_no],input_tensors_[neigh_no]));
        factor_input_tensors_[input_tensor_no][0]=TensorT(tensorA_indices);

        tensor_factor(input_tensors_[input_tensor_no],factor_input_tensors_[input_tensor_no][0],factor_input_tensors_[input_tensor_no][1],{factor_args_,"IndexName",nameint("leg",input_tensor_no)});

        //Print(input_tensor_no);
        //Print(tensor_coord);
        //Print(neigh_tensors_no);
        //Print(factor_input_tensors_[input_tensor_no][0].indices());
        //Print(factor_input_tensors_[input_tensor_no][1].indices());
    }
}
template
void TensorT_RG<ITensor>::obtain_factor_input_tensor(int input_tensor_no);
template
void TensorT_RG<IQTensor>::obtain_factor_input_tensor(int input_tensor_no);

template <class TensorT>
void TensorT_RG<TensorT>::obtain_zeroth_layer_trg_tensor(int tensor_no)
{
    if (lattice_.name().find("square")!=std::string::npos)
    {
        layered_trg_tensors_[0][tensor_no]=input_tensors_[tensor_no];
    }
    if (lattice_.name().find("kagome normal")!=std::string::npos)
    {
        int ix=tensor_no%lattice_.n_uc()[0],
            iy=tensor_no/lattice_.n_uc()[0];
        std::vector<int> utensors_no, dtensors_no;
        for (int subi=0; subi<lattice_.n_sites_uc(); subi++) utensors_no.push_back(lattice_.site_coord_to_list(ix,iy,subi));
        dtensors_no.push_back(lattice_.site_coord_to_list(ix,iy,0));
        dtensors_no.push_back(lattice_.site_coord_to_list(ix,iy-1,1));
        dtensors_no.push_back(lattice_.site_coord_to_list(ix-1,iy,2));

        //Print(utensors_no);
        //Print(dtensors_no);

        TensorT utensor=factor_input_tensors_[utensors_no[0]][0],
                dtensor=factor_input_tensors_[dtensors_no[0]][1];
        for (int subi=1; subi<lattice_.n_sites_uc(); subi++)
        {
            utensor*=factor_input_tensors_[utensors_no[subi]][0];
            dtensor*=factor_input_tensors_[dtensors_no[subi]][1];
        }

        layered_trg_tensors_[0][tensor_no]=utensor*dtensor;
    }
    if (layered_trg_tensors_[0][tensor_no].norm()<std::numeric_limits<double>::min()) 
    {
        iszero_=true;
        //Print(tensor_no);
        //Print(layered_trg_tensors_[0][tensor_no].norm());
    }
    else
    {
        layered_trg_tensors_norm_[0][tensor_no]=layered_trg_tensors_[0][tensor_no].norm();
        layered_trg_tensors_[0][tensor_no]/=layered_trg_tensors_norm_[0][tensor_no];
    }
}
template
void TensorT_RG<ITensor>::obtain_zeroth_layer_trg_tensor(int tensor_no);
template
void TensorT_RG<IQTensor>::obtain_zeroth_layer_trg_tensor(int tensor_no);


//TODO: divide trg tensor by some factor to avoid too large or too small singular value
template <class TensorT>
void TensorT_RG<TensorT>::obtain_trg_tensor(int layer_no, int tensor_no)
{
    const auto &factor_tensors=layered_factor_tensors_[layer_no-1];
    std::vector<int> surr_tensor_inds;
    //we should treat even layer and odd layer separately
    //for odd layer, the lattice_dim is the same as last layer
    //tensors of this layer sit at even plaquette centers of the lattice
    if (layer_no%2==1)
    {
        for (int surri=0; surri<4; surri++) surr_tensor_inds.push_back(site_ind_surr_plaq(plaq_ind_from_odd_layer_tensor_ind(tensor_no,layer_no),surri,layer_no-1));
    }
    //for even dim, the lattice_dim[i] is 1/2 of the last layer
    //assuming the tensor_no coord of this layer is (x,y), then it sits at plaq coord of last layer (2*x,2*y+1)
    //tensors of last layer sits at plaquette center around
    else
    {
        std::vector<int> tensor_no_coord=coord_from_ind(tensor_no,layered_lattice_dim_[layer_no]);
        //plaquette coord of last layer
        int plaq_ind=ind_from_coord({tensor_no_coord[0]*2,tensor_no_coord[1]*2+1},layered_lattice_dim_[layer_no-1]);
        for (int surri=0; surri<4; surri++) surr_tensor_inds.push_back(plaq_ind_surr_plaq(plaq_ind,surri,layer_no-1)/2);
    }
    layered_trg_tensors_[layer_no][tensor_no]=factor_tensors[surr_tensor_inds[0]][1]*factor_tensors[surr_tensor_inds[1]][0]*factor_tensors[surr_tensor_inds[2]][0]*factor_tensors[surr_tensor_inds[3]][1];
     
    if (layered_trg_tensors_[layer_no][tensor_no].norm()<std::numeric_limits<double>::min()) iszero_=true;
    else
    {
        layered_trg_tensors_norm_[layer_no][tensor_no]=layered_trg_tensors_[layer_no][tensor_no].norm();
        layered_trg_tensors_[layer_no][tensor_no]/=layered_trg_tensors_norm_[layer_no][tensor_no];
    }

    //Print(layer_no);
    //Print(tensor_no);
    //Print(layered_trg_tensors_[layer_no][tensor_no].norm());
    //Print(surr_tensor_inds);
    //Print(factor_tensors[surr_tensor_inds[0]][1]);
    //Print(factor_tensors[surr_tensor_inds[1]][0]);
    //Print(factor_tensors[surr_tensor_inds[2]][0]);
    //Print(factor_tensors[surr_tensor_inds[3]][1]);
    //Print(layered_trg_tensors_[layer_no][tensor_no]);
}
template
void TensorT_RG<ITensor>::obtain_trg_tensor(int layer_no, int tensor_no);
template
void TensorT_RG<IQTensor>::obtain_trg_tensor(int layer_no, int tensor_no);


template <class TensorT>
void TensorT_RG<TensorT>::obtain_factor_tensors(int layer_no, int tensor_no)
{
    //we do not need factor tensor of last layer
    if (layer_no==N_layer_-1) return;
    //get factor inds of first factor tensor
    std::vector<IndexT> factor_inds;
    std::vector<int> neighs(2);
    if (layer_no%2==0) 
    {
        neighs[0]=((tensor_no/layered_lattice_dim_[layer_no][0]+tensor_no%layered_lattice_dim_[layer_no][1])%2==0)?0:2;
        neighs[1]=3;
    }
    else
    {
        neighs[0]=0;
        neighs[1]=(((tensor_no*2)/layered_lattice_dim_[layer_no][0])%2==0)? 3:1;
    }

    factor_inds.push_back(commonIndex(layered_trg_tensors_[layer_no][tensor_no],layered_trg_tensors_[layer_no][neigh_tensor_ind(tensor_no,neighs[0],layer_no)]));
    factor_inds.push_back(commonIndex(layered_trg_tensors_[layer_no][tensor_no],layered_trg_tensors_[layer_no][neigh_tensor_ind(tensor_no,neighs[1],layer_no)]));

    layered_factor_tensors_[layer_no][tensor_no][0]=TensorT(factor_inds[0],factor_inds[1]);
    tensor_factor(layered_trg_tensors_[layer_no][tensor_no],layered_factor_tensors_[layer_no][tensor_no][0],layered_factor_tensors_[layer_no][tensor_no][1],{factor_args_,"IndexName",nameint("leg",tensor_no)});

    //Print(layer_no);
    //Print(tensor_no);
    //Print(neighs);
    //Print(factor_inds);
    //Print(neigh_tensor_ind(tensor_no,0,layer_no));
    //Print(neigh_tensor_ind(tensor_no,1,layer_no));
    //Print(neigh_tensor_ind(tensor_no,2,layer_no));
    //Print(neigh_tensor_ind(tensor_no,3,layer_no));
    //Print((layered_trg_tensors_[layer_no][tensor_no]-layered_factor_tensors_[layer_no][tensor_no][0]*layered_factor_tensors_[layer_no][tensor_no][1]).norm());
}
template
void TensorT_RG<ITensor>::obtain_factor_tensors(int layer_no, int tensor_no);
template
void TensorT_RG<IQTensor>::obtain_factor_tensors(int layer_no, int tensor_no);

template <class TensorT>
void TensorT_RG<TensorT>::obtain_trg_result()
{
    //Print(iszero_);
    if (iszero_==true) 
    {
        trg_result_=0;
        return;
    }
    const auto &trg_tensors=layered_trg_tensors_.back();
    auto result_tensor=(trg_tensors[0]*trg_tensors[1])*(trg_tensors[2]*trg_tensors[3]);
    trg_result_=result_tensor.toComplex();
    for (int layer_no=0; layer_no<layered_trg_tensors_.size(); layer_no++)
    {
        for (double tensor_norm: layered_trg_tensors_norm_[layer_no]) trg_result_*=tensor_norm;
    }
    //the case where trg_result_ is too small, we need to renormalize tensor
    if (std::abs(trg_result_)<std::numeric_limits<double>::min()) 
    {
        almost_zero_trg_=true;
    }
    //Print(trg_result_);
}
template
void TensorT_RG<ITensor>::obtain_trg_result();
template
void TensorT_RG<IQTensor>::obtain_trg_result();


template <class TensorT>
void TensorT_RG<TensorT>::update_trg_network(const std::vector<int> &input_inds, const std::vector<TensorT> &update_input_tensors)
{
    //update input tensors as well as zeroth layer
    for (int inputi=0; inputi<input_inds.size(); inputi++)
    {
        input_tensors_[input_inds[inputi]]=update_input_tensors[inputi];
        obtain_factor_input_tensor(input_inds[inputi]);
    }
    std::vector<int> update_inds_curr;
    obtain_zeroth_layer_update_inds(update_inds_curr,input_inds);
    //Print(input_inds);
    //Print(update_inds_curr);
    for (int tensor_no: update_inds_curr) obtain_zeroth_layer_trg_tensor(tensor_no);
    //update other layers
    std::vector<int> update_inds_last=update_inds_curr;
    for (int layeri=1; layeri<N_layer_; layeri++)
    {
        if (iszero_==true) break;
        //get updated factor tensors of last layer
        obtain_update_inds(update_inds_curr,update_inds_last,layeri);
        //Print(update_inds_last);
        //Print(update_inds_curr);
        for (int tensor_no: update_inds_last) obtain_factor_tensors(layeri-1,tensor_no);
        for (int tensor_no: update_inds_curr) obtain_trg_tensor(layeri,tensor_no);
        update_inds_last=update_inds_curr;
    }

    obtain_trg_result();
}
template
void TensorT_RG<ITensor>::update_trg_network(const std::vector<int> &input_inds, const std::vector<ITensor> &update_input_tensors);
template
void TensorT_RG<IQTensor>::update_trg_network(const std::vector<int> &input_inds, const std::vector<IQTensor> &update_input_tensors);


template <class TensorT>
void TensorT_RG<TensorT>::obtain_zeroth_layer_update_inds(std::vector<int> &update_inds, const std::vector<int> &input_inds)
{
    if (lattice_.name().find("square")!=std::string::npos) update_inds=input_inds;
    if (lattice_.name().find("kagome normal")!=std::string::npos) 
    {
        std::vector<bool> update_status(layered_trg_tensors_[0].size(),false);
        for (int ind: input_inds)
        {
            auto input_coord=lattice_.site_list_to_coord(ind);
            std::vector<int> update_sites;
            update_sites.push_back(ind_from_coord({input_coord[0],input_coord[1]},layered_lattice_dim_[0]));
            if (input_coord[2]==1) update_sites.push_back(ind_from_coord({input_coord[0],input_coord[1]+1},layered_lattice_dim_[0]));
            if (input_coord[2]==2) update_sites.push_back(ind_from_coord({input_coord[0]+1,input_coord[1]},layered_lattice_dim_[0]));
            for (int sitei: update_sites) update_status[sitei]=true;
        }
        for (int ind=0; ind<update_status.size(); ind++)
        {
            if (update_status[ind]==true) update_inds.push_back(ind);
        }
    }
    //Print(input_inds);
    //Print(update_inds);
}
template
void TensorT_RG<ITensor>::obtain_zeroth_layer_update_inds(std::vector<int> &update_inds, const std::vector<int> &input_inds);
template
void TensorT_RG<IQTensor>::obtain_zeroth_layer_update_inds(std::vector<int> &update_inds, const std::vector<int> &input_inds);

template <class TensorT>
void TensorT_RG<TensorT>::obtain_update_inds(std::vector<int> &update_inds_curr, const std::vector<int> &update_inds_last, int layer_no)
{
    update_inds_curr.clear();
    std::vector<bool> update_status(layered_trg_tensors_[layer_no].size(),false);
    if (layer_no%2==1)
    {
        for (int ind: update_inds_last)
        {
            std::vector<int> plaq_coord=coord_from_ind(ind,layered_lattice_dim_[layer_no-1]);
            update_status[plaq_ind_surr_site(ind,(plaq_coord[0]+plaq_coord[1])%2,layer_no-1)/2]=true;
            update_status[plaq_ind_surr_site(ind,(plaq_coord[0]+plaq_coord[1])%2+2,layer_no-1)/2]=true;
        }
    }
    else
    {
        for (int ind : update_inds_last)
        {
            std::vector<int> plaq_coord=coord_from_ind(plaq_ind_from_odd_layer_tensor_ind(ind,layer_no-1),layered_lattice_dim_[layer_no-1]);
            std::vector<int> update_tensor_inds;
            if (plaq_coord[1]%2==0)
            {
                update_tensor_inds.push_back(ind_from_coord({plaq_coord[0]/2,plaq_coord[1]/2},layered_lattice_dim_[layer_no]));
                update_tensor_inds.push_back(ind_from_coord({plaq_coord[0]/2,plaq_coord[1]/2-1},layered_lattice_dim_[layer_no]));
            }
            else
            {
                update_tensor_inds.push_back(ind_from_coord({plaq_coord[0]/2,plaq_coord[1]/2},layered_lattice_dim_[layer_no]));
                update_tensor_inds.push_back(ind_from_coord({plaq_coord[0]/2+1,plaq_coord[1]/2},layered_lattice_dim_[layer_no]));
            }
            update_status[update_tensor_inds[0]]=true;
            update_status[update_tensor_inds[1]]=true;
        }
    }
    for (int ind=0; ind<update_status.size(); ind++)
    {
        if (update_status[ind]==true) update_inds_curr.push_back(ind);
    }
}
template
void TensorT_RG<ITensor>::obtain_update_inds(std::vector<int> &update_inds_curr, const std::vector<int> &update_inds_last, int layer_no);
template
void TensorT_RG<IQTensor>::obtain_update_inds(std::vector<int> &update_inds_curr, const std::vector<int> &update_inds_last, int layer_no);

