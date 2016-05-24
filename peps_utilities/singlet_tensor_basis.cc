
#include "singlet_tensor_basis.h"

Singlet_Tensor_Basis::Singlet_Tensor_Basis(const std::vector<IQIndex> &iqinds):
    is_(iqinds)
{
    //Print(is_);
    init_spin_deg_and_basis();
    init_singlet_tensors();
}

Singlet_Tensor_Basis::Singlet_Tensor_Basis(const IndexSet<IQIndex> &iqinds_set):
    is_(iqinds_set)
{
    init_spin_deg_and_basis();
    init_singlet_tensors();
}


void Singlet_Tensor_Basis::init_spin_deg_and_basis()
{
    for (const auto &sz_leg : is_)
    {
        std::vector<int> flavor_deg;
        std::vector<Spin_Basis> spin_basis;

        bool spin_rep=iqind_spin_rep(sz_leg,flavor_deg);
        iqind_to_spin_basis(sz_leg,flavor_deg,spin_basis);
        assert(spin_rep);

        is_flavor_deg_.push_back(flavor_deg);
        is_spin_basis_.push_back(spin_basis);

        //Print(sz_leg);
        //Print(flavor_deg);
        //Print(spin_basis);
    }

}


void Singlet_Tensor_Basis::init_singlet_tensors()
{
    int total_leg_num=is_.r();

    //get direction for every leg
    std::vector<Arrow> dirs;

    for (const auto &iqind : is_)
    {
        dirs.push_back(iqind.dir());
    }

    
    int total_spin_sets_num=1;

    for (const auto &flavor_deg : is_flavor_deg_)
    {
        total_spin_sets_num*=flavor_deg.size();
        max_spins_.push_back(flavor_deg.size());
    }

    spin_flavor_list_to_num_=std::vector<std::vector<std::vector<int>>>(total_spin_sets_num);

    //iqinds is used to initialize IQTensors
    std::vector<IQIndex> iqinds;
    for (const auto &iqind : is_)
    {
        iqinds.push_back(iqind);
    }
    

    for (int spin_seti=0; spin_seti<total_spin_sets_num; spin_seti++)
    {
        std::vector<int> spin_list;

        spin_list=list_from_num(spin_seti,max_spins_);


        //spin_set_flavors stores degeneracy for this particualr spin set
        //max_sz_nums stores # of different sz's for this particular spin set
        //total_flavors=\prod(every index flavor_deg of this spin_list)
        //total_sz_sets_num=\prod(# of sz qn of every leg for this spin_list)
        //We should avoid the case where this spin set are not occur in the indexset
        std::vector<int> spin_set_flavors, max_sz_nums;
        int total_flavors=1, total_sz_sets_num=1;
        bool valid_spins=true;
        for (int i=0; i<total_leg_num; i++)
        {
            int flavor_deg=is_flavor_deg_[i][spin_list[i]];
            if (flavor_deg==0)
            {
                valid_spins=false;
                total_flavors=0;
                break;
            }
            spin_set_flavors.push_back(flavor_deg);
            max_sz_nums.push_back(spin_list[i]+1);

            total_flavors*=flavor_deg;
            total_sz_sets_num*=spin_list[i]+1;
        }
        if (!valid_spins) continue;

        spin_flavor_list_to_num_[spin_seti]=std::vector<std::vector<int>>(total_flavors);

        CGTensors cg_tensors(spin_list,dirs);

        if (!cg_tensors.valid()) continue;

        //cout << "spin_list:" << endl << spin_list << endl;

        for (const auto &cg_tensor : cg_tensors.K())
        {
            //PrintDat(cg_tensor);

            for (int degs_num=0; degs_num<total_flavors; degs_num++)
            {
                std::vector<int> flavor_list=list_from_num(degs_num,spin_set_flavors);
                //cout << "flavor_list: " << flavor_list << endl;

                IQTensor singlet_tensor(iqinds);

                for (int szs_num=0; szs_num<total_sz_sets_num; szs_num++)
                {
                    std::vector<int> cg_val_list=list_from_num(szs_num,max_sz_nums);
                    //NMAX==8
                    std::vector<IQIndexVal> cg_vals(NMAX,IQIndexVal::Null());
                    std::vector<int> sz_list;

                    for (int legi=0; legi<total_leg_num; legi++)
                    {
                        //Notice we count from 1 for IQIndexVal
                        cg_vals[legi]=cg_tensors.iqindice(legi)(cg_val_list[legi]+1);
                        sz_list.push_back(cg_vals[legi].qn().sz());
                    }

                    auto elem=cg_tensor(cg_vals[0],cg_vals[1],cg_vals[2],cg_vals[3],cg_vals[4],cg_vals[5],cg_vals[6],cg_vals[7]);
                    if (std::abs(elem) < EPSILON) continue;

                    //for (const auto &cg_val : cg_vals) cout << cg_val;
                    //cout << "sz_list: " << sz_list << endl << elem << endl;

                    //obtain the IQIndexVal from spin_list, cg_val_list and flavor_list
                    std::vector<IQIndexVal> singlet_vals(NMAX,IQIndexVal::Null());
                    for (int legi=0; legi<is_.r(); legi++)
                    {
                        int S=spin_list[legi], m=sz_list[legi], t=flavor_list[legi];
                        Spin_Basis spin_state(S,m,t);
                        auto iter=std::find(is_spin_basis_[legi].begin(),is_spin_basis_[legi].end(),spin_state);
                        int val=iter-is_spin_basis_[legi].begin()+1;
                        singlet_vals[legi]=is_[legi](val);
                    }

                    //NMAX==8
                    singlet_tensor(singlet_vals[0],singlet_vals[1],singlet_vals[2],singlet_vals[3],singlet_vals[4],singlet_vals[5],singlet_vals[6],singlet_vals[7])=elem;
                }//end of for loop to get one basis

                singlet_tensors_.push_back(singlet_tensor);
                spin_configs_.push_back(spin_list);
                flavor_configs_.push_back(flavor_list);
                spin_flavor_list_to_num_[spin_seti][degs_num].push_back(singlet_tensors_.size()-1);
                fusion_channel_.push_back(spin_flavor_list_to_num_[spin_seti][degs_num].size()-1);

                //PrintDat(singlet_tensor);

            }//end of for loop for flavors for the same singlet state

        }//end of for loop for different singlet for same spin_list
    }

    //set norm of singlet_tensors_ to be 1 and stores unnormalized tensor to origin_norm_singlet_tensors_;
    origin_norm_singlet_tensors_=singlet_tensors_;
    for (auto &tensor :singlet_tensors_)
        tensor/=tensor.norm();

    //Print(spin_flavor_list_to_num_.size());
    //PrintDat(singlet_tensors_);
}



IQTensor singlet_tensor_from_basis_params(const Singlet_Tensor_Basis &tensor_basis, const std::vector<double> &params)
{
    IQTensor singlet_tensor=params[0]*tensor_basis[0];
    int size=tensor_basis.dim();

    if (size==1) return singlet_tensor;

    for (int i=1; i<size; i++)
    {
        singlet_tensor+=params[i]*tensor_basis[i];
    }

    return singlet_tensor;
}

IQTensor singlet_tensor_from_basis_params(const Singlet_Tensor_Basis &tensor_basis, const std::vector<Complex> &params)
{
    IQTensor singlet_tensor=params[0]*tensor_basis[0];
    int size=tensor_basis.dim();

    if (size==1) return singlet_tensor;

    for (int i=1; i<size; i++)
    {
        singlet_tensor+=params[i]*tensor_basis[i];
    }

    return singlet_tensor;
}

IQTensor singlet_tensor_from_basis_params(const Singlet_Tensor_Basis &tensor_basis, const arma::Col<double> &params)
{
    IQTensor singlet_tensor=params[0]*tensor_basis[0];
    int size=tensor_basis.dim();

    if (size==1) return singlet_tensor;

    for (int i=1; i<size; i++)
    {
        singlet_tensor+=params[i]*tensor_basis[i];
    }

    return singlet_tensor;
}


void obtain_singlet_tensor_params(const IQTensor &singlet_tensor, const Singlet_Tensor_Basis &tensor_basis, std::vector<double> &params)
{
    params.clear();

    for (const auto &base : tensor_basis.tensors())
    {
        params.push_back((singlet_tensor*dag(base)).toReal());
    }

}

void obtain_singlet_tensor_params(const IQTensor &singlet_tensor, const Singlet_Tensor_Basis &tensor_basis, std::vector<Complex> &params)
{
    params.clear();

    for (const auto &base : tensor_basis.tensors())
    {
        params.push_back((singlet_tensor*dag(base)).toComplex());
    }

}


IQTensor singlet_tensor_from_projection(const IQTensor &tensor_project)
{
    Singlet_Tensor_Basis tensor_basis(tensor_project.indices());
    return singlet_tensor_from_projection(tensor_project,tensor_basis);
}

IQTensor singlet_tensor_from_projection(const IQTensor &tensor_project, const Singlet_Tensor_Basis &tensor_basis)
{
    if (tensor_project.isComplex())
    {
        std::vector<Complex> params;
        obtain_singlet_tensor_params(tensor_project,tensor_basis,params);
        return singlet_tensor_from_basis_params(tensor_basis,params);
    }
    else
    {
        std::vector<double> params;
        obtain_singlet_tensor_params(tensor_project,tensor_basis,params);
        return singlet_tensor_from_basis_params(tensor_basis,params);
    }
}

double diff_tensor_and_singlet_projection(const IQTensor &tensor)
{
    Singlet_Tensor_Basis tensor_basis(tensor.indices());
    return diff_tensor_and_singlet_projection(tensor,tensor_basis);
}

double diff_tensor_and_singlet_projection(const IQTensor &tensor, const Singlet_Tensor_Basis &tensor_basis)
{
    auto tensor_singlet=singlet_tensor_from_projection(tensor,tensor_basis);
    return (tensor-tensor_singlet).norm();
}


IQTensor spin_fusion_tensor(std::vector<IQIndex> input_spin_legs, Arrow spin_dir)
{
    //change the direction of input legs, s.t. same as spin_dir
    std::vector<IQTensor> tensors_leg_dir;
    for (const auto &spin_leg : input_spin_legs)
    {
        if (input_spin_legs.dir()!=spin_dir)
        {
            IQTensor temp_tensor(dag(spin_leg),prime(dag(spin_leg)));
            for (int val=0; val<spin_leg.m(); val++) temp_tensor(dag(spin_leg)(val+1),prime(dag(spin_leg))(val+1))=1-2*(val%2);
            tensors_leg_dir.push_back(temp_tensor);
            spin_leg.dag().prime();
        }
    }
    
    //get possible final spin qn's
    int Smax=0;
    std::vector<std::vector<int>> input_spin_reps(input_spin_legs.size());
    for (int inputi=0; inputi<input_spin_legs.size(); inputi++) 
    {
        iqind_spin_rep(input_spin_legs[inputi],input_spin_reps[inputi]);
        Smax+=input_spin_reps[inputi].size()-1;
    }

    //get all possible fusion channel
    std::vector<IQIndex> legs_form_singlet;
    std::vector<IQIndex> fixed_spin_legs;
    std::vector<IQTensor> fusion_tensor_basis;
    for (const auto &spin_leg: input_spin_legs) legs_form_singlet.push_back(dag(spin_leg));
    std::vector<int> output_spin_rep;
    for (int spini=0; spini<Smax+1; spini++)
    {
        std::vector<int> spin_qn(spini+1,0);
        spin_qn[spini]=1;
        fixed_spin_legs.push_back(Spin_leg(spin_qn,nameint("fixed_spin",spini)+"_legs",spin_dir));
        legs_form_singlet.push_back(fixed_spin_legs);
        fusion_tensor_basis.push_back(Singlet_Tensor_Basis(legs_form_singlet));
        output_spin_rep.push_back(fusion_tensor_basis.dim());
        legs_form_singlet.pop_back();
    }


    //obtain the fusion_tensor
    IQTensor output_spin_leg=Spin_leg(output_spin_rep,"output_spin_leg",spin_dir);
    std::vector<Spin_Basis> output_spin_basis;
    iqind_to_spin_basis(output_spin_leg,flavor_deg,output_spin_basis);

    legs_form_singlet.push_back(output_spin_leg);
    IQTensor fusion_tensor(legs_form_singlet);
    std::vector<int> input_legs_dim;
    int input_total_dim=1;
    for (const auto &spin_leg: input_spin_legs) 
    {
        input_legs_dim.push_back(spin_leg.m());
        input_total_dim*=spin_leg.m();
    }

    for (int ival=0; ival<input_total_dim; ival++)
    {
        std::vector<int> ival_list=list_from_num(ival,input_legs_dim);
        for (int oval=0; oval<output_spin_leg.m(); oval++)
        {
            std::array<IQIndexVal,8> fusion_legs_val, basis_legs_val;
            for (int legi=0; legi<ival_list.size(); legi++) 
            {
                fusion_legs_val[legi]=input_spin_legs[legi](ival_list+1);
                basis_legs_val[legi]=input_spin_legs[legi](ival_list+1);
            }
            fusion_legs_val[ival_list.size()+1]=output_spin_leg(oval+1);
            int oS=output_spin_basis[oval].S(),
                om=output_spin_basis[oval].m(),
                ot=output_spin_basis[oval].t();
            basis_legs_val[ival_list.size()+1]=fixed_spin_legs[oS][(oS-om)/2];

            fusion_tensor(fusion_legs_val[0],fusion_legs_val[1],fusion_legs_val[2],fusion_legs_val[3],fusion_legs_val[4],fusion_legs_val[5],fusion_legs_val[6],fusion_legs_val[7])=fusion_tensor_basis[oS][ot](basis_legs_val[0],basis_legs_val[1],basis_legs_val[2],basis_legs_val[3],basis_legs_val[4],basis_legs_val[5],basis_legs_val[6],basis_legs_val[7]);
        }
    }

    for (const auto &dir_tensor: tensors_leg_dir) fusion_tensor*=dir_tensor;

    return fusion_tensor;
}
