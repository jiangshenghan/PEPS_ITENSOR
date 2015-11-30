
#ifndef _SINGLET_TENSOR_BASIS_
#define _SINGLET_TENSOR_BASIS_

#include "cgtensor.h"

//Singlet_Tensor_Basis generates all spin singlet basis for a given set of legs which accommodate rep for SU(2) symmetry
class Singlet_Tensor_Basis
{
    public:
        //
        //Type Aliens
        //
        using const_iterator=typename std::vector<IQTensor>::const_iterator;

        //
        //Constructor
        //
        Singlet_Tensor_Basis() {}

        explicit Singlet_Tensor_Basis(const std::vector<IQIndex> &iqinds);

        explicit Singlet_Tensor_Basis(const IndexSet<IQIndex> &iqinds_set);

        //
        //Access Methods
        //
        const IndexSet<IQIndex> &indices() const { return is_; }

        const_iterator begin() const { return singlet_tensors_.begin(); }

        const_iterator end() const { return singlet_tensors_.end(); }

        //number of basis
        int dim() const { return singlet_tensors_.size(); }

        const std::vector<IQTensor> &tensors() const { return singlet_tensors_; }
        const IQTensor &tensor(int i) const { return singlet_tensors_[i]; }

        const std::vector<int> &spin_configs(int i) const { return spin_configs_[i]; }
        const std::vector<int> &flavor_configs(int i) const { return flavor_configs_[i]; }
        int fusion_channel(int i) const { return fusion_channel_[i]; }

        const std::vector<int> &flavor_deg(int legi) const { return is_flavor_deg_[legi]; }

        const IQTensor &operator()(int i) const { return singlet_tensors_[i]; }
        const IQTensor &operator[](int i) const { return singlet_tensors_[i]; }

        //Given a spin_list as well as flavor_list, get the correpsonding No. of base
        int spin_flavor_list_to_basis_no(const std::vector<int> &spin_list, const std::vector<int> &flavor_list, int fusion_channel=0) const
        {
            std::vector<int> spin_set_degs;
            int total_leg_num=is_.r();

            for (int i=0; i<total_leg_num; i++)
            {
                int deg=is_flavor_deg_[i][spin_list[i]];
                spin_set_degs.push_back(deg);
            }

            int spin_list_num=num_from_list(spin_list,max_spins_),
                flavor_list_num=num_from_list(flavor_list,spin_set_degs);

            return spin_flavor_list_to_num_[spin_list_num][flavor_list_num][fusion_channel];
        }

        //
        //Constructor Helper
        //
        void init_spin_deg_and_basis();
        void init_singlet_tensors();

    private:
        IndexSet<IQIndex> is_;
        //is_flavor_deg_[i] stores spin deg for IQIndex is_[i]
        std::vector<std::vector<int>> is_flavor_deg_;
        //max_spins_ stores max spin for each leg
        std::vector<int> max_spins_;
        //is_spin_basis_[i][j] stores spin basis |S,m,t\rangle for is_[i][j]
        std::vector<std::vector<Spin_Basis>> is_spin_basis_;

        //spin_configs_ stores the spin_list for every singlet_tensors_
        std::vector<std::vector<int>> spin_configs_;
        //flavor_configs stores the flavor_list for every singlet_tensors_
        std::vector<std::vector<int>> flavor_configs_;
        //for identical spin and deg configs, there may still be choice of different fusion_channels
        std::vector<int> fusion_channel_;
        //spin_flavor_list_to_num_ translate a particular spin_list, flavor_list and fusion channel to the no. of tensor basis
        //spin_list, flavor_list and fusion channel are encoded as three numbers
        std::vector<std::vector<std::vector<int>>> spin_flavor_list_to_num_;


        std::vector<IQTensor> singlet_tensors_;

};

inline std::ostream &operator<<(std::ostream &s, const Singlet_Tensor_Basis &tensor_basis)
{
    s << "Total dimension: " << tensor_basis.dim() << endl << endl;
    for (int i=0; i<tensor_basis.dim(); i++)
    {
        s << "Basis no: " << i << endl;
        s << "Spin config: " << tensor_basis.spin_configs(i) << endl;
        s << "Deg config: " << tensor_basis.flavor_configs(i) << endl;
        s << "Fusion channel: " << tensor_basis.fusion_channel(i) << endl;
        s << "Check basis no: " << tensor_basis.spin_flavor_list_to_basis_no(tensor_basis.spin_configs(i),tensor_basis.flavor_configs(i),tensor_basis.fusion_channel(i)) << endl;
        s << tensor_basis(i);
    }

    return s;
}


//obtain singlet tensors from params and tensor basis
//This function works for real tensor
IQTensor singlet_tensor_from_basis_params(const Singlet_Tensor_Basis &tensor_basis, const std::vector<double> &params);
IQTensor singlet_tensor_from_basis_params(const Singlet_Tensor_Basis &tensor_basis, const std::vector<Complex> &params);
IQTensor singlet_tensor_from_basis_params(const Singlet_Tensor_Basis &tensor_basis, const arma::Col<double> &params);

//obtain the parameters of singlet tensors from singlet tensors basis
//this function works for real tensor
void obtain_singlet_tensor_params(const IQTensor &singlet_tensor, const Singlet_Tensor_Basis &tensor_basis, std::vector<double> &params);
//this function work for complex tensor
void obtain_singlet_tensor_params(const IQTensor &singlet_tensor, const Singlet_Tensor_Basis &tensor_basis, std::vector<Complex> &params);

//obtain tensors project to spin singlet sector
IQTensor singlet_tensor_from_projection(const IQTensor &tensor_project);
IQTensor singlet_tensor_from_projection(const IQTensor &tensor_project, const Singlet_Tensor_Basis &tensor_basis);

//obtain differ of tensor and its projection to singlet tensors
double diff_tensor_and_singlet_projection(const IQTensor &tensor);
double diff_tensor_and_singlet_projection(const IQTensor &tensor, const Singlet_Tensor_Basis &tensor_basis);


#endif
