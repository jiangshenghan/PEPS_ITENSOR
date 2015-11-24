
#ifndef _TROTTER_H_
#define _TROTTER_H_

#include "singlet_tensor_basis.h"

//class for two-site Heisenberg Hamiltonian
class NN_Heisenberg_Hamiltonian
{
    public:
        //
        //Constructor
        //
        NN_Heisenberg_Hamiltonian(std::array<IQIndex,2> nn_sites, double J=1.) : J_(J)
        {
            for (int sitei=0; sitei<2; sitei++)
            {
                phys_legs_[sitei][0]=dag(nn_sites[sitei]);
                phys_legs_[sitei][1]=prime(nn_sites[sitei]);

                virt_legs_[sitei]=Spin_leg({0,0,1},nameint("virt_leg",sitei),Out);

                site_tensors_[sitei]=IQTensor(phys_legs_[sitei][0],phys_legs_[sitei][1],virt_legs_[sitei]);
            }
            bond_tensor_=IQTensor(dag(virt_legs_[0]),dag(virt_legs_[1]));

            //cout << "Site tensors: " << endl << site_tensors_[0] << site_tensors_[1];
            //cout << "Bond tensors: " << endl << bond_tensors_[0];

            init_site_tensors();
            init_bond_tensor();
        }

        //
        //Access Method
        //
        double J() const { return J_; }
        const IQTensor &site_tensors(int i) const { return site_tensors_[i]; }
        const IQTensor &bond_tensor() const { return bond_tensor_; }

        //
        //Constructor Helper
        //
        void init_site_tensors()
        {
            for (auto &tensor : site_tensors_)
            {
                Singlet_Tensor_Basis tensor_basis(tensor.indices());
                tensor=std::sqrt(6.)/2.*tensor_basis[0];
                //PrintDat(tensor);
            }
        }

        void init_bond_tensor()
        {
            Singlet_Tensor_Basis bond_tensor_basis(bond_tensor_.indices());
            bond_tensor_=-J_*std::sqrt(3.)*bond_tensor_basis.tensor(0);
        }

    private:
        double J_;
        
        std::array<std::array<IQIndex,2>,2> phys_legs_;
        std::array<IQIndex,2> virt_legs_;
        
        std::array<IQTensor,2> site_tensors_;
        IQTensor bond_tensor_;

};

//This struct stores information about time evolution of simple update
struct Evolution_Params
{
    public:
        Evolution_Params() {}
        Evolution_Params(int n_iter): iter_nums(n_iter) {}
        Evolution_Params(int n_iter, std::vector<int> n_steps, std::vector<double> imag_t) : 
            iter_nums(n_iter), steps_nums(n_steps), ts(imag_t)
        {
            assert(n_iter==step_nums.size());
            assert(n_iter==ts.size());
        }

        int iter_nums;
        std::vector<int> steps_nums;
        std::vector<double> ts;
};


//class for general trotter_gate
class Trotter_Gate
{
    public:
        //
        //Constructor
        //
        Trotter_Gate() {}
        
        Trotter_Gate(int site_num, int bond_num, double t=1.): 
            t_(t),
            phys_legs_(site_num),
            virt_legs_(site_num),
            site_tensors_(site_num), 
            bond_tensors_(bond_num) 
        {}

        //
        //Access Method
        //
        double t() { return t_; }

        const std::vector<IQTensor> &site_tensors() const { return site_tensors_; }
        const IQTensor &site_tensors(int i) const { return site_tensors_[i]; }

        const std::vector<IQTensor> &bond_tensors() const { return bond_tensors_; }
        const IQTensor &bond_tensors(int i) const { return bond_tensors_[i]; }

        //
        //Other Method
        //
        //Change time t_, and modify corresponding gates
        virtual void change_time(double t) = 0;


    protected:
        //t_ is the imaginary time
        double t_;

        std::vector<std::array<IQIndex,2>> phys_legs_;
        std::vector<std::vector<IQIndex>> virt_legs_;

        std::vector<IQTensor> site_tensors_;
        std::vector<IQTensor> bond_tensors_;
};

//
//class for trotter gate for NN Heisenberg Hamiltonian exp(-tH)~1-t\vec{S}_i\cdot\vec{S}_j
//
class NN_Heisenberg_Trotter_Gate : public Trotter_Gate
{
    public:
        //
        //Constructors
        //
        NN_Heisenberg_Trotter_Gate() {}

        NN_Heisenberg_Trotter_Gate(std::array<IQIndex,2> nn_sites, double t=1.);

        //
        //Other Method
        //
        //change time
        virtual void change_time(double t)
        {
            t_=t;
            init_bond_tensor();
        }


        //
        //Constructor Helpers
        //
        void init_site_tensors();
        void init_bond_tensor();

    private:
        //
        //Data Member
        //
};

inline NN_Heisenberg_Trotter_Gate::NN_Heisenberg_Trotter_Gate(std::array<IQIndex,2> nn_sites, double t): 
    Trotter_Gate(2,1,t)
{
    for (int sitei=0; sitei<2; sitei++)
    {
        phys_legs_[sitei][0]=dag(nn_sites[sitei]);
        phys_legs_[sitei][1]=prime(nn_sites[sitei]);

        virt_legs_[sitei].push_back(Spin_leg({1,0,1},nameint("virt_leg",sitei),Out));

        site_tensors_[sitei]=IQTensor(phys_legs_[sitei][0],phys_legs_[sitei][1],virt_legs_[sitei][0]);
    }
    bond_tensors_[0]=IQTensor(dag(virt_legs_[0][0]),dag(virt_legs_[1][0]));

    //cout << "Site tensors: " << endl << site_tensors_[0] << site_tensors_[1];
    //cout << "Bond tensors: " << endl << bond_tensors_[0];

    init_site_tensors();
    init_bond_tensor();
}


inline void NN_Heisenberg_Trotter_Gate::init_site_tensors()
{
    for (auto &tensor : site_tensors_)
    {
        Singlet_Tensor_Basis tensor_basis(tensor.indices());

        tensor+=tensor_basis[0]+std::sqrt(6.)/2.*tensor_basis[1];

        //PrintDat(tensor);
    }
}


inline void NN_Heisenberg_Trotter_Gate::init_bond_tensor()
{
    Singlet_Tensor_Basis bond_tensor_basis(bond_tensors_[0].indices());
    bond_tensors_[0]=bond_tensor_basis.tensor(0)+t_*std::sqrt(3.)*bond_tensor_basis.tensor(1);
}


#endif
