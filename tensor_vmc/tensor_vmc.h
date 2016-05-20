
#ifndef _TENSOR_VMC_H_
#define _TENSOR_VMC_H_

#include "peps.h"
#include "tensor_rg.h"
#include "mpi.h"

//using variational MC combined with tensor rg to get expectation value
//
//class TensorT_VMC_WF
//provides tensor wavefunctions for vmc. the weight for current spin config is obtained by trg
//Here we focus on physical leg d=2
//
template <class TensorT>
class TensorT_VMC_WF
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
        TensorT_VMC_WF(const std::vector<int> init_spin_config, const PEPSt<TensorT> &peps, int maxm=100);

        //
        //Acess Methods
        //
        const Lattice_Base &lattice() const { return tensor_rg_.lattice(); }
        int n_sites() const { return spin_config_.size(); }
        int n_bonds() const { return this->lattice().n_bonds_total(); }
        int spin_config(int sitei) const { return spin_config_[sitei]; }
        Complex wf_weight() const { return tensor_rg_.trg_result(); }
        bool is_zero() const { return tensor_rg_.is_zero(); }

        //Update methods
        //we reverse the spins of flip inds and update tensor_wf
        void update_wf(std::vector<int> flip_inds);

    private:
        std::vector<int> spin_config_;
        std::vector<TensorT> spin_prod_wf_;
        //the combined tensors should have one to one correspondance to lattice
        std::vector<TensorT> combined_tensors_;
        TensorT_RG<TensorT> tensor_rg_;
};


//option for measure_args:
//getInt: ThermalSteps, MeasureSteps, Bin_no, Maxm
//getString: Operator(SzSz, Heisenberg), InitSpins(random,antiferro)
template <class TensorT> 
void tensor_vmc(const PEPSt<TensorT> &peps, const Args &measure_args);
//vmc for parallel computing
template <class TensorT>
void tensor_vmc_parallel(const PEPSt<TensorT> &peps, const Args &measure_args);



//flip spin once and update corresponding wf
template <class TensorT>
void vmc_one_step(TensorT_VMC_WF<TensorT> &tensor_vmc_wf)
{
    std::vector<int> flip_inds(2);
    TensorT_VMC_WF<TensorT> tensor_vmc_wf_update=tensor_vmc_wf;

    //get flipped spin until two spins antiparallel
    do
    {
        flip_inds[0]=floor((rand_gen()+1)/2.*tensor_vmc_wf.n_sites());
        flip_inds[1]=floor((rand_gen()+1)/2.*tensor_vmc_wf.n_sites());
        //Print(flip_inds);
        //Print(tensor_vmc_wf.spin_config(flip_inds[0]));
        //Print(tensor_vmc_wf.spin_config(flip_inds[1]));
    }
    while (tensor_vmc_wf.spin_config(flip_inds[0])==tensor_vmc_wf.spin_config(flip_inds[1]));
    //Print(flip_inds);
    //Print(tensor_vmc_wf.spin_config(flip_inds[0]));
    //Print(tensor_vmc_wf.spin_config(flip_inds[1]));
    tensor_vmc_wf_update.update_wf(flip_inds);

    double flip_prob=std::pow(std::abs(tensor_vmc_wf_update.wf_weight()/tensor_vmc_wf.wf_weight()),2.);
    //Print(tensor_vmc_wf.wf_weight());
    //Print(tensor_vmc_wf_update.wf_weight());
    //Print(flip_prob);
    if (flip_prob<((rand_gen()+1)/2.)) return;
    tensor_vmc_wf=tensor_vmc_wf_update;
    //cout << "Spin Flipped!" << endl;
}
//parallel update, using rand_gen_parallel
template <class TensorT, class RandGen>
void vmc_one_step_parallel(TensorT_VMC_WF<TensorT> &tensor_vmc_wf, RandGen &generator_parallel)
{
    std::vector<int> flip_inds(2);
    TensorT_VMC_WF<TensorT> tensor_vmc_wf_update=tensor_vmc_wf;
    //auto rand_gen_parallel=std::bind(distribution,generator_parallel);

    //get flipped spin until two spins antiparallel
    do
    {
        flip_inds[0]=floor((distribution(generator_parallel)+1)/2.*tensor_vmc_wf.n_sites());
        flip_inds[1]=floor((distribution(generator_parallel)+1)/2.*tensor_vmc_wf.n_sites());
        //Print(flip_inds[0]);
        //Print(flip_inds[1]);
    }
    while (tensor_vmc_wf.spin_config(flip_inds[0])==tensor_vmc_wf.spin_config(flip_inds[1]));
    tensor_vmc_wf_update.update_wf(flip_inds);

    double flip_prob=std::pow(std::abs(tensor_vmc_wf_update.wf_weight()/tensor_vmc_wf.wf_weight()),2.);

    //Print(tensor_vmc_wf_update.wf_weight());
    //Print(flip_prob);
    if (flip_prob<((distribution(generator_parallel)+1)/2.)) return;
    tensor_vmc_wf=tensor_vmc_wf_update;
}



//energy average over all bonds
template <class TensorT>
Complex vmc_Heisenberg_energy(const TensorT_VMC_WF<TensorT> &tensor_vmc_wf)
{
    Complex heisenberg_energy_total=0;
    for (int bondi=0; bondi<tensor_vmc_wf.n_bonds(); bondi++)
    {
        std::vector<int> site_inds=tensor_vmc_wf.lattice().bond_neighbour_sites(bondi);
        //two site with the same spin
        if (tensor_vmc_wf.spin_config(site_inds[0])==tensor_vmc_wf.spin_config(site_inds[1])) heisenberg_energy_total+=0.25;
        else
        {
            TensorT_VMC_WF<TensorT> tensor_vmc_wf_temp=tensor_vmc_wf;
            tensor_vmc_wf_temp.update_wf(site_inds);
            heisenberg_energy_total+=-0.25+0.5*tensor_vmc_wf_temp.wf_weight()/tensor_vmc_wf.wf_weight();
        }
    }
    //Print(heisenberg_energy_total/(tensor_vmc_wf.n_bonds()*1.));
    return heisenberg_energy_total/(tensor_vmc_wf.n_bonds()*1.);
}
template <class TensorT>
Complex vmc_SzSz_bonds_energy(const TensorT_VMC_WF<TensorT> &tensor_vmc_wf)
{
    Complex szsz_energy_total=0;
    for (int bondi=0; bondi<tensor_vmc_wf.n_bonds(); bondi++)
    {
        std::vector<int> site_inds=tensor_vmc_wf.lattice().bond_neighbour_sites(bondi);
        if (tensor_vmc_wf.spin_config(site_inds[0])==tensor_vmc_wf.spin_config(site_inds[1])) szsz_energy_total+=0.25;
        else szsz_energy_total-=0.25;
    }
    return szsz_energy_total/(tensor_vmc_wf.n_bonds()*1.);
}

#endif
