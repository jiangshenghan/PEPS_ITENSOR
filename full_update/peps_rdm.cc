
#include "peps_rdm.h"

template <class TensorT>
PEPSt_RDM<TensorT>::PEPSt_RDM(const std::string &name, const std::vector<int> &cutting_sites, const std::vector<int> &cutting_bonds, const std::vector<TensorT> &env_tensors, const PEPSt<TensorT> &peps, std::vector<int> env_contract_seq): 
    name_(name), 
    cutting_sites_(cutting_sites), 
    cutting_bonds_(cutting_bonds),
    env_tensors_(env_tensors),
    peps_(peps),
    env_contract_seq_(env_contract_seq)
{
    if (env_contract_seq_.empty())
        env_contract_seq_=std::vector<int>(env_tensors_.size(),0);

    obtain_RDM_and_wf_norm();
}
template
PEPSt_RDM<ITensor>::PEPSt_RDM(const std::string &name, const std::vector<int> &cutting_sites, const std::vector<int> &cutting_bonds, const std::vector<ITensor> &env_tensors, const PEPSt<ITensor> &peps, std::vector<int> env_contract_seq);
template
PEPSt_RDM<IQTensor>::PEPSt_RDM(const std::string &name, const std::vector<int> &cutting_sites, const std::vector<int> &cutting_bonds, const std::vector<IQTensor> &env_tensors, const PEPSt<IQTensor> &peps, std::vector<int> env_contract_seq);


template <class TensorT>
void PEPSt_RDM<TensorT>::obtain_RDM_and_wf_norm()
{
    //get reduced density matrix
    std::vector<TensorT> combined_tensors;
    combine_site_bond_tensors(cutting_site_tensors(),cutting_bond_tensors(),combined_tensors);
    RDM_=tensor_from_contract_patch(combined_tensors,false);
    //Print(RDM_);
    clean(RDM_);
    //Print(RDM_.norm());

    //get wf_norm_sq_
    auto RDM_trace=RDM_;
    for (int cuti=0; cuti<cutting_sites_.size(); cuti++)
    {
        const auto &phys_leg=peps_.phys_legs(cutting_sites_[cuti]);
        RDM_trace.trace(phys_leg,dag(prime(phys_leg)));
        //Print(RDM_trace.norm());
    }
    //Print(RDM_trace);
    wf_norm_sq_=RDM_trace.toComplex();
}
template
void PEPSt_RDM<ITensor>::obtain_RDM_and_wf_norm();
template
void PEPSt_RDM<IQTensor>::obtain_RDM_and_wf_norm();


template <class TensorT>
Complex PEPSt_RDM<TensorT>::expect_val_from_replaced_tensors(std::array<std::vector<TensorT>,2> replaced_site_tensors_ket_bra, std::array<std::vector<TensorT>,2> replaced_bond_tensors_ket_bra)
{
    //calculate time
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start=std::chrono::system_clock::now();

    //init bra tensors
    for (int sitei=0; sitei<replaced_site_tensors_ket_bra[0].size(); sitei++)
    {
        if (replaced_site_tensors_ket_bra[0][sitei].indices()[0].dir()==replaced_site_tensors_ket_bra[1][sitei].indices()[0].dir())
        {
            replaced_site_tensors_ket_bra[1][sitei].prime(Link).dag();
        }
    }

    for (int bondi=0; bondi<replaced_bond_tensors_ket_bra[0].size(); bondi++)
    {
        if (replaced_bond_tensors_ket_bra[0][bondi].indices()[0].dir()==replaced_bond_tensors_ket_bra[1][bondi].indices()[0].dir())
        {
            replaced_bond_tensors_ket_bra[1][bondi].prime(Link).dag();
        }
    }

    //get combined tensor, and do contraction
    std::array<std::vector<TensorT>,2> replaced_combined_tensors_ket_bra;
    combine_site_bond_tensors(replaced_site_tensors_ket_bra[0],replaced_bond_tensors_ket_bra[0],replaced_combined_tensors_ket_bra[0]);
    combine_site_bond_tensors(replaced_site_tensors_ket_bra[1],replaced_bond_tensors_ket_bra[1],replaced_combined_tensors_ket_bra[1]);
    TensorT scalar_tensor=tensor_from_contract_patch(replaced_combined_tensors_ket_bra);

    end=std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed=end-start;
    //cout << "time to get expect val: " << endl;
    //Print(elapsed.count());

    return scalar_tensor.toComplex();
}
template
Complex PEPSt_RDM<ITensor>::expect_val_from_replaced_tensors(std::array<std::vector<ITensor>,2> replaced_site_tensors_ket_bra, std::array<std::vector<ITensor>,2> replaced_bond_tensors_ket_bra);
template
Complex PEPSt_RDM<IQTensor>::expect_val_from_replaced_tensors(std::array<std::vector<IQTensor>,2> replaced_site_tensors_ket_bra, std::array<std::vector<IQTensor>,2> replaced_bond_tensors_ket_bra);


template <class TensorT>
TensorT PEPSt_RDM<TensorT>::tensor_from_contract_patch(const std::array<std::vector<TensorT>,2> &tensors_ket_bra)
{
    //TODO:check the modified contraction sequence
    TensorT patch_tensor;
    std::vector<bool> env_contracted(env_tensors_.size(),false);
    for (int i=0; i<tensors_ket_bra[0].size(); i++)
    {
        TensorT double_layer_tensor=tensors_ket_bra[0][i];
        //contract some env tensors with single layer tensors
        for (int envi=0; envi<env_tensors_.size(); envi++)
        {
            //Print(envi);
            //Print(env_contracted[envi]);
            //Print(env_contract_seq_[envi]);
            //Print(env_tensors_[envi].indices());
            //Print(double_layer_tensor.indices());

            const auto &env_tensor=env_tensors_[envi];
            if (env_contracted[envi]) continue;
            if (env_contract_seq_[envi]==1) continue;
            if (commonIndex(env_tensor,tensors_ket_bra[0][i])!=IndexT::Null())
            {
                //Print(envi);
                //Print(commonIndex(env_tensor,tensors_ket_bra[0][i]));
                double_layer_tensor*=env_tensor;
                env_contracted[envi]=true;
            }
        }

        double_layer_tensor*=tensors_ket_bra[1][i];
        //Print(double_layer_tensor.indices());

        //contract other env tensors with double layer tensors
        for (int envi=0; envi<env_tensors_.size(); envi++)
        {
            const auto &env_tensor=env_tensors_[envi];
            if (env_contracted[envi]) continue;
            if (env_contract_seq_[envi]==0) continue;
            if (commonIndex(env_tensor,tensors_ket_bra[0][i])!=IndexT::Null())
            {
                double_layer_tensor*=env_tensor;
                env_contracted[envi]=true;
            }
        }
        //Print(double_layer_tensor.indices());

        if (i==0)
        {
            patch_tensor=double_layer_tensor;
        }
        else
        {
            patch_tensor*=double_layer_tensor;
        }
    }
    //Print(patch_tensor);
    return patch_tensor;
}
template
ITensor PEPSt_RDM<ITensor>::tensor_from_contract_patch(const std::array<std::vector<ITensor>,2> &tensors_ket_bra);
template
IQTensor PEPSt_RDM<IQTensor>::tensor_from_contract_patch(const std::array<std::vector<IQTensor>,2> &tensors_ket_bra);


template <class TensorT>
void PEPSt_RDM<TensorT>::combine_site_bond_tensors(const std::vector<TensorT> &site_tensors, const std::vector<TensorT> &bond_tensors, std::vector<TensorT> &combined_tensors)
{
    combined_tensors.clear();
    //TODO:square geometry

    //kagome normal geometry
    if (peps_.name().find("kagome normal")!=std::string::npos)
    {
        if (name_.find("two sites shape")!=std::string::npos)
        {
            combined_tensors.push_back(site_tensors[0]*bond_tensors[0]*bond_tensors[1]);
            combined_tensors.push_back(site_tensors[1]*bond_tensors[2]);
        }
    }

    //kagome cirac geometry
    if (peps_.name().find("kagome cirac")!=std::string::npos)
    {
        if (name_.find("tree shape I")!=std::string::npos)
        {
            combined_tensors.push_back(site_tensors[0]*bond_tensors[0]);
            combined_tensors.push_back(site_tensors[1]);
            combined_tensors.push_back(site_tensors[2]);
        }
    }

    //PrintDat(combined_tensors);

}
template
void PEPSt_RDM<ITensor>::combine_site_bond_tensors(const std::vector<ITensor> &site_tensors, const std::vector<ITensor> &bond_tensors, std::vector<ITensor> &combined_tensors);
template
void PEPSt_RDM<IQTensor>::combine_site_bond_tensors(const std::vector<IQTensor> &site_tensors, const std::vector<IQTensor> &bond_tensors, std::vector<IQTensor> &combined_tensors);




double heisenberg_energy_from_RDM(PEPSt_RDM<IQTensor> peps_rdm)
{
    std::vector<IQIndex> phys_legs;
    for (int cuti=0; cuti<peps_rdm.cutting_sites().size(); cuti++)
    {
        phys_legs.push_back(peps_rdm.cutting_phys_legs(cuti));
    }

    Complex energy=0;

    if (peps_rdm.peps_name().find("square")!=std::string::npos || peps_rdm.peps_name().find("kagome normal")!=std::string::npos)
    {
        NN_Heisenberg_Hamiltonian heisenberg_gate({phys_legs[0],phys_legs[1]});
        energy=(peps_rdm.RDM()*(heisenberg_gate.site_tensors(0)*heisenberg_gate.bond_tensor()*heisenberg_gate.site_tensors(1))).toComplex();
    }

    if (peps_rdm.peps_name().find("kagome cirac")!=std::string::npos)
    {
        IQTPO heisenberg_gate=SpinSpin_kagome_cirac(phys_legs);
        energy=(peps_rdm.RDM()*heisenberg_gate.site_tensors(0)*heisenberg_gate.site_tensors(1)*heisenberg_gate.bond_tensors(0)*heisenberg_gate.site_tensors(2)).toComplex();
    }


    return (energy/peps_rdm.wf_norm_sq()).real();
}

