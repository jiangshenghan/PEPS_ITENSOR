
#include "tensor_vmc.h"

template<class TensorT>
TensorT_VMC_WF<TensorT>::TensorT_VMC_WF(const std::vector<int> init_spin_config, const PEPSt<TensorT> &peps, int maxm):
    spin_config_(init_spin_config),
    combined_tensors_(peps.combined_site_tensors()),
    tensor_rg_(peps.lattice())
{
    if (peps.d()!=2)
    {
        cout << "Have not implemented d!=2 case!" << endl;
        exit(1);
    }

    //init tensor_rg_
    std::vector<TensorT> input_tensors;
    for (int sitei=0; sitei<spin_config_.size(); sitei++) 
    {
        spin_prod_wf_.push_back(TensorT(dag(peps.phys_legs(sitei)(spin_config_[sitei]+1))));
        input_tensors.push_back(combined_tensors_[sitei]*spin_prod_wf_[sitei]);
    }
    tensor_rg_=TensorT_RG<TensorT>(peps.lattice(),input_tensors,maxm);
    Print(tensor_rg_.trg_result());
}
template
TensorT_VMC_WF<ITensor>::TensorT_VMC_WF(const std::vector<int> init_spin_config, const PEPSt<ITensor> &peps, int maxm);
template
TensorT_VMC_WF<IQTensor>::TensorT_VMC_WF(const std::vector<int> init_spin_config, const PEPSt<IQTensor> &peps, int maxm);


template <class TensorT>
void TensorT_VMC_WF<TensorT>::update_wf(std::vector<int> flip_inds)
{
    std::vector<TensorT> update_input_tensors;
    for (int ind: flip_inds)
    {
        spin_config_[ind]=1-spin_config_[ind];
        IndexT phys_indice=spin_prod_wf_[ind].indices()[0];
        spin_prod_wf_[ind]=TensorT(phys_indice(spin_config_[ind]+1));
        update_input_tensors.push_back(combined_tensors_[ind]*spin_prod_wf_[ind]);
    }
    //Print(flip_inds);
    tensor_rg_.update_trg_network(flip_inds,update_input_tensors);
}
template
void TensorT_VMC_WF<ITensor>::update_wf(std::vector<int> flip_inds);
template
void TensorT_VMC_WF<IQTensor>::update_wf(std::vector<int> flip_inds);


template <class TensorT>
void tensor_vmc(const PEPSt<TensorT> &peps, const Args &measure_args)
{
    int thermal_steps=measure_args.getInt("ThermalSteps",50),
        measure_steps=measure_args.getInt("MeasureSteps",1000),
        bin_no=measure_args.getInt("BinNo",10),
        bin_steps=measure_steps/bin_no,
        sweep_no=peps.n_sites_total();
    int maxm=measure_args.getInt("Maxm",100);
    std::string ope=measure_args.getString("Operator","SzSz"),
                init_cond=measure_args.getString("InitSpins","antiferro");

    Print(thermal_steps);
    Print(measure_steps);
    Print(bin_no);
    Print(maxm);
    Print(ope);

    //init tensor_vmc_wf
    //init spin configs
    std::vector<int> init_spin_config(peps.n_sites_total());
    if (init_cond.find("antiferro")!=std::string::npos)
    {
        for (int sitei=0; sitei<init_spin_config.size(); sitei++)
            init_spin_config[sitei]=sitei%2;
    }
    if (init_cond.find("random")!=std::string::npos)
    {
        std::vector<int> spin_no(2,0);
        for (int sitei=0; sitei<init_spin_config.size(); sitei++)
        {
            int spin=round((rand_gen()+1)/2.);
            init_spin_config[sitei]=spin;
            spin_no[spin]++;
            if (spin_no[0]==init_spin_config.size()/2)
            {
                for (int sitej=sitei+1; sitej<init_spin_config.size(); sitej++) init_spin_config[sitej]=1;
                break;
            }
            if (spin_no[1]==init_spin_config.size()/2)
            {
                for (int sitej=sitei+1; sitej<init_spin_config.size(); sitej++) init_spin_config[sitej]=0;
                break;
            }
        }
    }
    Print(init_spin_config);

    TensorT_VMC_WF<TensorT> tensor_vmc_wf(init_spin_config,peps,maxm);
    if (tensor_vmc_wf.is_zero()==true)
    {
        cout << "Not good init state!" << endl;
        exit(1);
    }

    std::vector<Complex> bins_energy(bin_no,0);
    for (int sweepi=-thermal_steps; sweepi<measure_steps; sweepi++)
    {
        Print(sweepi);
        //random flip spins
        for (int flipi=0; flipi<sweep_no; flipi++) 
        {
            //Print(sweepi);
            //Print(flipi);
            vmc_one_step<TensorT>(tensor_vmc_wf);
        }

        //sweepi<0 means thermalization process
        if (sweepi<0) continue;

        Complex energy_temp;
        if (ope.find("SzSz")!=std::string::npos) energy_temp=vmc_SzSz_bonds_energy(tensor_vmc_wf);
        if (ope.find("Heisenberg")!=std::string::npos) energy_temp=vmc_Heisenberg_energy(tensor_vmc_wf);
        bins_energy[sweepi%bin_no]+=energy_temp;

        Print(energy_temp);
    }
    
    Complex energy=0;
    double energysq=0;
    for (int bini=0; bini<bin_no; bini++)
    {
        bins_energy[bini]/=bin_steps*1.;
        energy+=bins_energy[bini];
        energysq+=std::pow(std::abs(bins_energy[bini]),2.);
    }
    energy/=bin_no*1.;
    energysq/=bin_no*1.;
    double energy_err=sqrt((energysq-std::pow(std::abs(energy),2.))/(bin_no-1));

    Print(energy);
    Print(energy_err);
}
template
void tensor_vmc(const PEPSt<ITensor> &peps, const Args &measure_args);
template
void tensor_vmc(const PEPSt<IQTensor> &peps, const Args &measure_args);


template <class TensorT>
void tensor_vmc_parallel(const PEPSt<TensorT> &peps, const Args &measure_args)
{
    int mpi_id, bin_no;
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_id);
    MPI_Comm_size(MPI_COMM_WORLD,&bin_no);
    
    int thermal_steps=measure_args.getInt("ThermalSteps",50),
        measure_steps=measure_args.getInt("MeasureSteps",1000),
        sweep_no=peps.n_sites_total();
    int maxm=measure_args.getInt("Maxm",100);
    std::string ope=measure_args.getString("Operator","SzSz"),
                init_cond=measure_args.getString("InitSpins","antiferro");

    if (mpi_id==0)
    {
        Print(thermal_steps);
        Print(measure_steps);
        Print(bin_no);
        Print(maxm);
        Print(ope);
    }

    //init random number generator
    static std::default_random_engine generator_parallel(std::time(0)+mpi_id*std::time(0));


    //init tensor_vmc_wf
    //init spin configs
    std::vector<int> init_spin_config(peps.n_sites_total());
    if (init_cond.find("antiferro")!=std::string::npos)
    {
        for (int sitei=0; sitei<init_spin_config.size(); sitei++)
            init_spin_config[sitei]=sitei%2;
    }
    if (init_cond.find("random")!=std::string::npos)
    {
        std::vector<int> spin_no(2,0);
        for (int sitei=0; sitei<init_spin_config.size(); sitei++)
        {
            int spin=round((distribution(generator_parallel)+1)/2.);
            init_spin_config[sitei]=spin;
            spin_no[spin]++;
            if (spin_no[0]==init_spin_config.size()/2)
            {
                for (int sitej=sitei+1; sitej<init_spin_config.size(); sitej++) init_spin_config[sitej]=1;
                break;
            }
            if (spin_no[1]==init_spin_config.size()/2)
            {
                for (int sitej=sitei+1; sitej<init_spin_config.size(); sitej++) init_spin_config[sitej]=0;
                break;
            }
        }
    }
    //Print(init_spin_config);

    TensorT_VMC_WF<TensorT> tensor_vmc_wf(init_spin_config,peps,maxm);
    if (tensor_vmc_wf.is_zero()==true)
    {
        cout << "Not good init state!" << endl;
        exit(1);
    }

    Complex bin_energy=0;
    for (int sweepi=-thermal_steps; sweepi<measure_steps; sweepi++)
    {
        if (mpi_id==0) Print(sweepi);
        //random flip spins
        for (int flipi=0; flipi<sweep_no; flipi++) 
        {
            vmc_one_step_parallel(tensor_vmc_wf,generator_parallel);
            if (mpi_id==0) 
            {
                Print(flipi);
                Print(tensor_vmc_wf.wf_weight());
            }
        }

        //sweepi<0 means thermalization process
        if (sweepi<0) continue;

        Complex energy_temp;
        if (ope.find("SzSz")!=std::string::npos) energy_temp=vmc_SzSz_bonds_energy(tensor_vmc_wf);
        if (ope.find("Heisenberg")!=std::string::npos) energy_temp=vmc_Heisenberg_energy(tensor_vmc_wf);
        bin_energy+=energy_temp;
        Print(energy_temp);
    }
    bin_energy/=measure_steps*1.;
    
    //collect energy data and analysis
    if (mpi_id!=0) MPI_Send(&bin_energy,1,MPI::DOUBLE_COMPLEX,0,1,MPI_COMM_WORLD);
    if (mpi_id==0)
    {
        Complex energy=bin_energy;
        double energysq=pow(std::abs(bin_energy),2.);

        MPI_Status mpi_stat;
        for (int bini=1; bini<bin_no; bini++)
        {
            Complex energy_recv;
            MPI_Recv(&energy_recv,1,MPI::DOUBLE_COMPLEX,bini,1,MPI_COMM_WORLD,&mpi_stat);
            energy+=energy_recv;
            energysq+=pow(std::abs(energy_recv),2.);
        }
        energy/=bin_no*1.;
        energysq/=bin_no*1.;
        double energy_err=sqrt((energysq-std::pow(std::abs(energy),2.))/(bin_no-1));

        Print(energy);
        Print(energy_err);
    }
}
template
void tensor_vmc_parallel(const PEPSt<ITensor> &peps, const Args &measure_args);
template
void tensor_vmc_parallel(const PEPSt<IQTensor> &peps, const Args &measure_args);




