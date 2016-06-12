
#include "spin_sym_decomp.h"
#include "kagome_rvb.h"
#include "tensor_rg.h"

int main()
{
    int Lx=4, Ly=4;

    Kagome_Normal_Lattice_Torus kagome_torus({Lx,Ly});
    IQPEPS_IndexSet_SpinHalf indexset(3,kagome_torus);

    std::vector<IQPEPS> kagome_srvbs(2,IQPEPS(kagome_torus,indexset));

    kagome_psg::mu_12=1;
    kagome_psg::mu_c6=1;
    kagome_normal_srvb_peps(kagome_srvbs[0]);

    kagome_psg::mu_12=1;
    kagome_psg::mu_c6=1;
    kagome_normal_srvb_peps(kagome_srvbs[1]);

    std::vector<std::vector<IQTensor>> kagome_srvbs_tensors;
    for (const auto &kagome_srvb: kagome_srvbs) kagome_srvbs_tensors.push_back(kagome_srvb.combined_site_tensors());
    for (auto &tensor: kagome_srvbs_tensors[1]) tensor.prime(Link).dag();

    std::vector<IQCombiner> inds_combiners;
    for (const auto &leg: indexset.virt_legs()) inds_combiners.push_back(IQCombiner(leg,dag(prime(leg))));
    std::vector<IQTensor> overlap_tensors;
    for (int sitei=0; sitei<kagome_torus.n_sites_total(); sitei++) overlap_tensors.push_back(kagome_srvbs_tensors[0][sitei]*kagome_srvbs_tensors[1][sitei]);
    for (int sitei=0; sitei<kagome_torus.n_sites_total(); sitei++)
    {
        for (int legi=0; legi<indexset.virt_legs().size(); legi++)
        {
            if (hasindex(overlap_tensors[sitei],indexset.virt_legs(legi))) 
            {
                if (hasindex(kagome_srvbs[0].site_tensors(sitei),indexset.virt_legs(legi))) overlap_tensors[sitei]=overlap_tensors[sitei]*inds_combiners[legi];
                else overlap_tensors[sitei]=overlap_tensors[sitei]*dag(inds_combiners[legi]);
            }
        }
        Print(sitei);
        Print(overlap_tensors[sitei].indices());
    }
    TensorT_RG<IQTensor> wf_overlap(kagome_torus,overlap_tensors,20);
    Print(wf_overlap.trg_result());
}

