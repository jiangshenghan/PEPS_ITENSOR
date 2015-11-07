
#ifndef _PEPS_H_
#define _PEPS_H_

#include "TPO.h"
#include "peps_indexset.h"

//
//class PEPSt
//the lower t is for template
//
//using PEPS for ITensor while IQPEPS for IQTensor
//
template <class TensorT>
class PEPSt
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
        //PEPSt() {}
        PEPSt(const Lattice_Base &lattice);
        PEPSt(const Lattice_Base &lattice, PEPSt_IndexSet_Base<IndexT> &index_set);

        //Initialize site tensors by tensor data in one unit cell.
        //Thus, the PEPS is translationally invariant.
        PEPSt(const Lattice_Base &lattice, PEPSt_IndexSet_Base<IndexT> &index_set, std::vector<TensorT> site_tensors_uc, std::vector<TensorT> bond_tensors_uc);

        //
        //Access Methods
        //
        int d() const
        {
            return d_;
        }

        int D() const
        {
            return D_;
        }

        const std::array<int,2> &n_uc() const
        {
            return lattice_.n_uc();
        }

        inline const int &n_sites_uc() const
        {
            return lattice_.n_sites_uc();
        }

        inline const int &n_bonds_uc() const
        {
            return lattice_.n_bonds_uc();
        }

        inline const int &n_bonds_to_one_site() const
        {
            return lattice_.n_bonds_to_one_site();
        }

        inline const int &n_sites_total() const
        {
            return lattice_.n_sites_total();
        }

        inline const int &n_bonds_total() const
        {
            return lattice_.n_bonds_total();
        }

        int n_boundary_legs() const
        {
            return lattice_.n_boundary_legs();
        }

        inline const Lattice_Base &lattice() const
        {
            return lattice_;
        }

        inline const IndexT &phys_legs(const int &site_i) const
        {
            return indexset_ptr_->phys_legs(site_i);
        }

        inline const IndexT &virt_legs(const int &leg_i) const
        {
            return indexset_ptr_->virt_legs_(leg_i);
        }

        inline const TensorT &site_tensors(int i) const
        {
            return site_tensors_[i];
        }

        TensorT &site_tensors(int i)
        {
            return site_tensors_[i];
        }

        inline const std::vector<TensorT> &site_tensors() const
        {
            return site_tensors_;
        }

        inline const TensorT &bond_tensors(int i) const
        {
            return bond_tensors_[i];
        }

        TensorT &bond_tensors(int i)
        {
            return bond_tensors_[i];
        }

        inline const std::vector<TensorT> &bond_tensors() const
        {
            return bond_tensors_;
        }

        const TensorT &boundary_tensors(int i) const
        {
            return boundary_tensors_[i];
        }

        TensorT &boundary_tensors(int i)
        {
            return boundary_tensors_[i];
        }

        const std::vector<TensorT> &boundary_tensors() const
        {
            return boundary_tensors_;
        }

        std::vector<TensorT> &boundary_tensors()
        {
            return boundary_tensors_;
        }

        const std::string &name() const
        {
            return name_;
        }

        //
        //Methods
        //
        //given tensors in one uc, generate all translation related tensors
        //this function requires order for indices of given tensors must match those in peps
        void generate_site_tensors(std::vector<TensorT> site_tensors_uc);
        void generate_bond_tensors(std::vector<TensorT> bond_tensors_uc, double mu_12=1);
        //void generate_boundary_tensors(TensorT single_boundary_tensor);

        //this function returns site tensors that absorb the neighbour bond tensors and boundary tensors
        std::vector<TensorT> combined_site_tensors() const;


        //
        //Constructor Helpers
        //

        //According to the itensor library, ITensor is constructed by IndexSet<Index>(indexset), while IQTensor is constructed by indexset directly
        void construct_tensor(TensorT &tensor, std::vector<IndexT> &indexset);
        //construct new tensors
        void new_site_tensors();
        void new_bond_tensors();
        void new_boundary_tensors();
        //using random initial site tensors
        //for IQTensor we need to make sure at least one block is created
        //void random_site_tensors();

        //read/write from/to file
        void read(std::istream &s);
        void write(std::ostream &s) const;

    private:
        //
        //Data Member
        //
        //d_ is physical index dimension while D_ is virtual leg dim.
        int d_, D_;

        const Lattice_Base &lattice_;

        std::shared_ptr<PEPSt_IndexSet_Base<IndexT>> indexset_ptr_;

        std::vector<TensorT> site_tensors_, bond_tensors_; 

        //using for boundary condition
        //TODO:replaced by vector of MPS?
        std::vector<TensorT> boundary_tensors_;


        //name stores information about lattice, spin symmetry, extra degenerate
        std::string name_;
};
using PEPS=PEPSt<ITensor>;
using IQPEPS=PEPSt<IQTensor>;


//void randomize_spin_sym_square_peps(IQPEPS &square_peps);

template <class TensorT>
inline Tnetwork_Storage<TensorT> peps_to_tnetwork_storage(const PEPSt<TensorT> &peps)
{
    Tnetwork_Storage<TensorT> tnetwork_storage;
    if (peps.name().find("square")!=std::string::npos) 
    {
        tnetwork_storage._tnetwork_type=1;
        tnetwork_storage._n_subl=1;
    }
    if (peps.name().find("torus")!=std::string::npos) tnetwork_storage._boundary_condition=1;
    int Lx=peps.n_uc()[0], Ly=peps.n_uc()[1];
    tnetwork_storage._Lx=Lx;
    tnetwork_storage._Ly=Ly;

    tnetwork_storage._tensor_list.set_size(peps.n_sites_total());
    auto tensor_list=peps.combined_site_tensors();
    for (int sitei=0; sitei<peps.n_sites_total(); sitei++) 
        tnetwork_storage._tensor_list(sitei)=tensor_list[sitei];

    tnetwork_storage._coor_to_siteind.set_size(Lx,Ly);
    for (int x=0; x<Lx; x++)
    {
        for (int y=0; y<Ly; y++)
        {
            tnetwork_storage._coor_to_siteind(x,y).set_size(tnetwork_storage._n_subl);
            for (int subli=0; subli<tnetwork_storage._n_subl; subli++)
                tnetwork_storage._coor_to_siteind(x,y)(subli)=peps.lattice().site_coord_to_list(x,y,subli);
        }
    }

    return tnetwork_storage;
}

#endif
