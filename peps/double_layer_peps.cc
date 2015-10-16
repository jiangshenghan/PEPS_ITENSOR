
#include "double_layer_peps.h"

template <class TensorT>
Double_Layer_PEPSt<TensorT>::Double_Layer_PEPSt(const Lattice_Base &lattice):
    single_layer_peps_(lattice),
    single_layer_tensors_(lattice.n_sites_total()),
    double_layer_tensors_(lattice.n_sites_total()),
    virt_leg_combiners_(lattice.n_sites_total())
{}
template
Double_Layer_PEPSt<ITensor>::Double_Layer_PEPSt(const Lattice_Base &lattice);
template
Double_Layer_PEPSt<IQTensor>::Double_Layer_PEPSt(const Lattice_Base &lattice);

template <class TensorT>
Double_Layer_PEPSt<TensorT>::Double_Layer_PEPSt(const PEPSt<TensorT> &peps):
    single_layer_peps_(peps),
    single_layer_tensors_(peps.n_sites_total()),
    double_layer_tensors_(peps.n_sites_total()),
    virt_leg_combiners_(peps.n_sites_total())
{
    obtain_single_layer_tensors();
    
    obtain_layered_tensors_with_combined_legs();
}
template
Double_Layer_PEPSt<ITensor>::Double_Layer_PEPSt(const PEPSt<ITensor> &peps);
template
Double_Layer_PEPSt<IQTensor>::Double_Layer_PEPSt(const PEPSt<IQTensor> &peps);


template <class TensorT>
void Double_Layer_PEPSt<TensorT>::obtain_tensor_sandwich_single_site_operators(std::vector<TensorT> single_site_operators, const std::vector<int> &acting_sites, std::vector<TensorT>& sandwiched_tensors)
{
    //auto sandwiched_tensors=double_layer_tensors_;

    for (auto site_no : acting_sites)
    {
        //get operator act on site_no
        IndexT phys_leg;
        for (const auto &leg : single_layer_tensors_[site_no].indices())
        {
            if (leg.type()==Site)
            {
                phys_leg=leg;
                break;
            }
        }
        for (int legi=0; legi<2; legi++)
        {
            auto old_leg=single_site_operators[site_no].indices()[legi];
            if (old_leg.primeLevel()==0) single_site_operators[site_no].replaceIndex(old_leg,dag(phys_leg));
            if (old_leg.primeLevel()==1) single_site_operators[site_no].replaceIndex(old_leg,prime(phys_leg));
        }

        //acting operator and replace the corresponding sandwiched_tensors
        TensorT lower_tensor=single_layer_tensors_[site_no],
                upper_tensor=prime(dag(lower_tensor));
        sandwiched_tensors[site_no]=lower_tensor*single_site_operators[site_no]*upper_tensor;
        for (const auto &combiner : virt_leg_combiners_[site_no])
        {
            sandwiched_tensors[site_no]=sandwiched_tensors[site_no]*combiner;
        }
    }
}
template
void Double_Layer_PEPSt<ITensor>::obtain_tensor_sandwich_single_site_operators(std::vector<ITensor> single_site_operators, const std::vector<int> &acting_sites, std::vector<ITensor>& sandwiched_tensors);
template
void Double_Layer_PEPSt<IQTensor>::obtain_tensor_sandwich_single_site_operators(std::vector<IQTensor> single_site_operators, const std::vector<int> &acting_sites, std::vector<IQTensor>& sandwiched_tensors);


template <class TensorT>
void Double_Layer_PEPSt<TensorT>::obtain_single_layer_tensors()
{
    for (int sitei=0; sitei<single_layer_peps_.n_sites_total(); sitei++)
    {
        auto site_tensor=single_layer_peps_.site_tensors(sitei);

        //multiply neighbouring bond tensors 
        for (int neighi=0; neighi<this->lattice().n_bonds_to_one_site(); neighi++)
        {
            int bond_no=this->lattice().site_neighbour_bonds(sitei,neighi);
            if (bond_no<0) continue;
            //For bulk bond, we only multiply those start from the site to avoid double counting. Namely bond_end_sites(bond_no,0)==sitei
            //For boundary bond, we will always absorb to the bond
            if (this->lattice().bond_end_sites(bond_no,0)==sitei ||
                this->lattice().bond_end_sites(bond_no,0)<0)
            {
                site_tensor*=single_layer_peps_.bond_tensors(bond_no);
            }
        }
        //multiply neighbouring boundary tensors
        for (const auto &boundary_no : this->lattice().site_neighbour_boundary(sitei))
        {
            site_tensor*=single_layer_peps_.boundary_tensors(boundary_no);
        }
        clean(site_tensor);

        single_layer_tensors_[sitei]=site_tensor;

        //cout << "Single layer tensor " << sitei << endl << site_tensor << endl;
    }
}
template
void Double_Layer_PEPSt<ITensor>::obtain_single_layer_tensors();
template
void Double_Layer_PEPSt<IQTensor>::obtain_single_layer_tensors();


template <class TensorT>
void Double_Layer_PEPSt<TensorT>::obtain_layered_tensors_with_combined_legs()
{
    //create combiners for virt_leg and prime(dag(virt_leg)), where virt_leg connects sitei and neighbour_sites[neighbour_i]
    //every combiner should appear twice, so we stores combined indiceswhich only appear once, and delete it if it already appear twice
    std::vector<IndexT> combined_indices;
    std::vector<CombinerT> indices_combiners;

    int combine_i=0;
    for (int sitei=0; sitei<this->lattice().n_sites_total(); sitei++)
    {
        for (const auto &virt_leg : single_layer_tensors_[sitei].indices())
        {
            //for physical leg, we do not combine
            if (virt_leg.type()==Site) continue;

            auto virt_leg_iter=std::find(combined_indices.begin(),combined_indices.end(),virt_leg);

            //if the leg has already been combined, add the combiner to virt_leg_combiners_[sitei], and delete the legs and combiners in combined_indices(_combiners) since each leg only appear twice
            if (virt_leg_iter!=combined_indices.end())
            {
                std::swap(*virt_leg_iter,*(combined_indices.end()-1));
                std::swap(indices_combiners[virt_leg_iter-combined_indices.begin()],*(indices_combiners.end()-1));

                virt_leg_combiners_[sitei].push_back(dag(*(indices_combiners.end()-1)));

                combined_indices.pop_back();
                indices_combiners.pop_back();

                continue;
            }

            //for the case where the virtual leg appears the first time
            auto leg_combiner=CombinerT(virt_leg,dag(virt_leg).prime());
            leg_combiner.init(nameint("leg_combiner ",combine_i),Link,virt_leg.dir());
            virt_leg_combiners_[sitei].push_back(leg_combiner);
            combine_i++;

            combined_indices.push_back(virt_leg);
            indices_combiners.push_back(leg_combiner);
        }
    }

    //for (int sitei=0; sitei<this->lattice().n_sites_total(); sitei++)
    //{
    //    cout << "Combiners for site " << sitei << endl;
    //    for (const auto &leg_combiner : virt_leg_combiners_[sitei])
    //    {
    //        cout << leg_combiner;
    //    }
    //    cout << endl;
    //}


    //create layered tensors with virtual legs combined
    std::vector<TensorT> upper_tensors(single_layer_tensors_);
    for (auto &tensor : upper_tensors)
    {
        tensor.dag();
        tensor.prime(Link);
    }

    for (int sitei=0; sitei<this->lattice().n_sites_total(); sitei++)
    {
        double_layer_tensors_[sitei]=single_layer_tensors_[sitei]*upper_tensors[sitei];
        for (const auto &combiner : virt_leg_combiners_[sitei])
        {
            double_layer_tensors_[sitei]=double_layer_tensors_[sitei]*combiner;
        }
        clean(double_layer_tensors_[sitei]);
    }

}
template
void Double_Layer_PEPSt<ITensor>::obtain_layered_tensors_with_combined_legs();
template
void Double_Layer_PEPSt<IQTensor>::obtain_layered_tensors_with_combined_legs();


template <class TensorT>
void Double_Layer_PEPSt<TensorT>::read(std::istream &s)
{
    single_layer_peps_.read(s);
    
    for (auto &tensor : single_layer_tensors_) tensor.read(s);

    //reconstruct double_layer_tensors_ and virt_leg_combiners_
    obtain_layered_tensors_with_combined_legs();
}
template
void Double_Layer_PEPSt<ITensor>::read(std::istream &s);
template
void Double_Layer_PEPSt<IQTensor>::read(std::istream &s);

template <class TensorT>
void Double_Layer_PEPSt<TensorT>::write(std::ostream &s) const
{
    single_layer_peps_.write(s);

    for (const auto &tensor : single_layer_tensors_) tensor.write(s);
}
template
void Double_Layer_PEPSt<ITensor>::write(std::ostream &s) const;
template
void Double_Layer_PEPSt<IQTensor>::write(std::ostream &s) const;



//
//class Cylinder_Square_Double_Layer_PEPSt
//
template <class TensorT>
Cylinder_Square_Double_Layer_PEPSt<TensorT>::Cylinder_Square_Double_Layer_PEPSt(const Lattice_Base &square_cylinder):
    Double_Layer_PEPSt<TensorT>(square_cylinder),
    col_lr_{{0,square_cylinder.n_uc()[0]-1}},
    iterative_combiners_{{std::vector<CombinerT>(square_cylinder.n_uc()[1]),std::vector<CombinerT>(square_cylinder.n_uc()[1])}}
{
    //assert(square_cylinder.name().find("cylinder")!=std::string::npos);

    cutting_col_=this->lattice().n_uc()[0]/2;
}
template
Cylinder_Square_Double_Layer_PEPSt<ITensor>::Cylinder_Square_Double_Layer_PEPSt(const Lattice_Base &square_cylinder);
template
Cylinder_Square_Double_Layer_PEPSt<IQTensor>::Cylinder_Square_Double_Layer_PEPSt(const Lattice_Base &square_cylinder);

template <class TensorT>
Cylinder_Square_Double_Layer_PEPSt<TensorT>::Cylinder_Square_Double_Layer_PEPSt(const PEPSt<TensorT> &square_peps, int cutting_col):
    Double_Layer_PEPSt<TensorT>(square_peps),
    cutting_col_(cutting_col),
    col_lr_{{0,square_peps.lattice().n_uc()[0]-1}},
    iterative_combiners_{{std::vector<CombinerT>(square_peps.lattice().n_uc()[1]),std::vector<CombinerT>(square_peps.lattice().n_uc()[1])}}
{
    if (cutting_col_==-1) 
        cutting_col_=this->lattice().n_uc()[0]/2;
}
template
Cylinder_Square_Double_Layer_PEPSt<ITensor>::Cylinder_Square_Double_Layer_PEPSt(const PEPSt<ITensor> &square_peps, int cutting_col);
template
Cylinder_Square_Double_Layer_PEPSt<IQTensor>::Cylinder_Square_Double_Layer_PEPSt(const PEPSt<IQTensor> &square_peps, int cutting_col);

template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::obtain_boundary_theory_iterative()
{
    obtain_sigma_lr_iterative(cutting_col_-1,cutting_col_);

    from_sigma_vec_to_mat();
    from_sigma_lr_to_sigma_b();
    obtain_density_matrix_spectrum();
}
template
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::obtain_boundary_theory_iterative();
template 
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::obtain_boundary_theory_iterative();

template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::obtain_sigma_lr_iterative(int left_end_col, int right_end_col)
{
    snake_walking_boundary_col();

    for (col_lr_[0]=1; col_lr_[0]<=left_end_col; )
    {
        snake_walking_bulk_col(col_lr_[0],1,-1);
        col_lr_[0]++;

        if (col_lr_[0]>left_end_col) break;

        snake_walking_bulk_col(col_lr_[0],1,1);
        col_lr_[0]++;
    }
    col_lr_[0]=left_end_col;

    for (col_lr_[1]=this->lattice().n_uc()[0]-2; col_lr_[1]>=right_end_col; )
    {
        snake_walking_bulk_col(col_lr_[1],-1,1);
        col_lr_[1]--;

        if (col_lr_[1]<right_end_col) break;

        snake_walking_bulk_col(col_lr_[1],-1,-1);
        col_lr_[1]--;
    }
    col_lr_[1]=right_end_col;

    //cout << "\n========================================\n" << endl;
    //cout << "sigma_left for contraction to " << col_lr_[0] << " cols:" << endl;
    //PrintDat(sigma_lr_[0]);
    //cout << "sigma_right for contraction to " << col_lr_[1] << " cols:" << endl;
    //PrintDat(sigma_lr_[1]);
    //cout << "\n========================================\n" << endl;

}
template
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::obtain_sigma_lr_iterative(int left_end_col, int right_end_col);
template
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::obtain_sigma_lr_iterative(int left_end_col, int right_end_col);


template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::snake_walking_boundary_col()
{
    //Initialize the first col (left boundary) from down to up
    int sitei=0;
    iterative_combiners_[0][0]=CombinerT(commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei+1])));
    sigma_lr_[0]=this->double_layer_tensors_[sitei]*iterative_combiners_[0][0];
    clean(sigma_lr_[0]);

    for (int rowi=1; rowi<this->lattice().n_uc()[1]; rowi++)
    {
        sitei+=this->lattice().n_uc()[0];
        auto combining_indice=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei+1]));
        iterative_combiners_[0][rowi]=CombinerT(iterative_combiners_[0][rowi-1].right(),combining_indice);
        sigma_lr_[0]*=this->double_layer_tensors_[sitei];
        sigma_lr_[0]=sigma_lr_[0]*iterative_combiners_[0][rowi];
    }


    //Initialize the last col (right boundary) from up to down
    sitei=this->lattice().n_sites_total()-1;
    iterative_combiners_[1][this->lattice().n_uc()[1]-1]=CombinerT(commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei-1])));
    sigma_lr_[1]=this->double_layer_tensors_[sitei]*iterative_combiners_[1][this->lattice().n_uc()[1]-1];
    clean(sigma_lr_[1]);

    for (int rowi=this->lattice().n_uc()[1]-2; rowi>=0; rowi--)
    {
        sitei-=this->lattice().n_uc()[0];
        auto combining_indice=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei-1]));
        iterative_combiners_[1][rowi]=CombinerT(iterative_combiners_[1][rowi+1].right(),combining_indice);
        sigma_lr_[1]*=this->double_layer_tensors_[sitei];
        sigma_lr_[1]=sigma_lr_[1]*iterative_combiners_[1][rowi];
    }

    for (auto &sigma : sigma_lr_)
    {
        sigma/=sigma.norm();
        clean(sigma);
    }

    cout << "\n========================================\n" << endl;
    cout << "Contraction for boundary cols:" << endl;
    //cout << "iterative combiners for left part:" << endl;
    //for (const auto &combiner : iterative_combiners_[0]) cout << combiner;
    //cout << "iterative combiners for right part:" << endl;
    //for (const auto &combiner : iterative_combiners_[1]) cout << combiner;
    //cout << "sigma left:" << endl;
    //PrintDat(sigma_lr_[0]);
    //cout << "sigma right:" << endl;
    //PrintDat(sigma_lr_[1]);
    //cout << "first ten elems of sigma left:" << endl;
    //print_vector_nonzero_elem(sigma_lr_[0],10);
    //cout << "first ten elems of sigma right:" << endl;
    //print_vector_nonzero_elem(sigma_lr_[1],10);
    cout << "Finishing contracting for boundary_cols!" << endl;
    cout << "\n========================================\n" << endl;

}
template
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::snake_walking_boundary_col();
template
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::snake_walking_boundary_col();


template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::snake_walking_bulk_col(int coli, int horizontal_dir, int vertical_dir)
{
    int lr_no=(1-horizontal_dir)/2;
    int start_row=(1-vertical_dir)/2*(this->lattice().n_uc()[1]-1);
    int rowi=start_row;
    while (rowi>=0 && rowi<this->lattice().n_uc()[1])
    {
        int sitei=this->lattice().site_coord_to_list(coli,rowi,0);
        //decombine indice to be multiplied, and then mutiply a new tensor
        sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(iterative_combiners_[lr_no][rowi]);
        sigma_lr_[lr_no]*=this->double_layer_tensors_[sitei];

        //combine the indice of new tensor which will be multiplied next time
        auto combining_indice=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei+horizontal_dir]));
        if (rowi==start_row)
        {
            iterative_combiners_[lr_no][rowi]=CombinerT(combining_indice);
        }
        else
        {
            iterative_combiners_[lr_no][rowi]=CombinerT(iterative_combiners_[lr_no][rowi-vertical_dir].right(),combining_indice);
        }
        sigma_lr_[lr_no]=sigma_lr_[lr_no]*iterative_combiners_[lr_no][rowi];

        rowi+=vertical_dir;
    }

    sigma_lr_[lr_no]/=sigma_lr_[lr_no].norm();
    clean(sigma_lr_[lr_no]);

    cout << "\n========================================\n" << endl;
    cout << "Iterative contraction for bulk col:" << endl;
    cout << "Horizontal Direction: " << horizontal_dir << endl
         << "Vertical Direction: " << vertical_dir << endl
         << "Col: " << coli << endl << endl;
    //for (const auto &combiner : iterative_combiners_[lr_no]) cout << combiner;
    //PrintDat(sigma_lr_[lr_no]);
    //print the first ten nonzero elems
    //cout << "First ten elems of sigma_lr[" << lr_no << "] are" << endl;
    //print_vector_nonzero_elem(sigma_lr_[lr_no],10);
    cout << "\n========================================\n" << endl;

}
template
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::snake_walking_bulk_col(int coli, int horizontal_dir, int vertical_dir);
template
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::snake_walking_bulk_col(int coli, int horizontal_dir, int vertical_dir);

//TODO: make index combiner match to each other such that the multiplication is always in right order
template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::from_sigma_vec_to_mat()
{
    std::array<int,2> vertical_dirs{-2*(cutting_col_%2)+1,2*((this->lattice().n_uc()[0]-cutting_col_)%2)-1};

    std::array<std::vector<CombinerT>,2> upper_combiners{std::vector<CombinerT>(this->lattice().n_uc()[1]),std::vector<CombinerT>(this->lattice().n_uc()[1])}, 
                                         lower_combiners{std::vector<CombinerT>(this->lattice().n_uc()[1]),std::vector<CombinerT>(this->lattice().n_uc()[1])};
    for (int lr_no=0; lr_no<2; lr_no++)
    {
        int start_row=(1-vertical_dirs[lr_no])/2*(this->lattice().n_uc()[1]-1);
        int rowi=start_row;

        while (rowi>=0 && rowi<this->lattice().n_uc()[1])
        {
            sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(iterative_combiners_[lr_no][rowi]);

            //decombine upper leg and lower leg
            int sitei=this->lattice().site_coord_to_list(cutting_col_-1+lr_no,rowi,0);
            for (const auto &leg_combiner : this->virt_leg_combiners_[sitei])
            {
                if (hasindex(sigma_lr_[lr_no],leg_combiner.right()))
                {
                    sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(leg_combiner);

                    if (rowi==start_row)
                    {
                        //lower leg has no prime while upper leg has prime
                        if (left_leg_of_combiners(leg_combiner,0).primeLevel()==0)
                        {
                            lower_combiners[lr_no][rowi]=CombinerT(dag(left_leg_of_combiners(leg_combiner,0)));
                            //upper_combiners[lr_no][rowi]=CombinerT(dag(left_leg_of_combiners(leg_combiner,1)));
                        }
                        else
                        {
                            lower_combiners[lr_no][rowi]=CombinerT(dag(left_leg_of_combiners(leg_combiner,1)));
                            //upper_combiners[lr_no][rowi]=CombinerT(dag(left_leg_of_combiners(leg_combiner,0)));
                        }
                    }

                    if (rowi!=start_row)
                    {
                        if (left_leg_of_combiners(leg_combiner,0).primeLevel()==0)
                        {
                            lower_combiners[lr_no][rowi]=CombinerT(dag(left_leg_of_combiners(leg_combiner,0)),lower_combiners[lr_no][rowi-vertical_dirs[lr_no]].right());
                            //upper_combiners[lr_no][rowi]=CombinerT(dag(left_leg_of_combiners(leg_combiner,1)),upper_combiners[lr_no][rowi-vertical_dirs[lr_no]].right());
                        }
                        else
                        {
                            //upper_combiners[lr_no][rowi]=CombinerT(dag(left_leg_of_combiners(leg_combiner,0)),upper_combiners[lr_no][rowi-vertical_dirs[lr_no]].right());
                            lower_combiners[lr_no][rowi]=CombinerT(dag(left_leg_of_combiners(leg_combiner,1)),lower_combiners[lr_no][rowi-vertical_dirs[lr_no]].right());
                        }
                    }

                    upper_combiners[lr_no][rowi]=lower_combiners[lr_no][rowi];
                    upper_combiners[lr_no][rowi].prime();
                    upper_combiners[lr_no][rowi].dag();

                    sigma_lr_[lr_no]=sigma_lr_[lr_no]*lower_combiners[lr_no][rowi]*upper_combiners[lr_no][rowi];

                    break;
                }
            }

            rowi+=vertical_dirs[lr_no];
        }
    }

    for (auto &sigma : sigma_lr_)
    {
        sigma/=sigma.norm();
        clean(sigma);
    }

    //cout << "\n========================================\n" << endl;
    //cout << "Output for sigma_l and sigma_r:" << endl;
    //cout << "left part:" << endl;
    //for (int rowi=0; rowi<this->lattice().n_uc()[1]; rowi++)
    //{
    //    cout << lower_combiners[0][rowi];
    //    cout << upper_combiners[0][rowi];
    //}
    //PrintDat(sigma_lr_[0]);
    //cout << "right part:" << endl;
    //for (int rowi=0; rowi<this->lattice().n_uc()[1]; rowi++)
    //{
    //    cout << lower_combiners[1][rowi];
    //    cout << upper_combiners[1][rowi];
    //}
    //PrintDat(sigma_lr_[1]);
    //cout << "\n========================================\n" << endl;
}
template
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::from_sigma_vec_to_mat();
template
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::from_sigma_vec_to_mat();


//TODO: sigma_lr_'s may not have a good quantum number due to boundary
//condition? In this case diagHermitian does not work for IQTensor?
template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::from_sigma_lr_to_sigma_b()
{
    //PrintDat(sigma_lr_[0]);
    //get eigenvalues of sigma_l
    TensorT U_l, D_l;
    diagHermitian(sigma_lr_[0],U_l,D_l);

    //cout << "sigma_l singular value:" << endl;
    //PrintDat(D_l);

    auto Dl_legs=D_l.indices();
    
    for (int val=1; val<=Dl_legs[0].m(); val++)
    {
        double elem=D_l(Dl_legs[0](val),Dl_legs[1](val));
        if (elem>EPSILON)
        {
            D_l(Dl_legs[0](val),Dl_legs[1](val))=std::sqrt(elem);
        }
    }
    clean(D_l);

    auto sqrt_sigma_l=dag(U_l)*D_l*prime(U_l);

    //cout << "sqrt of sigma_l singular value:" << endl;
    //PrintDat(D_l);
    //cout << "sqrt of sigma_l:" << endl;
    //PrintDat(sqrt_sigma_l);

    //obtain sigma_b_, sigma_b_ shares the same index as sigma_l
    sigma_b_=sigma_lr_[1];
    for (const auto &right_leg : sigma_lr_[1].indices())
    {
        for (const auto &left_leg : sigma_lr_[0].indices())
        {
            if (right_leg.primeLevel()==left_leg.primeLevel())
            {
                auto temp_tensor=sqrt_sigma_l;
                temp_tensor.replaceIndex(left_leg,dag(right_leg));
                sigma_b_*=temp_tensor;
                break;
            }
        }
    }

    //cout << "\n==============================\n" << endl;
    //PrintDat(sigma_b_);
    //cout << "\n==============================\n" << endl;
}
template
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::from_sigma_lr_to_sigma_b();
template
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::from_sigma_lr_to_sigma_b();


template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::obtain_density_matrix_spectrum()
{
    TensorT U_b, D_b;

    diagHermitian(sigma_b_,U_b,D_b);

    auto Db_legs=D_b.indices();
    for (int val=1; val<=Db_legs[0].m(); val++)
    {
        density_mat_spectrum_.push_back(D_b(Db_legs[0](val),Db_legs[1](val)));
    }

    //normalize and sort density_mat_spectrum
    double spectrum_sum=0;
    for (const auto &eigval : density_mat_spectrum_) spectrum_sum+=eigval;
    for (auto &eigval : density_mat_spectrum_)  eigval/=spectrum_sum; 
    std::sort(density_mat_spectrum_.begin(),density_mat_spectrum_.end());
    auto max_eigval=*(density_mat_spectrum_.end()-1);
    for (auto &eigval : density_mat_spectrum_)
    {
        if (eigval/max_eigval<EPSILON)
            eigval=0;
    }

    //cout << "Spectrum of density matrix: " << endl;
    //for (double eigval : density_mat_spectrum_) cout << eigval << " ";
    //cout << endl;
}
template
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::obtain_density_matrix_spectrum();
template
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::obtain_density_matrix_spectrum();


template <class TensorT>
std::vector<double> Cylinder_Square_Double_Layer_PEPSt<TensorT>::entanglement_spectrum()
{
    std::vector<double> entanglement_spectrum;
    std::transform(density_mat_spectrum_.begin(),density_mat_spectrum_.end(),entanglement_spectrum.begin(),[](double eigval){ return std::exp(-eigval); });
    return entanglement_spectrum;
}
template
std::vector<double> Cylinder_Square_Double_Layer_PEPSt<ITensor>::entanglement_spectrum();
template
std::vector<double> Cylinder_Square_Double_Layer_PEPSt<IQTensor>::entanglement_spectrum();


template <class TensorT>
double Cylinder_Square_Double_Layer_PEPSt<TensorT>::entanglement_entropy_vN()
{
    double entanglement_entropy_vN=0;
    for (const auto &eigval : density_mat_spectrum_)
    {
        if (eigval<EPSILON) continue;
        entanglement_entropy_vN+=-eigval*std::log(eigval);
    }
    return entanglement_entropy_vN;
}
template
double Cylinder_Square_Double_Layer_PEPSt<ITensor>::entanglement_entropy_vN();
template
double Cylinder_Square_Double_Layer_PEPSt<IQTensor>::entanglement_entropy_vN();


template <class TensorT>
double Cylinder_Square_Double_Layer_PEPSt<TensorT>::entanglement_entropy_Renyi(double renyi_n)
{
    double entanglement_entropy_Renyi=0;
    for (const auto &eigval : density_mat_spectrum_)
    {
        if (eigval<EPSILON) continue;
        entanglement_entropy_Renyi+=renyi_n/(1-renyi_n)*std::log(std::pow(eigval,renyi_n));
    }
    return entanglement_entropy_Renyi;
}
template
double Cylinder_Square_Double_Layer_PEPSt<ITensor>::entanglement_entropy_Renyi(double renyi_n);
template
double Cylinder_Square_Double_Layer_PEPSt<IQTensor>::entanglement_entropy_Renyi(double renyi_n);



//template <class TensorT>
//void Cylinder_Square_Double_Layer_PEPSt<TensorT>::obtain_transfer_matrix(int coli)
//{
//    int sitei=this->lattice().site_coord_to_list(coli,0,0);
//    //leg_combiners[0/1] is for left/right legs of transfer matrix
//    std::array<std::vector<CombinerT>,2> leg_combiners;
//    std::array<IndexT,2> lr_legs;
//
//    lr_legs[0]=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei-1]));
//    lr_legs[1]=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei+1]));
//
//    leg_combiners[0].push_back(CombinerT(lr_legs[0]));
//    leg_combiners[1].push_back(CombinerT(lr_legs[1]));
//    
//    transfer_mat_=this->double_layer_tensors_[sitei]*leg_combiners[0][0]*leg_combiners[1][0];
//
//    for (int rowi=1; rowi<this->lattice().n_uc()[1]; rowi++)
//    {
//        sitei+=this->lattice().n_uc()[0];
//
//        lr_legs[0]=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei-1]));
//        lr_legs[1]=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei+1]));
//
//        leg_combiners[0].push_back(CombinerT(lr_legs[0],leg_combiners[0][rowi-1].right()));
//        leg_combiners[1].push_back(CombinerT(lr_legs[1],leg_combiners[1][rowi-1].right()));
//
//        transfer_mat_*=this->double_layer_tensors_[sitei];
//        transfer_mat_=transfer_mat_*leg_combiners[0][rowi]*leg_combiners[1][rowi];
//    }
//
//    //We make transfer_mat_ as matrix --I--mat--I'--
//    auto oind=leg_combiners[1][this->lattice().n_uc()[1]-1].right(),
//         nind=prime(leg_combiners[0][this->lattice().n_uc()[1]-1].right());
//    if ((oind.dir()-nind.dir())!=0) nind.dag();
//    transfer_mat_.replaceIndex(oind,nind);
//
//    transfer_mat_/=transfer_mat_.norm();
//    clean(transfer_mat_);
//    
//
//    //cout << "\n----------------------------\n" << endl;
//    //cout << "transfer matrix for col " << coli << ":" << endl;
//    //PrintDat(transfer_mat_);
//    //cout << "\n----------------------------\n" << endl;
//}
//template
//void Cylinder_Square_Double_Layer_PEPSt<ITensor>::obtain_transfer_matrix(int coli);
//template
//void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::obtain_transfer_matrix(int coli);


template <class TensorT>
double Cylinder_Square_Double_Layer_PEPSt<TensorT>::sandwiched_peps_norm(const std::vector<TensorT> &sandwiched_tensors, int horizontal_dir)
{
    //store the original data
    auto origin_col_lr=col_lr_;
    auto origin_sigma_lr=sigma_lr_;
    auto origin_iterative_combiners=iterative_combiners_;
    auto origin_double_layer_tensors=this->double_layer_tensors_;
    this->double_layer_tensors_=sandwiched_tensors;

    //do contraction along horizontal_dir until sigma_l and sigma_r meet
    int n_cols=this->lattice().n_uc()[0];
    int lr_no=(1-horizontal_dir)/2;
    int vertical_dir=1-2*std::abs((lr_no*n_cols-col_lr_[lr_no])%2);

    col_lr_[lr_no]+=horizontal_dir;
    vertical_dir*=-1;
    while (col_lr_[0]<col_lr_[1])
    {
        snake_walking_bulk_col(col_lr_[lr_no],horizontal_dir,vertical_dir);

        col_lr_[lr_no]+=horizontal_dir;
        vertical_dir*=-1;
        if (col_lr_[0]==col_lr_[1]) break;

        snake_walking_bulk_col(col_lr_[lr_no],horizontal_dir,vertical_dir);
        col_lr_[lr_no]+=horizontal_dir;
        vertical_dir*=-1;
    }

    //decombine sigma_lr and then multiply them
    col_lr_[lr_no]-=horizontal_dir;
    decombine_sigma_lr();

    double expect_value=(sigma_lr_[0]*sigma_lr_[1]).toReal();
    

    //get the original data
    col_lr_=origin_col_lr;
    sigma_lr_=origin_sigma_lr;
    iterative_combiners_=origin_iterative_combiners;
    this->double_layer_tensors_=origin_double_layer_tensors;

    return expect_value;
}
template
double Cylinder_Square_Double_Layer_PEPSt<ITensor>::sandwiched_peps_norm(const std::vector<ITensor> &sandwiched_tensors, int horizontal_dir);
template
double Cylinder_Square_Double_Layer_PEPSt<IQTensor>::sandwiched_peps_norm(const std::vector<IQTensor> &sandwiched_tensors, int horizontal_dir);


template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::decombine_sigma_lr()
{
    //only do decombining when sigma_lr_ are vectors
    if (sigma_lr_[0].r()!=1 || sigma_lr_[1].r()!=1)
    {
        cout << "Decombine failed: sigma_lr are not vectors." << endl;
        return;
    }

    int n_rows=this->lattice().n_uc()[1];
    for (int lr_no=0; lr_no<2; lr_no++)
    {
        if (iterative_combiners_[lr_no][0].numLeft()==1) //decombine from top to down
        {
            for (int rowi=n_rows-1; rowi>=0; rowi--)
            {
                sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(iterative_combiners_[lr_no][rowi]);
            }
            continue;
        }

        if (iterative_combiners_[lr_no][n_rows-1].numLeft()==1)//decombine from down to top
        {
            for (int rowi=0; rowi<n_rows; rowi++)
            {
                sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(iterative_combiners_[lr_no][rowi]);
            }
            continue;
        }
        assert(cerr << "Invalid iterative combiners!" << endl);
    }
}
template
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::decombine_sigma_lr();
template
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::decombine_sigma_lr();

template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::recombine_sigma_lr()
{
    int n_cols=this->lattice().n_uc()[0];
    int n_rows=this->lattice().n_uc()[1];
    for (int lr_no=0; lr_no<2; lr_no++)
    {
        int horizontal_dir=1-2*lr_no,
            vertical_dir=1-2*std::abs((lr_no*n_cols-col_lr_[lr_no])%2);
        int start_row=(1-vertical_dir)/2*(n_rows-1);
        int rowi=start_row;
        while (rowi>=0 && rowi<n_rows)
        {
            int sitei=this->lattice().site_coord_to_list(col_lr_[lr_no],rowi,0);

            auto combining_indice=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei+horizontal_dir]));
            if (rowi==start_row)
            {
                iterative_combiners_[lr_no][rowi]=CombinerT(combining_indice);
            }
            else
            {
                iterative_combiners_[lr_no][rowi]=CombinerT(iterative_combiners_[lr_no][rowi-vertical_dir].right(),combining_indice);
            }
            sigma_lr_[lr_no]=sigma_lr_[lr_no]*iterative_combiners_[lr_no][rowi];

            rowi+=vertical_dir;
        }

        clean(sigma_lr_[lr_no]);
    }

}
template
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::recombine_sigma_lr();
template
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::recombine_sigma_lr();

template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::match_indices_sigma_lr()
{
    for (int lr_no=0; lr_no<2; lr_no++)
    {
        int horizontal_dir=1-2*lr_no;
        IndexSet<IndexT> new_indices;
        for (int rowi=0; rowi<this->lattice().n_uc()[1]; rowi++)
        {
            int sitei=this->lattice().site_coord_to_list(col_lr_[lr_no],rowi,0);
            auto leg=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei+horizontal_dir]));
            new_indices.addindex(leg);
            //Print(leg);
        }

        for (auto oind : sigma_lr_[lr_no].indices())
        {
            //Print(oind);
            for (const auto &nind : new_indices)
            {
                if (oind.name()==nind.name())
                {
                    //Print(nind);
                    sigma_lr_[lr_no].replaceIndex(oind,nind);
                    break;
                }
            }
        }
    }
}
template
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::match_indices_sigma_lr();
template
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::match_indices_sigma_lr();

template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::move_sigma_lr(const std::array<int,2> &new_col_lr)
{
    int n_rows=this->lattice().n_uc()[1];
        

    //decombine sigma_lr_ to tensor, then replace the indices by indices in new cols, and finally, recombine it
    decombine_sigma_lr();
    for (int lr_no=0; lr_no<2; lr_no++)
    {
        int horizontal_dir=1-2*lr_no;

        //get old_indices and new_indices in the order of rows, so that they have one to one correspondance
        for (int rowi=0; rowi<n_rows; rowi++)
        {
            int old_site=this->lattice().site_coord_to_list(col_lr_[lr_no],rowi,0),
                new_site=this->lattice().site_coord_to_list(new_col_lr[lr_no],rowi,0);
            auto old_leg=commonIndex(this->double_layer_tensors_[old_site],dag(this->double_layer_tensors_[old_site+horizontal_dir]));
            auto new_leg=commonIndex(this->double_layer_tensors_[new_site],dag(this->double_layer_tensors_[new_site+horizontal_dir]));

            sigma_lr_[lr_no].replaceIndex(old_leg,new_leg);
        }

    }

    col_lr_=new_col_lr;
}
template 
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::move_sigma_lr(const std::array<int,2> &new_col_lr);
template
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::move_sigma_lr(const std::array<int,2> &new_col_lr);


template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::read(std::istream &s)
{
    Double_Layer_PEPSt<TensorT>::read(s);
    //for (const auto &tensor : this->double_layer_tensors()) Print(tensor);

    s.read((char*)&cutting_col_,sizeof(cutting_col_));
    for (auto &col : col_lr_) s.read((char*)&col,sizeof(col));
    //cout << "col_lr=" << col_lr_[0] << " " << col_lr_[1] << endl;

    for (auto &sigma : sigma_lr_) sigma.read(s);
    //Print(sigma_lr_[0]);
    //Print(sigma_lr_[1]);
    //replace indices of sigma_lr_ such that they share the indices of double_layer_tensors_. Then recombine sigma_lr_ to a vector
    match_indices_sigma_lr();
    recombine_sigma_lr();

    sigma_b_.read(s);

    for (auto eigval : density_mat_spectrum_) s.read((char*)&eigval,sizeof(eigval));
}
template 
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::read(std::istream &s);
template 
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::read(std::istream &s);


template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::write(std::ostream &s) const
{
    Double_Layer_PEPSt<TensorT>::write(s);

    s.write((char*)&cutting_col_,sizeof(cutting_col_));
    for (auto col : col_lr_) s.write((char*)&col,sizeof(col));

    //we should only write decombined sigma_lr_
    for (const auto &sigma : sigma_lr_) sigma.write(s);

    sigma_b_.write(s);

    for (auto eigval : density_mat_spectrum_) s.write((char*)&eigval,sizeof(eigval));
}
template
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::write(std::ostream &s) const;
template
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::write(std::ostream &s) const;



double honeycomb_cylinder_SzSz_correlator(const std::array<Coordinate,2> &acting_sites_coord, Cylinder_Square_Double_Layer_PEPSt<IQTensor> &square_cylinder_double_peps)
{
    //we create Sz operator for honeycomb sites, and then transfer it to square lattice
    std::vector<IQTensor> square_Sz_operators;
    std::vector<int> square_acting_sites;
    //col_lr stores the closest cols to leftmost and rightmost cols
    std::array<int,2> col_lr={square_cylinder_double_peps.lattice().n_uc()[0]-1,0};

    for (auto site_coord : acting_sites_coord)
    {
        auto site_no=square_cylinder_double_peps.lattice().site_coord_to_list(site_coord[0],site_coord[1],0);
        square_acting_sites.push_back(site_no);

        //update col_lr
        if ((site_coord[0]-1)<col_lr[0]) col_lr[0]=site_coord[0]-1;
        if ((site_coord[0]+1)>col_lr[1]) col_lr[1]=site_coord[0]+1;

        //create Sz operator and Id operator for one uc of honeycomb
        std::vector<IQIndex> Sz_legs={Spin_leg(std::vector<int>{0,1},"Sz_leg 0",Out,Site), Spin_leg(std::vector<int>{0,1},"Sz_leg 1",Out,Site)};
        IQTensor Sz_operator(dag(Sz_legs[0])(1),prime(Sz_legs[0])(1)),
                 Id_operator(dag(Sz_legs[1])(1),prime(Sz_legs[1])(1));

        Sz_operator(dag(Sz_legs[0])(2),prime(Sz_legs[0])(2))=-1;
        Id_operator(dag(Sz_legs[1])(2),prime(Sz_legs[1])(2))=1;
        
        //transfer to square Sz operator in one site
        auto combined_SzId_operator=Sz_operator*Id_operator;
        IQCombiner Sz_legs_combiner;

        if (site_coord[2]==0) 
        {
            Sz_legs_combiner.addleft(Sz_legs[0]);
            Sz_legs_combiner.addleft(Sz_legs[1]);
        }
        if (site_coord[2]==1) 
        {
            Sz_legs_combiner.addleft(Sz_legs[1]);
            Sz_legs_combiner.addleft(Sz_legs[0]);
        }

        //TODO: check the order of combined operator
        Sz_legs_combiner.init("square_Sz_leg",Site,Out);
        combined_SzId_operator=combined_SzId_operator*dag(Sz_legs_combiner)*prime(Sz_legs_combiner);
        square_Sz_operators.push_back(combined_SzId_operator);
    }

    //move sigma_lr_ to col_lr_ to reduce contraction time, and then do contraction
    square_cylinder_double_peps.move_sigma_lr(col_lr);

    //obtain the sandwiched double peps, calculate wavefunction norm
    auto square_sandwiched_tensors=square_cylinder_double_peps.double_layer_tensors();
    double wf_norm=square_cylinder_double_peps.sandwiched_peps_norm(square_sandwiched_tensors);

    square_cylinder_double_peps.obtain_tensor_sandwich_single_site_operators(square_Sz_operators,square_acting_sites,square_sandwiched_tensors);

    double unnormlized_expect_val=square_cylinder_double_peps.sandwiched_peps_norm(square_sandwiched_tensors);

    return unnormlized_expect_val/wf_norm;
}
