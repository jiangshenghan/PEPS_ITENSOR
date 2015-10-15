
#include "double_layer_peps.h"

template <class TensorT>
Double_Layer_PEPSt<TensorT>::Double_Layer_PEPSt(const PEPSt<TensorT> &peps):
    single_layer_peps_(peps),
    layered_site_tensors_(peps.n_sites_total()),
    virt_leg_combiners_(peps.n_sites_total())
{
    std::vector<TensorT> combined_site_tensors;
    obtain_combined_site_tensors(peps,combined_site_tensors);
    obtain_layered_tensors_with_combined_legs(combined_site_tensors);
}
template
Double_Layer_PEPSt<ITensor>::Double_Layer_PEPSt(const PEPSt<ITensor> &peps);
template
Double_Layer_PEPSt<IQTensor>::Double_Layer_PEPSt(const PEPSt<IQTensor> &peps);


template <class TensorT>
void Double_Layer_PEPSt<TensorT>::obtain_combined_site_tensors(const PEPSt<TensorT> &peps, std::vector<TensorT> &combined_site_tensors)
{
    for (int sitei=0; sitei<peps.n_sites_total(); sitei++)
    {
        auto site_tensor=peps.site_tensors(sitei);

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
                site_tensor*=peps.bond_tensors(bond_no);
            }
        }
        //multiply neighbouring boundary tensors
        for (const auto &boundary_no : this->lattice().site_neighbour_boundary(sitei))
        {
            site_tensor*=peps.boundary_tensors(boundary_no);
        }
        clean(site_tensor);

        combined_site_tensors.push_back(site_tensor);

        //cout << "Combined tensor " << sitei << endl << site_tensor << endl;
    }
}
template
void Double_Layer_PEPSt<ITensor>::obtain_combined_site_tensors(const PEPSt<ITensor> &peps, std::vector<ITensor> &combined_site_tensors);
template
void Double_Layer_PEPSt<IQTensor>::obtain_combined_site_tensors(const PEPSt<IQTensor> &peps, std::vector<IQTensor> &combined_site_tensors);


template <class TensorT>
void Double_Layer_PEPSt<TensorT>::obtain_layered_tensors_with_combined_legs(const std::vector<TensorT> &lower_tensors)
{
    //create combiners for virt_leg and prime(dag(virt_leg)), where virt_leg connects sitei and neighbour_sites[neighbour_i]
    //every combiner should appear twice, so we stores combined indiceswhich only appear once, and delete it if it already appear twice
    std::vector<IndexT> combined_indices;
    std::vector<CombinerT> indices_combiners;

    int combine_i=0;
    for (int sitei=0; sitei<this->lattice().n_sites_total(); sitei++)
    {
        for (const auto &virt_leg : lower_tensors[sitei].indices())
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
    std::vector<TensorT> upper_tensors(lower_tensors);
    for (auto &tensor : upper_tensors)
    {
        tensor.dag();
        tensor.prime(Link);
    }

    for (int sitei=0; sitei<this->lattice().n_sites_total(); sitei++)
    {
        layered_site_tensors_[sitei]=lower_tensors[sitei]*upper_tensors[sitei];
        for (const auto &combiner : virt_leg_combiners_[sitei])
        {
            layered_site_tensors_[sitei]=layered_site_tensors_[sitei]*combiner;
        }
        clean(layered_site_tensors_[sitei]);
    }

}
template
void Double_Layer_PEPSt<ITensor>::obtain_layered_tensors_with_combined_legs(const std::vector<ITensor> &lower_tensors);
template
void Double_Layer_PEPSt<IQTensor>::obtain_layered_tensors_with_combined_legs(const std::vector<IQTensor> &lower_tensors);



//
//class Cylinder_Square_Double_Layer_PEPSt
//
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
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::obtain_sigma_lr_iterative(int left_end_col, int right_end_col, bool do_decombine)
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

    //To do decombine, we need to make sure Ly<=8
    if (do_decombine)
    {
        for (int lr_no=0; lr_no<2; lr_no++)
        {
            int n_rows=this->lattice().n_uc()[1];
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

    cout << "\n========================================\n" << endl;
    cout << "sigma_left for contraction to " << col_lr_[0] << " cols:" << endl;
    PrintDat(sigma_lr_[0]);
    cout << "sigma_right for contraction to " << col_lr_[1] << " cols:" << endl;
    PrintDat(sigma_lr_[1]);
    cout << "\n========================================\n" << endl;

}
template
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::obtain_sigma_lr_iterative(int left_end_col, int right_end_col, bool do_decombine);
template
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::obtain_sigma_lr_iterative(int left_end_col, int right_end_col, bool do_decombine);


template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::snake_walking_boundary_col()
{
    //Initialize the first col (left boundary) from down to up
    int sitei=0;
    iterative_combiners_[0][0]=CombinerT(commonIndex(this->layered_site_tensors_[sitei],dag(this->layered_site_tensors_[sitei+1])));
    sigma_lr_[0]=this->layered_site_tensors_[sitei]*iterative_combiners_[0][0];
    clean(sigma_lr_[0]);

    for (int rowi=1; rowi<this->lattice().n_uc()[1]; rowi++)
    {
        sitei+=this->lattice().n_uc()[0];
        auto combining_indice=commonIndex(this->layered_site_tensors_[sitei],dag(this->layered_site_tensors_[sitei+1]));
        iterative_combiners_[0][rowi]=CombinerT(iterative_combiners_[0][rowi-1].right(),combining_indice);
        sigma_lr_[0]*=this->layered_site_tensors_[sitei];
        sigma_lr_[0]=sigma_lr_[0]*iterative_combiners_[0][rowi];
    }


    //Initialize the last col (right boundary) from up to down
    sitei=this->lattice().n_sites_total()-1;
    iterative_combiners_[1][this->lattice().n_uc()[1]-1]=CombinerT(commonIndex(this->layered_site_tensors_[sitei],dag(this->layered_site_tensors_[sitei-1])));
    sigma_lr_[1]=this->layered_site_tensors_[sitei]*iterative_combiners_[1][this->lattice().n_uc()[1]-1];
    clean(sigma_lr_[1]);

    for (int rowi=this->lattice().n_uc()[1]-2; rowi>=0; rowi--)
    {
        sitei-=this->lattice().n_uc()[0];
        auto combining_indice=commonIndex(this->layered_site_tensors_[sitei],dag(this->layered_site_tensors_[sitei-1]));
        iterative_combiners_[1][rowi]=CombinerT(iterative_combiners_[1][rowi+1].right(),combining_indice);
        sigma_lr_[1]*=this->layered_site_tensors_[sitei];
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
    cout << "first ten elems of sigma left:" << endl;
    print_vector_nonzero_elem(sigma_lr_[0],10);
    cout << "first ten elems of sigma right:" << endl;
    print_vector_nonzero_elem(sigma_lr_[1],10);
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
        sigma_lr_[lr_no]*=this->layered_site_tensors_[sitei];

        //combine the indice of new tensor which will be multiplied next time
        auto combining_indice=commonIndex(this->layered_site_tensors_[sitei],dag(this->layered_site_tensors_[sitei+horizontal_dir]));
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
    cout << "First ten elems of sigma_lr[" << lr_no << "] are" << endl;
    print_vector_nonzero_elem(sigma_lr_[lr_no],10);
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



template <class TensorT>
void Cylinder_Square_Double_Layer_PEPSt<TensorT>::obtain_transfer_matrix(int coli)
{
    int sitei=this->lattice().site_coord_to_list(coli,0,0);
    //leg_combiners[0/1] is for left/right legs of transfer matrix
    std::array<std::vector<CombinerT>,2> leg_combiners;
    std::array<IndexT,2> lr_legs;

    lr_legs[0]=commonIndex(this->layered_site_tensors_[sitei],dag(this->layered_site_tensors_[sitei-1]));
    lr_legs[1]=commonIndex(this->layered_site_tensors_[sitei],dag(this->layered_site_tensors_[sitei+1]));

    leg_combiners[0].push_back(CombinerT(lr_legs[0]));
    leg_combiners[1].push_back(CombinerT(lr_legs[1]));
    
    transfer_mat_=this->layered_site_tensors_[sitei]*leg_combiners[0][0]*leg_combiners[1][0];

    for (int rowi=1; rowi<this->lattice().n_uc()[1]; rowi++)
    {
        sitei+=this->lattice().n_uc()[0];

        lr_legs[0]=commonIndex(this->layered_site_tensors_[sitei],dag(this->layered_site_tensors_[sitei-1]));
        lr_legs[1]=commonIndex(this->layered_site_tensors_[sitei],dag(this->layered_site_tensors_[sitei+1]));

        leg_combiners[0].push_back(CombinerT(lr_legs[0],leg_combiners[0][rowi-1].right()));
        leg_combiners[1].push_back(CombinerT(lr_legs[1],leg_combiners[1][rowi-1].right()));

        transfer_mat_*=this->layered_site_tensors_[sitei];
        transfer_mat_=transfer_mat_*leg_combiners[0][rowi]*leg_combiners[1][rowi];
    }

    //We make transfer_mat_ as matrix --I--mat--I'--
    auto oind=leg_combiners[1][this->lattice().n_uc()[1]-1].right(),
         nind=prime(leg_combiners[0][this->lattice().n_uc()[1]-1].right());
    if ((oind.dir()-nind.dir())!=0) nind.dag();
    transfer_mat_.replaceIndex(oind,nind);

    transfer_mat_/=transfer_mat_.norm();
    clean(transfer_mat_);
    

    //cout << "\n----------------------------\n" << endl;
    //cout << "transfer matrix for col " << coli << ":" << endl;
    //PrintDat(transfer_mat_);
    //cout << "\n----------------------------\n" << endl;
}
template
void Cylinder_Square_Double_Layer_PEPSt<ITensor>::obtain_transfer_matrix(int coli);
template
void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::obtain_transfer_matrix(int coli);


//template <class TensorT>
//void Cylinder_Square_Double_Layer_PEPSt<TensorT>::read(std::istream &s)
//{
//}
//template 
//void Cylinder_Square_Double_Layer_PEPSt<ITensor>::read(std::istream &s);
//template 
//void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::read(std::istream &s);
//
//
//template <class TensorT>
//void Cylinder_Square_Double_Layer_PEPSt<TensorT>::write(std::ostream &s) const
//{
//}
//template
//void Cylinder_Square_Double_Layer_PEPSt<ITensor>::write(std::ostream &s) const;
//template
//void Cylinder_Square_Double_Layer_PEPSt<IQTensor>::write(std::ostream &s) const;
