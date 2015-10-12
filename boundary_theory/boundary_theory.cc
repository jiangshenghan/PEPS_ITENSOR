
#include "boundary_theory.h"

template <class TensorT>
Boundary_Theory<TensorT>::Boundary_Theory(const PEPSt<TensorT> &square_peps, int cutting_col, const std::string &method):
    n_rows_(square_peps.n_uc()[1]),
    n_cols_(square_peps.n_uc()[0]),
    cutting_col_(cutting_col),
    square_layered_peps_(square_peps),
    iterative_combiners_{{std::vector<CombinerT>(n_rows_),std::vector<CombinerT>(n_rows_)}}
{
    if (cutting_col_==-1) 
        cutting_col_=n_cols_/2;

    //obtain_transfer_matrix();

    if (method.find("iterative")!=std::string::npos) 
        obtain_boundary_theory_iterative();

}
template
Boundary_Theory<ITensor>::Boundary_Theory(const PEPSt<ITensor> &square_peps, int cutting_col, const std::string &method);
template
Boundary_Theory<IQTensor>::Boundary_Theory(const PEPSt<IQTensor> &square_peps, int cutting_col, const std::string &method);

template <class TensorT>
void Boundary_Theory<TensorT>::obtain_boundary_theory_iterative()
{
    snake_walking_boundary_col();

    for (int lcoli=1; lcoli<cutting_col_; )
    {
        snake_walking_bulk_col(lcoli,1,-1);
        lcoli++;

        if (lcoli==cutting_col_)
        {
            break;
        }

        snake_walking_bulk_col(lcoli,1,1);
        lcoli++;
    }

    for (int rcoli=n_cols_-2; rcoli>=cutting_col_; )
    {
        snake_walking_bulk_col(rcoli,-1,1);
        rcoli--;

        if (rcoli<cutting_col_)
        {
            break;
        }

        snake_walking_bulk_col(rcoli,-1,-1);
        rcoli--;
    }

    from_sigma_vec_to_mat();
    from_sigma_lr_to_sigma_b();
    obtain_density_matrix_spectrum();
}
template
void Boundary_Theory<ITensor>::obtain_boundary_theory_iterative();
template 
void Boundary_Theory<IQTensor>::obtain_boundary_theory_iterative();


template <class TensorT>
void Boundary_Theory<TensorT>::snake_walking_boundary_col()
{
    //Initialize the first col (left boundary) from down to up
    int sitei=0;
    iterative_combiners_[0][0]=CombinerT(commonIndex(square_layered_peps_.layered_site_tensors(sitei),dag(square_layered_peps_.layered_site_tensors(sitei+1))));
    sigma_lr_[0]=square_layered_peps_.layered_site_tensors(sitei)*iterative_combiners_[0][0];
    clean(sigma_lr_[0]);

    for (int rowi=1; rowi<n_rows_; rowi++)
    {
        sitei+=n_cols_;
        auto combining_indice=commonIndex(square_layered_peps_.layered_site_tensors(sitei),dag(square_layered_peps_.layered_site_tensors(sitei+1)));
        iterative_combiners_[0][rowi]=CombinerT(iterative_combiners_[0][rowi-1].right(),combining_indice);
        sigma_lr_[0]*=square_layered_peps_.layered_site_tensors(sitei);
        sigma_lr_[0]=sigma_lr_[0]*iterative_combiners_[0][rowi];
    }


    //Initialize the last col (right boundary) from up to down
    sitei=n_cols_*n_rows_-1;
    iterative_combiners_[1][n_rows_-1]=CombinerT(commonIndex(square_layered_peps_.layered_site_tensors(sitei),dag(square_layered_peps_.layered_site_tensors(sitei-1))));
    sigma_lr_[1]=square_layered_peps_.layered_site_tensors(sitei)*iterative_combiners_[1][n_rows_-1];
    clean(sigma_lr_[1]);

    for (int rowi=n_rows_-2; rowi>=0; rowi--)
    {
        sitei-=n_cols_;
        auto combining_indice=commonIndex(square_layered_peps_.layered_site_tensors(sitei),dag(square_layered_peps_.layered_site_tensors(sitei-1)));
        iterative_combiners_[1][rowi]=CombinerT(iterative_combiners_[1][rowi+1].right(),combining_indice);
        sigma_lr_[1]*=square_layered_peps_.layered_site_tensors(sitei);
        sigma_lr_[1]=sigma_lr_[1]*iterative_combiners_[1][rowi];
    }

    for (auto &sigma : sigma_lr_)
    {
        sigma/=sigma.norm();
        clean(sigma);
    }

    cout << "\n----------------------------\n" << endl;
    cout << "Contraction for boundary cols:" << endl;
    //cout << "iterative combiners for left part:" << endl;
    //for (const auto &combiner : iterative_combiners_[0]) cout << combiner;
    //cout << "iterative combiners for right part:" << endl;
    //for (const auto &combiner : iterative_combiners_[1]) cout << combiner;
    cout << "sigma left:" << endl;
    PrintDat(sigma_lr_[0]);
    cout << "sigma right:" << endl;
    PrintDat(sigma_lr_[1]);
    cout << "\n----------------------------\n" << endl;

}
template
void Boundary_Theory<ITensor>::snake_walking_boundary_col();
template
void Boundary_Theory<IQTensor>::snake_walking_boundary_col();


template <class TensorT>
void Boundary_Theory<TensorT>::snake_walking_bulk_col(int coli, int horizontal_dir, int vertical_dir)
{
    int lr_no=(1-horizontal_dir)/2;
    int start_row=(1-vertical_dir)/2*(n_rows_-1);
    int rowi=start_row;
    while (rowi>=0 && rowi<n_rows_)
    {
        int sitei=square_layered_peps_.lattice().site_coord_to_list(coli,rowi,0);
        //decombine indice to be multiplied, and then mutiply a new tensor
        sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(iterative_combiners_[lr_no][rowi]);
        sigma_lr_[lr_no]*=square_layered_peps_.layered_site_tensors(sitei);

        //combine the indice of new tensor which will be multiplied next time
        auto combining_indice=commonIndex(square_layered_peps_.layered_site_tensors(sitei),dag(square_layered_peps_.layered_site_tensors(sitei+horizontal_dir)));
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

    cout << "\n----------------------------\n" << endl;
    cout << "Iterative contraction for bulk col:" << endl;
    cout << "Horizontal Direction: " << horizontal_dir << endl
         << "Vertical Direction: " << vertical_dir << endl
         << "Col: " << coli << endl << endl;
    //for (const auto &combiner : iterative_combiners_[lr_no]) cout << combiner;
    PrintDat(sigma_lr_[lr_no]);
    cout << "\n----------------------------\n" << endl;

}
template
void Boundary_Theory<ITensor>::snake_walking_bulk_col(int coli, int horizontal_dir, int vertical_dir);
template
void Boundary_Theory<IQTensor>::snake_walking_bulk_col(int coli, int horizontal_dir, int vertical_dir);


template <class TensorT>
void Boundary_Theory<TensorT>::from_sigma_vec_to_mat()
{
    std::array<int,2> vertical_dirs{-2*(cutting_col_%2)+1,2*((n_cols_-cutting_col_)%2)-1};

    std::array<std::vector<CombinerT>,2> upper_combiners{std::vector<CombinerT>(n_rows_),std::vector<CombinerT>(n_rows_)}, 
                                         lower_combiners{std::vector<CombinerT>(n_rows_),std::vector<CombinerT>(n_rows_)};
    for (int lr_no=0; lr_no<2; lr_no++)
    {
        int start_row=(1-vertical_dirs[lr_no])/2*(n_rows_-1);
        int rowi=start_row;

        while (rowi>=0 && rowi<n_rows_)
        {
            sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(iterative_combiners_[lr_no][rowi]);

            //decombine upper leg and lower leg
            int sitei=square_layered_peps_.lattice().site_coord_to_list(cutting_col_-1+lr_no,rowi,0);
            for (const auto &leg_combiner : square_layered_peps_.virt_leg_combiners(sitei))
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

    cout << "\n----------------------------\n" << endl;
    cout << "Output for sigma_l and sigma_r:" << endl;
    //cout << "left part:" << endl;
    //for (int rowi=0; rowi<n_rows_; rowi++)
    //{
    //    cout << lower_combiners[0][rowi];
    //    cout << upper_combiners[0][rowi];
    //}
    PrintDat(sigma_lr_[0]);
    //cout << "right part:" << endl;
    //for (int rowi=0; rowi<n_rows_; rowi++)
    //{
    //    cout << lower_combiners[1][rowi];
    //    cout << upper_combiners[1][rowi];
    //}
    PrintDat(sigma_lr_[1]);
    cout << "\n----------------------------\n" << endl;
}
template
void Boundary_Theory<ITensor>::from_sigma_vec_to_mat();
template
void Boundary_Theory<IQTensor>::from_sigma_vec_to_mat();


//TODO: sigma_lr_'s may not have a good quantum number due to boundary
//condition? In this case diagHermitian does not work for IQTensor?
template <class TensorT>
void Boundary_Theory<TensorT>::from_sigma_lr_to_sigma_b()
{
    //PrintDat(sigma_lr_[0]);
    //get eigenvalues of sigma_l
    TensorT U_l, D_l;
    diagHermitian(sigma_lr_[0],U_l,D_l);

    cout << "sigma_l singular value:" << endl;
    PrintDat(D_l);

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

    cout << "sqrt of sigma_l singular value:" << endl;
    PrintDat(D_l);
    cout << "sqrt of sigma_l:" << endl;
    PrintDat(sqrt_sigma_l);

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

    cout << "\n---------------------------------\n" << endl;
    PrintDat(sigma_b_);
    cout << "\n---------------------------------\n" << endl;
}
template
void Boundary_Theory<ITensor>::from_sigma_lr_to_sigma_b();
template
void Boundary_Theory<IQTensor>::from_sigma_lr_to_sigma_b();


template <class TensorT>
void Boundary_Theory<TensorT>::obtain_density_matrix_spectrum()
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
void Boundary_Theory<ITensor>::obtain_density_matrix_spectrum();
template
void Boundary_Theory<IQTensor>::obtain_density_matrix_spectrum();


template <class TensorT>
std::vector<double> Boundary_Theory<TensorT>::entanglement_spectrum()
{
    std::vector<double> entanglement_spectrum;
    std::transform(density_mat_spectrum_.begin(),density_mat_spectrum_.end(),entanglement_spectrum.begin(),[](double eigval){ return std::exp(-eigval); });
    return entanglement_spectrum;
}
template
std::vector<double> Boundary_Theory<ITensor>::entanglement_spectrum();
template
std::vector<double> Boundary_Theory<IQTensor>::entanglement_spectrum();


template <class TensorT>
double Boundary_Theory<TensorT>::entanglement_entropy_vN()
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
double Boundary_Theory<ITensor>::entanglement_entropy_vN();
template
double Boundary_Theory<IQTensor>::entanglement_entropy_vN();


template <class TensorT>
double Boundary_Theory<TensorT>::entanglement_entropy_Renyi(double renyi_n)
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
double Boundary_Theory<ITensor>::entanglement_entropy_Renyi(double renyi_n);
template
double Boundary_Theory<IQTensor>::entanglement_entropy_Renyi(double renyi_n);



template <class TensorT>
void Boundary_Theory<TensorT>::obtain_transfer_matrix(int coli)
{
    int sitei=square_layered_peps_.lattice().site_coord_to_list(coli,0,0);
    //leg_combiners[0/1] is for left/right legs of transfer matrix
    std::array<std::vector<CombinerT>,2> leg_combiners;
    std::array<IndexT,2> lr_legs;

    lr_legs[0]=commonIndex(square_layered_peps_.layered_site_tensors(sitei),dag(square_layered_peps_.layered_site_tensors(sitei-1)));
    lr_legs[1]=commonIndex(square_layered_peps_.layered_site_tensors(sitei),dag(square_layered_peps_.layered_site_tensors(sitei+1)));

    leg_combiners[0].push_back(CombinerT(lr_legs[0]));
    leg_combiners[1].push_back(CombinerT(lr_legs[1]));
    
    transfer_mat_=square_layered_peps_.layered_site_tensors(sitei)*leg_combiners[0][0]*leg_combiners[1][0];

    for (int rowi=1; rowi<n_rows_; rowi++)
    {
        sitei+=n_cols_;

        lr_legs[0]=commonIndex(square_layered_peps_.layered_site_tensors(sitei),dag(square_layered_peps_.layered_site_tensors(sitei-1)));
        lr_legs[1]=commonIndex(square_layered_peps_.layered_site_tensors(sitei),dag(square_layered_peps_.layered_site_tensors(sitei+1)));

        leg_combiners[0].push_back(CombinerT(lr_legs[0],leg_combiners[0][rowi-1].right()));
        leg_combiners[1].push_back(CombinerT(lr_legs[1],leg_combiners[1][rowi-1].right()));

        transfer_mat_*=square_layered_peps_.layered_site_tensors(sitei);
        transfer_mat_=transfer_mat_*leg_combiners[0][rowi]*leg_combiners[1][rowi];
    }

    //We make transfer_mat_ as matrix --I--mat--I'--
    auto oind=leg_combiners[1][n_rows_-1].right(),
         nind=prime(leg_combiners[0][n_rows_-1].right());
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
void Boundary_Theory<ITensor>::obtain_transfer_matrix(int coli);
template
void Boundary_Theory<IQTensor>::obtain_transfer_matrix(int coli);


//template <class TensorT>
//void Boundary_Theory<TensorT>::apply_transfer_matrix_to_sigma_lr_vec(int lr_no)
//{
//    auto transfer_mat_temp=transfer_mat_;
//    auto transfer_mat_indices=transfer_mat_.indices();
//    IndexT oind, nind=dag(sigma_lr_[lr_no].indices()[0]);
//
//    if(transfer_mat_indices[0].primeLevel()==lr_no)
//    {
//        oind=transfer_mat_indices[0];
//    }
//    else
//    {
//        oind=transfer_mat_indices[1];
//    }
//    transfer_mat_.replaceIndex(oind,nind);
//
//    auto applied_vec=transfer_mat_*sigma_lr_[lr_no];
//
//    cout << "\n----------------------------\n" << endl;
//    cout << "transfer_mat*sigma/sigma:" << endl;
//    for (int val=1; val<=nind.m(); val++)
//    {
//        auto elem=sigma_lr_[lr_no](dag(nind)(val));
//        if (elem<EPSILON) continue;
//        cout << applied_vec(applied_vec.indices()[0](val))/elem << endl;
//    }
//    cout << "\n----------------------------\n" << endl;
//}
//template 
//void Boundary_Theory<ITensor>::apply_transfer_matrix_to_sigma_lr_vec(int lr_no);
//template 
//void Boundary_Theory<IQTensor>::apply_transfer_matrix_to_sigma_lr_vec(int lr_no);
