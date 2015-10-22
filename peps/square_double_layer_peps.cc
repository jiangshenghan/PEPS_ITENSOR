
#include "square_double_layer_peps.h"


//
//class Square_Double_Layer_PEPSt
//
template <class TensorT>
Square_Double_Layer_PEPSt<TensorT>::Square_Double_Layer_PEPSt(const Lattice_Base &square_lattice):
    Double_Layer_PEPSt<TensorT>(square_lattice),
    col_lr_{{0,square_lattice.n_uc()[0]-1}},
    lower_combiners_{{std::vector<CombinerT>(square_lattice.n_uc()[1]),std::vector<CombinerT>(square_lattice.n_uc()[1])}},
    upper_combiners_{{std::vector<CombinerT>(square_lattice.n_uc()[1]),std::vector<CombinerT>(square_lattice.n_uc()[1])}}
{
    //we always focus on long cylinder or ribbon geometry, where Ly is the short direction
    if (square_lattice.name().find("torus")!=std::string::npos)
    {
        cout << "We only focus on cylinder or ribbon!" << endl;
        exit(EXIT_FAILURE);
    }
    cutting_col_=this->lattice().n_uc()[0]/2;
}
template
Square_Double_Layer_PEPSt<ITensor>::Square_Double_Layer_PEPSt(const Lattice_Base &square_lattice);
template
Square_Double_Layer_PEPSt<IQTensor>::Square_Double_Layer_PEPSt(const Lattice_Base &square_lattice);

template <class TensorT>
Square_Double_Layer_PEPSt<TensorT>::Square_Double_Layer_PEPSt(const PEPSt<TensorT> &square_peps, int cutting_col):
    Double_Layer_PEPSt<TensorT>(square_peps),
    cutting_col_(cutting_col),
    col_lr_{{0,square_peps.lattice().n_uc()[0]-1}},
    lower_combiners_{{std::vector<CombinerT>(square_peps.lattice().n_uc()[1]),std::vector<CombinerT>(square_peps.lattice().n_uc()[1])}},
    upper_combiners_{{std::vector<CombinerT>(square_peps.lattice().n_uc()[1]),std::vector<CombinerT>(square_peps.lattice().n_uc()[1])}}
{
    //we always focus on long cylinder or ribbon geometry, where Ly is the short direction
    if (square_peps.name().find("torus")!=std::string::npos)
    {
        cout << "We only focus on cylinder or ribbon!" << endl;
        exit(EXIT_FAILURE);
    }
    if (cutting_col_==-1) 
        cutting_col_=this->lattice().n_uc()[0]/2;
}
template
Square_Double_Layer_PEPSt<ITensor>::Square_Double_Layer_PEPSt(const PEPSt<ITensor> &square_peps, int cutting_col);
template
Square_Double_Layer_PEPSt<IQTensor>::Square_Double_Layer_PEPSt(const PEPSt<IQTensor> &square_peps, int cutting_col);

template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::obtain_boundary_theory_iterative()
{
    obtain_sigma_lr_iterative(cutting_col_-1,cutting_col_);
    from_sigma_lr_to_sigma_b();
    obtain_density_matrix_spectrum();
}
template
void Square_Double_Layer_PEPSt<ITensor>::obtain_boundary_theory_iterative();
template 
void Square_Double_Layer_PEPSt<IQTensor>::obtain_boundary_theory_iterative();

template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::obtain_sigma_lr_iterative(int left_end_col, int right_end_col)
{
    snake_walking_boundary_col();

    for (col_lr_[0]=1; col_lr_[0]<=left_end_col; )
    {
        snake_walking_bulk_col(col_lr_[0],1,-1,true);
        col_lr_[0]++;

        if (col_lr_[0]>left_end_col) break;

        snake_walking_bulk_col(col_lr_[0],1,1,true);
        col_lr_[0]++;
    }
    col_lr_[0]=left_end_col;

    for (col_lr_[1]=this->lattice().n_uc()[0]-2; col_lr_[1]>=right_end_col; )
    {
        snake_walking_bulk_col(col_lr_[1],-1,1,true);
        col_lr_[1]--;

        if (col_lr_[1]<right_end_col) break;

        snake_walking_bulk_col(col_lr_[1],-1,-1,true);
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
void Square_Double_Layer_PEPSt<ITensor>::obtain_sigma_lr_iterative(int left_end_col, int right_end_col);
template
void Square_Double_Layer_PEPSt<IQTensor>::obtain_sigma_lr_iterative(int left_end_col, int right_end_col);


template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::snake_walking_boundary_col()
{
    int n_rows=this->lattice().n_uc()[1], 
        n_cols=this->lattice().n_uc()[0];
    //Initialize the first col (left boundary) from down to up
    int sitei=0;
    auto double_virt_ind=commonIndex(this->double_layer_tensors(sitei),dag(this->double_layer_tensors(sitei+1)));
    IndexT lower_ind;
    auto double_virt_leg_combiner=this->decombine_double_virt_indice(sitei,double_virt_ind,lower_ind);
    lower_combiners_[0][0]=CombinerT(lower_ind);
    upper_combiners_[0][0]=prime(dag(lower_combiners_[0][0]));

    //Print(lower_ind);
    //Print(lower_combiners_[0][0]);
    //Print(upper_combiners_[0][0]);
    //Print(double_virt_leg_combiner);

    sigma_lr_[0]=this->double_layer_tensors_[sitei]*dag(double_virt_leg_combiner)*lower_combiners_[0][0]*upper_combiners_[0][0];
    clean(sigma_lr_[0]);

    for (int rowi=1; rowi<n_rows; rowi++)
    {
        sitei+=n_cols;
        double_virt_ind=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei+1]));
        double_virt_leg_combiner=this->decombine_double_virt_indice(sitei,double_virt_ind,lower_ind);
        lower_combiners_[0][rowi]=CombinerT(lower_combiners_[0][rowi-1].right(),lower_ind);
        upper_combiners_[0][rowi]=prime(dag(lower_combiners_[0][rowi]));
        sigma_lr_[0]*=this->double_layer_tensors_[sitei]*dag(double_virt_leg_combiner);
        sigma_lr_[0]=sigma_lr_[0]*lower_combiners_[0][rowi]*upper_combiners_[0][rowi];
    }


    //Initialize the last col (right boundary) from up to down
    sitei=this->lattice().n_sites_total()-1;
    double_virt_ind=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei-1]));
    double_virt_leg_combiner=this->decombine_double_virt_indice(sitei,double_virt_ind,lower_ind);
    lower_combiners_[1][n_rows-1]=CombinerT(lower_ind);
    upper_combiners_[1][n_rows-1]=prime(dag(lower_combiners_[1][n_rows-1]));
    sigma_lr_[1]=this->double_layer_tensors_[sitei]*dag(double_virt_leg_combiner)*lower_combiners_[1][n_rows-1]*upper_combiners_[1][n_rows-1];
    clean(sigma_lr_[1]);

    for (int rowi=n_rows-2; rowi>=0; rowi--)
    {
        sitei-=n_cols;
        double_virt_ind=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei-1]));
        double_virt_leg_combiner=this->decombine_double_virt_indice(sitei,double_virt_ind,lower_ind);
        lower_combiners_[1][rowi]=CombinerT(lower_combiners_[1][rowi+1].right(),lower_ind);
        upper_combiners_[1][rowi]=prime(dag(lower_combiners_[1][rowi]));
        sigma_lr_[1]*=this->double_layer_tensors_[sitei]*dag(double_virt_leg_combiner);
        sigma_lr_[1]=sigma_lr_[1]*lower_combiners_[1][rowi]*upper_combiners_[1][rowi];
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
void Square_Double_Layer_PEPSt<ITensor>::snake_walking_boundary_col();
template
void Square_Double_Layer_PEPSt<IQTensor>::snake_walking_boundary_col();


template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::snake_walking_bulk_col(int coli, int horizontal_dir, int vertical_dir, bool do_normalize)
{
    int n_rows=this->lattice().n_uc()[1],
        n_cols=this->lattice().n_uc()[0];
    int lr_no=(1-horizontal_dir)/2;
    int start_row=(1-vertical_dir)/2*(n_rows-1);
    int rowi=start_row;
    while (rowi>=0 && rowi<this->lattice().n_uc()[1])
    {
        int sitei=this->lattice().site_coord_to_list(coli,rowi,0);

        //cout << "Col " << coli << " Row " << rowi << endl;
        //Print(sigma_lr_[lr_no]);
        //Print(lower_combiners_[lr_no][rowi]);
        //Print(upper_combiners_[lr_no][rowi]);
        //Print(sigma_lr_[lr_no]*dag(lower_combiners_[lr_no][rowi])*dag(upper_combiners_[lr_no][rowi]));

        //decombine indice to be multiplied, and then mutiply a new tensor
        sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(lower_combiners_[lr_no][rowi])*dag(upper_combiners_[lr_no][rowi]);
        sigma_lr_[lr_no]*=this->double_layer_tensors_[sitei]*dag(this->virt_leg_combiners({sitei,sitei-horizontal_dir}));


        //combine the indice of new tensor which will be multiplied next time
        IndexT lower_ind;
        auto double_virt_leg_combiner=decombine_double_virt_indice({sitei,sitei+horizontal_dir},lower_ind);
        if (rowi==start_row) 
        {
            lower_combiners_[lr_no][rowi]=CombinerT(lower_ind);
            upper_combiners_[lr_no][rowi]=prime(dag(lower_combiners_[lr_no][rowi]));
        }
        else
        {
            lower_combiners_[lr_no][rowi]=CombinerT(lower_combiners_[lr_no][rowi-vertical_dir].right(),lower_ind);
            upper_combiners_[lr_no][rowi]=prime(dag(lower_combiners_[lr_no][rowi]));
        }
        //PrintDat(sigma_lr_[lr_no]);
        sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(double_virt_leg_combiner)*lower_combiners_[lr_no][rowi]*upper_combiners_[lr_no][rowi];

        rowi+=vertical_dir;
    }

    if (do_normalize)
    {
        sigma_lr_[lr_no]/=sigma_lr_[lr_no].norm();
    }

    clean(sigma_lr_[lr_no]);

    //cout << "\n========================================\n" << endl;
    //cout << "Iterative contraction for bulk col:" << endl;
    //cout << "Horizontal Direction: " << horizontal_dir << endl
    //     << "Vertical Direction: " << vertical_dir << endl
    //     << "Col: " << coli << endl << endl;
    //for (const auto &combiner : iterative_combiners_[lr_no]) cout << combiner;
    //PrintDat(sigma_lr_[lr_no]);
    //print the first ten nonzero elems
    //cout << "First ten elems of sigma_lr[" << lr_no << "] are" << endl;
    //print_vector_nonzero_elem(sigma_lr_[lr_no],10);
    //cout << "\n========================================\n" << endl;

}
template
void Square_Double_Layer_PEPSt<ITensor>::snake_walking_bulk_col(int coli, int horizontal_dir, int vertical_dir, bool do_normalize);
template
void Square_Double_Layer_PEPSt<IQTensor>::snake_walking_bulk_col(int coli, int horizontal_dir, int vertical_dir, bool do_normalize);


template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::match_sigma_left_right(int lr_no)
{
    if (col_lr_[0]!=(col_lr_[1]-1))
    {
        cout << "sigma_left and sigma_right are not meeting!" << endl;
        exit(EXIT_FAILURE);
    }

    //decombine sigma_lr_(lr_no) and then recombine it use combiners of (1-lr_no)
    //we can handle up to Ly=6 case
    decombine_sigma_lr(lr_no);
    //Print(sigma_lr_[lr_no]);

    int n_rows=this->lattice().n_uc()[1];
    int start_row;
    int horizontal_dir=1-2*lr_no,
        vertical_dir;
    if (lower_combiners_[1-lr_no][0].numLeft()==1) 
    {
        start_row=0;
        vertical_dir=1;
    }
    if (lower_combiners_[1-lr_no][n_rows-1].numLeft()==1)
    {
        start_row=n_rows-1;
        vertical_dir=-1;
    }

    int rowi=start_row;
    while (rowi>=0 && rowi<n_rows)
    {
        int sitei=this->lattice().site_coord_to_list(col_lr_[lr_no],rowi,0);
        auto double_virt_ind=commonIndex(this->double_layer_tensors_[sitei],this->double_layer_tensors_[sitei+horizontal_dir]);
        sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(this->virt_leg_combiners(sitei,double_virt_ind))*dag(lower_combiners_[1-lr_no][rowi])*dag(upper_combiners_[1-lr_no][rowi]);
        rowi+=vertical_dir;
    }

    //cout << "\n========================================\n" << endl;
    //cout << "Output for sigma_l and sigma_r:" << endl;
    //cout << "left part:" << endl;
    //PrintDat(sigma_lr_[0]);
    //cout << "right part:" << endl;
    //PrintDat(sigma_lr_[1]);
    //cout << "sigma_left and sigma_right matched!" << endl;
    //PrintDat(sigma_lr_[0]*sigma_lr_[1]);
    //cout << "\n========================================\n" << endl;
}
template
void Square_Double_Layer_PEPSt<ITensor>::match_sigma_left_right(int lr_no);
template
void Square_Double_Layer_PEPSt<IQTensor>::match_sigma_left_right(int lr_no);


//TODO: sigma_lr_'s may not have a good quantum number due to boundary
//condition? In this case diagHermitian does not work for IQTensor?
template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::from_sigma_lr_to_sigma_b()
{
    if (col_lr_[0]!=(cutting_col_-1) || col_lr_[1]!=cutting_col_)
    {
        cout << "sigma_lr have not reached cutting column!" << endl;
        exit(EXIT_FAILURE);
    }

    //matches the indices of sigma_lr to each other
    match_sigma_left_right();
    
    //PrintDat(sigma_lr_[0]);
    //get eigenvalues of sigma_l
    std::array<ITensor,2> sigma_lr_itensor={toITensor(sigma_lr_[0]),toITensor(sigma_lr_[1])};
    ITensor U_l, D_l;
    diagHermitian(sigma_lr_itensor[0],U_l,D_l);

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

    //obtain sigma_b_, sigma_b_ shares the same index as sigma_lr_
    int protect_plevel=3;
    sigma_b_=sigma_lr_itensor[1];
    for (int plevel=0; plevel<2; plevel++)
    {
        sqrt_sigma_l.mapprime(plevel,protect_plevel);
        sigma_b_*=sqrt_sigma_l;
        sqrt_sigma_l.mapprime(protect_plevel,plevel);
        sigma_b_.mapprime(protect_plevel,1-plevel);
        //Print(sigma_b_);
    }

    //cout << "\n==============================\n" << endl;
    //PrintDat(sigma_b_);
    //cout << "\n==============================\n" << endl;
}
template
void Square_Double_Layer_PEPSt<ITensor>::from_sigma_lr_to_sigma_b();
template
void Square_Double_Layer_PEPSt<IQTensor>::from_sigma_lr_to_sigma_b();


template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::obtain_density_matrix_spectrum()
{
    //TensorT U_b, D_b;
    //diagHermitian(sigma_b_,U_b,D_b);

    //We use ITensor since sigma_b_ may not have good quantum number
    ITensor U_b, D_b;
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
    //auto max_eigval=*(density_mat_spectrum_.end()-1);
    //for (auto &eigval : density_mat_spectrum_)
    //{
    //    if (eigval/max_eigval<EPSILON)
    //        eigval=0;
    //}

    //cout << "Spectrum of density matrix: " << endl;
    //for (double eigval : density_mat_spectrum_) cout << eigval << " ";
    //cout << endl;
}
template
void Square_Double_Layer_PEPSt<ITensor>::obtain_density_matrix_spectrum();
template
void Square_Double_Layer_PEPSt<IQTensor>::obtain_density_matrix_spectrum();


template <class TensorT>
std::vector<double> Square_Double_Layer_PEPSt<TensorT>::entanglement_spectrum()
{
    std::vector<double> entanglement_spectrum;
    std::transform(density_mat_spectrum_.begin(),density_mat_spectrum_.end(),entanglement_spectrum.begin(),[](double eigval){ return std::exp(-eigval); });
    return entanglement_spectrum;
}
template
std::vector<double> Square_Double_Layer_PEPSt<ITensor>::entanglement_spectrum();
template
std::vector<double> Square_Double_Layer_PEPSt<IQTensor>::entanglement_spectrum();


template <class TensorT>
double Square_Double_Layer_PEPSt<TensorT>::entanglement_entropy_vN()
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
double Square_Double_Layer_PEPSt<ITensor>::entanglement_entropy_vN();
template
double Square_Double_Layer_PEPSt<IQTensor>::entanglement_entropy_vN();


template <class TensorT>
double Square_Double_Layer_PEPSt<TensorT>::entanglement_entropy_Renyi(double renyi_n)
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
double Square_Double_Layer_PEPSt<ITensor>::entanglement_entropy_Renyi(double renyi_n);
template
double Square_Double_Layer_PEPSt<IQTensor>::entanglement_entropy_Renyi(double renyi_n);



template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::obtain_transfer_matrix(int coli)
{
    int sitei=this->lattice().site_coord_to_list(coli,0,0);
    //leg_combiners[0/1] is for left/right legs of transfer matrix
    std::array<std::vector<CombinerT>,2> leg_combiners;
    std::array<IndexT,2> lr_legs;

    lr_legs[0]=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei-1]));
    lr_legs[1]=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei+1]));

    leg_combiners[0].push_back(CombinerT(lr_legs[0]));
    leg_combiners[1].push_back(CombinerT(lr_legs[1]));
    
    transfer_mat_=this->double_layer_tensors_[sitei]*leg_combiners[0][0]*leg_combiners[1][0];

    for (int rowi=1; rowi<this->lattice().n_uc()[1]; rowi++)
    {
        sitei+=this->lattice().n_uc()[0];

        lr_legs[0]=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei-1]));
        lr_legs[1]=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei+1]));

        leg_combiners[0].push_back(CombinerT(lr_legs[0],leg_combiners[0][rowi-1].right()));
        leg_combiners[1].push_back(CombinerT(lr_legs[1],leg_combiners[1][rowi-1].right()));

        transfer_mat_*=this->double_layer_tensors_[sitei];
        transfer_mat_=transfer_mat_*leg_combiners[0][rowi]*leg_combiners[1][rowi];
    }

    //We make transfer_mat_ as matrix --I--mat--I'--
    IndexT oind=leg_combiners[1][this->lattice().n_uc()[1]-1].right(),
           nind=prime(leg_combiners[0][this->lattice().n_uc()[1]-1].right());
    if ((oind.dir()-nind.dir())!=0) nind.dag();
    transfer_mat_.replaceIndex(oind,nind);

    //transfer_mat_/=transfer_mat_.norm();
    clean(transfer_mat_);
    

    //cout << "\n========================================\n" << endl;
    //cout << "transfer matrix for col " << coli << ":" << endl;
    //PrintDat(transfer_mat_);
    //cout << "\n========================================\n" << endl;
}
template
void Square_Double_Layer_PEPSt<ITensor>::obtain_transfer_matrix(int coli);
template
void Square_Double_Layer_PEPSt<IQTensor>::obtain_transfer_matrix(int coli);

template <class TensorT>
std::vector<Complex> Square_Double_Layer_PEPSt<TensorT>::transfer_matrix_eigvals()
{
    //transfer matrix may not be hermitian
    auto transfer_mat_legs=transfer_mat_.indices();
    arma::mat transfer_mat_temp(transfer_mat_legs[0].m(),transfer_mat_legs[1].m());
    for (int val0=1; val0<=transfer_mat_legs[0].m(); val0++)
    {
        for (int val1=1; val1<=transfer_mat_legs[1].m(); val1++)
        {
            transfer_mat_temp(val0-1,val1-1)=transfer_mat_(transfer_mat_legs[0](val0),transfer_mat_legs[1](val1));
        }
    }
    //arma::cx_vec eigvals_transfer_cx=arma::eig_gen(transfer_mat_temp);
    arma::cx_vec eigvals_transfer_temp;
    arma::cx_mat eigvecs_transfer_temp;
    arma::eig_gen(eigvals_transfer_temp,eigvecs_transfer_temp,transfer_mat_temp);
    eigvals_transfer_temp=arma::sort(eigvals_transfer_temp);
    //cout << eigvals_transfer_temp << endl;

    std::vector<Complex> eigvals_transfer=arma::conv_to<std::vector<Complex>>::from(eigvals_transfer_temp);
    //std::sort(eigvals_transfer.begin(),eigvals_transfer.end(),[](Complex eigval1, Complex eigval2){ return std::norm(eigval1)<std::norm(eigval2); });
    return eigvals_transfer;

    //TensorT U_transfer, D_transfer;
    //diagHermitian(transfer_mat_,U_transfer,D_transfer);
    ////PrintDat(D_transfer);

    //auto D_legs=D_transfer.indices();
    //std::vector<double> eigvals_transfer;
    //for (int val=1; val<=D_legs[0].m(); val++)
    //{
    //    eigvals_transfer.push_back(D_transfer(D_legs[0](val),D_legs[1](val)));
    //}

    ////sort eigvals_transfer
    //std::sort(eigvals_transfer.begin(),eigvals_transfer.end());
   //// for (auto &eigval : eigvals_transfer)
   ////     cout << eigval << endl;

    //return eigvals_transfer;
}
template
std::vector<Complex> Square_Double_Layer_PEPSt<ITensor>::transfer_matrix_eigvals();
template
std::vector<Complex> Square_Double_Layer_PEPSt<IQTensor>::transfer_matrix_eigvals();



template <class TensorT>
double Square_Double_Layer_PEPSt<TensorT>::obtain_correlators(const std::vector<int> &acting_sites_list, const std::vector<TensorT> &direct_prod_operators)
{
    //new_col stores the closest cols to leftmost and rightmost operators
    std::array<int,2>  new_col_lr={this->lattice().n_uc()[0]-1,0};
    for (auto sitei : acting_sites_list)
    {
        auto sitei_coord=this->lattice().site_list_to_coord(sitei);
        if ((sitei_coord[0]-1)<new_col_lr[0]) new_col_lr[0]=sitei_coord[0]-1;
        if ((sitei_coord[0]+1)>new_col_lr[1]) new_col_lr[1]=sitei_coord[0]+1;
    }
    //cout << "new cols: " << new_col_lr[0] << " " << new_col_lr[1] << endl;
    //move sigma_lr_ to new_col_lr to reduce contraction time
    move_sigma_lr(new_col_lr);

    //calculate wavefunction norm
    auto sandwiched_tensors=this->double_layer_tensors_;
    double wf_norm=sandwiched_peps_norm(sandwiched_tensors);

    //cout << "wavefunction norm: " << wf_norm << endl;
    
    //calculate the unnormalized expect_val
    this->obtain_peps_sandwich_single_site_operators(direct_prod_operators,acting_sites_list,sandwiched_tensors);

    //for (int sitei=0; sitei<this->lattice().n_sites_total(); sitei++)
    //{
    //    auto diff_tensor=sandwiched_tensors[sitei]-this->double_layer_tensors_[sitei];
    //    if (diff_tensor.norm()>EPSILON)
    //    {
    //        cout << "Tensor at site " << sitei << " has changed!" << endl;
    //        PrintDat(this->double_layer_tensors_[sitei]);
    //        PrintDat(sandwiched_tensors[sitei]);
    //    }
    //}

    double expect_val=this->sandwiched_peps_norm(sandwiched_tensors);

    //cout << "unnormalized correlator: " << expect_val << endl;

    return expect_val/wf_norm;
}
template
double Square_Double_Layer_PEPSt<ITensor>::obtain_correlators(const std::vector<int> &acting_sites_list, const std::vector<ITensor> &direct_prod_operators);
template
double Square_Double_Layer_PEPSt<IQTensor>::obtain_correlators(const std::vector<int> &acting_sites_list, const std::vector<IQTensor> &direct_prod_operators);

template <class TensorT>
double Square_Double_Layer_PEPSt<TensorT>::obtain_correlators(const std::vector<std::vector<int>> &acting_sites_list, const std::vector<TPOt<TensorT>> &tensor_prod_operators)
{
    //new_col stores the closest cols to leftmost and rightmost operators
    std::array<int,2> new_col_lr={this->lattice().n_uc()[0]-1,0};
    for (const auto &acting_sites : acting_sites_list)
    {
        for (auto sitei : acting_sites)
        {
            auto sitei_coord=this->lattice().site_list_to_coord(sitei);
            if ((sitei_coord[0]-1)<new_col_lr[0]) new_col_lr[0]=sitei_coord[0]-1;
            if ((sitei_coord[0]+1)>new_col_lr[1]) new_col_lr[1]=sitei_coord[0]+1;
        }
    }
    //cout << "new cols: " << new_col_lr[0] << " " << new_col_lr[1] << endl;
    //move sigma_lr_ to new_col_lr to reduce contraction time
    move_sigma_lr(new_col_lr);

    //calculate wavefunction norm
    auto sandwiched_tensors=this->double_layer_tensors_;
    double wf_norm=sandwiched_peps_norm(sandwiched_tensors);

    //cout << "wavefunction norm: " << wf_norm << endl;
    
    //calculate the unnormalized expect_val
    this->obtain_peps_sandwich_tensor_prod_operators(tensor_prod_operators,acting_sites_list,sandwiched_tensors);

    double expect_val=this->sandwiched_peps_norm(sandwiched_tensors);

    //cout << "unnormalized correlator: " << expect_val << endl;

    return expect_val/wf_norm;
}
template
double Square_Double_Layer_PEPSt<ITensor>::obtain_correlators(const std::vector<std::vector<int>> &acting_sites_list, const std::vector<TPOt<ITensor>> &tensor_prod_operators);
template
double Square_Double_Layer_PEPSt<IQTensor>::obtain_correlators(const std::vector<std::vector<int>> &acting_sites_list, const std::vector<TPOt<IQTensor>> &tensor_prod_operators);



template <class TensorT>
double Square_Double_Layer_PEPSt<TensorT>::sandwiched_peps_norm(const std::vector<TensorT> &sandwiched_tensors, int horizontal_dir)
{
    //store the original data
    auto origin_col_lr=col_lr_;
    auto origin_sigma_lr=sigma_lr_;
    auto origin_lower_combiners=lower_combiners_;
    auto origin_upper_combiners=upper_combiners_;
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
    lower_combiners_=origin_lower_combiners;
    upper_combiners_=origin_upper_combiners;
    this->double_layer_tensors_=origin_double_layer_tensors;

    return expect_value;
}
template
double Square_Double_Layer_PEPSt<ITensor>::sandwiched_peps_norm(const std::vector<ITensor> &sandwiched_tensors, int horizontal_dir);
template
double Square_Double_Layer_PEPSt<IQTensor>::sandwiched_peps_norm(const std::vector<IQTensor> &sandwiched_tensors, int horizontal_dir);


//TODO:during the decombing/recombing process, we can reach at most 6 rows since NMAX=8. We can modify algrithm to handle 8 rows
template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::decombine_sigma_lr(int lr_no)
{
    int n_rows=this->lattice().n_uc()[1],
        n_cols=this->lattice().n_uc()[0];

    int horizontal_dir=1-2*lr_no;
    int vertical_dir;
    int start_row;

    //decombine from top to down
    if (lower_combiners_[lr_no][0].numLeft()==1) 
    {
        vertical_dir=-1;
        start_row=n_rows-1;
    }
    //decombine from down to top
    if (lower_combiners_[lr_no][n_rows-1].numLeft()==1) 
    {
        vertical_dir=1;
        start_row=0;
    }

    int rowi=start_row;
    while (rowi>=0 && rowi<n_rows)
    {
        int sitei=this->lattice().site_coord_to_list(col_lr_[lr_no],rowi,0);
        sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(lower_combiners_[lr_no][rowi])*dag(upper_combiners_[lr_no][rowi]);
        //combine a pair of upper and lower indice
        auto double_virt_ind=commonIndex(this->double_layer_tensors(sitei),dag(this->double_layer_tensors(sitei+horizontal_dir)));
        sigma_lr_[lr_no]=sigma_lr_[lr_no]*this->virt_leg_combiners(sitei,double_virt_ind);
        rowi+=vertical_dir;
    }

    //if (lower_combiners_[lr_no][0].numLeft()==1) //decombine from top to down
    //{
    //    int sitei=this->lattice().site_coord_to_list(col_lr_[lr_no],n_rows-1,0);
    //    for (int rowi=n_rows-1; rowi>=0; rowi--)
    //    {
    //        sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(lower_combiners_[lr_no][rowi])*dag(upper_combiners_[lr_no][rowi]);

    //        //combine pair of upper and lower indice
    //        auto double_virt_ind=commonIndex(this->double_layer_tensors(sitei),dag(this->double_layer_tensors(sitei+horizontal_dir)));
    //        sigma_lr_[lr_no]=sigma_lr_[lr_no]*this->virt_leg_combiners(sitei,double_virt_ind);

    //        sitei-=n_cols;
    //    }
    //}

    //if (lower_combiners_[lr_no][n_rows-1].numLeft()==1)//decombine from down to top
    //{
    //    int sitei=this->lattice().site_coord_to_list(col_lr_[lr_no],0,0);
    //    for (int rowi=0; rowi<n_rows; rowi++)
    //    {
    //        sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(lower_combiners_[lr_no][rowi])*dag(upper_combiners_[lr_no][rowi]);

    //        //combine pair of upper and lower indice
    //        auto double_virt_ind=commonIndex(this->double_layer_tensors(sitei),dag(this->double_layer_tensors(sitei+horizontal_dir)));
    //        sigma_lr_[lr_no]=sigma_lr_[lr_no]*this->virt_leg_combiners(sitei,double_virt_ind);

    //        sitei+=n_cols;
    //    }
    //}
}
template
void Square_Double_Layer_PEPSt<ITensor>::decombine_sigma_lr(int lr_no);
template
void Square_Double_Layer_PEPSt<IQTensor>::decombine_sigma_lr(int lr_no);


template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::decombine_sigma_lr()
{
    decombine_sigma_lr(0);
    decombine_sigma_lr(1);
}
template
void Square_Double_Layer_PEPSt<ITensor>::decombine_sigma_lr();
template
void Square_Double_Layer_PEPSt<IQTensor>::decombine_sigma_lr();

template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::recombine_sigma_lr()
{
    int n_cols=this->lattice().n_uc()[0],
        n_rows=this->lattice().n_uc()[1];
    for (int lr_no=0; lr_no<2; lr_no++)
    {
        int horizontal_dir=1-2*lr_no,
            vertical_dir=1-2*std::abs((lr_no*n_cols-col_lr_[lr_no])%2);
        int start_row=(1-vertical_dir)/2*(n_rows-1);
        int rowi=start_row;
        while (rowi>=0 && rowi<n_rows)
        {
            int sitei=this->lattice().site_coord_to_list(col_lr_[lr_no],rowi,0);

            auto double_virt_ind=commonIndex(this->double_layer_tensors_[sitei],dag(this->double_layer_tensors_[sitei+horizontal_dir]));
            IndexT lower_ind;
            auto double_virt_leg_combiner=this->decombine_double_virt_indice(sitei,double_virt_ind,lower_ind);
            if (rowi==start_row)
            {
                lower_combiners_[lr_no][rowi]=CombinerT(lower_ind);
                upper_combiners_[lr_no][rowi]=prime(dag(lower_combiners_[lr_no][rowi]));
            }
            else
            {
                lower_combiners_[lr_no][rowi]=CombinerT(lower_combiners_[lr_no][rowi-vertical_dir].right(),lower_ind);
                upper_combiners_[lr_no][rowi]=prime(dag(lower_combiners_[lr_no][rowi]));
            }
            sigma_lr_[lr_no]=sigma_lr_[lr_no]*dag(double_virt_leg_combiner)*lower_combiners_[lr_no][rowi]*upper_combiners_[lr_no][rowi];

            rowi+=vertical_dir;
        }

        clean(sigma_lr_[lr_no]);
    }

}
template
void Square_Double_Layer_PEPSt<ITensor>::recombine_sigma_lr();
template
void Square_Double_Layer_PEPSt<IQTensor>::recombine_sigma_lr();

template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::match_sigma_lr_double_peps()
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
void Square_Double_Layer_PEPSt<ITensor>::match_sigma_lr_double_peps();
template
void Square_Double_Layer_PEPSt<IQTensor>::match_sigma_lr_double_peps();

template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::move_sigma_lr(const std::array<int,2> &new_col_lr)
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
    recombine_sigma_lr();
}
template 
void Square_Double_Layer_PEPSt<ITensor>::move_sigma_lr(const std::array<int,2> &new_col_lr);
template
void Square_Double_Layer_PEPSt<IQTensor>::move_sigma_lr(const std::array<int,2> &new_col_lr);


template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::read(std::istream &s)
{
    Double_Layer_PEPSt<TensorT>::read(s);
    //for (const auto &tensor : this->double_layer_tensors()) Print(tensor);

    s.read((char*)&cutting_col_,sizeof(cutting_col_));
    for (auto &col : col_lr_) s.read((char*)&col,sizeof(col));
    //cout << "col_lr=" << col_lr_[0] << " " << col_lr_[1] << endl;

    //what we read is decombined sigma_lr_ whose indices do not match the new double layer peps
    for (auto &sigma : sigma_lr_) sigma.read(s);
    //Print(sigma_lr_[0]);
    //Print(sigma_lr_[1]);
    //replace indices of sigma_lr_ such that they share the indices of double_layer_tensors_. Then recombine sigma_lr_ to a vector
    match_sigma_lr_double_peps();
    recombine_sigma_lr();

    sigma_b_.read(s);

    for (auto eigval : density_mat_spectrum_) s.read((char*)&eigval,sizeof(eigval));
}
template 
void Square_Double_Layer_PEPSt<ITensor>::read(std::istream &s);
template 
void Square_Double_Layer_PEPSt<IQTensor>::read(std::istream &s);


template <class TensorT>
void Square_Double_Layer_PEPSt<TensorT>::write(std::ostream &s) const
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
void Square_Double_Layer_PEPSt<ITensor>::write(std::ostream &s) const;
template
void Square_Double_Layer_PEPSt<IQTensor>::write(std::ostream &s) const;



//double honeycomb_cylinder_multiple_Sz_correlator(const std::array<Coordinate,2> &acting_sites_coord, Square_Double_Layer_PEPSt<IQTensor> &square_cylinder_double_peps)
//{
//    //we create Sz operator for honeycomb sites, and then transfer it to square lattice
//    std::vector<IQTensor> square_Sz_operators;
//    std::vector<int> square_acting_sites;
//    //col_lr stores the closest cols to leftmost and rightmost cols
//    std::array<int,2> col_lr={square_cylinder_double_peps.lattice().n_uc()[0]-1,0};
//
//    for (auto site_coord : acting_sites_coord)
//    {
//        auto site_no=square_cylinder_double_peps.lattice().site_coord_to_list(site_coord[0],site_coord[1],0);
//        square_acting_sites.push_back(site_no);
//
//        //update col_lr
//        if ((site_coord[0]-1)<col_lr[0]) col_lr[0]=site_coord[0]-1;
//        if ((site_coord[0]+1)>col_lr[1]) col_lr[1]=site_coord[0]+1;
//
//        //create Sz operator and Id operator for one uc of honeycomb
//        std::vector<IQIndex> Sz_legs={Spin_leg(std::vector<int>{0,1},"Sz_leg 0",Out,Site), Spin_leg(std::vector<int>{0,1},"Sz_leg 1",Out,Site)};
//        IQTensor Sz_operator(dag(Sz_legs[0])(1),prime(Sz_legs[0])(1)),
//                 Id_operator(dag(Sz_legs[1])(1),prime(Sz_legs[1])(1));
//
//        Sz_operator(dag(Sz_legs[0])(2),prime(Sz_legs[0])(2))=-1;
//        Id_operator(dag(Sz_legs[1])(2),prime(Sz_legs[1])(2))=1;
//        
//        //transfer to square Sz operator in one site
//        auto combined_SzId_operator=Sz_operator*Id_operator;
//        IQCombiner Sz_legs_combiner;
//
//        if (site_coord[2]==0) 
//        {
//            Sz_legs_combiner.addleft(Sz_legs[0]);
//            Sz_legs_combiner.addleft(Sz_legs[1]);
//        }
//        if (site_coord[2]==1) 
//        {
//            Sz_legs_combiner.addleft(Sz_legs[1]);
//            Sz_legs_combiner.addleft(Sz_legs[0]);
//        }
//
//        //TODO: check the order of combined operator
//        Sz_legs_combiner.init("square_Sz_leg",Site,Out);
//        combined_SzId_operator=combined_SzId_operator*dag(Sz_legs_combiner)*prime(Sz_legs_combiner);
//        square_Sz_operators.push_back(combined_SzId_operator);
//    }
//
//    //move sigma_lr_ to col_lr_ to reduce contraction time, and then do contraction
//    square_cylinder_double_peps.move_sigma_lr(col_lr);
//
//    //obtain the sandwiched double peps, calculate wavefunction norm
//    auto square_sandwiched_tensors=square_cylinder_double_peps.double_layer_tensors();
//    double wf_norm=square_cylinder_double_peps.sandwiched_peps_norm(square_sandwiched_tensors);
//
//    square_cylinder_double_peps.obtain_peps_sandwich_single_site_operators(square_Sz_operators,square_acting_sites,square_sandwiched_tensors);
//
//    double unnormlized_expect_val=square_cylinder_double_peps.sandwiched_peps_norm(square_sandwiched_tensors);
//
//    return unnormlized_expect_val/wf_norm;
//}

