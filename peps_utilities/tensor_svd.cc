
#include "tensor_svd.h"
#include "svdalgs.h"
#include "lapack_wrap.h"

bool vmc_isZero(const ITensor&T)
    {
    if(T.scale().isZero())
    {
      return true;
    }
    if(std::abs(T.normNoScale())<1.E-17) {
      return true;
    }
    return false;
    }


bool vmc_isZero(IQTensor& T)
    {
      
    if(T.empty()) return true;
    for(const ITensor& t : T.blocks())
        {
        if(!vmc_isZero(t)) return false;
        }
    return true;
    }



struct EigQN
    {
    Real eig;
    QN qn;

    EigQN(Real eg,QN q)
        : eig(eg),qn(q) 
        { }

    bool
    operator<(const EigQN& other) const 
        { 
        return eig < other.eig;
        }
    };


Real static
truncate(Vector& D,
         int maxm,
         int minm,
         Real cutoff,
         bool absoluteCutoff,
         bool doRelCutoff)
    {
    int m = D.Length();
    if(m == 1) return 0;

    Real truncerr = 0;

    //Zero out any negative weight
    for(int zerom = m; zerom > 0; --zerom)
        {
        if(D(zerom) >= 0) break;
        D(zerom) = 0;
        }

    if(absoluteCutoff)
        {
        for(;m > maxm || (D(m) < cutoff && m > minm); --m)
            {
            truncerr += D(m);
            }
        }
    else
        {
        const Real scale = doRelCutoff ? D(1) : 1.0;
        for(;m > maxm || (truncerr+D(m) < cutoff*scale && m > minm); --m)
            {
            truncerr += D(m);
            }
        truncerr = (D(1) == 0 ? 0 : truncerr/scale);
        }

    D.ReduceDimension(m); 

    return truncerr;
    }

Real static
truncate(std::vector<EigQN>& alleig, 
         int& m, 
         Real& docut, 
         int maxm,
         int minm,
         Real cutoff,
         bool absoluteCutoff,
         bool doRelCutoff)
    {
    m = (int)alleig.size();
    if(m == 1)
        {
        docut = alleig.front().eig/2.;
        return 0;
        }
    long mdisc = 0;

    Real truncerr = 0;

    if(absoluteCutoff)
        {
        while(m > maxm || ( (alleig.at(mdisc).eig < cutoff && m > minm)
            && mdisc < (long)alleig.size() ) )
            {
            if(alleig.at(mdisc).eig > 0)
                truncerr += alleig.at(mdisc).eig;
            else
                alleig.at(mdisc).eig = 0;

            ++mdisc;
            --m;
            }
        docut = (mdisc > 0 
                ? (alleig.at(mdisc-1).eig + alleig.at(mdisc).eig)*0.5 - 1E-5*alleig.at(mdisc-1).eig
                : -1);
        }
    else
        {
        Real scale = doRelCutoff ? alleig.back().eig : 1.0;
        while(   m > maxm 
             || ( (mdisc < (int)alleig.size()) && (truncerr+alleig.at(mdisc).eig < cutoff*scale && m > minm))
             )
            {
            if(alleig.at(mdisc).eig > 0)
                truncerr += alleig.at(mdisc).eig;
            else
                alleig.at(mdisc).eig = 0;

            ++mdisc;
            --m;
            }
        if(mdisc >= int(alleig.size())) mdisc = alleig.size() - 1;
        docut = (mdisc > 0 
                ? (alleig.at(mdisc-1).eig + alleig.at(mdisc).eig)*0.5 - 1E-5*alleig.at(mdisc-1).eig
                : -1);
        truncerr = (alleig.back().eig == 0 ? 0 : truncerr/scale);
        }


    return truncerr;
    }



template<class Tensor>
bool tensor_svd(Tensor AA, Tensor& U, Tensor& D, Tensor& V, const Args& args)
    {
    using IndexT = typename Tensor::IndexT;
    using CombinerT = typename Tensor::CombinerT;

    const Real noise = args.getReal("Noise",0.);
    const bool useOrigM = args.getBool("UseOrigM",false);
    const Args* args_ = &args;
    if(AA.valid()==false){  return false;}
    if(vmc_isZero(AA)) {
      return false;
    }

    if(noise > 0)
        Error("Noise term not implemented for svd");

    //Combiners which transform AA
    //into a rank 2 tensor
    CombinerT Ucomb, Vcomb;

    //Divide up indices based on U
    //If U is null, use V instead
    const Tensor &L = (U ? U : V);
    CombinerT &Lcomb = (U ? Ucomb : Vcomb),
              &Rcomb = (U ? Vcomb : Ucomb);
    for(const IndexT& I : AA.indices())
        { 
        if(hasindex(L,I))
            Lcomb.addleft(I);
        else
            Rcomb.addleft(I);
        }

    AA = Ucomb * AA * Vcomb;

    Args newArgs(args);
    if(useOrigM)
        {
        //Try to determine current m,
        //then set minm_ and maxm_ to this.
        newArgs.add("Cutoff",-1);
        int minm = 1,
            maxm = MAX_M;
        if(D.r() == 0)
            {
            IndexT mid = commonIndex(U,V,Link);
            if(mid) minm = maxm = mid.m();
            else    minm = maxm = 1;
            }
        else
            {
            minm = maxm = D.indices().front().m();
            }
        newArgs.add("Minm",minm);
        newArgs.add("Maxm",maxm);
        args_ = &newArgs;
        }

    bool svd_ok_q = 
    tensor_svdRank2(AA,Ucomb.right(),Vcomb.right(),U,D,V,*args_);
    if(svd_ok_q==false){
      return false;
    }

    U = dag(Ucomb) * U;
    V = V * dag(Vcomb);

    return true;

    }
template bool tensor_svd(ITensor AA, ITensor& U, ITensor& D, ITensor& V, const Args& args);
template bool tensor_svd(IQTensor AA, IQTensor& U, IQTensor& D, IQTensor& V, const Args& args);


bool 
tensor_svdRank2(ITensor A, const Index& ui, const Index& vi,
         ITensor& U, ITensor& D, ITensor& V,
         const Args& args)
    {
    const Real thresh = args.getReal("SVDThreshold",1E-4);
    const Real cutoff = args.getReal("Cutoff",MIN_CUT);
    const int maxm = args.getInt("Maxm",MAX_M);
    const int minm = args.getInt("Minm",1);
    const bool do_truncate = args.getBool("Truncate",true);
    const bool doRelCutoff = args.getBool("DoRelCutoff",false);
    const bool absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    const bool cplx = A.isComplex();

    auto lname = args.getString("LeftIndexName","ul");
    auto rname = args.getString("RightIndexName","vl");
    auto itype = getIndexType(args,"IndexType",Link);
    auto litype = getIndexType(args,"LeftIndexType",itype);
    auto ritype = getIndexType(args,"RightIndexType",itype);

    if(A.r() != 2)
        {
        Error("A must be matrix-like");
        }

    Matrix UU,VV,
           iUU,iVV;
    Vector DD;

    if(!cplx)
        {
        Matrix M;
        A.toMatrix11NoScale(ui,vi,M);

        //SVD(M,UU,DD,VV,thresh);
        SVD_bf(M,UU,DD,VV);
        }
    else
        {
        ITensor Are = realPart(A),
                Aim = imagPart(A);
        Are.scaleTo(A.scale());
        Aim.scaleTo(A.scale());
        Matrix Mre,Mim;
        Are.toMatrix11NoScale(ui,vi,Mre);
        Aim.toMatrix11NoScale(ui,vi,Mim);

        //SVDComplex(Mre,Mim,UU,iUU,DD,VV,iVV);
        //SVD(Mre,Mim,UU,iUU,DD,VV,iVV,thresh);
        SVD_bf(Mre,Mim,UU,iUU,DD,VV,iVV);
        }

    //Truncate

    Spectrum spec;

    int m = DD.Length();
    Real terr = 0;
    if(do_truncate)
        {
        Vector sqrD(DD);
        for(int j = 1; j <= sqrD.Length(); ++j)
            sqrD(j) = sqr(DD(j));
        terr = truncate(sqrD,maxm,minm,cutoff,absoluteCutoff,doRelCutoff);
        m = sqrD.Length();
        DD.ReduceDimension(m);
        }


    if(args.getBool("ShowEigs",false))
        {
        cout << endl;
        printfln("minm = %d, maxm = %d, cutoff = %.3E",minm,maxm,cutoff);
        printfln("truncate = %s",(do_truncate?"true":"false"));
        printfln("doRelCutoff = %s",(doRelCutoff?"true":"false"));
        printfln("absoluteCutoff = %s",(absoluteCutoff?"true":"false"));
        printfln("Kept m=%d states in svdRank2 line 169", m);
        printfln("svdtruncerr = %.3E",spec.truncerr());

        int stop = min(10,DD.Length());
        Vector Ds = DD.SubVector(1,stop);

        Real orderMag = log(fabs(DD(1))) + A.scale().logNum();
        if(fabs(orderMag) < 5 && A.scale().isFiniteReal())
            {
            Ds *= fabs(A.scale().real0());
            cout << "Singular values: ";
            }
        else
            {
            cout << "Singular values (not including scale = " << A.scale() << "):";
            }

        for(int j = 1; j <= stop; ++j)
            {
            const Real sval = Ds(j);
            printf(( sval > 1E-3 && sval < 1000) ? ("%.3f") : ("%.3E") , sval); 
            print((j != stop) ? ", " : "\n");
            }
        cout << endl;
        }
    
    Index uL(lname,m,Link),vL(rname,m,Link);

    D = ITensor(uL,vL,DD);
    D *= A.scale();
    U = ITensor(ui,uL,UU.Columns(1,m));
    V = ITensor(vL,vi,VV.Rows(1,m));

    //Fix for cases where A.scale() may be negative
    if(A.scale().sign() == -1)
        {
        D *= -1;
        U *= -1;
        }

    if(cplx)
        {
        ITensor iU(ui,uL,iUU.Columns(1,m)),
                iV(vL,vi,iVV.Rows(1,m));
        if(iU.norm() > 1E-14)
            U = U + iU*Complex_i;
        if(iV.norm() > 1E-14)
            V = V + iV*Complex_i;
        }

    //Square all singular values
    //since convention is to report
    //density matrix eigs
    for(int j = 1; j <= m; ++j)
        {
        DD(j) *= DD(j);
        }

    if(A.scale().isFiniteReal())
        {
        DD *= sqr(A.scale().real0());
        }
    else
        {
        cout << "Warning: scale not finite real" << endl;
        }

    //Global::lastd() = DD;

    //Include A's scale to get the actual eigenvalues kept
    //as long as the leading eigenvalue is within a few orders
    //of magnitude of 1.0. Otherwise just report the scaled eigs.
    //Real orderMag = log(fabs(DD(1))) + A.scale().logNum();
    //if(fabs(orderMag) < 5 && A.scale().isFiniteReal())
    //    {
    //    Global::lastd() *= A.scale().real();
    //    }

    return true;

    } // void svdRank2


bool
tensor_svdRank2(IQTensor A, const IQIndex& uI, const IQIndex& vI,
         IQTensor& U, IQTensor& D, IQTensor& V,
         const Args& args)
    {
    auto cplx = A.isComplex();
    auto thresh = args.getReal("SVDThreshold",1E-4);
    auto cutoff = args.getReal("Cutoff",MIN_CUT);
    auto maxm = args.getInt("Maxm",MAX_M);
    auto minm = args.getInt("Minm",1);
    auto do_truncate = args.getBool("Truncate",true);
    auto doRelCutoff = args.getBool("DoRelCutoff",false);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    auto logrefNorm = args.getReal("LogRefNorm",0.);

    auto lname = args.getString("LeftIndexName","ul");
    auto rname = args.getString("RightIndexName","vl");
    auto itype = getIndexType(args,"IndexType",Link);
    auto litype = getIndexType(args,"LeftIndexType",itype);
    auto ritype = getIndexType(args,"RightIndexType",itype);

    if(A.r() != 2)
        {
        Error("A must be matrix-like");
        }

    const int Nblock = A.blocks().size();
    if(Nblock == 0){
        return false;
    }

    std::vector<Matrix> Umatrix(Nblock),
                   Vmatrix(Nblock),
                   iUmatrix,
                   iVmatrix;
    if(cplx)
        {
        iUmatrix.resize(Nblock);
        iVmatrix.resize(Nblock);
        }

    std::vector<Vector> dvector(Nblock);

    std::vector<EigQN> alleig;
    alleig.reserve(min(uI.m(),vI.m()));

    if(uI.m() == 0){
        return false;
    }
    if(vI.m() == 0){
        return false;
    }

    LogNumber refNorm(logrefNorm,1);
    if(doRelCutoff)
        {
        Real maxLogNum = -200;
        A.scaleOutNorm();
        for(const ITensor& t : A.blocks())
            {
            maxLogNum = max(maxLogNum,t.scale().logNum());
            }
        refNorm = LogNumber(maxLogNum,1);
        }
    A.scaleTo(refNorm);

    //1. SVD each ITensor within A.
    //   Store results in mmatrix and mvector.
    int itenind = 0;
    for(const ITensor& t : A.blocks())
        {
        Matrix &UU = Umatrix.at(itenind);
        Matrix &VV = Vmatrix.at(itenind);
        Vector &d =  dvector.at(itenind);

        const Index *ui=0,*vi=0;
        bool gotui = false;
        for(const Index& I : t.indices())
            {
            if(!gotui) 
                {
                ui = &I;
                gotui = true;
                }
            else       
                {
                vi = &I;
                break;
                }
            }

        if(!hasindex(uI,*ui))
            std::swap(ui,vi);

        if(!cplx)
            {
            Matrix M(ui->m(),vi->m());
	    if(t.r()!=2){return false;}
            t.toMatrix11NoScale(*ui,*vi,M);

            //SVD(M,UU,d,VV,thresh);
	    SVD_bf(M,UU,d,VV);
            }
        else
            {
            ITensor ret = realPart(t),
                    imt = imagPart(t);
            ret.scaleTo(refNorm);
            imt.scaleTo(refNorm);
            Matrix Mre(ui->m(),vi->m()),
                   Mim(ui->m(),vi->m());
            ret.toMatrix11NoScale(*ui,*vi,Mre);
            imt.toMatrix11NoScale(*ui,*vi,Mim);

            //SVDComplex(Mre,Mim,
            //           UU,iUmatrix.at(itenind),
            //           d,
            //           VV,iVmatrix.at(itenind));
	
            //SVD(Mre,Mim, UU,iUmatrix.at(itenind),d,VV,iVmatrix.at(itenind),thresh);
	    SVD_bf(Mre,Mim,UU,iUmatrix.at(itenind),d,VV,iVmatrix.at(itenind));
            }

        //Store the squared singular values
        //(denmat eigenvalues) in alleig
        QN q = qn(uI,*ui);
        for(int j = 1; j <= d.Length(); ++j) 
            {
            alleig.push_back(EigQN(sqr(d(j)),q));
            }

        ++itenind;
        }

    //2. Truncate eigenvalues

    //Determine number of states to keep m
    int m = (int)alleig.size();
    Real svdtruncerr = 0;
    Real docut = -1;

    //Sort all eigenvalues from smallest to largest
    //irrespective of quantum numbers
    sort(alleig.begin(),alleig.end());

    if(do_truncate)
        {
        svdtruncerr = truncate(alleig,m,docut,maxm,minm,cutoff,
                               absoluteCutoff,doRelCutoff);
        }

    if(args.getBool("ShowEigs",false))
        {
        cout << endl;
        println("svdRank2 (IQTensor):");
        printfln("    minm = %d, maxm = %d, cutoff = %.3E",minm,maxm,cutoff);
        printfln("    Kept m = %d states in svdRank2",m);
        printfln("    svdtruncerr = %.2E",svdtruncerr);
        printfln("    docut = %.2E",docut);
        println("    doRelCutoff is ",(doRelCutoff ? "true" : "false"));
        println("    absoluteCutoff is ",(absoluteCutoff ? "true" : "false"));
        println("    refNorm is ",refNorm);

        const int s = alleig.size();
        const int max_show = 20;
        int stop = s-min(s,max_show);

        //Include refNorm in printed eigs as long as
        //the leading eig is within a few orders of magnitude
        //of 1.0. Otherwise just print the scaled eigs.
        Real orderMag = log(fabs(alleig.at(s-1).eig)) + refNorm.logNum();
        Real real_fac = 1;
        if(fabs(orderMag) < 5 && refNorm.isFiniteReal())
            {
            real_fac = refNorm.real();
            cout << "    Singular values: ";
            }
        else
            {
            cout << "    Singular values [omitting scale factor " << refNorm << "]: \n";
            if(alleig.at(s-1).eig > 1.e10)
                {
                Error("bad alleig");
                }
            cout << "    ";
            }

        for(int j = s-1; j >= stop; --j)
            {
            const Real sval = std::sqrt(alleig.at(j).eig)*real_fac;
            printf( (sval >= 1E-3 && sval < 1E3) ? ("%.3f") : ("%.3E"), sval);
            print((j != stop) ? ", " : "\n");
            }
        cout << endl;
        } //end if(showeigs_)

    //Truncate denmat eigenvalue vectors
    //Also form new Link index with appropriate m's for each block
    IQIndex::Storage Liq, Riq;
    Liq.reserve(Nblock);
    Riq.reserve(Nblock);

    std::vector<ITensor> Ublock,
                    Vblock,
                    iUblock,
                    iVblock;
    Ublock.reserve(Nblock);
    Vblock.reserve(Nblock);
    if(cplx)
        {
        iUblock.reserve(Nblock);
        iVblock.reserve(Nblock);
        }

    std::vector<ITensor> Dblock;
    Dblock.reserve(Nblock);

    itenind = 0;
    int total_m = 0;
    for(const ITensor& t : A.blocks())
        {
        const Matrix& UU = Umatrix.at(itenind);
        const Matrix& VV = Vmatrix.at(itenind);
        Vector& thisD = dvector.at(itenind);

        int this_m = 1;
        while(this_m <= thisD.Length() && sqr(thisD(this_m)) > docut) 
            {
            ++total_m;
            if(thisD(this_m) < 0) thisD(this_m) = 0;
            ++this_m;
            }
        --this_m; //since the loop overshoots by 1

        if(m == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
            { this_m = 1; m = 1; docut = 1; }

        if(this_m == 0) { ++itenind; continue; }

        const Index *ui=0,*vi=0;
        bool gotui = false;
        for(const Index& I : t.indices())
            {
            if(!gotui) 
                {
                ui = &I;
                gotui = true;
                }
            else       
                {
                vi = &I;
                break;
                }
            }

        if(!hasindex(uI,*ui))
            std::swap(ui,vi);

        Index l("l",this_m);
        Liq.push_back(IndexQN(l,qn(uI,*ui)));

        Index r("r",this_m);
        Riq.push_back(IndexQN(r,qn(vI,*vi)));

        Dblock.push_back(ITensor(l,r,thisD.SubVector(1,this_m)));

        Ublock.push_back(ITensor(*ui,l,UU.Columns(1,this_m)));
        Vblock.push_back(ITensor(r,*vi,VV.Rows(1,this_m)));

        if(cplx)
            {
            iUblock.push_back(ITensor(*ui,l,iUmatrix.at(itenind).Columns(1,this_m)));
            iVblock.push_back(ITensor(r,*vi,iVmatrix.at(itenind).Rows(1,this_m)));
            }

        ++itenind;
        }

    if(Liq.size() == 0){
      return false;
    }

    IQIndex L(lname,Liq,uI.dir()), R(rname,Riq,vI.dir());

    D = IQTensor(L,R);
    U = IQTensor(uI,dag(L));
    V = IQTensor(dag(R),vI);

    //Load blocks into D,U, and V
    for(size_t j = 0; j < Dblock.size(); ++j)
        {
        D += Dblock.at(j);
        U += Ublock.at(j);
        V += Vblock.at(j);
        }

    if(cplx)
        {
        IQTensor iU(uI,dag(L));
        IQTensor iV(dag(R),vI);
        for(size_t j = 0; j < Dblock.size(); ++j)
            {
            if(iUblock.at(j).norm() > 1E-14)
                {
                iU += iUblock.at(j);
                }
            if(iVblock.at(j).norm() > 1E-14)
                {
                iV += iVblock.at(j);
                }
            }
        if(!iU.blocks().empty())
            {
            U = U + iU*Complex_i;
            }
        if(!iV.blocks().empty())
            {
            V = V + iV*Complex_i;
            }
        }

    //Originally eigs were found by calling
    //toMatrix11NoScale, so put the scale back in
    D *= refNorm;

    int aesize = int(alleig.size());
    int neig = std::min(aesize,L.m());
    Vector DD(neig);
    std::vector<QN> qns(neig);
    for(int i = 0; i < neig; ++i) 
        {
        DD[i] = alleig.at(aesize-1-i).eig;
        qns[i] = alleig.at(aesize-1-i).qn;
        }

    return true;

    } //void svdRank2


void 
SVD_bf(const MatrixRef& A, Matrix& U, Vector& d, Matrix& V){
      LAPACK_INT m = A.Nrows(), 
               n = A.Ncols(); 

    if(m < n)
        {
        Matrix mret = A.t(), 
               UUre,
               VVre;
        SVD_bf(mret, UUre, d, VVre);
        V = UUre.t();
        U = VVre.t();
        return;
        }

    char jobz = 'S';

    Matrix AA(n,m);
    for(int i = 1; i <= n; i++)
	for(int j = 1; j <= m; j++)
        {
	    AA(i,j) = A(j,i); 
        }

    Matrix UU(n,m), VV(n,n);
    d.ReDimension(n);
    LAPACK_INT info = 0;

    dgesdd_wrapper(&jobz,&m,&n,
                   (LAPACK_REAL*)AA.Store(),
                   d.Store(), 
                   (LAPACK_REAL*)UU.Store(),
                   (LAPACK_REAL*)VV.Store(),
                   &info);

    if(info != 0) 
        {
        cout <<"dgesdd failed, "<< "info = " << info<<". Now try dgesvd." << endl;
        char jobu = 'S', jobvt='S';
        for(int i = 1; i <= n; i++)
	for(int j = 1; j <= m; j++)
        {
	    AA(i,j) = A(j,i); 
        }
        info=0;
        dgesvd_wrapper(&jobu,&jobvt,&m,&n,
                   (LAPACK_REAL*)AA.Store(),
                   d.Store(), 
                   (LAPACK_REAL*)UU.Store(),
                   (LAPACK_REAL*)VV.Store(),
                   &info);
	if(info != 0) {
	    cout <<"dgesvd also failed, "<< "info = " << info << endl;
            Error("Error condition in dgesvd");
	  }
        }

    U.ReDimension(m,n);
    V.ReDimension(n,n);

    for(int i = 1; i <= n; ++i)
    for(int j = 1; j <= m; ++j)
        {
        U(j,i) = UU(i,j); 
        }

    for(int i = 1; i <= n; ++i)
    for(int j = 1; j <= n; ++j)
        {
        V(j,i) = VV(i,j); 
        }

  return;
}

void 
SVD_bf(const MatrixRef& Are, const MatrixRef& Aim, 
    Matrix& Ure, Matrix& Uim, 
    Vector& d, 
    Matrix& Vre, Matrix& Vim){
      LAPACK_INT m = Are.Nrows(), 
               n = Are.Ncols(); 
#ifdef DEBUG
    if(Aim.Nrows() != m || Aim.Ncols() != n)
        {
        Error("Aim must have same dimensions as Are");
        }
#endif

    if(m < n)
        {
        Matrix mret = Are.t(), 
               mimt = -Aim.t(),
               UUre,UUim,
               VVre,VVim;
        SVD_bf(mret,mimt, UUre, UUim, d, VVre, VVim);
        Vre = UUre.t();
        Vim = -UUim.t();
        Ure = VVre.t();
        Uim = -VVim.t();
        return;
        }

    char jobz = 'S';

    Matrix AA(n,2*m);
    for(int i = 1; i <= n; i++)
	for(int j = 1; j <= m; j++)
        {
	    AA(i,2*j-1) = Are(j,i); 
            AA(i,2*j) = Aim(j,i);
        }

    Matrix UU(n,2*m), VV(n,2*n);
    d.ReDimension(n);
    LAPACK_INT info = 0;

    zgesdd_wrapper(&jobz,&m,&n,
                   (LAPACK_COMPLEX*)AA.Store(),
                   d.Store(), 
                   (LAPACK_COMPLEX*)UU.Store(),
                   (LAPACK_COMPLEX*)VV.Store(),
                   &info);

    if(info != 0) 
        {
        cout <<"zgesdd failed, "<< "info = " << info<<". Now try zgesvd." << endl;
        char jobu = 'S', jobvt='S';
        for(int i = 1; i <= n; i++)
	for(int j = 1; j <= m; j++)
        {
	    AA(i,2*j-1) = Are(j,i); 
            AA(i,2*j) = Aim(j,i);
        }
        info=0;
        zgesvd_wrapper(&jobu,&jobvt,&m,&n,
                   (LAPACK_COMPLEX*)AA.Store(),
                   d.Store(), 
                   (LAPACK_COMPLEX*)UU.Store(),
                   (LAPACK_COMPLEX*)VV.Store(),
                   &info);
	if(info != 0) {
	    cout <<"zgesvd also failed, "<< "info = " << info << endl;
            Error("Error condition in zgesvd");
	  }
        }

    Ure.ReDimension(m,n);
    Uim.ReDimension(m,n);
    Vre.ReDimension(n,n);
    Vim.ReDimension(n,n);

    for(int i = 1; i <= n; ++i)
    for(int j = 1; j <= m; ++j)
        {
        Ure(j,i) = UU(i,2*j-1); 
        Uim(j,i) = UU(i,2*j);
        }

    for(int i = 1; i <= n; ++i)
    for(int j = 1; j <= n; ++j)
        {
        Vre(j,i) = VV(i,2*j-1); 
        Vim(j,i) = VV(i,2*j);
        }
    
  return;
}


template <typename Tensor>
void tensor_factor(Tensor const& T, Tensor& A, Tensor& B, Args const& args)
{
    auto name = args.getString("IndexName","c");
    Tensor D;
    tensor_svd(T,A,D,B,{args,"LeftIndexName",name});
    auto dl = commonIndex(A,D);
    auto dr = commonIndex(B,D);

    D.mapElems([](Real x){ return std::sqrt(std::fabs(x)); });
    A *= D;
    B *= D;
    //Replace index dl with dr
    Tensor delta=dag(D);
    delta.mapElems([](double x){ return ((x>0) ? 1:0); });
    A*=delta;

    //Print(dl);
    //Print(dr);
    //Print(A);
    //Print(B);
    //Print(D);
}
template void
tensor_factor(ITensor const& T,ITensor& A, ITensor& B, Args const& args);
template void
tensor_factor(IQTensor const& T, IQTensor& A, IQTensor& B, Args const& args);
