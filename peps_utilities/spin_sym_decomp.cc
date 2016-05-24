
#include "spin_sym_decomp.h"

void spin_sym_svdRank2(IQTensor A, const IQIndex &ui, const IQIndex &vi, IQTensor &U, IQTensor &D, IQTensor &V, const Args &args)
{
    if (A.r()!=2) Error("A must be matrix-like");
    //we fix the convention that index of ui and vi have different direction, so that CG coeff is just identity
    if (ui.dir()==vi.dir()) Error("ui and vi should have different direction");

    //get the spin rep of ui,vi as well as the quantum number for each indexval of ui and vi
    std::vector<int> ui_spin_rep, vi_spin_rep;
    std::vector<Spin_Basis> ui_spin_basis, vi_spin_basis;
    iqind_spin_rep(ui,ui_spin_rep);
    iqind_spin_rep(vi,vi_spin_rep);
    iqind_to_spin_basis(ui,ui_spin_rep,ui_spin_basis);
    iqind_to_spin_basis(vi,vi_spin_rep,vi_spin_basis);

    //we can get the indval of a particular state |S,m,t\rangle of ui and vi
    std::vector<std::vector<std::vector<int>>> ui_spin_basis_to_indval, vi_spin_basis_to_indval;
    indval_from_spin_rep_basis(ui_spin_rep,ui_spin_basis,ui_spin_basis_to_indval);
    indval_from_spin_rep_basis(vi_spin_rep,vi_spin_basis,vi_spin_basis_to_indval);

    //perform SVD on the fixed S, and for each S, we only need to focus on one paticular m due to Wigner-Eckart theorem
    int S_max=std::min(ui_spin_rep.size(),vi_spin_rep.size());
    std::vector<Index> ui_SS(S_max), vi_SS(S_max), Lind_SS(S_max), Rind_SS(S_max);
    std::vector<ITensor> A_SS(S_max),U_SS(S_max), V_SS(S_max), D_SS(S_max);
    std::vector<int> Dflavor_deg(S_max,0);
    for (int spini=0; spini<S_max; spini++)
    {
        if (ui_spin_rep[spini]==0 || vi_spin_rep[spini]==0) continue;
        //create ITensor A_SS with fixed S and m=S
        ui_SS[spini]=Index(nameint("uspin_",spini),ui_spin_rep[spini]);
        vi_SS[spini]=Index(nameint("vspin_",spini),vi_spin_rep[spini]);
        A_SS[spini]=ITensor(ui_SS[spini],vi_SS[spini]);
        U_SS[spini]=ITensor(ui_SS[spini]);
        for (int ut=0; ut<ui_SS[spini].m(); ut++)
        {
            for (int vt=0; vt<vi_SS[spini].m(); vt++)
            {
                A_SS[spini](ui_SS[spini](ut+1),vi_SS[spini](vt+1))=A(ui(ui_spin_basis_to_indval[spini][spini][ut]+1),vi(vi_spin_basis_to_indval[spini][spini][vt]+1));
            }
        }
        //perform SVD on this A_SS
        tensor_svdRank2(A_SS[spini],ui_SS[spini],vi_SS[spini],U_SS[spini],D_SS[spini],V_SS[spini],args);
        Lind_SS[spini]=commonIndex(D_SS[spini],U_SS[spini]);
        Rind_SS[spini]=commonIndex(D_SS[spini],V_SS[spini]);
        Dflavor_deg[spini]=Lind_SS[spini].m();
    }
    
    //get U,D,V
    IQIndex Lind=Spin_leg(Dflavor_deg,"Lind",ui.dir()),
            Rind=Spin_leg(Dflavor_deg,"Rind",vi.dir());
    //the left leg and right leg of D should be isomorphic (except direction)
    std::vector<Spin_Basis> Dspin_basis;
    std::vector<std::vector<std::vector<int>>> Dspin_basis_to_indval;
    iqind_to_spin_basis(Lind,Dspin_basis);
    indval_from_spin_rep_basis(Dflavor_deg,Dspin_basis,Dspin_basis_to_indval);

    D=IQTensor(Lind,Rind);
    U=IQTensor(ui,dag(Lind));
    V=IQTensor(vi,dag(Rind));

    for (int spini=0; spini<S_max; spini++)
    {
        for (int mi=-spini; mi<=spini; mi+=2)
        {
            //get U
            for (int ut=0; ut<ui_spin_rep[spini]; ut++)
            {
                for (int lt=0; lt<Dflavor_deg[spini]; lt++)
                {
                    int ui_val=ui_spin_basis_to_indval[spini][(spini+mi)/2][ut],
                        Lval=Dspin_basis_to_indval[spini][(spini+mi)/2][lt];
                    U(ui(ui_val+1),dag(Lind)(Lval+1))=U_SS[spini](ui_SS[spini](ut+1),dag(Lind_SS[spini](lt+1)));
                }
            }
            //get V
            for (int vt=0; vt<vi_spin_rep[spini]; vt++)
            {
                for (int rt=0; rt<Dflavor_deg[spini]; rt++)
                {
                    int vi_val=vi_spin_basis_to_indval[spini][(spini+mi)/2][vt],
                        Rval=Dspin_basis_to_indval[spini][(spini+mi)/2][rt];
                    V(vi(vi_val+1),dag(Rind)(Rval+1))=V_SS[spini](vi_SS[spini](vt+1),dag(Rind_SS[spini](rt+1)));
                }
            }
            //get D
            for (int dt=0; dt<Dflavor_deg[spini]; dt++)
            {
                int dval=Dspin_basis_to_indval[spini][(spini+mi)/2][dt];
                D(Lind(dval+1),Rind(dval+1))=D_SS[spini](Lind_SS[spini](dt+1),Rind_SS[spini](dt+1));
            }
        }
    }
}
