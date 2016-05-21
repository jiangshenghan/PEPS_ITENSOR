
#include "spin_sym_decomp.h"

void spin_sym_svdRank2(IQTensor A, const IQIndex &ui, const Index &vi, IQTensor &U, IQTensor &D, IQTensor &V, const Args &args=Global::args())
{
    if (A.r()!=2) Error("A must be matrix-like");
    //we fix the convention that index of ui and vi have different direction, so that CG coeff is just identity
    if (ui.dir()==vi.dir()) Error("ui and vi should have different direction");

    //get the spin rep of ui,vi as well as the quantum number for each indexval of ui and vi
    std::vector<int> ui_spin_rep, vi_spin_rep;
    std::vector<Spin_Basis> ui_spin_basis, vi_spin_basis;
    iqind_spin_rep(ui,ui_spin_rep);
    iqind_spin_rep(vi,vi_spin_rep);
    iqind_to_spin_basis(ui_spin_basis,vi_spin_basis);

    //perform SVD on the fixed S, and for each S, we only need to focus on one paticular m due to Wigner-Eckart theorem
    for (int spini=0; spini<std::min(ui_spin_rep.size(),vi_spin_rep.size()); spini++)
    {
        if (ui_spin_rep[spini]==0 || vi_spin_rep[spini]==0) continue;
        //create ITensor with fixed S and m=S
        Index sm_ui(nameint("uspin_",spini),ui_spin_rep[spini]),
              sm_vi(nameint("vspin_",spini),vi_spin_rep[spini]);
        ITensor A_S(sm_ui,sm_vi);
    }
}
