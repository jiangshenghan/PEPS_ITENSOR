
//This file declares decomposition method (QR, SVD...) for spin symmetric tensor

#ifndef _SPIN_SYM_DECOMP_
#define _SPIN_SYM_DECOMP_

#include "singlet_tensor_basis.h"
#include "tensor_svd.h"

void spin_sym_svdRank2(IQTensor A, const IQIndex &ui, const IQIndex &vi, IQTensor &U, IQTensor &D, IQTensor &V, const Args &args=Global::args());

#endif
