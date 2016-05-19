
#ifndef _TENSOR_SVD_H_
#define _TENSOR_SVD_H_

#include "utilities.h"

template<class Tensor> 
bool tensor_svd(Tensor AA, Tensor& U, Tensor& D, Tensor& V, const Args& args=Global::args());

bool
tensor_svdRank2(ITensor A, const Index& ui, const Index& vi,
         ITensor& U, ITensor& D, ITensor& V,
         const Args& args = Global::args());

bool 
tensor_svdRank2(IQTensor A, const IQIndex& uI, const IQIndex& vI,
         IQTensor& U, IQTensor& D, IQTensor& V,
         const Args& args = Global::args());

void 
SVD_bf(const MatrixRef& A, Matrix& U, Vector& D, Matrix& V);

void 
SVD_bf(const MatrixRef& Are, const MatrixRef& Aim, 
    Matrix& Ure, Matrix& Uim, 
    Vector& D, 
    Matrix& Vre, Matrix& Vim);

//
// The "factor" decomposition is based on the SVD, but factorizes a tensor T into only two tensors T=A*B where A and B share a single common index.
//
// If the SVD of T is T=U*S*V where S is a diagonal matrix of singular values, then A and B are schematically A=U*sqrt(S) and B=sqrt(S)*V.
//
// In addition to the named Args recognized by the svd routine, factor accepts an Arg "IndexName" which wil be the name of the common index connecting A and B.
template <typename Tensor>
void tensor_factor(Tensor const& T, Tensor &A, Tensor &B, Args const &args = Args::global());
#endif
