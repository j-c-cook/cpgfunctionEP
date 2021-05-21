//
// Created by jackcook on 5/21/21.
//

#include <LinearAlgebra/lapack.h>

namespace jcc {
namespace lapack {

// gesv solves for x in Ax=b
// double gesv (d_gesv)
void gesv(int &n, int &nrhs, vector<double> &a, int &lda, vector<int> &ipiv,
          vector<double> &b, int &lbd, int &info) {
    dgesv_(&n, &nrhs, &*a.begin(), &lda, &*ipiv.begin(), &*b.begin(), &lbd,
           &info);
}  // d_gesv();

}  // namespace lapack
}  // namespace jcc