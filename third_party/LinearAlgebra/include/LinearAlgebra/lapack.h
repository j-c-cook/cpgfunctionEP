//
// Created by jackcook on 5/21/21.
//

#include <vector>

using namespace std;

#ifndef LINEARALGEBRA_LAPACK_H
#define LINEARALGEBRA_LAPACK_H

namespace jcc {
namespace lapack {

// gesv solves for x in Ax=b
// double gesv (d_gesv)
extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv,
                       double *b, int *lbd, int *info);
void gesv(int &n, int &nrhs, vector<double> &a, int &lda, vector<int> &ipiv,
          vector<double> &b, int &lbd, int &info);

}  // namespace lapack
}  // namespace jcc

#endif //LINEARALGEBRA_LAPACK_H
