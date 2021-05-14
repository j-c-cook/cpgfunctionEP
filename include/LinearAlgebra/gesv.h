//
// Created by jackcook on 5/14/21.
//

#include <iostream>
#include <vector>

using namespace std;

#ifndef LINEARALGEBRA_DGESV_H
#define LINEARALGEBRA_DGESV_H

namespace jcc {
    namespace la {

        extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *lbd, int *info);

        void gesv(int &n, int &nrhs, vector<double> &a, int &lda, vector<int> &ipiv, vector<double> &b, int &lbd, int &info);

    }  // namespace la
}  // namespace jcc

#endif //LINEARALGEBRA_DGESV_H
