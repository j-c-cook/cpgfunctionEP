//
// Created by jackcook on 5/14/21.
//

#ifndef LINEARALGEBRA_DGESV_H
#define LINEARALGEBRA_DGESV_H

#include <iostream>
#include <vector>

using namespace std;

namespace la {
    namespace _gesv {

        extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *lbd, int *info);

        void gesv(int &n, int &nrhs, vector<double> &a, int &lda, vector<int> &ipiv, vector<double> &b, int &lbd, int &info);

    }  // namespace _gesv
}  // namespace la

#endif //LINEARALGEBRA_DGESV_H
