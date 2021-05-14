//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/gesv.h>
#include <iostream>
#include <vector>

using namespace std;

namespace jcc {
    namespace la {

        void gesv(int &n, int &nrhs, vector<double> &a, int &lda, vector<int> &ipiv, vector<double> &b, int &lbd, int &info) {
            dgesv_(&n, &nrhs, &*a.begin(), &lda, &*ipiv.begin(), &*b.begin(), &lbd, &info);
        }

    }  // namespace la
}  // namespace jcc
