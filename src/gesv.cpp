//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/gesv.h>

namespace la {
    namespace _gesv {

        void gesv(int &n, int &nrhs, vector<double> &a, int &lda, vector<int> &ipiv, vector<double> &b, int &lbd, int &info) {
            dgesv_(&n, &nrhs, &*a.begin(), &lda, &*ipiv.begin(), &*b.begin(), &lbd, &info);
        }

    }  // namespace _gesv
}  // namespace la
