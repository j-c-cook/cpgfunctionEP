//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/gemv.h>

namespace la {
    namespace _gemv {

        void gemv(char &trans, int &m, int &n, double &alpha, vector<double> &A, int &lda, vector<double> &x,
                  int &incx, double &beta, vector<double> &y, int &incy) {
            dgemv_(&trans, &m, &n, &alpha, &*A.begin(), &lda, &*x.begin(), &incx, &beta, &*y.begin(), &incy);
        }

    }  // namespace _gemv
}  // namespace la