//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/axpy.h>

namespace jcc {
    namespace la {

        void axpy(int &n, double &a, vector<double> &x, int &incx, vector<double> &y, int &incy) {
            daxpy_(&n, &a, &*x.begin(), &incx, &*y.begin(), &incy);
        }

    }  // namespace la
}  // namespace jcc