//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/copy.h>

namespace jcc {
    namespace la {

        void copy(int &n, vector<double> &x, int &incx, vector<double> &y, int &incy) {
            dcopy_(&n, &*x.begin(), &incx, &*y.begin(), &incy);
        }

    }  // namespace la
}  // namespace jcc