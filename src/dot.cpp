//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/dot.h>

#include <iostream>
#include <vector>

namespace jcc {
    namespace la {

        double dot(int &n, std::vector<double> &x, int &incx, std::vector<double> &y, int &incy) {
            double result = ddot_(&n, &*x.begin(), &incx, &*y.begin(), &incy);
            return result;
        }

    }  // namespace la
}  // namespace jcc