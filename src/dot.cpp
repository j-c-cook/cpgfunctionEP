//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/dot.h>

#include <iostream>
#include <vector>

namespace la {
    namespace _dot {

        double dot(int &n, std::vector<double> &x, std::vector<double> &y, int incx=1, int incy=1) {
            double result = ddot_(&n, &*x.begin(), &incx, &*y.begin(), &incy);
            return result;
        }

    }  // namespace _dot
}  // namespace la