//
// Created by jackcook on 5/14/21.
//

#include <iostream>
#include <vector>

#ifndef LINEARALGEBRA__DOT_H
#define LINEARALGEBRA__DOT_H

namespace la {
    namespace _dot {

        extern "C" double ddot_(int *n, double *x, int *incx, double *y, int * incy);
        double dot(int &n, std::vector<double> &x, std::vector<double> &y, int incx, int incy);

    }  // namespace _dot
}  // namespace la

#endif //LINEARALGEBRA__DOT_H
