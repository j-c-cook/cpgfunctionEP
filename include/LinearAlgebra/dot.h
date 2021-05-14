//
// Created by jackcook on 5/14/21.
//

#include <iostream>
#include <vector>

#ifndef LINEARALGEBRA__DOT_H
#define LINEARALGEBRA__DOT_H

namespace jcc {
    namespace la {

        extern "C" double ddot_(int *n, double *x, int *incx, double *y, int * incy);
        double dot(int &n, std::vector<double> &x, int &incx, std::vector<double> &y, int &incy);

    }  // namespace la
}  // namespace jcc

#endif //LINEARALGEBRA__DOT_H
