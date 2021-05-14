//
// Created by jackcook on 5/14/21.
//

#include <iostream>
#include <vector>

using namespace std;

#ifndef LINEARALGEBRA_AXPY_H
#define LINEARALGEBRA_AXPY_H

namespace jcc {
    namespace la {

        extern "C" void daxpy_(int *n, double *a, double *x, int *incx, double *y, int *incy);
        void axpy(int &n, double &a, vector<double> &x, int &incx, vector<double> &y, int &incy);

    }  // namespace la
}  // namespace jcc

#endif //LINEARALGEBRA_AXPY_H
