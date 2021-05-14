//
// Created by jackcook on 5/14/21.
//

#include <vector>

using namespace std;

#ifndef LINEARALGEBRA_COPY_H
#define LINEARALGEBRA_COPY_H

namespace jcc {
    namespace la {

        extern "C" void dcopy_(int *n, double *x, int *incx, double *y, int *incy);
        void copy(int &n, vector<double> &x, int &incx, vector<double> &y, int &incy);

    }  // namespace la
}  // namespace jcc

#endif //LINEARALGEBRA_COPY_H
