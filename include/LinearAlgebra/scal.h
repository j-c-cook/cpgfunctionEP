//
// Created by jackcook on 5/14/21.
//

#include <vector>

using namespace std;

#ifndef LINEARALGEBRA_SCAL_H
#define LINEARALGEBRA_SCAL_H

namespace jcc {
    namespace la {

        extern "C" void dscal_(int *n, double *a, double *x, int *incx);
        void scal(int &n, double &a, vector<double> &x, int &incx);

    }  // namespace la
}  // namespace jcc

#endif //LINEARALGEBRA_SCAL_H
