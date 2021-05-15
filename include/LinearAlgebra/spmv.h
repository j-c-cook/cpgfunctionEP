//
// Created by jackcook on 5/14/21.
//

#include <vector>

using namespace std;

#ifndef LINEARALGEBRA_SPMV_H
#define LINEARALGEBRA_SPMV_H

namespace jcc {
    namespace la {

        extern "C" void dspmv_(char *uplo, int *n, double *alpha, double *A,
                               double *x, int *incx, double *beta, double *y, int *incy);
        void spmv(char &uplo, int &n, double &alpha, vector<double> &A,
                  vector<double> &x, int &incx, double &beta, vector<double> &y, int &incy);

    }  // namespace la
}  // namespace jcc

#endif //LINEARALGEBRA_SPMV_H
