//
// Created by jackcook on 5/14/21.
//

#include <iostream>
#include <vector>

using namespace std;

#ifndef LINEARALGEBRA_GEMV_H
#define LINEARALGEBRA_GEMV_H

namespace jcc {
    namespace la {

        extern "C" void dgemv_(char *Trans, int *m, int *n, double *alpha, double *A, int * lda,
                               double *x, int * incx, double *beta, double *y, int *incy);

        void gemv(char &trans, int &m, int &n, double &alpha, vector<double> &A, int &lda, vector<double> &x,
                  int &incx, double &beta, vector<double> &y, int &incy);

    }  // namespace la
}  // namespace jcc

#endif //LINEARALGEBRA_GEMV_H
