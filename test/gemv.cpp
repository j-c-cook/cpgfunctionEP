//
// Created by jackcook on 5/14/21.
//

//extern "C" void dgemv_(char *Trans, int *m, int *n, double *alpha, double *A, int * lda,
//                       double *x, int * incx, double *beta, double *y, int *incy);

#include <LinearAlgebra/gemv.h>
#include <iostream>
#include <vector>


int main() {
    char trans = 'T';  // do a transpose to make it row-wise order

    // matrix A is a 2x5
    std::vector<double> A = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int m = 5;  // number of columns in matrix A (when trans = T)
    int n = 2;  // number of rows in matrix A (when trans = T)

    double alpha = 1.;  // constant to be multiplied by alpha * A * x

    int lda = m;  // specifies the first dimension of A, must be at least max(1, m)

    std::vector<double> x = {4, 8, 10, 6, 8};
    int incx = 1;

    double beta = 0.;
    std::vector<double> y(n, 0);  // y is a vector of size n full of zeros
    int incy = 1;

    la::_gemv::gemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy);

    for (int i = 0; i < y.size(); i++) {
        std::cout << y[i] << std::endl;
    }

    return 0;
}