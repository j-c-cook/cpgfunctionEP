//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/blas.h>

#include <iostream>
#include <vector>


int main(){

    char uplo = 'U';  // upper triangular in packed form

    int n = 2;
    double alpha = 1.;
    std::vector<double> A = {1, 2, 3};

    std::vector<double> x = {4, 8};
    int incx = 1;

    double beta = 1.;
    std::vector<double> y(n);
    int incy = 1;

    jcc::blas::spmv(uplo, n, alpha, A, x, incx, beta, y, incy);

    for (double v : y) {
        std::cout << v << std::endl;
    }

    return 0;
}