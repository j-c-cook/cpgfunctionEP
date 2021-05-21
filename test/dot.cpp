//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/blas.h>

#include <iostream>
#include <vector>

int main() {

    int n = 3;
    std::vector<double> x = {1.2, 2.4, 3.8};
    std::vector<double> y = {4.8, 5.5, 6.2};
    int incx = 1;
    int incy = 1;
    double result;

    result = jcc::blas::dot(n, x, incx, y, incy);

    std::cout.precision(32);
    std::cout << result << std::endl;

    return 0;
}