//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/blas.h>

#include <algorithm>
#include <iterator>

#include <iostream>
#include <vector>

int main() {

    std::vector<double> x = {2, 8, 9, 10, 8, 4, 3};
    std::vector<double> y(x.size(), 0);

    int n = x.size();
    int incx = 1;
    int incy = 1;

//    jcc::blas::copy(n, x, incx, y, incy);

    std::vector<double>::iterator begin;
    std::vector<double>::iterator end;

    int a = 0;

    begin = x.begin() + a;
    end = x.begin() + n;

    std::copy(begin, end, y.begin());

    for (double v : y) {
        std::cout << v << std::endl;
    }

    return 0;
}

