//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/copy.h>

#include <iostream>
#include <vector>

int main() {

    std::vector<double> x = {2, 8, 9};
    std::vector<double> y(x.size(), 0);

    int n = x.size();
    int incx = 1;
    int incy = 1;

    jcc::la::copy(n, x, incx, y, incy);

    for (double v : y) {
        std::cout << v << std::endl;
    }

    return 0;
}

