//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/scal.h>

#include <iostream>
#include <vector>

int main(){

    double a = -1;
    std::vector<double> x = {2, 5, 8};
    int n = x.size();
    int incx = 1;

    jcc::la::scal(n, a, x, incx);

    for (double v : x) {
        std::cout << v << std::endl;
    }

    return 0;
}