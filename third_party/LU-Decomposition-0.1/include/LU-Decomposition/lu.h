//
// Created by jackcook on 8/17/21.
//

#include <vector>
#include <cmath>
#include <stdexcept>

#ifndef LU_DECOMPOSITION_LU_H
#define LU_DECOMPOSITION_LU_H

namespace jcc {

    void CheckSingularity(std::vector<double> &A, int &n);

    void CroutDecomposition(std::vector<double> &A, int &n, std::vector<int> &indx);

    void CroutSolve(std::vector<double> &LU, std::vector<double> &b, int &n,
                    std::vector<int> &indx);

}

#endif //LU_DECOMPOSITION_LU_H
