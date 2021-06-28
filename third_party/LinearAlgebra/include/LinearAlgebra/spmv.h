//
// Created by jackcook on 6/28/21.
//

#include <vector>

using namespace std;

#ifndef LINEARALGEBRA_SPMV_H
#define LINEARALGEBRA_SPMV_H

namespace jcc { namespace blas {

    void spmv(int &n, double &alpha, vector<double> &A, vector<double> &x,
              double &beta, vector<double> &y, int &start, int &n_threads);

}  // namespace blas
}  // namespace jcc

#endif //LINEARALGEBRA_SPMV_H
