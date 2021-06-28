//
// Created by jackcook on 6/27/21.
//

#include <vector>

using namespace std;

#ifndef LINEARALGEBRA_AXPY_H
#define LINEARALGEBRA_AXPY_H

namespace jcc  { namespace blas {

    void axpy(int &n, double &a, vector<double> &x, vector<double> &y,
              int &start, int &n_threads);

}  // namespace blas
}  // namespace jcc



#endif //LINEARALGEBRA_AXPY_H
