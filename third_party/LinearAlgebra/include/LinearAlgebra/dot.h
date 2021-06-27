//
// Created by jackcook on 6/27/21.
//

#include <vector>

using namespace std;

#ifndef LINEARALGEBRA_DOT_H
#define LINEARALGEBRA_DOT_H

namespace jcc { namespace blas {

    double dot(int &n, vector<double> &x, vector<double> &y, int &start,
               int &n_threads);

}  // namespace blas
}  // namespace jcc

#endif //LINEARALGEBRA_DOT_H
