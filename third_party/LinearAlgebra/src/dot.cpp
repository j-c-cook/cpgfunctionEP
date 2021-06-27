//
// Created by jackcook on 6/27/21.
//

#include <LinearAlgebra/dot.h>

namespace jcc { namespace blas {

    double dot(int &n, vector<double> &x, vector<double> &y, int &start,
               int &n_threads) {
        double dot_product = 0;

        #pragma omp parallel for num_threads(n_threads)
        for (int i = start; i<n; i++) {
            dot_product += x[i] * y[i];
        }
        return dot_product;
    }  // dot();
    }  // namespace blas
}  // namespace jcc