//
// Created by jackcook on 6/27/21.
//

#include <LinearAlgebra/axpy.h>

namespace jcc  { namespace blas {

    void axpy(int &n, double &a, vector<double> &x, vector<double> &y,
              int &start, int &n_threads){
        // y = a*x + y
        #pragma omp parallel for num_threads(n_threads)
        for (int i=0; i<n; i++) {
            y[i] = a * x[i+start] + y[i];
        }

    }  // axpy();

}  // namespace blas
}  // namespace jcc