//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/spmv.h>

namespace jcc {
    namespace la {

        void spmv(char &uplo, int &n, double &alpha, vector<double> &A,
                  vector<double> &x, int &incx, double &beta, vector<double> &y, int &incy){
            dspmv_(&uplo, &n, &alpha, &*A.begin(), &*x.begin(), &incx, &beta, &*y.begin(), &incy);
        }

    }  // namespace la
}  // namespace jcc