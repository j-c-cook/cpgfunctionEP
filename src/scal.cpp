//
// Created by jackcook on 5/14/21.
//

#include <LinearAlgebra/scal.h>

namespace jcc {
    namespace la {

        void scal(int &n, double &a, vector<double> &x, int &incx){
            dscal_(&n, &a, &*x.begin(), &incx);
        }

    }  // namespace la
}  // namespace jcc