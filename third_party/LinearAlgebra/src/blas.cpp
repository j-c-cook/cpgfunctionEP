//
// Created by jackcook on 5/21/21.
//

#include <LinearAlgebra/blas.h>

namespace jcc {
namespace blas {

// Level 1 routines

// - scal x = a*x
// double scal (dscal)
void scal(int &n, double &a, vector<double> &x, int &incx) {
    dscal_(&n, &a, &*x.begin(), &incx);
}  // d_scal();

// - copy x into y
// double copy (dcopy)
void copy(int &n, vector<double> &x, int &incx, vector<double> &y, int &incy) {
    dcopy_(&n, &*x.begin(), &incx, &*y.begin(), &incy);
}  // d_copy();

// - axpy y = a*x + y
// double axpy (daxpy)
void axpy(int &n, double &a, vector<double> &x, int &incx, vector<double> &y,
          int &incy) {
    daxpy_(&n, &a, &*x.begin(), &incx, &*y.begin(), &incy);
}  // d_axpy();

// -- Level 2 routines --

// - gemv matrix vector multiply y = alpha*A*x + beta*y
// double gemv (d_gemv)
void gemv(char &trans, int &m, int &n, double &alpha, vector<double> &A,
          int &lda, vector<double> &x, int &incx, double &beta,
          vector<double> &y, int &incy) {
    dgemv_(&trans, &m, &n, &alpha, &*A.begin(), &lda, &*x.begin(), &incx,
           &beta, &*y.begin(), &incy);
}  // d_gemv();

// - spmv symmetric packed matrix vector multiply y = alpha*A*x + beta*y
// double spmv (d_spmv)
void spmv(char &uplo, int &n, double &alpha, vector<double> &A,
          vector<double> &x, int &incx, double &beta, vector<double> &y,
          int &incy){
    dspmv_(&uplo, &n, &alpha, &*A.begin(), &*x.begin(), &incx, &beta,
           &*y.begin(), &incy);
}  // d_spmv();

// - dot computes the dot product of two vectors
// double dot (ddot)
double dot(int &n, std::vector<double> &x, int &incx, std::vector<double> &y,
           int &incy) {
    double result = ddot_(&n, &*x.begin(), &incx, &*y.begin(), &incy);
    return result;
}  // d_dot();

}  // namespace blas
}  // namespace jcc

