//
// Created by jackcook on 5/21/21.
//

#include <vector>

using namespace std;

#ifndef LINEARALGEBRA_BLAS_H
#define LINEARALGEBRA_BLAS_H

namespace jcc {
namespace blas {

// -- Level 1 routines --

// - scal x = a*x
// double scal (dscal)
extern "C" void dscal_(int *n, double *a, double *x, int *incx);
void scal(int &n, double &a, vector<double> &x, int &incx);

// - copy x into y
// double copy (dcopy)
extern "C" void dcopy_(int *n, double *x, int *incx, double *y, int *incy);
void copy(int &n, vector<double> &x, int &incx, vector<double> &y, int &incy);

// - axpy y = a*x + y
// double axpy (daxpy)
extern "C" void daxpy_(int *n, double *a, double *x, int *incx, double *y,
                       int *incy);
void axpy(int &n, double &a, vector<double> &x, int &incx, vector<double> &y,
          int &incy);

// - dot computes the dot product of two vectors
// double dot (ddot)
extern "C" double ddot_(int *n, double *x, int *incx, double *y, int * incy);
double dot(int &n, std::vector<double> &x, int &incx, std::vector<double> &y,
           int &incy);

// -- Level 2 routines --

// - gemv matrix vector multiply y = alpha*A*x + beta*y
// double gemv (d_gemv)
extern "C" void dgemv_(char *Trans, int *m, int *n, double *alpha, double *A,
                       int * lda, double *x, int * incx, double *beta,
                       double *y, int *incy);

void gemv(char &trans, int &m, int &n, double &alpha, vector<double> &A,
          int &lda, vector<double> &x, int &incx, double &beta,
          vector<double> &y, int &incy);

// - spmv symmetric packed matrix vector multiply y = alpha*A*x + beta*y
// double spmv (d_spmv)
extern "C" void dspmv_(char *uplo, int *n, double *alpha, double *A, double *x,
                       int *incx, double *beta, double *y, int *incy);
void spmv(char &uplo, int &n, double &alpha, vector<double> &A,
          vector<double> &x, int &incx, double &beta, vector<double> &y,
          int &incy);

}  // namespace blas
}  // namespace jcc

#endif //LINEARALGEBRA_BLAS_H
