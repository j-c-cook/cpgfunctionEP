# LinearAlgebra
Linear algebra in C++ using BLAS and LAPACK libraries

## BLAS functions

- Level 1
  - [scal](test/scal.cpp) - Scale vector `x` by `a`, `x = a * x`
  - [copy](test/copy.cpp) - Copy vector `x` into vector `y`
  - [axpy](test/axpy.cpp) - Vector product `y = a*x + y`
  - [dot](test/dot.cpp) - Computes the dot product of two 1D vectors
- Level 2 
  - [gemv](test/gemv.cpp) - Performs matrix multiplication of a matrix by a vector `y = alpha*A*x + beta*y`
  - [spmv](test/spmv.cpp) - Matrix vector mutiplication when `A` is a packed symmetric matrix, `y = alpha*A*x + beta*y`

## LAPACK functions
- [gesv](test/gesv.cpp) - LU decomposition to solve system of linear equations for x in `Ax=b` 
