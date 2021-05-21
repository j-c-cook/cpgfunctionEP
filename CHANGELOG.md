# History of changes

## Current version

### Enhancements

- [Issue 8](https://github.com/j-c-cook/LinearAlgebra/issues/8) - 
  The BLAS and LAPACK functions are each put into their own namespace in side of 
  their own modules.  

## Version 1.0.0 (2021-05-21)

### New features

- [Issue 7](https://github.com/j-c-cook/LinearAlgebra/issues/7) - Matrix vector multiplication
  using upper triangular symmetric matrix (`spmv`) 

- [Issue 6](https://github.com/j-c-cook/LinearAlgebra/issues/6) - Add scalar ability `x = a*x` 

- [Issue 5](https://github.com/j-c-cook/LinearAlgebra/issues/5) - Adds LAPACK's copy function for 
  `std::vector`'s

- [Issue 4](https://github.com/j-c-cook/LinearAlgebra/issues/4) - Add ability to do products of vectors 

- [Issue 3](https://github.com/j-c-cook/LinearAlgebra/issues/3) - Add ability to compute matrix 
  multiplication of a matrix by a vector using `gemv` function 

- [Issue 2](https://github.com/j-c-cook/LinearAlgebra/issues/2) - Add ability to compute
  the dot product of two 1D vectors 

- [Issue 1](https://github.com/j-c-cook/LinearAlgebra/issues/1) - 
  Adds `LU` decomposition of `Ax=b` using `std::vector` by calling `gesv` function
  (same function is used in numpy.linalg.solve())
