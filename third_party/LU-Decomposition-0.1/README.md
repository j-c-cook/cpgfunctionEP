# LU-Decomposition
The purpose is to provide a fast, standalone, lightweight native C++ library that can solve a system of equations faster than Eigen and OpenBLAS (Fortran). Currently, neither of those are true.

The C++ standard library techniques for performing dot products and vector products appear to be slow due to the use of only a single thread. The construction and then added dependency of a native C++ BLAS library with a heavy dependency on OpenMP should provide a speed increase. 

The current method implemented is Crout's algorithm with implicit pivoting, inspired by Press et al. (2000). Any future work should be put into developing a multi-threaded BLAS library. The recursive tile technique of Dongarra et al. (2014) seems promising.  

## References

Dongarra, J., Faverge, M., Ltaief, H., & Luszczek, P. (2014). Achieving numerical accuracy and high performance using recursive tile LU factorization with partial pivoting. Concurrency and Computation, 26(7), 1408–1431. https://doi.org/10.1002/cpe.3110

Press, William H. Numerical Recipes in C++ : the Art of Scientific Computing . 2nd ed. Cambridge, UK ;: Cambridge University Press, 2002. Print.
