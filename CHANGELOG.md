# History of changes (cpgfunctionEP)

## Version 0.1.0 (2021-xx-xx)

### Maintenance

* [Issue 35](https://github.com/j-c-cook/cpgfunction/issues/35) - The C++ API
  library for an interface to BLAS and LAPACK, 
  [LinearAlgebra](https://github.com/j-c-cook/LinearAlgebra), is made into a 
  `subtree` rather than `submodule`.

* [Issue 50](https://github.com/j-c-cook/cpgfunction/issues/50) - Include boost
  files in `third_party/` rather than CMake finding the library in path. Boost
  is stripped to only the required libraries so that the size of the dependency
  is reduced.
  
### Fixes
  
* [Issue 10](https://github.com/j-c-cook/cpgfunctionEP/issues/10) - Vector over flow error associated with number of 
  time steps being less than 6 is fixed. 
  Commit [558c63c](https://github.com/j-c-cook/cpgfunctionEP/commit/558c63c9b5c4041282fdb473327463b8fe183fbd)
  fixes a vector overflow error that occured when the number of time steps described were less than 6 (`nt<6`). Commit
  [94e222b](https://github.com/j-c-cook/cpgfunctionEP/commit/94e222bb0cf670f1a4f3058bd1b7b6c7f8bf8233) now allows the 
  function time_geometric to have less than 6 time steps. 


# History of changes (cpgfunction)

## Version 2.0.0 (2021-05-23)

### Enhancements

* [Issue 25](https://github.com/j-c-cook/cpgfunction/issues/25) - Removes all references to the 3D `h_ij`
  segment response matrix. See [PR 30](https://github.com/j-c-cook/cpgfunction/pull/30).

* [Issue 32](https://github.com/j-c-cook/cpgfunction/issues/32) - The multi-dimensional matrices, 
  `q_reconstructed` and `h_ij`, are made one dimensional prior to passage into the temporal superposition
  function so that `BLAS` routines can be heavily depended on and the loops completely unraveled.  
  See [PR 30](https://github.com/j-c-cook/cpgfunction/pull/30).

* [Issue 33](https://github.com/j-c-cook/cpgfunction/issues/33) - It is found that the 
  packed segment resopnse matrix can be directly made us of in `BLAS spmv`, and that addition
  greatly optimizes the temporal superposition function. For now the assumption is made that all
  segments in the field are of equivalent length, which is true and fine, but at some point in the
  future unequal segment lengths should be made possible again. 
  See [PR 30](https://github.com/j-c-cook/cpgfunction/pull/30).

### New features

* [Issue 28](https://github.com/j-c-cook/cpgfunction/issues/28) -
  The third party library LinearAlgebra (`jcc:la`) is included and made use of for `LU`
  factorization in `gfunction.cpp`

* [Issue 12](https://github.com/j-c-cook/cpgfunction/issues/12) -
  A boolean toggle option is added for multi-threading for computing the 
  uniform borehole wall temperature (UBHWT) g-function

### API Changes

* [Issue 16](https://github.com/j-c-cook/cpgfunction/issues/16) - The `uniform borehole wall temperature` 
  g-function definition is defined for planned use in EnergyPlus with all arguments. Not all the arguments
  currently have a purpose, the adaptive discretization and number of thread arguments are place holders.

## Version 1.0.0 (2021-05-12)

### New features

* [Issue 20](https://github.com/j-c-cook/cpgfunction/issues/20) - 
  Added OpenBlas as the basic linear algebra subprogram (BLAS) vendor to CMakeLists.txt

* [Issue 18](https://github.com/j-c-cook/cpgfunction/issues/18) - 
  Added new borefield interface with API access to typical borehole configurations

* [Issue 13](https://github.com/j-c-cook/cpgfunction/issues/13) - 
  Implemented g-function accuracy tests via CMakeLists.txt for a Rectangle, Open Rectangle, U shape, 
  L shape and a custom (Poisson disk) configuration

* [Commit 45141fa](https://github.com/j-c-cook/cpgfunction/pull/14/commits/45141fa745d92ac8a08eea2a06801d7a01fac367) - 
  Create new uniform borehole wall temperature API to consider the new borefield and time API's

* [Commit f8863ad](https://github.com/j-c-cook/cpgfunction/pull/14/commits/f8863ad6879bdcb43d8bbed48ab1be1701eb56f5) - 
  Added time vector API and associated test

* [Commit 654160f](https://github.com/j-c-cook/cpgfunction/pull/14/commits/654160f9b508f57b917fc0630437cff726dc8440) - 
  Modify API for creating a vector of boreholes (borefield)




