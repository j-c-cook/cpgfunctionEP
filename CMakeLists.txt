cmake_minimum_required(VERSION 3.6)
project(LinearAlgebra)

include_directories(include)

# Create a DLL with an API
add_library(LinearAlgebra SHARED src/_gesv.cpp)

# Depends on LAPACK and BLAS
find_package(LAPACK)
target_link_libraries(LinearAlgebra ${LAPACK_LIBRARIES})

# Open blas is specified because the default (Netlib BLAS) is much slower on my intel
set(BLA_VENDOR OpenBLAS)
find_package(BLAS)
if (BLAS_FOUND)
    target_link_libraries(LinearAlgebra ${BLAS_LIBRARIES})
else()
    # ERROR
endif()

# Create a test binary that depends on the DLL
add_executable(test1 test/gesv.cpp)

target_link_libraries(test1 LinearAlgebra)

# Turn on testing
enable_testing()

# Test the binary itself
add_test(NAME RunTest1 COMMAND ${CMAKE_BINARY_DIR}/test1)
