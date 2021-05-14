cmake_minimum_required(VERSION 3.6)
project(LinearAlgebra)

include_directories(include)

# Create a DLL with an API
add_library(LinearAlgebra SHARED
        src/gesv.cpp
        src/dot.cpp
        src/gemv.cpp
        src/axpy.cpp
        src/copy.cpp
        src/scal.cpp)

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
add_executable(test2 test/dot.cpp)
add_executable(test3 test/gemv.cpp)
add_executable(test4 test/axpy.cpp)
add_executable(test5 test/copy.cpp)
add_executable(test6 test/scal.cpp)

target_link_libraries(test1 LinearAlgebra)
target_link_libraries(test2 LinearAlgebra)
target_link_libraries(test3 LinearAlgebra)
target_link_libraries(test4 LinearAlgebra)
target_link_libraries(test5 LinearAlgebra)
target_link_libraries(test6 LinearAlgebra)

# Turn on testing
enable_testing()

# Test the binary itself
add_test(NAME RunTest1 COMMAND ${CMAKE_BINARY_DIR}/test1)
add_test(NAME RunTest2 COMMAND ${CMAKE_BINARY_DIR}/test2)
add_test(NAME RunTest3 COMMAND ${CMAKE_BINARY_DIR}/test3)
add_test(NAME RunTest4 COMMAND ${CMAKE_BINARY_DIR}/test4)
add_test(NAME RunTest5 COMMAND ${CMAKE_BINARY_DIR}/test5)
add_test(NAME RunTest6 COMMAND ${CMAKE_BINARY_DIR}/test6)
