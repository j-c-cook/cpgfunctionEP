cmake_minimum_required(VERSION 3.6)
project(CpGFunctionEP)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++17")

# Everything will need to access this
include_directories(include)
# Third party includes
include_directories(third_party/boost-1.71.0)  # It depends on boost (stripped)
include_directories(third_party/nlohmann-3.9.1)
include_directories(third_party/eigen-3.3.9)  # Eigen for LU decomposition

add_compile_definitions(_USE_MATH_DEFINES)

# Create a DLL with an API
add_library(cpgfunction SHARED
        src/boreholes.cpp
        src/gfunction.cpp
        src/heat_transfer.cpp
        src/interpolation.cpp
        src/utilities.cpp
        src/coordinates.cpp
        src/statistics.cpp
        src/blas.cpp)


# It depends on threads
find_package(Threads REQUIRED)
target_link_libraries(cpgfunction ${CMAKE_THREAD_LIBS_INIT})

OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
    FIND_PACKAGE(OpenMP)
    IF(OPENMP_FOUND)
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    ENDIF()
ENDIF()

# Create a test binary that depends on the DLL
add_executable(test1 test/gFunction_minimal.cpp)
add_executable(test2 test/interpolation.cpp)
# add_executable(test3 test/configurations.cpp)
add_executable(test4 test/borefield_definition.cpp)
add_executable(test5 test/time_definition.cpp)
add_executable(test6 test/compute_UBHWT_gFunction.cpp)
add_executable(test7 test/packed.cpp)

target_link_libraries(test1 cpgfunction)
target_link_libraries(test2 cpgfunction)
# target_link_libraries(test3 cpgfunction)
target_link_libraries(test4 cpgfunction)
target_link_libraries(test5 cpgfunction)
target_link_libraries(test6 cpgfunction)
target_link_libraries(test7 cpgfunction)

# Copy validation files to build directory so tests can open
file(GLOB JsonValidation test/validation/*.json)
file(GLOB CustomLayout test/custom_configurations/*.json)

file(COPY ${JsonValidation} DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CustomLayout} DESTINATION ${CMAKE_BINARY_DIR})
message(STATUS "The binary directory: ${CMAKE_BINARY_DIR}")

# Turn on testing
enable_testing()

# Test the test binary itself
add_test(NAME RunTest1 COMMAND "${CMAKE_BINARY_DIR}/test1")
add_test(NAME RunTest2 COMMAND "${CMAKE_BINARY_DIR}/test2")
# add_test(NAME RunTest3 COMMAND "${CMAKE_BINARY_DIR}/test3")
add_test(NAME RunTest4 COMMAND "${CMAKE_BINARY_DIR}/test4")
add_test(NAME RunTest5 COMMAND "${CMAKE_BINARY_DIR}/test5")
# Pass variable path into test 6 for json files
add_test(NAME RunTest6 COMMAND ${CMAKE_BINARY_DIR}/test6)
add_test(NAME RunTest7 COMMAND ${CMAKE_BINARY_DIR}/test7)
