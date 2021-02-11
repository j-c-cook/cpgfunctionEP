cmake_minimum_required(VERSION 2.8)
project(CpGFunction)

set(CMAKE_CXX_STANDARD 11)
# Compiler requests -std=c++11 be set for compiler option
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Everything will need to access this
include_directories(include)

# Create a DLL with an API
add_library(cpgfunction SHARED src/boreholes.cpp src/gfunction.cpp src/heat_transfer.cpp src/interpolation.cpp src/SegmentResponse.cpp src/utilities.cpp)

# It depends on boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
# If the boost installation already on your computer is garbage or if you don't have one, uncomment this and
# manually set the path to BOOST_ROOT (https://stackoverflow.com/a/58447515/11637415)
# set(BOOST_ROOT /path/to/boost)
find_package(Boost)
message(STATUS "Boost_FOUND=${Boost_FOUND}")
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    # ERROR
endif()
target_link_libraries(cpgfunction ${Boost_LIBRARIES})

# It depends on threads
find_package(Threads REQUIRED)
target_link_libraries(cpgfunction ${CMAKE_THREAD_LIBS_INIT})

find_package(LAPACK)
target_link_libraries(cpgfunction ${LAPACK_LIBRARIES})

find_package(BLAS)
target_link_libraries(cpgfunction ${BLAS_LIBRARIES})

# Create a test binary that depends on the DLL
add_executable(test1 test/main.cpp)
target_link_libraries(test1 cpgfunction)


# Turn on testing
enable_testing()

# Test the test binary itself
add_test(NAME RunTest1 COMMAND "${CMAKE_BINARY_DIR}/test1")