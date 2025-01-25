cmake_minimum_required(VERSION 3.28)
project(double-pendulum LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
# SET(CMAKE_CXX_COMPILER "/usr/local/bin/g++-14" CACHE STRING "C++ compiler" FORCE)

add_executable(
    dp
    src/main.cc
    src/graph.cc
    src/double_pendulum.cc
)
 target_compile_features(dp PRIVATE cxx_std_17)
