cmake_minimum_required(VERSION 3.16)
project(dp)

set(CMAKE_CXX_STANDARD 14)

SET(CMAKE_CXX_COMPILER "/usr/local/bin/g++-14" CACHE STRING "C++ compiler" FORCE)

add_compile_options(-Wall -Wextra -Werror=all -g)

add_executable(
    dp
    main.cc
    graph.cc
    double_pendulum.cc
)