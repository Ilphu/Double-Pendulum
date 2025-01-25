# cmake_minimum_required(VERSION 3.16)
# project(dp)

# set(CMAKE_CXX_STANDARD 14)

# SET(CMAKE_CXX_COMPILER "/usr/local/bin/g++-14" CACHE STRING "C++ compiler" FORCE)

# add_compile_options(-g)

# add_executable(
#     dp
#     main.cc
#     graph.cc
#     double_pendulum.cc
# )

cmake_minimum_required(VERSION 3.28)
project(double-pendulum LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_CXX_COMPILER "/usr/local/bin/g++-14" CACHE STRING "C++ compiler" FORCE)

# include(FetchContent)
# FetchContent_Declare(SFML
#     GIT_REPOSITORY https://github.com/SFML/SFML.git
#     GIT_TAG 3.0.0
#     GIT_SHALLOW ON
#     EXCLUDE_FROM_ALL
#     SYSTEM)
# FetchContent_MakeAvailable(SFML)

add_executable(
    dp
    src/main.cc
    src/graph.cc
    src/double_pendulum.cc
)
 target_compile_features(dp PRIVATE cxx_std_17)
# target_link_libraries(main PRIVATE SFML::Graphics)
