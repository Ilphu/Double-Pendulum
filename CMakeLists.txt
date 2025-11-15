cmake_minimum_required(VERSION 3.28)
project(double-pendulum LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(SDL2 REQUIRED)

add_executable(
    dp
    src/main.cc
    src/graph.cc
    src/double_pendulum.cc
)
target_link_libraries(dp ${SDL2_LIBRARIES})
target_include_directories(dp PRIVATE ${SDL2_INCLUDE_DIRS})
target_compile_features(dp PRIVATE cxx_std_17)
