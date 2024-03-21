cmake_minimum_required(VERSION 3.6)

SET(CMAKE_CXX_COMPILER "g++")
set(CMAKE_C_COMPILER "gcc" )

include_directories(BCE)
include_directories(biGraph)
include_directories(tools)

find_package(pybind11 REQUIRED)
set(SOURCE_FILE "run/run.cpp" "BCE/BCE.cpp")
pybind11_add_module(mbiclique ${SOURCE_FILE})

target_compile_options(mbiclique PUBLIC -std=c++11)
target_compile_options(mbiclique PUBLIC -march=native)
target_compile_options(mbiclique PUBLIC -mavx)
target_compile_options(mbiclique PUBLIC -O3)

