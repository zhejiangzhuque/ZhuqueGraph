cmake_minimum_required(VERSION 3.6)

SET(CMAKE_CXX_COMPILER "g++")
set( CMAKE_C_COMPILER "gcc" )

find_package(pybind11 REQUIRED)
include_directories(biGraph)
include_directories(plex)
include_directories(tools)

set(SOURCE_FILE "run/run.cpp" "plex/biplexv2.cpp")
pybind11_add_module(multiPivot ${SOURCE_FILE})

target_compile_options(multiPivot PUBLIC -std=c++11)
target_compile_options(multiPivot PUBLIC -march=native)
target_compile_options(multiPivot PUBLIC -mavx)
target_compile_options(multiPivot PUBLIC -O3)

