cmake_minimum_required(VERSION 3.6)

SET(CMAKE_CXX_COMPILER "g++")

find_package(pybind11 REQUIRED)
set(SOURCE_FILE "gps.cpp")

pybind11_add_module(gpsbiclique ${SOURCE_FILE})

target_compile_options(gpsbiclique PUBLIC -std=c++11)
target_compile_options(gpsbiclique PUBLIC -march=native)
target_compile_options(gpsbiclique PUBLIC -mavx)
target_compile_options(gpsbiclique PUBLIC -O3)