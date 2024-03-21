cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_C_COMPILER "gcc")

find_package(pybind11 REQUIRED)
set(SOURCE_FILES main.cc methods.cc graph.cc)
pybind11_add_module(landmark ${SOURCE_FILES})

target_compile_options(landmark PUBLIC -O2)
target_compile_options(landmark PUBLIC -g)
target_compile_options(landmark PUBLIC -lgflags_nothreads)
