cmake_minimum_required(VERSION 3.6)

SET(CMAKE_CXX_COMPILER "g++")
set( CMAKE_C_COMPILER "gcc" )

find_package(pybind11 REQUIRED)
include_directories(biClique)
include_directories(biGraph)
include_directories(tools)
include_directories(densestSubgraph)

add_subdirectory(VLDB17)

set(SOURCE_FILE "runner/run.cpp" "biClique/colorPathSpecificPQ.cpp" "biClique/bcAndPath.cpp" "biClique/colorPath.cpp"
        "biClique/colorPathAV5.cpp" "biClique/colorPathPequalsQ.cpp" "biClique/colorPathAV5PequalsQ.cpp"
         "biClique/rawEdgePivot.cpp" "biClique/BCListPlusPlus.cpp" "biClique/bcAndPathV5.cpp"
        "biClique/pivotAndPath.cpp" "biClique/edgePivotSpecificPQ.cpp" "biClique/pivotAndPathSamplev5.cpp"
        "densestSubgraph/exactFlowAlgorithm.cpp" "densestSubgraph/rawEdgePivot2.cpp" "biClique/turan.cpp" densestSubgraph/pqBicliqeDensest.cpp
        "biClique/pivotAndPathPequalsQ.cpp" "biClique/pivotAndPathSamplev5PequalsQ.cpp" "biClique/BK.cpp")
pybind11_add_module(biclique ${SOURCE_FILE})

target_compile_options(biclique PUBLIC -std=c++11)
target_compile_options(biclique PUBLIC -march=native)
target_compile_options(biclique PUBLIC -mavx)
target_compile_options(biclique PUBLIC -O3)



