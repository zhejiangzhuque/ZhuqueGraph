#include "../tools/getArgs.hpp"
#include "../densestSubgraph/pqBicliqeDensest.h"
#include <pybind11/pybind11.h>
#include <cassert>
#include <string>
#include <iostream>
#include <ctime>
#include <chrono>
using std::string;
//int p = 4; int q = 4;
int densest(string filePath,string outFilePath,int p ,int q){
    if(filePath==""){
        std::cout << "filePath" << std::endl;
        exit(-1);
    }
    std::cout << filePath << ' ' << outFilePath << std::endl;

    bicliqueDensest * runner = new bicliqueDensest(filePath, outFilePath, p, q);
    auto t1 = std::chrono::steady_clock::now();
    runner->run();
    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "time:" << duration.count() << "ms" << std::endl;
    delete runner;
    return 0;
}

//PYBIND11_MODULE(biclique, m){
//m.def("densest", &densest,
//pybind11::arg("filePath")="",pybind11::arg("outFilePath")="",
//pybind11::arg("p")=4,pybind11::arg("q")=4);
//}