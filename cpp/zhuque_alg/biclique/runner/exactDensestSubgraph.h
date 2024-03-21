#include "../tools/getArgs.hpp"
#include "../densestSubgraph/exactFlowAlgorithm.h"
#include <pybind11/pybind11.h>
#include <cassert>
#include <string>
#include <iostream>
#include <ctime>
#include <chrono>
using std::string;
//int p = 4; int q = 4; double initialDensity = 0.0;
int exdensest(string filePath,string outFilePath,int p,int q,double initialDensity){
    if(filePath==""){
        std::cout << "filePath" << std::endl;
        exit(-1);
    }

    std::cout << filePath << ' ' << outFilePath << std::endl;
    std::cout << "p:" << p << std::endl;
    std::cout << "q:" << q << std::endl;
    std::cout << "initialDensity:" << initialDensity << std::endl;
    exactFlow * runner = new exactFlow(filePath, outFilePath, p, q, initialDensity);
    auto t1 = std::chrono::steady_clock::now();
    runner->run();
    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "time:" << duration.count() << "ms" << std::endl;
    delete runner;
    return 0;
}
//PYBIND11_MODULE(biclique, m1){
//m1.def("exdensest", &exdensest,
//pybind11::arg("filePath")="",pybind11::arg("outFilePath")="",
//pybind11::arg("p")=4,pybind11::arg("q")=4,pybind11::arg("initialDensity")=0.0);
//}