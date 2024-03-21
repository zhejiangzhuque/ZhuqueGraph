#include "../tools/getArgs.hpp"
// #include "../plex/biplex.h"
#include "../plex/biplexv2.h"
#include <string>
#include <chrono>
#include <iostream>
#include <pybind11/pybind11.h>
using namespace std;
//uint32_t ls = 1, uint32_t rs = 1; int k = 2;uint64_t outPutT = INT64_MAX;
int run(string fPath,string order,uint32_t ls, uint32_t rs,int graphMode,int k,uint64_t outPutT) {
    if(fPath=="") {
        printf("file path: -f\n");
        return 0;
    }

    if(order=="") {
        printf("graph order: -d(core, two)\n");
        return 0;
    }

    std::cout << "file path: " << fPath << " "
              << "k: " << k << " "
              << "order: " << order << " "
              << "l : " << ls << " "
              << "r : " << rs << " "
              << "mode: " << graphMode << std::endl;

    // if(ac.exist("-v2")) {
        std::cout << "-v2" << std::endl;
        biplexv2 bcev2(fPath, graphMode, k, order, ls, rs, outPutT);
        auto t1 = std::chrono::steady_clock::now();
        bcev2.run();
        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
    // }
    // else {
//         biplex bce(fPath, graphMode, k);
//
//         auto t1 = std::chrono::steady_clock::now();
//
//         bce.run();
//
//         auto t2 = std::chrono::steady_clock::now();
//         auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
//         std::cout << "time:" << duration.count() << "ms" << std::endl;
    // }
    return 0;
}


PYBIND11_MODULE(multiPivot, m){
    m.def("run", &run,
        pybind11::arg("fPath")="",pybind11::arg("order")="core",
        pybind11::arg("ls")=1,pybind11::arg("rs")=1,pybind11::arg("graphMode")=0,pybind11::arg("k")=2,
        pybind11::arg("outPutT")=INT64_MAX);
}