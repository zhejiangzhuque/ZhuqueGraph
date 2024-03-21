#include "../tools/getArgs.hpp"
#include "../BCE/BCE.h"
#include <string>
#include <chrono>
#include <iostream>
#include <pybind11/pybind11.h>
using namespace std;

//uint32_t ls = 1, uint32_t rs = 1;
int run(string fPath,string order,uint32_t ls, uint32_t rs,int noUVM) {
    if(fPath=="") {
        printf("file path: -f\n");
        return 0;
    }

    if(order=="") {
        printf("graph order: -d(core, two)\n");
        return 0;
    }
    std::cout << fPath << ' ' << order << std::endl;
    std::cout << "l " << ls << std::endl;
    std::cout << "r " << rs << std::endl;

    if(noUVM==1) {
        BCE bce(fPath, 0, order, ls, rs);
        auto t1 = std::chrono::steady_clock::now();
        bce.run();
        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
    }
    else {
        BCE bce(fPath, 1, order, ls, rs);
        auto t1 = std::chrono::steady_clock::now();
        bce.run();
        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
    }
    return 0;
}

PYBIND11_MODULE(mbiclique, m){
    m.def("run", &run,
          pybind11::arg("fPath")="",pybind11::arg("order")="core",
          pybind11::arg("ls")=1,pybind11::arg("rs")=1,pybind11::arg("noUVM")=1);
}