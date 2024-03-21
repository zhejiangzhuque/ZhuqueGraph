#include "../biClique/rawEdgePivot.h"
#include "../tools/getArgs.hpp"
#include "../biClique/BCListPlusPlus.h"
#include "../biClique/BK.h"
#include "../biClique/edgePivotSpecificPQ.h"
#include "../biClique/colorPath.h"
#include "../biClique/pivotAndPath.h"
#include "../biClique/turan.h"
#include "../biClique/bcAndPath.h"
#include "../biClique/colorPathPequalsQ.h"
#include "../biClique/pivotAndPathPequalsQ.h"
#include "../biClique/colorPathSpecificPQ.h"
#include "densestSubgraph.h"
#include "exactDensestSubgraph.h"
#include <pybind11/pybind11.h>
#include <cassert>
#include <string>
#include <iostream>
#include <ctime>
#include <chrono>
using std::string;

int run(string filePath, string outFilePath,int p,int q,string algo,string v5,double r,string bar,int H,uint64_t T,double realV){
    if(filePath==""){
        std::cout << "filePath" << std::endl;
        exit(-1);
    }
    std::cout << filePath << ' ' << outFilePath << std::endl;
    if(algo=="m") {
        BK * counter = new BK(filePath, outFilePath);
        auto t1 = std::chrono::steady_clock::now();
        counter->exactCountMaximal(p, q);
        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }
    else if(algo=="pm") {
        rawEdgePivot * counter = new rawEdgePivot(filePath, outFilePath);

        auto t1 = std::chrono::steady_clock::now();

        if(v5=="v5") counter->exactCountMaximalPivotV2();
        else counter->exactCountMaximalPivot();

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }
    else if(algo=="apm") {//approximate pm
        rawEdgePivot * counter = new rawEdgePivot(filePath, outFilePath);
        auto t1 = std::chrono::steady_clock::now();
        counter->approximateCountMaximalPivot(r);
        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }
    // else if(aC->exist("-fpm")) {//no use delete it
    //     fastEdgePivot * counter = new fastEdgePivot(filePath, outFilePath);
    //     auto t1 = std::chrono::steady_clock::now();

    //     counter->exactCountMaximalPivotFast();

    //     auto t2 = std::chrono::steady_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    //     std::cout << "time:" << duration.count() << "ms" << std::endl;
    //     delete counter;
    // }
    else if(algo=="fpmPQ") {
        edgePivotSpecificPQ * counter = new edgePivotSpecificPQ(filePath, outFilePath, p, q);
        auto t1 = std::chrono::steady_clock::now();

        if(bar!="") {
            uint32_t b = atoi(bar.c_str());
            counter->sparseCounting(b);
        }
        else counter->exactCountMaximalPivotFast();

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }
    else if(algo=="cp") {//pure color path
        colorPath * counter = new colorPath(filePath, outFilePath, H, p, q);
        auto t1 = std::chrono::steady_clock::now();
        // counter->approximateCounting();
        // counter->testSubgraphSize();

        // counter->approximateCountingAll(T);

        if(v5=="v5")
            counter->approximateCountingAllVersion5(T);
        else counter->approximateCountingAllVersion2(T);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }
    else if(algo=="cppq") {//pure color path, specific pq
        colorPathSpecificPQ * counter = new colorPathSpecificPQ(filePath, outFilePath,  p, q);
        auto t1 = std::chrono::steady_clock::now();
        // counter->approximateCounting();
        // counter->testSubgraphSize();
        // counter->approximateCountingAll(T);

        // if(aC->exist("-v5"))
        //     counter->approximateCountingAllVersion5(T);
        // else
        counter->approximateCountingAllVersion2(T);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }
    else if(algo=="peqq") {//p equals q sampling, p=q<H
        colorPathPequalsQ * counter = new colorPathPequalsQ(filePath, outFilePath, H, p, q);
        auto t1 = std::chrono::steady_clock::now();
        // counter->approximateCounting();
        // counter->testSubgraphSize();
        // counter->approximateCountingAll(T);

        if(v5=="v5")
            counter->approximateCountingAllVersion5(T);
        else counter->approximateCountingAllVersion2(T);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }
    else if(algo=="pp") {//exact+sampling p<H,q<H
        pivotAndPath * counter = new pivotAndPath(filePath, outFilePath, H);

        auto t1 = std::chrono::steady_clock::now();

        if(v5=="v5") counter->countingV5(T);
        else counter->counting(T);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }
    else if(algo=="pppeqq") {//exact+sampling p=q<H
        pivotAndPathPequalsQ * counter = new pivotAndPathPequalsQ(filePath, outFilePath, H);

        auto t1 = std::chrono::steady_clock::now();
        if(v5=="v5") counter->countingV5(T);
        else
        counter->counting(T);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }
    else if(algo=="tu") {
        turan * counter = new turan(filePath, outFilePath, p, q);
        auto t1 = std::chrono::steady_clock::now();
        counter->sample(T);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }
    else if(algo=="bcpath") {
        // int H = 11;
        // if(aC->exist("-H")) H = atoi(aC->get("-H").c_str());
        uint32_t b = atoi(bar.c_str());

        bcAndPath * counter = new bcAndPath(filePath, outFilePath, p, q);

        auto t1 = std::chrono::steady_clock::now();

        if(v5=="v5") counter->countingV5(T, realV, b);
        else counter->counting(T, realV, b);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }
    else {//BCList++
        BCListPlusPlus * counter = new BCListPlusPlus(filePath, outFilePath, p, q);
        auto t1 = std::chrono::steady_clock::now();

        counter->exactCount();
        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }
    return 0;
}

PYBIND11_MODULE(biclique, m){
    m.def("run", &run,
        pybind11::arg("filePath")="",pybind11::arg("outFilePath")="",
        pybind11::arg("p")=4,pybind11::arg("q")=4,pybind11::arg("algo")="",pybind11::arg("v5")="",pybind11::arg("r")=0.1,
        pybind11::arg("bar")="1000",pybind11::arg("H")=11,pybind11::arg("T")=100000,pybind11::arg("realV")=1.0);
    m.def("exdensest", &exdensest,
        pybind11::arg("filePath")="",pybind11::arg("outFilePath")="",
        pybind11::arg("p")=4,pybind11::arg("q")=4,pybind11::arg("initialDensity")=0.0);
    m.def("densest", &densest,
        pybind11::arg("filePath")="",pybind11::arg("outFilePath")="",
        pybind11::arg("p")=4,pybind11::arg("q")=4);
}