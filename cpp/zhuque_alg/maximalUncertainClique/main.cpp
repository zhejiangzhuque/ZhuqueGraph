#include "maximalClique.h"
#include <string>
#include <pybind11/pybind11.h>

int algorithm(string filename, int alg, int k, double eta){
    Algorithm *mc = new maximalClique();
    mc->read_graph(filename.c_str());
    mc->setparemeters(eta, alg, k);
    mc->run();
    printf(" Alg=%d\n" , alg);
    delete mc;
    return 0;
}

PYBIND11_MODULE(maximalCliques, m){
m.def("algorithm", &algorithm,
pybind11::arg("filename")="./example/example.txt",
pybind11::arg("alg")=2,pybind11::arg("k")=3,pybind11::arg("eta")=0.5);
}