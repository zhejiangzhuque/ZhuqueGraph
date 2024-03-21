#include "Graph.h"
#include "Utility.h"
#include <pybind11/pybind11.h>
using namespace std;

int algorithm(string filename, string attribute, string algorithm, string threshold, string order, string delta){
    Graph* graph = new Graph(filename.c_str());
    graph->ReadGraph(filename.c_str(), attribute.c_str(), atoi(threshold.c_str()), "r");//read and store graph
    if(strcmp(algorithm.c_str(),"FairClique")==0) graph->FindClique(order.c_str());
    if(strcmp(algorithm.c_str(),"StrongClique")==0 ) graph->FindStrongClique(order.c_str());
    //based on alternative selection
    if(strcmp(algorithm.c_str(),"RelativeStrong") ==0 ) graph->FindRelatedFairClique_S(order.c_str(), atoi(delta.c_str()));
    //base on weak fair clique enumeration
    if(strcmp(algorithm.c_str(),"RelativeWeak") ==0 ) graph->FindRelatedFairClique_W("colorful", atoi(delta.c_str()));
    if(strcmp(algorithm.c_str(),"Baseline")==0 ) graph->Baseline(order.c_str(),atoi(delta.c_str())); //order strong weak
    printf("algorithm done\n");
    printf("max_mem=%lld\n", graph->get_max());
    return 0;
}


PYBIND11_MODULE(fairnessclique, m){
m.def("algorithm", &algorithm,
pybind11::arg("filename")="./fairnessclique/example.txt",
pybind11::arg("attribute")="./fairnessclique/example2.txt",pybind11::arg("algorithm")="FairClique",pybind11::arg("threshold")="2",
pybind11::arg("order")="sorted",pybind11::arg("delta")="1");
}