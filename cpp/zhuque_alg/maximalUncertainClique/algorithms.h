#pragma once
#include "structs.h"
#include <cstdlib>
#include <ctime>
#include <queue>
#include <algorithm>
#include <cstring>
#include <vector>

using namespace std;
#define pairs pair<int, double>

class Algorithm
{
protected:
	int V, E;
    int maxdeg;
    int *deg;
    pairs **adj;
    pairs *data;
    
    double eta;
    int mincliquesize;
    int algorithm;

    vector<int> core, topcore, colors;
public:
    Algorithm(/* args */);
    virtual ~Algorithm();
    virtual void run() {}

    void setparemeters(double eta, int alg, int k) {
        this->eta = eta; this->algorithm = alg;
		    mincliquesize = k;
        printf("eta=%.2e, alg=%d, mincliquesize=%d\n", eta, alg, k);
    }

    void read_graph(const char *str);
    void scalability(bool randomv, float scal);
    void testprintGraph();
    //void setSortedAdj();

    int core_decompsition(int *nodeset, int nodesize);
    int topKEtaCore(int *nodeset);
    int topKEtaCoreDecompsition(int *nodeset);
    int topKEtatriangle(int *nodeset);
    int coloring(int *nodeset, int nodesize);
    void resetGraph(int *nodeset, int nodesize);
};