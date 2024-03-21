#ifndef BCLISTPLUSPLUS_H
#define BCLISTPLUSPLUS_H

#include <string>
#include <algorithm>
#include <cmath>
#include <random>
#include <unordered_set>
#include <chrono>

#include "biGraph.hpp"
#include "linearSet.hpp"

struct Edge {
    ui v;
    double p;

    bool operator < (const Edge & e) const { return v < e.v; }
    bool operator == (const Edge & e) const { return v == e.v; }
    bool operator > (const Edge & e) const { return v > e.v; }
};

class BCListPlusPlus {
private:
    const int maxD = 100000;
    int p, q;

    // biGraph * g;
    std::vector<ui> * pIdx;
    std::vector<Edge> * pEdge;
    ui n[2];
    

    double ** C, *bf3;
    void computeC() {
        int maxPQ = std::max(p, q) + 1;
        C = new double*[maxD];
        bf3 = new double[maxPQ * maxD];
        for(int i = 0; i < maxD; i++) {
            C[i] = bf3 + i * maxPQ;
        }
        C[0][0] = 1;
        C[1][0] = 1;
        C[1][1] = 1;
        for(int i = 2; i < maxD; i++) {
            C[i][0] = 1;
            if(i < maxPQ) C[i][i] = 1;
            for(int j = 1; j < i && j < maxPQ; j++) {
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
            }
        }
    }
    double ans;
    std::vector<double> ansTmp;

private:
    struct twoHopGraph {
        std::vector<std::vector<uint32_t>> lists;
        std::vector<std::vector<int>> d;
        // std::vector<LinearSet> nodes;
        LinearSet nodes;
        std::vector<bool> contained;
        std::vector<std::vector<uint32_t>> containNodes;

        void print() {
            for(uint32_t i = 0; i < lists.size(); i++) {
                printf("%u:", i);
                for(auto v:lists[i]) printf("%u ", v);
                printf("\n");
            }
        }

        bool contain(uint32_t u, uint32_t v) {
            return std::binary_search(containNodes[u].begin(), containNodes[u].end(), v);
        }
    } H;
    void collect2HopNeighbors();
    void collect2HopNeighborsContained();
    void collect2HopNeighborsContainedCDAG();

    LinearSet S;
    std::vector< std::vector<uint32_t> > tmpNodes;


public:
    ~BCListPlusPlus() {
        delete [] bf3;
        delete [] C;
    }
    BCListPlusPlus(ui p, ui q, std::vector<ui> *pIdx, std::vector<Edge> *pEdge):p(p),q(q),pIdx(pIdx),pEdge(pEdge) {
        computeC();
        n[0] = pIdx[0].size() - 1;
        n[1] = pIdx[1].size() - 1;
    }
    
    double exactCount();
    void layerBasedListing(int l, int pH, int pS);
    void listAllSubsetsOfS(int l, int i, int pS, double w);

    bool costEstimate();
};

#endif