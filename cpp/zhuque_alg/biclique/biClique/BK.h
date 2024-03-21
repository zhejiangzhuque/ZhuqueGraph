#ifndef BK_H
#define BK_H

#include <string>
#include <algorithm>
#include <cmath>
#include <random>
#include <unordered_set>
#include <chrono>

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"


class BK {
private:
    const int maxD = 10000;

    biGraph * g;

    double ** C, *bf3;
    void computeC() {
        C = new double*[maxD];
        bf3 = new double[maxD / 50 * maxD];
        for(int i = 0; i < maxD; i++) {
            C[i] = bf3 + i * maxD / 50;
        }
        C[0][0] = 1;
        C[1][0] = 1;
        C[1][1] = 1;
        for(int i = 2; i < maxD; i++) {
            C[i][0] = 1;
            if(i < maxD / 50) C[i][i] = 1;
            for(int j = 1; j < i && j < maxD / 50; j++) {
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
            }
        }
    }
    double ans;
private:
    struct twoHopGraph {
        std::vector<std::vector<uint32_t>> lists;
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
    void collect2HopNeighbors(uint32_t p, uint32_t q);
    void collect2HopNeighborsContained(int p, int q);
    void collect2HopNeighborsContainedCDAG(int p, int q);

    LinearSet S;
    LinearSet Q;
    std::vector< std::vector<uint32_t> > tmpNodes;

public:
    ~BK() {
        delete [] bf3;
        delete [] C;
    }
    BK(const std::string & filePath, const std::string & outFilePath) {
        computeC();
        g = new biGraph(filePath);
        printf("load graph\n");fflush(stdout);
    }

    void exactCountMaximal(int p, int q);
    void countFormMaximalClique(int l, int L, int pH, int pS, int pQ, int p, int q);

};

#endif