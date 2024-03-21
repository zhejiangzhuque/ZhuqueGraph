#ifndef RAWEDGEPIVOT2_H
#define RAWEDGEPIVOT2_H
#include <string>
#include <algorithm>
#include <cmath>
#include <random>
#include <unordered_set>
#include <chrono>

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"
#include "dinic.h"

class rawEdgePivot2 {
private:
    const int maxD = 1100000;
    const int maxD2 = 20;

    biGraph * g;
    uint32_t p, q;

    double ** C, *bf3;
    void computeC() {
        C = new double*[maxD];
        bf3 = new double[maxD2 * maxD];
        for(int i = 0; i < maxD; i++) {
            C[i] = bf3 + i * maxD2;
        }
        C[0][0] = 1;
        C[1][0] = 1;
        C[1][1] = 1;
        for(int i = 2; i < maxD; i++) {
            C[i][0] = 1;
            if(i < maxD2) C[i][i] = 1;
            for(int j = 1; j < i && j < maxD2; j++) {
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
            }
        }
    }

    double ans;
    std::vector<double> * localAns = nullptr;
    double maxLocalCount;

    Dinic * dinic = nullptr;
    uint32_t cliqueId = 0;
private:
    //pivot
    struct treePath {
        int p1, h1, p2, h2;
    };
    std::vector<std::vector<double>> ansAll;
    LinearSet candL, candR;
    std::vector< std::vector<uint32_t> > tmpNodesL, tmpNodesR;
    std::vector<uint32_t> PL, PR, HL, HR;


private:
    void pivotCount(int l, int pL, int pR, treePath t);

    void pivotForDensestSubgraph(int l, int pL, int pR, treePath t);
    void listAllPminus1QBicliquesL(int i, int l, treePath t, uint32_t idxl);
    void listAllPminus1QBicliquesR(int i, int r, treePath t, uint32_t idxl);
    void listAllPQminus1BicliquesL(int i, int l, treePath t, uint32_t idxr);
    void listAllPQminus1BicliquesR(int i, int r, treePath t, uint32_t idxr);
    std::vector<uint32_t> PLtmp, PRtmp;
    LinearSet comNeiL, comNeiR;
    //comNeiL is the vertices on the L side.
    uint32_t updateComNeiL(uint32_t idxl, uint32_t v) {
        uint32_t newIdx = 0;
        for(uint32_t i = g->pV[v]; i < g->pV[v + 1]; i++) {
            uint32_t u = g->e2[i];
            if(comNeiL.idx(u) < idxl) comNeiL.changeTo(u, newIdx++);
        }
        return newIdx;
    };
    uint32_t updateComNeiR(uint32_t idxr, uint32_t u) {
        uint32_t newIdx = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            if(comNeiR.idx(v) < idxr) comNeiR.changeTo(v, newIdx++);
        }
        return newIdx;
    };
    // std::vector<uint32_t> idxComNeiL, idxComNeiR;
public:
    ~rawEdgePivot2() {
        delete [] bf3;
        delete [] C;
    }
    rawEdgePivot2(const std::string & filePath, const std::string & outFilePath) {
        computeC();
        g = new biGraph(filePath);
        // g->coreReduction(2, 2);
        printf("load graph\n");fflush(stdout);
    }
    
    rawEdgePivot2(uint32_t p, uint32_t q):p(p), q(q) {
        computeC();
    }

    void init(biGraph * g, std::vector<double> * localAns, uint32_t p, uint32_t q) {
        this->g = g;
        this->localAns = localAns;
        this->p = p;
        this->q = q;
        ans = 0.0;
        maxLocalCount = -1.0;

        for(uint32_t i = 0; i < g->n1 + g->n2; i++) (*localAns)[i] = 0.0;

        if(candL.size() < g->n1) candL.resize(g->n1);
        if(candR.size() < g->n2) candR.resize(g->n2);
        if(tmpNodesL.size() == 0) tmpNodesL.resize(std::max(g->maxDu, g->maxDv) + 5);
        if(tmpNodesR.size() == 0) tmpNodesR.resize(std::max(g->maxDu, g->maxDv) + 5);

        // ansAll.resize(std::max(g->maxDu, g->maxDv) + 5);

        if(PL.size() == 0) {
            PL.resize(g->maxDv);
            HL.resize(g->maxDv);
            PR.resize(g->maxDu);
            HR.resize(g->maxDu);
        }

        PL.clear();
        PR.clear();
        HL.clear();
        HR.clear();
    }
    
    void localCount();
    double totalCount() { return ans; }
    double getMaxLocalCount() {
        if(maxLocalCount >= 0.0) return maxLocalCount;
        for(auto w : (*localAns)) maxLocalCount = std::max(maxLocalCount, w);
        return maxLocalCount;
    }

    void addEdgesType2and3(Dinic * dinic);
};

#endif