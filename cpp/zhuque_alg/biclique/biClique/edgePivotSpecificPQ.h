#ifndef EDGEPIVOTSPECIFICPQ_H
#define EDGEPIVOTSPECIFICPQ_H

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"
#include <iostream>

class edgePivotSpecificPQ {
private:
    const int maxD = 500000;
    biGraph * g;
    int p, q;

    double ** C, *bf3;
    void computeC() {
        int maxPQ = std::max(p, q) + 1;

        C = new double*[maxD];
        bf3 = new double[maxD * maxPQ];
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
private:
    //pivot
    struct treePath {
        int p1, h1, p2, h2;
    };
    std::vector<std::vector<double>> ansAll;
    LinearSet candL, candR;
    std::vector< std::vector<uint32_t> > tmpNodesL, tmpNodesR;

private:
    void pivotCountFast(int l, int pL, int pR, treePath t);

public:
    ~edgePivotSpecificPQ() {
        delete [] bf3;
        delete [] C;
    }
    edgePivotSpecificPQ(const std::string & filePath, const std::string & outFilePath, int p_, int q_) {
        p = p_;
        q = q_;

        computeC();
        g = new biGraph(filePath);
        
        // g->coreReduction(p, q);
        // std::cout << "n1 n2 m:" << g->n1 << ' '<< g->n2 << ' ' << g->m<< std::endl;

        g->createHashTables();
        printf("load graph\n");fflush(stdout);

        ansAll.resize(p + 1);
        for(int i = 0; i <= p; i++) {
            ansAll[i].resize(q + 1);
        }
    }

    void exactCountMaximalPivotFast();
    void sparseCounting(uint32_t bar);
};

#endif
