#ifndef PIVOTANDPATHPEQUALSQ_H
#define PIVOTANDPATHPEQUALSQ_H

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"

class pivotAndPathPequalsQ {
private:
    biGraph * g;
    int p, q;
    LinearSet candL, candR;
    // std::vector<std::vector<double>> ansAll;
    std::vector<double> ansAll;
    int minPQ = 100;

    double ** C, *bf3;
    void computeC() {
        // int maxPQ = std::max(p, q) + 1;
        int maxPQ = minPQ + 1;
        int maxD = std::max(g->maxDu, g->maxDv) + 1;

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

private:
    //pivot
    struct treePath {
        int p1, h1, p2, h2;
    };
    // std::vector<std::vector<double>> ansAll;
    // LinearSet candL, candR;
    std::vector< std::vector<uint32_t> > tmpNodesL, tmpNodesR;

    void pivotCount(int l, int pL, int pR, treePath t);

private://partition
    std::vector<uint32_t> parts[2];

private:
    void pivot(std::vector<uint32_t> & nodes);
    void sample(std::vector<uint32_t> & nodes, uint64_t T);
    void sample2(std::vector<uint32_t> & nodes, uint64_t T);
    void samplev5(std::vector<uint32_t> & nodes, uint64_t T);

public:
    ~pivotAndPathPequalsQ() {
        delete g;
        delete [] C;
        delete [] bf3;
    }
    pivotAndPathPequalsQ(const std::string & filePath, const std::string & outFilePath, int H) {
        minPQ = H;
        g = new biGraph(filePath);
        printf("load graph\n");fflush(stdout);

        computeC();

        ansAll.resize(minPQ + 1);
        // for(uint32_t i = 0; i <= minPQ; i++) {
        //     ansAll[i].resize(minPQ + 1);
        // }

        candL.resize(g->n1);
        candR.resize(g->n2);
    }

    void counting(uint64_t T);
    void countingV5(uint64_t T);
};

#endif