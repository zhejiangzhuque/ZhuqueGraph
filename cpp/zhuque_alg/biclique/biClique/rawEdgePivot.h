#ifndef RAWEDGEPIVOT_H
#define RAWEDGEPIVOT_H
#include <string>
#include <algorithm>
#include <cmath>
#include <random>
#include <unordered_set>
#include <chrono>

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"

// #define PQWEDGE
// #define COEI_BASELINE

class rawEdgePivot {
private:
    const int maxD = 1100000;
    const int maxD2 = 100;

    biGraph * g;

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
private:
    //pivot
    struct treePath {
        int p1, h1, p2, h2;
    };
    std::vector<std::vector<double>> ansAll;
    LinearSet candL, candR;
    std::vector< std::vector<uint32_t> > tmpNodesL, tmpNodesR;

#ifdef PQWEDGE
std::vector<std::vector<double>> ansWedge;
// std::vector<uint32_t> Pl, Pr, Hl, Hr;
double sumDegPl, sumDegPr, sumDegHl, sumDegHr;
#endif

private:
    void pivotCount(int l, int pL, int pR, treePath t);
    
public:
    ~rawEdgePivot() {
        delete [] bf3;
        delete [] C;
    }
    rawEdgePivot(const std::string & filePath, const std::string & outFilePath) {
        computeC();
        g = new biGraph(filePath);
        // g->coreReduction(2, 2);
        printf("load graph\n");fflush(stdout);
    }
    
    void exactCountMaximalPivot();
    void approximateCountMaximalPivot(double rate);

    void exactCountMaximalPivotV2();
};

#endif