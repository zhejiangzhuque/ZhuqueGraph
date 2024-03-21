#ifndef FASTEDGEPIVOT_H
#define FASTEDGEPIVOT_H
#include <string>
#include <algorithm>
#include <cmath>
#include <random>
#include <unordered_set>
#include <chrono>

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"

class fastEdgePivot {
private:
    const int maxD = 50000;
    const int maxD2 = 500;

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

    double Cnm(uint32_t n, uint32_t m) {
        if(n > maxD || m > maxD / 50) {
            double ans = 1;
            for(int i = 1; i <= m; i++) {
                ans *= (n-i+1) / i;
            }
            return ans;
        }
        else return C[n][m];
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
    ~fastEdgePivot() {
        delete [] bf3;
        delete [] C;
    }
    fastEdgePivot(const std::string & filePath, const std::string & outFilePath) {
        computeC();
        g = new biGraph(filePath);
        g->createHashTables();
        printf("load graph\n");fflush(stdout);
    }

    void exactCountMaximalPivotFast();
};

#endif