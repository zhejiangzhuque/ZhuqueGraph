#ifndef TURAN
#define TURAN

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"

class turan{
private:
    biGraph * g;
    int p, q;

    const int maxD = 1100000;
    const int maxD2 = 10;
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

    LinearSet candL, candR;

public:
    turan(const std::string & filePath, const std::string & outFilePath, 
        int p = 3, int q = 3):p(p), q(q) {
        g = new biGraph(filePath);
        printf("load graph\n");fflush(stdout);

        computeC();

        candL.resize(g->n1);
        candR.resize(g->n2);
    }

    void sample(uint64_t T);
};

#endif