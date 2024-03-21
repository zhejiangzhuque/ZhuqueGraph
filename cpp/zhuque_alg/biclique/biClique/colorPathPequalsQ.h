#ifndef COLORPATHPEQUALSQ
#define COLORPATHPEQUALSQ

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"

class colorPathPequalsQ {
private:
    biGraph * g;
    int p, q;
    LinearSet candL, candR;
    std::vector<double> ansAll;
    // std::vector<std::vector<double> > ansAll;
    int minPQ;

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

public:
    ~colorPathPequalsQ() {
        delete g;
        delete [] C;
        delete [] bf3;
    }
    colorPathPequalsQ(const std::string & filePath, const std::string & outFilePath, 
        int H=11, int p_=0, int q_=0) {
        p = p_;
        q = q_;
        minPQ = H;

        g = new biGraph(filePath);
        printf("load graph\n");fflush(stdout);

        computeC();
    }

    
//差分优化
    void approximateCountingAllVersion2(uint64_t T);
    void approximateCountingAllVersion5(uint64_t T);
};

#endif