#ifndef BCANDPATH
#define BCANDPATH

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"

#include <chrono>

class bcAndPath {
private:
    const int maxD = 100000;
    int minPQ;
    int p, q;
    LinearSet candL, candR;
    biGraph * g;

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
    double ans = 0;

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
    ~bcAndPath() {
        delete [] bf3;
        delete [] C;
    }
    bcAndPath(const std::string & filePath, const std::string & outFilePath, int p_, int q_) {
        auto t1 = std::chrono::steady_clock::now();
        p = p_;
        q = q_;
        minPQ = std::min(p, q);
        computeC();

        g = new biGraph(filePath);
        
        // g->coreReduction(p, q);
        // auto t2 = std::chrono::steady_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        // std::cout << "core reduction time:" << duration.count() << "ms" << std::endl;
        // std::cout << "n1 n2 m:" << g->n1 << ' '<< g->n2 << ' ' << g->m<< std::endl;
        // t2 = t1;

        if(costEstimate()) {
            g->swapUV();
            std::swap(p, q);
            printf("swap\n");
        }

        // t2 = std::chrono::steady_clock::now();
        // duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        // printf("costEstimate %lld ms\n", (long long)duration.count());
        printf("load graph\n");
        fflush(stdout);
    }
    
    void counting(uint64_t T, double realV = 0.0, uint32_t bar = 1000);

    bool costEstimate();
    void exactCount(std::vector<uint32_t> & nodes);
    void layerBasedListing(int l, int pH, int pS);

    void sample(std::vector<uint32_t> & nodes, uint64_t T);

    void countingV5(uint64_t T, double realV = 0.0, uint32_t bar = 1000);
    void sampleV5(std::vector<uint32_t> & nodes, uint64_t T);
};

#endif