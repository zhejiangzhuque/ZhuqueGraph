#ifndef BCE_H
#define BCE_H

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSetThree.hpp"
#include <string>


class BCE {
private:
    biGraph * g;
    uint32_t ls, rs, lrs[2];
    std::string od;
    
    uint64_t totalMaximalCount = 0;

private:
    linearSetThree S[2];
    std::vector<std::vector<uint32_t>> ws;
    std::vector<uint32_t> realR[2];
    uint32_t okDegreeSide = 0;
    std::vector<uint32_t> deg;
    uint32_t realRSize[2];
    void bbranch(uint32_t);
    void bbranch2(uint32_t deep, uint32_t x[2], uint32_t c[2], uint32_t r[2]);

public:
    BCE(const std::string & fPath, int mode, const std::string & order, uint32_t ls, uint32_t rs) {
        g = new biGraph(fPath, mode, order);
        printf("load graph\n");fflush(stdout);  
        od = order; 

        this->ls = ls;
        this->rs = rs;
        lrs[0] = ls;
        lrs[1] = rs;
        
        S[0].resize(g->n[0]);
        S[1].resize(g->n[1]);
        deg.resize(std::max(g->n[0], g->n[1]));
        
        ws.resize(g->maxDu + g->maxDv);
        realRSize[0] = realRSize[1] = 0;
    }

    void run();
};

#endif