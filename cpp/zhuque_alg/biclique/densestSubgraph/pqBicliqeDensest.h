#ifndef PQBICLIQUEDENSEST_H
#define PQBICLIQUEDENSEST_H

#include "../biGraph/biGraph.hpp"
#include "rawEdgePivot2.h"

class bicliqueDensest {
private:
    uint32_t p, q;
    biGraph * g;
    rawEdgePivot2 * counter;

public:
    bicliqueDensest(const std::string & filePath, const std::string & outFilePath,
        uint32_t p, uint32_t q):p(p), q(q) {
        g = new biGraph(filePath);
        counter = new rawEdgePivot2(p, q);
    }
    ~bicliqueDensest() {
        delete g;
        delete counter;
    }

    void run();
};

#endif