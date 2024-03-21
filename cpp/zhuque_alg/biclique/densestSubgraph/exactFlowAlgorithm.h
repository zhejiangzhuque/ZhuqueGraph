#ifndef EXACTFLOWALGORITHM_H
#define EXACTFLOWALGORITHM_H

#include "../biGraph/biGraph.hpp"
#include "rawEdgePivot2.h"

class exactFlow {
private:
    uint32_t p, q;
    biGraph * g;
    double initialMaxDensity;
    double maxDensity = 0.0;
public:
    exactFlow(const std::string & filePath, const std::string & outFilePath,
        uint32_t p, uint32_t q, double initial):p(p), q(q), initialMaxDensity(initial) {
        g = new biGraph(filePath);
    }
    ~exactFlow() { delete g; }

    void run();
};

#endif