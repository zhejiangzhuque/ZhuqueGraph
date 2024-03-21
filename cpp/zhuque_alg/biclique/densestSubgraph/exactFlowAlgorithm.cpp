#include "exactFlowAlgorithm.h"
#include "dinic.h"
#include "rawEdgePivot2.h"

#include <cassert>

void exactFlow::run() {
// g->print();
    //get local count
    rawEdgePivot2 * counter = new rawEdgePivot2(p, q);
    std::vector<double> localBiclique(g->n1 + g->n2);
    counter->init(g, &localBiclique, p, q);
    counter->localCount();
    double totalCnt = counter->totalCount();

    printf("total count:%.0f\n", totalCnt);

    uint32_t noZeroN = 0;
    for(uint32_t u = 0; u < g->n1 + g->n2; u++) {
        if(localBiclique[u] > 0) noZeroN++;
    }
    printf("noZeroCnt:%u\n", noZeroN);

    //nodes:|V| + 2 + totalCnt
    //edges:|V|*2 + totalCnt*(p+q)*2 
    uint32_t n = g->n1 + g->n2 + 2 + totalCnt*(p + q);
    uint32_t m = (g->n1 + g->n2) * 2;//type1 & 4 edges
    m += totalCnt * (p + q); //upper bound of type 2 edges
    m += totalCnt * (p + q) * (p + q - 1);//upper bound of type 3 edges

    uint32_t s = n - 2, t = n - 1;
    printf("n:%u %u\n", n, g->n1 +g->n2);
    printf("m:%u\n", m);

    Dinic dinic(n, m, s, t);

    //type 1 edge. from source s to v, with weight Cpq
    for(uint32_t u = 0; u < g->n1 + g->n2; u++) if(localBiclique[u] > 0) dinic.addEdge(s, u, localBiclique[u]);
    uint32_t tp1Edges = dinic.getEdges();
    printf("type1 edges, %u, %u toatl %u\n", g->n1 + g->n2, (g->n1 + g->n2) * 2, tp1Edges);
    
    //type 2 edge. from v to biclique, with weight 1
    //type 3 edge. from biclique to v, with weight INF
    counter->addEdgesType2and3(&dinic);
    delete counter;
    
    uint32_t nowEdges = dinic.getEdges();
    printf("type23 edges should %.0f, %.0f, toatl %u\n", 
        totalCnt * (p + q)*2, tp1Edges + totalCnt * (p + q) * 2 * 2, nowEdges);
    std::vector<double> capacity(nowEdges);
    dinic.copyCapacity(capacity);

    auto floatEqual = [](double a, double b, double eps = 1e-8) { return std::abs(a - b) < eps; };
// {
//     double mid = 5.0 / 6.0;
//     for(uint32_t u = 0; u < g->n1 + g->n2; u++) {
//         dinic.addEdge(u, t, (p + q) * mid);
//     }
// // dinic.tmp();
//     double flow = dinic.maxFlow();

//     printf("flow %f, mid %f, flow-(p+q)*cpq %f\n", flow, mid, flow - (p + q)*totalCnt);
// // dinic.tmp();

// // std::vector<bool> vis(n, false);
// // dinic.print(n, vis);
// // dinic.updateEdgesFrom(nowEdges, (p + q) * mid, capacity);
// // dinic.getCut(n, vis);
// }

    double mid, ans = -1;
    uint32_t ansSize = 0;
    double l = std::max(initialMaxDensity, totalCnt / (g->n1+g->n2));
//     double l = 0.0;
    double r = totalCnt / (p + q);
    bool first = true;
    while(!floatEqual(l, r) && l < r) {
        mid = (l + r) / 2;

        if(first) {
            //type 4 edge. from v to sink t.
            for(uint32_t u = 0; u < g->n1 + g->n2; u++) {
                dinic.addEdge(u, t, (p + q) * mid);
            }
// assert(dinic.getEdges() - nowEdges == (g->n1 + g->n2) * 2);
// dinic.printEdges();
            first = false;
        }
        else {
            dinic.updateEdgesFrom(nowEdges, (p + q) * mid, capacity);
// dinic.printEdges();
        }
        double flow = dinic.maxFlow();
        uint32_t sSize = dinic.SSize(g->n1 + g->n2);
printf("flow %f, mid %f, flow-(p+q)*cpq %f, sSize %u\n", 
    flow, mid, flow - (p + q)*totalCnt, sSize);fflush(stdout);
        
        if(sSize == 0) {
            r = mid;
printf("r = mid, l r: %f %f\n", l, r);
        }
        else {
            l = mid;
            ans = mid;
            ansSize = sSize;
printf("l = mid, l r: %f %f\n", l, r);
        }

        // if(flow < (p + q) * totalCnt) {
        //     // if(floatEqual(flow, (p + q) * totalCnt, 1e-5)) ans = mid;
        //     l = mid;
        // }
        // else r = mid;
    }
    
    dinic.updateEdgesFrom(nowEdges, (p + q) * ans, capacity);
    dinic.maxFlow();
    dinic.print(g->n1 + g->n2);
    printf("l:%f\n", l);
    printf("r:%f\n", r);
    printf("ans:%f\n", ans);
    printf("ansSize:%u\n", ansSize);

// std::vector<uint32_t> nodes(g->n1 + g->n2);
// nodes.clear();
// std::vector<bool> vis(g->n2, false);
// dinic.getS(g->n1 + g->n2, nodes);
// uint32_t edges = 0;
// uint32_t n1 = 0, n2 = 0;
// for(auto u : nodes) {
//     if(u < g->n1) n1++;
//     else {
//         n2++; 
//         vis[u - g->n1] = true;
//     }
// }
// for(auto u : nodes) {
//     if(u >= g->n1) break;
//     for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//         uint32_t v = g->e1[i];
//         if(vis[v]) edges++;
//     }
// }

// printf("%u %u %u\n", n1, n2, edges);
// for(auto u : nodes) {
//     if(u >= g->n1) break;
//     for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//         uint32_t v = g->e1[i];
//         if(vis[v]) printf("%u %u\n", u, v);
//     }
// }

    
    // }

    // delete counter;
}