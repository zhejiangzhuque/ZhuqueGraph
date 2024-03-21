#include "biGraph.hpp"
#include "getArgs.hpp"
#include "BCListPlusPlus.h"

#include <string>
#include <chrono>
#include <vector>
#include <utility>
#include <random>
#include <queue>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <fstream>

using ui = uint32_t;
using P = std::pair<ui, ui>;
using Pdu = std::pair<double, ui>;

// struct Edge {
//     ui v;
//     double p;

//     bool operator < (const Edge & e) const { return v < e.v; }
//     bool operator == (const Edge & e) const { return v == e.v; }
//     bool operator > (const Edge & e) const { return v > e.v; }
// };

double ** C, *bf3;
void computeC() {
    const int maxD = 1000000;
    int maxPQ = 10 + 1;
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

double bclistRunner(ui p, ui q, ui * n, std::vector<ui> *pIdx, std::vector<Edge> *pEdge) {
    BCListPlusPlus * runner = new BCListPlusPlus(p, q, pIdx, pEdge);
    double ans = runner->exactCount();
    delete runner;
    return ans;
}

void solve(biGraph * g, int p, int q, ui sampleEdges, double alpha, bool onePQ) {
    std::cout << p << ' ' << q << ' ' << sampleEdges << std::endl;

    auto t1 = std::chrono::steady_clock::now();

    std::vector<P> edges(sampleEdges);
    std::vector<double> w(sampleEdges);
    std::vector<double> r(sampleEdges);
    std::vector<double> por(sampleEdges);
    std::priority_queue<Pdu, std::vector<Pdu>, std::greater<Pdu>> que;
    std::vector<std::unordered_set<ui>> sg[2];
    sg[0].resize(g->n[0]);
    sg[1].resize(g->n[1]);
    double z = 0.0;

    std::vector<ui> wButtefly(g->m);
    std::vector<uint32_t> sum(g->n[0]);
    std::vector<uint32_t> tmp(g->n[0]);
    std::vector<ui> vis(g->n[1]);
    uint32_t l = 0;
    for(ui u = 0; u < g->n[0]; u++) {
        for(ui i = g->p[0][u]; i < g->p[0][u + 1]; i++) {
            ui v = g->e[0][i];
            vis[v] = i + 1;
        }
        for(ui i = g->p[0][u]; i < g->p[0][u + 1]; i++) {
            ui v = g->e[0][i];

            if(g->p[1][v + 1] > 0) 
            for(ui j = g->p[1][v + 1] - 1; j >= 0; j--) {
                ui w = g->e[1][j];
                if(w == u) break;

                if(sum[w] == 0) tmp[l++] = w;
                sum[w]++;
            }

            for(ui j = 0; j < l; j++) {
                ui w = tmp[j];
                sum[w] = 0;

                ui rj = g->p[0][w + 1] - 1;
                if(g->p[0][w + 1] > 0) {
                    while(g->e[0][rj] != v) --rj;
                }
                if(g->p[0][w + 1] > 0)
                for(ui k = g->p[0][w + 1] - 1; k > rj; k--) {
                    ui z = g->e[0][k];

                    if(vis[z]) {
                        wButtefly[i]++;
                        wButtefly[vis[z] - 1]++;
                        wButtefly[rj]++;
                        wButtefly[k]++;
                    }
                }
            }
            l = 0;
        }
        for(ui i = g->p[0][u]; i < g->p[0][u + 1]; i++) {
            ui v = g->e[0][i];
            vis[v] = 0;
        }
    }
auto tt = std::chrono::steady_clock::now();
auto durationtt = std::chrono::duration_cast<std::chrono::milliseconds>(tt - t1);
std::cout << "localButterflyTime:" << durationtt.count() << "ms" << std::endl;
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> uiDistribution(alpha, 1);
double maxUi = 0;
    for(ui u = 0; u < g->n[0]; u++) {
        for(ui i = g->p[0][u]; i < g->p[0][u + 1]; i++) {
            ui v = g->e[0][i];

            double rdu = uiDistribution(generator);
            // double rdu = 1.0;
            // if(rdu < 1e-12) rdu = 1.0;

            if(i < sampleEdges) {
                // w[i] = sg[0][u].size() + sg[1][v].size();
                // w[i] = 1;
                // w[i] = g->p[0][u+1]-g->p[0][u] + g->p[1][v+1]-g->p[1][v];
                w[i] = wButtefly[i];
                // w[i] = w[i]*w[i];
// std::cout<<i<<' '<<u<<' '<<v<<std::endl;
                // w[i] = 0;
                // if(!sg[0][u].empty() && !sg[1][v].empty())
                // for(auto a : sg[0][u]) {
                //     for(auto b : sg[1][v]) {
                //         if(sg[1][a].find(b) != sg[1][a].end()) w[i]++;
                //     }
                // }
maxUi = std::max(maxUi, w[i]);
                r[i] = w[i] / rdu;

                sg[0][u].insert(v);
                sg[1][v].insert(u);
                edges[i] = P(u, v);
                que.push(P(r[i], i));

                // z = std::max(z, r[i]);
            }
            else {

                // ui wi = sg[0][u].size() + sg[1][v].size();
                // ui  wi = 1;
                // ui wi = g->p[0][u+1]-g->p[0][u] + g->p[1][v+1]-g->p[1][v];
                ui wi = wButtefly[i];
                // ui wi = 0;
                // if(!sg[0][u].empty() && !sg[1][v].empty())
                // for(auto a : sg[0][u]) {
                //     for(auto b : sg[1][v]) {
                //         if(sg[1][a].find(b) != sg[1][a].end()) wi++;
                //     }
                // }
                double ri = wi / rdu;
maxUi = std::max(maxUi, w[i]);

                Pdu head = que.top();

                if(ri <= head.first) {
                    z = std::max(z, ri);
                    continue;
                }
                z = std::max(z, head.first);

                ui j = head.second;
                que.pop();
                edges[j] = P(u, v);
                sg[0][u].erase(v);
                sg[1][v].erase(u);
                w[j] = wi;
                r[j] = ri;
                que.push(P(r[j], j));
            }
        }

    }

    std::cout << "z:" << z << std::endl;
    std::cout << "que_size:" << que.size() << std::endl;
    // std::cout << "maxUi:" << maxUi << std::endl;
    for(ui i = 0; i < sampleEdges; i++) por[i] = std::min(1.0, w[i] / z);
    double sumP = 0;
    for(ui i = 0; i < sampleEdges; i++) sumP += por[i];
    std::cout << "average probability:" << sumP / sampleEdges << std::endl;


    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "sampleTime:" << duration.count() << "ms" << std::endl;

    double ans = 0.0;

    //build graph
    std::map<ui, ui> reId[2];
    ui newId[2] = {0, 0};
    for(ui i = 0; i < edges.size(); i++) {
        ui u = edges[i].first;
        ui v = edges[i].second;
        if(reId[0].find(u) == reId[0].end()) reId[0][u] = newId[0]++;
        if(reId[1].find(v) == reId[1].end()) reId[1][v] = newId[1]++;
    }
    for(ui i = 0; i < edges.size(); i++) {
        edges[i].first = reId[0][edges[i].first];
        edges[i].second = reId[1][edges[i].second];
    }
    std::cout << "reId:" << newId[0] << ' ' << newId[1] << std::endl;

    std::vector<ui> pIdx[2];
    std::vector<Edge> pEdge[2];
    pIdx[0].resize(newId[0] + 1);
    pEdge[0].resize(edges.size() * 2);
    pIdx[1].resize(newId[1] + 1);
    pEdge[1].resize(edges.size() * 2);
    for(ui i = 0; i < edges.size(); i++) {
        ui u = edges[i].first;
        ui v = edges[i].second;
        pIdx[0][u + 1]++;
        pIdx[1][v + 1]++;
    }
    for(ui u = 0; u < newId[0]; u++) pIdx[0][u + 1] += pIdx[0][u];
    for(ui u = 0; u < newId[1]; u++) pIdx[1][u + 1] += pIdx[1][u];

    for(ui i = 0; i < edges.size(); i++) {
        ui u = edges[i].first;
        ui v = edges[i].second;
        pEdge[0][pIdx[0][u]].v = v;
        pEdge[0][pIdx[0][u]++].p = por[i];
        pEdge[1][pIdx[1][v]].v = u;
        pEdge[1][pIdx[1][v]++].p = por[i];
    }

    for(ui u = newId[0]; u > 0; u--) pIdx[0][u] = pIdx[0][u - 1];
    for(ui u = newId[1]; u > 0; u--) pIdx[1][u] = pIdx[1][u - 1];
    pIdx[0][0] = pIdx[1][0] = 0;
    std::cout << "here" << std::endl;

    for(ui u = 0; u < newId[0]; u++) {
// std::cout << u << ' ' << pIdx[0][u + 1] - pIdx[0][u] << std::endl;
        std::sort(pEdge[0].begin() + pIdx[0][u], pEdge[0].begin() + pIdx[0][u + 1]);
    }
    std::cout << "here" << std::endl;
    for(ui u = 0; u < newId[1]; u++) {
        std::sort(pEdge[1].begin() + pIdx[1][u], pEdge[1].begin() + pIdx[1][u + 1]);
    }
    std::cout << "sorted" << std::endl;

    auto t3 = std::chrono::steady_clock::now();
    auto duration32 = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
    std::cout << "buildGraph:" << duration32.count() << "ms" << std::endl;


#define OUTPUTSAMPLEDGRAPH
#ifdef OUTPUTSAMPLEDGRAPH
std::fstream fout("dblp_sampled_2_5_1e6.txt", std::ios::out);
fout << newId[0] << ' ' << newId[1] << ' ' << edges.size() << '\n';
for(ui u = 0; u < newId[0]; u++) {
    for(ui i = pIdx[0][u]; i < pIdx[0][u + 1]; i++) {
        fout << u << ' ' << pEdge[0][i].v << '\n';
    }
}
return;
#endif

    if(!onePQ)
    for(ui p = 2; p < 10; p++) {
        // for(ui q = 2; q < 10; q++) {
        for(ui q = p; q < p + 1; q++) {
            auto t4 = std::chrono::steady_clock::now();

            double cnt = bclistRunner(p, q, newId, pIdx, pEdge);
            std::cout << p << '-' << q << ':' << cnt << std::endl;

            auto t5 = std::chrono::steady_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4);
            std::cout << "pqtime," <<  p << '-' << q << ' ' << duration.count() << "ms" << std::endl;
        }
    }
    else {
        auto t4 = std::chrono::steady_clock::now();

        double cnt = bclistRunner(p, q, newId, pIdx, pEdge);
        std::cout << p << '-' << q << ':' << cnt << std::endl;

        auto t5 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4);
        std::cout << "pqtime" << p << '-' << q << ' ' << duration.count() << "ms" << std::endl;
    }
}

int main(int argc, char *argv[])
{
    argsController ac(argc, argv);

    if(!ac.exist("-f")) {
        printf("file path: -f\n");
        return 0;
    }
    std::string fPath = ac["-f"];

    int graphMode = (ac.exist("noUVM") == false);

    std::string order = "core";
    if(ac.exist("-d")) order = ac["-d"];

    int p = 4, q = 4;
    if(ac.exist("-p")) p = std::stoi(ac["-p"]);
    if(ac.exist("-q")) q = std::stoi(ac["-q"]);

    int sampleEdges = 1000000;
    if(ac.exist("-e")) sampleEdges = std::stoi(ac["-e"]);

    bool onePQ = true;
    if(ac.exist("-all")) {
        onePQ = false;
        std::cout << "all" << std::endl;
    }

    double alpha = 0.0;
    if(ac.exist("-a")) {
        alpha = std::stof(ac["-a"]);
    }
    std::cout << "alpha:" << alpha << std::endl;

    std::cout << " -f " << fPath << 
                " -d " << order << 
                " -p " << p << " -q " << q << 
                " -e" << sampleEdges << std::endl;

    biGraph * g = new biGraph(fPath, graphMode, order);

    std::cout << "g: " << g->n[0] << ' ' << g->n[1] << ' ' << g->m << std::endl;

    computeC();

    auto t1 = std::chrono::steady_clock::now();

    solve(g, p, q, sampleEdges, alpha, onePQ);

    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "totaltime:" << duration.count() << "ms" << std::endl;

    delete g;
    
    return 0;
}