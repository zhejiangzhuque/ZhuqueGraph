#include "pqBicliqeDensest.h"

// #include "../tools/lazyQueue.hpp"

#include <fstream>
#include <unordered_set>
#include <assert.h>

void bicliqueDensest::run() {
    struct queNode
    {
        uint32_t v;
        double cliqueCnt;

        bool operator >(const queNode & t) const {
             return cliqueCnt < t.cliqueCnt;
        }

        bool operator <(const queNode & t) const {
            return cliqueCnt > t.cliqueCnt;
        }
    };
    // lazyQueue<queNode> heap;
    std::priority_queue<queNode, std::vector<queNode> > heap;

// queNode a = queNode{1, 1.0};
// queNode b = queNode{2, 2.0};
// heap.push(a);
// heap.push(b);
// heap.push(queNode{0, 3.0});
// heap.push(queNode{0, 2.0});
// heap.push(queNode{0, 1.0});
// int ttt = 5;
// while(ttt--) {
// printf("test %u %f\n", heap.top().v, heap.top().cliqueCnt);heap.pop();
// }
// printf("%d\n", int(a < b));
// return;

// g->print();
    uint32_t n = g->n1 + g->n2;
    
    std::vector<double> localBiclique(n), localBicliqueSg(n);
    printf("n %u\n", n);
    counter->init(g, &localBiclique, p, q);
    counter->localCount();

    printf("totalCount %u-%u:%.0f\n", p, q, counter->totalCount());
// for(uint32_t u = 0; u < g->n1; u++) printf("U_%u:%.0f\n", u, localBiclique[u]);
// for(uint32_t v = 0; v < g->n2; v++) printf("V_%u:%.0f\n", v, localBiclique[v+g->n1]);
// fflush(stdout);
    
    // std::vector<uint32_t> ids(n);
    // std::vector<uint32_t> keys(n);
    std::vector<uint32_t> cliqueCoreNumber(n);
    std::vector<double> leftCliques(n);
    std::vector<uint32_t> labelsL(g->n1);
    std::vector<uint32_t> labelsR(g->n2);

    //subgraph
    std::vector<uint32_t> stk(n);
    uint32_t lstk = 0;
    std::vector<bool> inStk(n);
    std::vector<uint32_t> reId(n);
    
    uint32_t labL = 0, labR = 0;
    std::vector<bool> vis(n);
    std::fill(vis.begin(), vis.end(), false);

    for(uint32_t i = 0; i < n; i++) {
        if(localBiclique[i] > 0) heap.push(queNode{i, localBiclique[i]});
        else vis[i] = true;
    }
    printf("initial heap size %u\n", heap.size());
    uint32_t mostLoops = heap.size();

    leftCliques[0] = counter->totalCount();
    double maxDensity = counter->totalCount() / mostLoops;
    uint32_t maxDensityIndex = 0;

    biGraph * sg = new biGraph();
    sg->pU.resize(g->n1 + 1);
    sg->pV.resize(g->n2 + 1);
    sg->e1.resize(g->m);
    sg->e2.resize(g->m);
    sg->pU[0] = sg->pV[0] = 0;

    auto buildSubgraph = [&](uint32_t u) {
        if(u < g->n1) {
assert(lstk == 0);
            for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
                uint32_t vv = g->e1[i];
                uint32_t v = vv + g->n1;
                if(vis[v]) continue;
                for(uint32_t j = g->pV[vv]; j < g->pV[vv + 1]; j++) {
                    uint32_t w = g->e2[j];
                    if(vis[w] || inStk[w]) continue;
                    stk[lstk++] = w;
                    inStk[w] = true;
                }
            }

            std::sort(stk.begin(), stk.begin() + lstk);
            sg->n1 = lstk;
            sg->m = 0;
            for(uint32_t i = 0; i < lstk; i++) {
                reId[stk[i]] = i;
            }

            for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
                uint32_t v = g->e1[i] + g->n1;
                if(vis[v]) continue;
                stk[lstk] = v;
                reId[v] = lstk++;
                inStk[v] = true;
            }
            sg->n2 = lstk - sg->n1;
            for(uint32_t i = 1; i <= sg->n1; i++) sg->pU[i] = 0;
            for(uint32_t i = 1; i <= sg->n2; i++) sg->pV[i] = 0;

            for(uint32_t i = sg->n1; i < lstk; i++) {
                uint32_t vv = stk[i] - g->n1;
                for(uint32_t j = g->pV[vv]; j < g->pV[vv + 1]; j++) {
                    uint32_t w = g->e2[j];
                    if(!inStk[w]) continue;
                    sg->pU[reId[w] + 1]++;
                    sg->pV[i - sg->n1 + 1]++;
                }
            }
            for(uint32_t i = 1; i <= sg->n1; i++) sg->pU[i] += sg->pU[i - 1];
            for(uint32_t i = 1; i <= sg->n2; i++) sg->pV[i] += sg->pV[i - 1];
            sg->m = sg->pU[sg->n1];
            for(uint32_t i = sg->n1; i < lstk; i++) {
                uint32_t vv = stk[i] - g->n1;
                for(uint32_t j = g->pV[vv]; j < g->pV[vv + 1]; j++) {
                    uint32_t w = g->e2[j];
                    if(!inStk[w]) continue;
                    sg->e1[sg->pU[reId[w]]++] = i - sg->n1;
                    sg->e2[sg->pV[i - sg->n1]++] = reId[w];
                }
            }
            for(uint32_t i = sg->n1; i >= 1; i--) sg->pU[i] = sg->pU[i - 1];
            sg->pU[0] = 0;
            for(uint32_t i = sg->n2; i >= 1; i--) sg->pV[i] = sg->pV[i - 1];
            sg->pV[0] = 0;
            for(uint32_t i = 0; i < sg->n2; i++) {
                std::sort(sg->e2.begin() + sg->pV[i],
                    sg->e2.begin() + sg->pV[i + 1]);
            }
            
            for(uint32_t i = 0; i < lstk; i++) inStk[stk[i]] = false;
            lstk = 0;
        }
        else {
            uint32_t vv = u - g->n1;
            for(uint32_t i = g->pV[vv]; i < g->pV[vv + 1]; i++) {
                uint32_t w = g->e2[i];
                if(vis[w]) continue;
                stk[lstk++] = w;
                inStk[w] = true;
// printf("push %u, lstk %u\n", w, lstk);fflush(stdout);
            }

            // std::sort(stk.begin(), stk.begin() + lstk);
            sg->n1 = lstk;
            sg->m = 0;
            for(uint32_t i = 0; i < lstk; i++) {
                reId[stk[i]] = i;
            }

            for(uint32_t j = 0; j < sg->n1; j++) {
                uint32_t w = stk[j];

                for(uint32_t i = g->pU[w]; i < g->pU[w + 1]; i++) {
                    uint32_t v = g->e1[i] + g->n1;
                    if(vis[v] || inStk[v]) continue;
                    stk[lstk] = v;
                    reId[v] = lstk++;
                    inStk[v] = true;
                }
            }
            
            sg->n2 = lstk - sg->n1;
            for(uint32_t i = 1; i <= sg->n1; i++) sg->pU[i] = 0;
            for(uint32_t i = 1; i <= sg->n2; i++) sg->pV[i] = 0;

            for(uint32_t j = 0; j < sg->n1; j++) {
                uint32_t w = stk[j];
                for(uint32_t i = g->pU[w]; i < g->pU[w + 1]; i++) {
                    uint32_t v = g->e1[i] + g->n1;
                    if(!inStk[v]) continue;
                    sg->pU[j + 1]++;
// printf("w %u, j %u, lstk %u, %u\n", v, j, lstk, reId[v]);fflush(stdout);
                    sg->pV[reId[v] - sg->n1 + 1]++;
                }
            }
            for(uint32_t i = 1; i <= sg->n1; i++) sg->pU[i] += sg->pU[i - 1];
            for(uint32_t i = 1; i <= sg->n2; i++) sg->pV[i] += sg->pV[i - 1];
            sg->m = sg->pU[sg->n1];
            for(uint32_t j = 0; j < sg->n1; j++) {
                uint32_t w = stk[j];
                for(uint32_t i = g->pU[w]; i < g->pU[w + 1]; i++) {
                    uint32_t v = g->e1[i] + g->n1;
                    if(!inStk[v]) continue;
                    sg->e1[sg->pU[j]++] = reId[v] - sg->n1;
                    sg->e2[sg->pV[reId[v] - sg->n1]++] = j;
                }
            }
            for(uint32_t i = sg->n1; i >= 1; i--) sg->pU[i] = sg->pU[i - 1];
            sg->pU[0] = 0;
            for(uint32_t i = sg->n2; i >= 1; i--) sg->pV[i] = sg->pV[i - 1];
            sg->pV[0] = 0;
            for(uint32_t i = 0; i < sg->n1; i++) {
                std::sort(sg->e1.begin() + sg->pU[i],
                    sg->e1.begin() + sg->pU[i + 1]);
            }
            
            for(uint32_t i = 0; i < lstk; i++) inStk[stk[i]] = false;
            lstk = 0;
        }
    };

// std::unordered_set<uint32_t> st;
// st.reserve(mostLoops);
// for(uint32_t i = 0; i < n; i++) if(!vis[i]) st.insert(i);
    if(heap.size() > 0)
    for(uint32_t i = 0; heap.size() > 0 && i < mostLoops - p - q; i++) {
        uint32_t u;
// printf("    i %u\n", i);
        // if(!heap.pop_min(u, degU)) printf("errorLheap\n");
        queNode h = heap.top(); heap.pop();

        while(vis[h.v] || std::abs(h.cliqueCnt - localBiclique[h.v]) > 1e-8) {
// if(i==57) {
// printf("%u %u tw(%.1f %.1f) %f %d %d\n", i, h.v, localBiclique[h.v],
//      h.cliqueCnt, std::abs(h.cliqueCnt - localBiclique[h.v]), 
//      int(std::abs(h.cliqueCnt - localBiclique[h.v]) > 1e-6), (int)vis[h.v]); fflush(stdout);
// }
            h = heap.top(); heap.pop();
            
        }
// printf("here\n");fflush(stdout);
        u = h.v;
        vis[u] = true;
// st.erase(u);
// for(auto x: st) {
//     if(h.cliqueCnt > localBiclique[x]) {
//         printf("%u %f, %u %f\n", x, localBiclique[x], u, h.cliqueCnt);
//         while(!heap.empty()) {
//             h = heap.top(); heap.pop();
//             if(h.v == x) printf("x:%f\n", h.cliqueCnt);
//         }
//         fflush(stdout);
//     }
//     assert(localBiclique[x] >= h.cliqueCnt);
// }
// printf("    hv %u\n", h.v);
        if(u < g->n1) labelsL[u] = labL++;
        else labelsR[u - g->n1] = labR++;
        cliqueCoreNumber[u] = i;

        if(h.cliqueCnt == 0) {
            printf("i %u, zero cnt\n", i);
            continue;
        }

        buildSubgraph(u);
// if(i == 4) {
//     std::fstream ot("dd.txt", std::ios::out);
//     ot << sg->n1 << ' ' << sg->n2 << ' ' << sg->m << std::endl;
//     for(uint32_t u = 0; u < sg->n1; u++) {
//         for(uint32_t i = sg->pU[u]; i < sg->pU[u + 1]; i++) {
//             ot << u << ' ' << sg->e1[i] << std::endl;
//         }
//     }
//     ot.close();
// }
// sg->print();
// printf("xxxx sgn1, %u sgn2, %u, lc %.0f, hv %u, hcc %.0f\n",
    //  sg->n1, sg->n2, localBiclique[h.v], h.v, h.cliqueCnt);
        if(u < g->n1) counter->init(sg, &localBicliqueSg, p - 1, q);
        else counter->init(sg, &localBicliqueSg, p, q - 1);
        counter->localCount();
// printf("total %.0f\n", counter->totalCount());


        for(uint32_t su = 0; su < sg->n1 + sg->n2; su++) {
            uint32_t reU = stk[su];
// assert(!vis[reU]);
// assert(reU != u);
            localBiclique[reU] -= localBicliqueSg[su];
            if(localBiclique[reU]> 0) {
// printf("larger lc %f lcsg %f , reU %u su %u\n", localBiclique[reU], localBicliqueSg[su], reU, su);
//     fflush(stdout);
                
                heap.push(queNode{reU, localBiclique[reU]});
            }
            else {
// if(localBiclique[reU] < localBicliqueSg[su]) {
//     printf("smaller lc %f lcsg %f , reU %u su %u\n", localBiclique[reU], localBicliqueSg[su], reU, su);
//     fflush(stdout);
// }
                assert(localBiclique[reU] == 0.0);
                mostLoops--;//remove reU
// st.erase(reU);
                vis[reU] = true;
            }
            
        }

        if(i > 0) leftCliques[i] = leftCliques[i - 1] - counter->totalCount();
        else leftCliques[0] -= counter->totalCount();

        // maxDensity = std::max(maxDensity, leftCliques[i] / (n - i - 1));
printf("i %u, u %u, hcc %.0f, totalCnt %.0f, heapsz %u, leftCliques %.0f, maxDNow %f\n", 
            i, u, counter->totalCount(), h.cliqueCnt, heap.size(), leftCliques[i], maxDensity);
// fflush(stdout);
        if(maxDensity < leftCliques[i] / (mostLoops - i - 1)) {
            maxDensity = leftCliques[i] / (mostLoops - i - 1);
            maxDensityIndex = i;
        }

        if(leftCliques[i] == 0 || leftCliques[i] < (p+q) * maxDensity) {
            printf("leftCliques/(p+q)<maxDensityNow: %.0f/(%u+%u)=%f<%f\n", 
                leftCliques[i], p, q, leftCliques[i]/(p+q), maxDensity);
            break;
        }
    }

    printf("maxDensity:%f\n", maxDensity);
    printf("The densest subgraph size:%u\n", mostLoops - maxDensityIndex);

    delete sg;
}