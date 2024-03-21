#include "BCListPlusPlus.h"
#include <chrono>

// #define DEBUG
double BCListPlusPlus::exactCount() {
    collect2HopNeighbors();
    printf("collect 2\n"); fflush(stdout);

    S.resize(n[1]);

    uint32_t maxDegree = 0;
    for(uint32_t u = 0; u < n[0]; u++) {
        maxDegree = std::max(maxDegree, (uint32_t)H.lists[u].size());
    }


    tmpNodes.resize(p + 1);
    tmpNodes[0].resize(n[0]);
    for(int i = 1; i <= p; i++) {
        tmpNodes[i].resize(maxDegree);
    }
    H.nodes.resize(n[0]);
    // H.nodes.resize(p + 1);
    // for(int i = 0; i <= p; i++) {
    //     H.nodes[i].resize(n[0]);
    // }

    H.d.resize(p + 1);
    for(int i = 0; i <= p; i++) {
        H.d[i].resize(n[0]);
    }
    for(uint32_t u = 0; u < n[0]; u++) {
        H.d[0][u] = H.lists[u].size();
    }

    printf("maxDegree of H %u\n", maxDegree);
    fflush(stdout);
    
    ansTmp.resize(n[1]);

    ans = 0;

    layerBasedListing(0, n[0], n[1]);

    return ans;
}

void BCListPlusPlus::listAllSubsetsOfS(int l, int i, int pS, double w) {
    if(l == q) {
        ans += 1.0 / w;
        return;
    }

    if(l + pS - i < q) return;

    for(int j = i; j < pS; j++) {
        listAllSubsetsOfS(l + 1, j + 1, pS, w * ansTmp[j]);
    }
}

void BCListPlusPlus::layerBasedListing(int l, int pH, int pS) {
#ifdef DEBUG
printf("l:%d pH:%d pS:%d ans:%.2f\n", l, pH, pS, ans);
#endif
    if(l == p) {
        // ans += C[pS][q];
        for(int i = 0; i < pS; i++) ansTmp[i] = 1.0;
        for(int i = 0; i < pH; i++) {
            uint32_t u = H.nodes[i];
            for(ui j = pIdx[0][u]; j < pIdx[0][u + 1]; j++) {
                ui v = pEdge[0][j].v;
                if(S.idx(v) < pS) {
                    ansTmp[S.idx(v)] *= pEdge[0][j].p;
                }
            }
        }
        
        listAllSubsetsOfS(0, 0, pS, 1.0);
        
        return;
    }

    H.nodes.copy(tmpNodes[l].data(), pH);
// for(int i = 0; i < pH; i++) {
//     printf("%u ", tmpNodes[l][i]);
// }printf("\n");fflush(stdout);

// auto t1 = std::chrono::steady_clock::now();
    for(int j = 0; j < pH; j++) {
        uint32_t u = tmpNodes[l][j];
        // uint32_t u = H.nodes[l][j];
        if(H.lists[u].size() < uint32_t(p - l - 1)) {
            continue;
        }

        int pSS = 0;
        if(l == 0) {
            for(uint32_t i = pIdx[0][u]; i < pIdx[0][u + 1]; i++) {
                S.changeTo(pEdge[0][i].v, pSS++);
            }
        }
        else {
            for(uint32_t i = pIdx[0][u]; i < pIdx[0][u + 1]; i++) {
                if(S.idx(pEdge[0][i].v) < (uint32_t)pS) {
                    S.changeTo(pEdge[0][i].v, pSS++);
                }
            }
        }

        if(pSS < q) continue;
        
        int pHH = 0;

        for(int i = 0; i < H.d[l][u]; i++) {
            auto v = H.lists[u][i];
            if(H.nodes.idx(v) < (uint32_t)pH) {
                H.nodes.changeTo(v, pHH++);
            }
            // if(H.nodes[l].idx(v) < (uint32_t)pH) {
            //     H.nodes[l + 1].changeTo(v, pHH++);
            // }
        }

        if(l + 1 < p)
        for(int i = 0; i < pHH; i++) {
            // uint32_t u = H.nodes[l + 1][i];
            uint32_t u = H.nodes[i];
            int d = H.d[l][u];
            for(int k = 0; k < d; k++) {
                auto v = H.lists[u][k];
                if(H.nodes.idx(v) >= pHH) {
                    std::swap(H.lists[u][k], H.lists[u][--d]);
                    --k;
                }
                // if(H.nodes[l + 1].idx(v) >= pHH) {
                //     std::swap(H.lists[u][k], H.lists[u][--d]);
                //     --k;
                // }
            }
            H.d[l + 1][u] = d;
        }

        layerBasedListing(l + 1, pHH, pSS);
    }
}

void BCListPlusPlus::collect2HopNeighbors() {
    H.lists.resize(n[0]);

    std::vector<uint32_t> sum(n[0]);
    std::vector<uint32_t> tmp(n[0]);
    uint32_t l = 0;

double twotwo = 0;
    for(uint32_t u = 0; u < n[0]; u++) {
        for(uint32_t i = pIdx[0][u]; i < pIdx[0][u + 1]; i++) {
            uint32_t v = pEdge[0][i].v;
            for(uint32_t j = pIdx[1][v]; j < pIdx[1][v + 1]; j++) {
                uint32_t w = pEdge[1][j].v;
                if(w > u) {
                    if(sum[w] == 0) tmp[l++] = w;
                    sum[w]++;
                }
            }
        }

        for(uint32_t i = 0; i < l; i++) {
            uint32_t w = tmp[i];
            if(sum[w] >= q) {
twotwo += C[sum[w]][2];
                H.lists[u].push_back(w);
            }
            sum[w] = 0;
        }
        l = 0;
    }

if(p == q && p == 2)
    printf("2-2 clique %.0f\n", twotwo);
}

