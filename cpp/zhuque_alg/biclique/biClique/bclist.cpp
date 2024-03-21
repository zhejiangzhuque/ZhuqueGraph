#include "BCListPlusPlus.h"

bool BCListPlusPlus::costEstimate(int p, int q) {
    // std::default_random_engine e(10007);
    // std::uniform_int_distribution<unsigned> chooseU(0, g->n1 - 1);
    // std::uniform_int_distribution<unsigned> chooseV(0, g->n2 - 1);

    std::vector<uint32_t> sum(std::max(g->n1, g->n2));
    std::vector<uint32_t> tmp(std::max(g->n1, g->n2));
    int l = 0;

    uint32_t rd = 10;
    // uint32_t t = rd;
    uint32_t Du = 0;
    // while(t--) {
        // uint32_t u = chooseU(e);
    for(uint32_t u = 0; u < g->n1; u += rd) {
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];

                if(w != u) {
                    if(sum[w] == 0) tmp[l++] = w;
                    sum[w]++;
                }
            }
        }
    }
    for(int i = 0; i < l; i++) {
        if(sum[tmp[i]] >= (uint32_t)q) {
            Du++;
        }
        sum[tmp[i]] = 0;
    }
    l = 0;
    
    // }
    
    // t = rd;
    uint32_t Dv = 0;
    // while(t--) {
        // uint32_t v = chooseV(e);
    for(uint32_t v = 0; v < g->n2; v += rd) {
        for(uint32_t i = g->pV[v]; i < g->pV[v + 1]; i++) {
            uint32_t u = g->e2[i];
    // if(u >= g->n1) {
    //     printf("n2 %u, %u, %u, i %u\n", g->n2, v, u, i);fflush(stdout);
    // }
            for(uint32_t j = g->pU[u]; j < g->pU[u + 1]; j++) {
                uint32_t w = g->e1[j];

                if(w != v) {
                    if(sum[w] == 0) tmp[l++] = w;
                    sum[w]++;
                }
            }
        }
    }

    for(int i = 0; i < l; i++) {
        if(sum[tmp[i]] >= (uint32_t)p) {
            Dv++;
        }
        sum[tmp[i]] = 0;
    }
    l = 0;
    
    // }

    double costU = pow(1.0 * Du / (g->n1 / 10), p - 2) * Du * g->n1 * g->maxDu;
    double costV = pow(1.0 * Dv / (g->n2 / 10), p - 2) * Dv * g->n2 * g->maxDv;
    printf("cost:%.2f %.2f, du %u, dv %u\n", costU, costV, Du, Dv);

    return costU > costV;
}

void BCListPlusPlus::exactCount(int p, int q) {
    auto t1 = std::chrono::steady_clock::now();

    if(costEstimate(p, q)) {
        g->swapUV();
        std::swap(p, q);
        printf("swap\n");
    }

    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    printf("costEstimate %lld ms\n", (long long)duration.count());
#ifdef DEBUG
    g->print();
#endif
    collect2HopNeighbors(p, q);
    printf("collect 2\n"); fflush(stdout);
#ifdef DEBUG
    H.print();
#endif
    S.resize(g->n2);

    uint32_t maxDegree = 0;
    for(uint32_t u = 0; u < g->n1; u++) {
        maxDegree = std::max(maxDegree, (uint32_t)H.lists[u].size());
    }
    tmpNodes.resize(p);
    tmpNodes[0].resize(g->n1);
    for(int i = 1; i < p; i++) {
        tmpNodes[i].resize(maxDegree);
    }
    printf("maxDegree of H %u\n", maxDegree);
    fflush(stdout);
    
#ifdef DEBUG
double tmppp = 0;
for(uint32_t v = 0; v < g->n2; v++) {
    if(g->deg2(v) >= 3) {
        tmppp += C[g->deg2(v)][3];
    }
}
printf("3-1 %.2f\n", tmppp);
#endif
    

    ans = 0;

    layerBasedListing(0, g->n1, g->n2, p, q);

    printf("ans %.2f\n", ans);
    fflush(stdout);
}

void BCListPlusPlus::layerBasedListing(int l, int pH, int pS, int p, int q) {
#ifdef DEBUG
printf("l:%d pH:%d pS:%d\n", l, pH, pS);
#endif
    if(l == p) {
        ans += C[pS][q];
        return;
    }

    H.nodes.copy(tmpNodes[l].data(), pH);
// for(int i = 0; i < pH; i++) {
//     printf("%u ", tmpNodes[l][i]);
// }printf("\n");fflush(stdout);

    for(int j = 0; j < pH; j++) {
        uint32_t u = tmpNodes[l][j];
        if(H.lists[u].size() < uint32_t(p - l - 1)) {
            continue;
        }

        int pSS = 0;
        if(l == 0) {
            for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
                S.changeTo(g->e1[i], pSS++);
            }
        }
        else {
            for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
                if(S.idx(g->e1[i]) < (uint32_t)pS) {
                    S.changeTo(g->e1[i], pSS++);
                }
            }
        }

        if(pSS < q) continue;
        
        int pHH = 0;
        for(auto v: H.lists[u]) {
            if(H.nodes.idx(v) < (uint32_t)pH) {
                H.nodes.changeTo(v, pHH++);
            }
        }

        layerBasedListing(l + 1, pHH, pSS, p, q);
    }
}

void BCListPlusPlus::collect2HopNeighbors(uint32_t p, uint32_t q) {
    H.lists.resize(g->n1);
    H.nodes.resize(g->n1);

    std::vector<uint32_t> sum(g->n1);
    std::vector<uint32_t> tmp(g->n1);
    uint32_t l = 0;

int twotwo = 0;
    for(uint32_t u = 0; u < g->n1; u++) {
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) {
                    if(sum[w] == 0) tmp[l++] = w;
                    sum[w]++;
                }
            }
        }

        for(uint32_t i = 0; i < l; i++) {
            uint32_t w = tmp[i];
            if(sum[w] >= q) {
twotwo += C[sum[w]][q];
                H.lists[u].push_back(w);
            }
            sum[w] = 0;
        }
        l = 0;
    }
printf("2-2 clique %d\n", twotwo);
}
