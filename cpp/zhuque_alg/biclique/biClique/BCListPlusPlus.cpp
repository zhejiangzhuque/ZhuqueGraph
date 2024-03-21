#include "BCListPlusPlus.h"
#include <chrono>

bool BCListPlusPlus::costEstimateRaw() {
    std::default_random_engine e(10007);
    std::uniform_int_distribution<unsigned> chooseU(0, g->n1 - 1);
    std::uniform_int_distribution<unsigned> chooseV(0, g->n2 - 1);

    std::vector<uint32_t> sum(std::max(g->n1, g->n2));
    std::vector<uint32_t> tmp(std::max(g->n1, g->n2));
    int l = 0;

    uint32_t rd = 10;
    // uint32_t t = rd;
    uint32_t Du = 0;
    // while(t--) {
    uint32_t u = chooseU(e);
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
    uint32_t v = chooseV(e);
    // for(uint32_t v = 0; v < g->n2; v += rd) {
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
    // }

    for(int i = 0; i < l; i++) {
        if(sum[tmp[i]] >= (uint32_t)p) {
            Dv++;
        }
        sum[tmp[i]] = 0;
    }
    l = 0;
    
    // }

    double costU = pow(1.0 * Du / (g->n1 / 10), p - 2) * Du * g->maxDu;
    double costV = pow(1.0 * Dv / (g->n2 / 10), q - 2) * Dv * g->maxDv;
    printf("cost:%.2f %.2f, du %u, dv %u\n", costU, costV, Du, Dv);

    return costU > costV;
}

bool BCListPlusPlus::costEstimate() {
    // std::default_random_engine e(10007);
    // std::uniform_int_distribution<unsigned> chooseU(0, g->n1 - 1);
    // std::uniform_int_distribution<unsigned> chooseV(0, g->n2 - 1);

    std::vector<uint32_t> sum(std::max(g->n1, g->n2));
    std::vector<uint32_t> tmp(std::max(g->n1, g->n2));
    int l = 0;

    uint32_t rd = 100;
    // uint32_t t = rd;
    uint32_t Du = 0, maxDu, rdu = 0;
    double sumDu = 0.0;
    // while(t--) {
        // uint32_t u = chooseU(e);
    for(uint32_t u = 0; u < g->n1; u += rd) {
        rdu++;
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

        for(int i = 0; i < l; i++) {
            if(sum[tmp[i]] >= (uint32_t)q) {
                Du++;
            }
            sum[tmp[i]] = 0;
        }
        l = 0;

        maxDu = std::max(maxDu, Du);
        sumDu += Du;
    }
    
    uint32_t Dv = 0, maxDv, rdv = 0;
    double sumDv = 0.0;
    // while(t--) {
        // uint32_t v = chooseV(e);
    for(uint32_t v = 0; v < g->n2; v += rd) {
        rdv++;
        for(uint32_t i = g->pV[v]; i < g->pV[v + 1]; i++) {
            uint32_t u = g->e2[i];
    // if(u >= g->n1) {
    //     printf("n2 %u, %u, %u, i %u\n", g->n2, v, u, i);fflush(stdout);
    // }
            for(uint32_t j = g->pU[u]; j < g->pU[u + 1]; j++) {
                uint32_t w = g->e1[j];

                if(w > v) {
                    if(sum[w] == 0) tmp[l++] = w;
                    sum[w]++;
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

        maxDv = std::max(maxDv, Dv);
        sumDv += Dv;
    }
    
    sumDu = sumDu / rdu * g->n1;
    sumDv = sumDv / rdv * g->n2;

    double avgDu = std::max(2.0, sumDu / g->n1);
    double avgDv = std::max(2.0, sumDv / g->n2);

    double costU = pow(avgDu, p - 2) * sumDu;
    double costV = pow(avgDv, q - 2) * sumDv;
    printf("cost:%.2f %.2f, du %u, dv %u\n", costU, costV, Du, Dv);

    return costU > costV;
}

void BCListPlusPlus::exactCount() {
    
#ifdef DEBUG
    g->print();
#endif
    collect2HopNeighbors();
    printf("collect 2\n"); fflush(stdout);
#ifdef DEBUG
    H.print();
#endif
    S.resize(g->n2);

    uint32_t maxDegree = 0;
    for(uint32_t u = 0; u < g->n1; u++) {
        maxDegree = std::max(maxDegree, (uint32_t)H.lists[u].size());
    }


    tmpNodes.resize(p + 1);
    tmpNodes[0].resize(g->n1);
    for(int i = 1; i <= p; i++) {
        tmpNodes[i].resize(maxDegree);
    }
    H.nodes.resize(g->n1);
    // H.nodes.resize(p + 1);
    // for(int i = 0; i <= p; i++) {
    //     H.nodes[i].resize(g->n1);
    // }

    H.d.resize(p + 1);
    for(int i = 0; i <= p; i++) {
        H.d[i].resize(g->n1);
    }
    for(uint32_t u = 0; u < g->n1; u++) {
        H.d[0][u] = H.lists[u].size();
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

    layerBasedListing(0, g->n1, g->n2);

    printf("ans %.2f\n", ans);
    fflush(stdout);
}

void BCListPlusPlus::layerBasedListing(int l, int pH, int pS) {
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

auto t1 = std::chrono::steady_clock::now();
    for(int j = 0; j < pH; j++) {
// if(l == 0 && j % 200 == 0) {
//     auto t2 = std::chrono::steady_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
//     printf("%u, %lld\n", j/10000, (long long)duration.count());fflush(stdout);
//     t1 = t2;
// }

        uint32_t u = tmpNodes[l][j];
        // uint32_t u = H.nodes[l][j];
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
    H.lists.resize(g->n1);

    std::vector<uint32_t> sum(g->n1);
    std::vector<uint32_t> tmp(g->n1);
    uint32_t l = 0;

double twotwo = 0;
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
printf("2-2 clique %.0f\n", twotwo);
}


void BCListPlusPlus::collect2HopNeighborsContainedCDAG() {
    H.contained.resize(g->n1);
    H.lists.resize(g->n1);
    H.nodes.resize(g->n1);
    H.containNodes.resize(g->n1);

    std::vector<uint32_t> sum(g->n1);
    std::vector<uint32_t> tmp(g->n1);
    uint32_t l = 0;

// int twotwo = 0;
    for(uint32_t u = g->n1; u >= 0; u--) {
        H.contained[u] = false;

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w >= u) {//w < u也可能dw == du
                    if(sum[w] == 0) tmp[l++] = w;
                    sum[w]++;
                }
            }
        }

        for(uint32_t i = 0; i < l; i++) {
            uint32_t w = tmp[i];
            
            if(sum[w] >= (uint32_t)p) {
                if(sum[w] == g->deg1(u)) {
                    H.containNodes[w].push_back(u);
                    H.contained[u] = true;
                }
                
                H.lists[w].push_back(u);
            }
            sum[w] = 0;
        }
        l = 0;
    }
}


// void BCListPlusPlus::exactCountMaximalPivot(int p, int q) {
//     memset(ansAll, 0, sizeof ansAll);

//     collect2HopNeighborsContainedCDAG(p, q);
//     printf("collect 2/contain\n"); fflush(stdout);

//     S.resize(g->n2);
//     // Q.resize(g->n1);

//     uint32_t maxDegree = 0;
//     for(uint32_t u = 0; u < g->n1; u++) {
//         maxDegree = std::max(maxDegree, 
//             (uint32_t)(H.lists[u].size()+H.containNodes[u].size()));
//     }
//     printf("maxDegree of H %u\n", maxDegree);

//     tmpNodes.resize(maxDegree + 1);
//     for(uint32_t i = 0; i < maxDegree + 1; i++) {
//         tmpNodes[i].resize(maxDegree + 1);
//     }

//     for(uint32_t u = g->n1; u >= p; u--) {
//         if(H.contained[u]) continue;
//         if(H.lists[u].size() < p) continue;
//         if(g->deg1(u) < q) continue;

//         int L = 0;
//         int pH = 0;
//         for(auto v: H.lists[u]) {
//             if(g->deg1(v) == g->deg1(u)) {
//                 L++;
//             }
//             else H.nodes.changeTo(v, pH++);
//         }

//         int pS = 0;
//         for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//             S.changeTo(g->e1[i], pS++);
//         }

//         pivotCount(0, L, pH, pS, 0, p, q);
//     }

//     printf("ans[%d][%d] %.2f\n", p, q, ansAll[p][q]);
//     fflush(stdout);
// }

// void BCListPlusPlus::pivotCount(int l, int L, int pH, int pS, int pQ, int p, int q) {
//     int pivotIdx = 0;
//     int pivotDeg = 0;

//     for(int i = 0; i < pH; i++) {
//         int deg = 0;
//         for(int j = 0; j < pH; j++) if(i != j) {
//             if(H.contain(H.nodes[i], H.nodes[j])) {
//                 deg++;
//             }
//         }

//         if(deg > pivotDeg) {
//             pivotIdx = i;
//             pivotDeg = deg;
//         }
//     }

//     int pivot = H.nodes[pivotIdx];
//     for(int j = 0; j < pH; j++) {
//         if(H.contain(pivot, H.nodes[j])) {
//             H.nodes.changeToByPos(j, --pH);
//         }
//     }

//   //  H.nodes.copy(tmpNodes[l].data(), pH);

//     while(pH > 0) {
//         uint32_t u = H.nodes[0];
//         H.nodes.changeTo(u, --pH);

//         int pSS = 0;
//         for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//             if(S.idx(g->e1[i]) < pS) {
//                 S.changeTo(g->e1[i], pSS++);
//             }
//         }

//         int addL = 1;
//         int pHH = 0;
//         for(int i = 0; i < pH; i++) {
//             uint32_t v = H.nodes[i];

//             int nV = 0;
//             for(uint32_t j = g->pU[v]; j < g->pU[v + 1]; j++) {
//                 if(S.idx(g->e1[j]) < pSS) {
//                     nV++;
//                 }
//             }

//             if(nV == pSS) {
//                 addL++;
//             }
//             else if(nV > 0) {
//                 H.nodes.changeToByPos(i, pHH++);
//             }
//         }

                
//     }
// }