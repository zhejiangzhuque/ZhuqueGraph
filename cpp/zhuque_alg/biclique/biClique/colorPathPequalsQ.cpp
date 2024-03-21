#include "colorPathPequalsQ.h"
#include <chrono>
#include <random>
#include <cassert>





void colorPathPequalsQ::approximateCountingAllVersion2(uint64_t T) {
// g->print();
    // g->coreReductionFast22();

    // printf("%u %u %llu\n", g->n1, g->n2, g->m);

    printf("maxDu:%u maxDv:%u\n", g->maxDu, g->maxDv);
    fflush(stdout);
    uint32_t maxDuv = std::max(g->maxDu, g->maxDv);

    // const int minPQ = 11;
    printf("H:%d\n", minPQ);
    printf("T:%llu\n", T);
    uint32_t maxE = std::min((1u<<21), g->maxDu * g->maxDv);
    // uint32_t maxE = 105000;

    candL.resize(g->n1);
    candR.resize(g->n2);

    double * bufferForDpU = new double[minPQ * maxE]();
    double * bufferForDpV = new double[minPQ * maxE]();

    double ** dpU = new double *[minPQ];
    for(uint32_t k = 0; k < minPQ; k++) {
        dpU[k] = bufferForDpU + k * maxE;
    }
    for(uint32_t e = 0; e < maxE; e++) {
        dpU[0][e] = 0;
        dpU[1][e] = 1;
    }

    double ** dpV = new double*[minPQ];
    for(uint32_t k = 0; k < minPQ; k++) {
        dpV[k] = bufferForDpV + k * maxE;
    }
    for(uint32_t e = 0; e < maxE; e++) {
        dpV[0][e] = 1;
    }

    ansAll.resize(minPQ + 1);
    // for(uint32_t i = 0; i <= minPQ; i++) {
    //     ansAll[i].resize(minPQ + 1);
    // }

    //subgraph
    std::vector<int> pU(g->maxDv+5), e1(maxE), pV(g->maxDu+5), e2(maxE);
    std::vector<int> mapUtoV(maxE), mapVtoU(maxE);
    std::vector<int> cnt(g->n1);
    std::vector<uint32_t> candLtemp(g->maxDv);
    //(minPQ-minPQ)-bipath
    std::vector<double> sumW(minPQ + 2);
    std::vector<double> ddp(maxE);

    //first compute dp
printf("start dp first\n");fflush(stdout);
    // uint32_t realMaxE = 0;

    auto computeDP = [&](int pL, int pR)->int {
        std::fill(pV.begin(), pV.begin() + pR + 1, 0);
        for(int j = 0; j < pL; j++) {
            uint32_t x = candL[j];
            pU[j + 1] = pU[j];
            
            if(pR < g->deg1(x))
            for(int k = 0; k < pR; k++) {
                uint32_t y = candR[k];
                if(g->connectUV(x, y)) {
                    e1[pU[j + 1]++] = k;
                    pV[k + 1]++;
                }
            }
            else {
                auto st = g->e1.begin() + g->pU[x];
                auto ed = g->e1.begin() + g->pU[x + 1];
                uint32_t i = std::upper_bound(st, ed, candR[0]) - g->e1.begin();
                for(; i < g->pU[x + 1]; i++) {
                    int k = candR.idx(g->e1[i]);
                    if(k < pR) {
                        e1[pU[j + 1]++] = k;
                        pV[k + 1]++;
                    }
                }
            }
        }
        assert(pU[pL] < maxE);
// printf(" %.0f\n", pU[pL]);
        if(pU[pL] == 0) return 0;

        for(int v = 0; v < pR; v++) {
            pV[v + 1] += pV[v];
        }
        assert(pU[pL] == pV[pR]);

        for(int u = 0; u < pL; u++) {
            for(int i = pU[u]; i < pU[u + 1]; i++) {
                int v = e1[i];

                e2[pV[v]] = u;

                mapUtoV[i] = pV[v];
                mapVtoU[pV[v]] = i;

                pV[v]++;
            }
        }

        for(int v = pR; v >= 1; v--) {
            pV[v] = pV[v - 1];
        }
        pV[0] = 0;
        //graph constructed done

        for(int u = 0; u < pL; u++) {
            for(int i = pU[u]; i < pU[u + 1]; i++) {
                // int v = e1[i];
                dpV[1][i] = pU[u + 1] - i - 1;
            }
        }

        int minLR = std::min(pL, pR);
        int k = 2;
        for(; k <= minLR && k < minPQ; k++) {
            // bool f = false;
            // memset(dpU[k], 0, sizeof(double) * pU[pL]);

            std::fill(ddp.begin(), ddp.begin() + pU[pL], 0);
            for(int v = 0; v < pR; v++) {
                for(int i = pV[v] + 1; i < pV[v + 1]; i++) {
                    int u = e2[i];
                    if(dpV[k - 1][mapVtoU[i]] < 0.5) continue;

                    ddp[pV[v]] += dpV[k - 1][mapVtoU[i]]; 
                    ddp[i] -= dpV[k - 1][mapVtoU[i]];
                    // ddp[pU[u]] += dpV[k - 1][i];
                    // int ed = std::lower_bound(e1.begin() + pU[u], 
                    //     e1.begin() + pU[u + 1], v) - e1.begin();
                    // ddp[ed] -= dpV[k - 1][i];
                    // for(int j = pU[u]; e1[j] != v; j++) {
                    //     dpU[k][j] += dpV[k - 1][i];
                    // }
                }
            }
            // dpU[k][0] = ddp[0];
            for(int v = 0; v < pR; v++) {
                dpU[k][pV[v]] = ddp[pV[v]] < 0 ? 0 : ddp[pV[v]];
                // if(dpU[k][pV[v]] < 0) dpU[k][pV[v]] = 0;
// assert(dpU[k][pV[v]] >= -1);
                for(int e = pV[v] + 1; e < pV[v + 1]; e++) {
                    dpU[k][e] = dpU[k][e-1] + ddp[e];
// assert(dpU[k][e] >= -1);
                    // if(dpU[k][e] < 0) dpU[k][e] = 0;
                }
            }
//             for(int e = 1; e < pV[pR]; e++) {
//                 dpU[k][e] = dpU[k][e-1] + ddp[e];
// assert(dpU[k][e] >= -1e-5);
//             }

            std::fill(ddp.begin(), ddp.begin() + pU[pL], 0);
            // memset(dpV[k], 0, sizeof(double) * pU[pL]);
            for(int u = 0; u < pL; u++) {
                for(int i = pU[u] + 1; i < pU[u + 1]; i++) {
                    int v = e1[i];

                    if(dpU[k][mapUtoV[i]] < 0.5) continue;

                    ddp[pU[u]] += dpU[k][mapUtoV[i]];
                    ddp[i] -= dpU[k][mapUtoV[i]];

                    // int j = std::lower_bound(e2.begin() + pV[v], 
                    //     e2.begin() + pV[v + 1], u) - e2.begin();

                    // ddp[pV[v]] += dpU[k][i];
                    // ddp[j] -= dpU[k][i];

                    // for(int j = pV[v]; e2[j] != u; j++) {
                    //     dpV[k][j] += dpU[k][i];
                    // }
                }
            }
//             dpV[k][0] = ddp[0];
//             for(int e = 1; e < pU[pL]; e++) {
//                 dpV[k][e] = dpV[k][e - 1] + ddp[e];
// assert(dpV[k][e] >= -1e-5);
//             }
            for(int u = 0; u < pL; u++) {
                dpV[k][pU[u]] = ddp[pU[u]] < 0? 0 : ddp[pU[u]];
// if(dpV[k][pU[u]] < -1) {
//     printf("dpV:%.2f\n", dpV[k][pU[u]]);
//     fflush(stdout);
// }
// assert(dpV[k][pU[u]] >= -1);
                for(int e = pU[u] + 1; e < pU[u + 1]; e++) {
                    dpV[k][e] = dpV[k][e-1] + ddp[e];
                    // if(dpV[k][e] < 0) dpV[k][e] = 0;
// assert(dpV[k][e] >= -1);
                }
            }
        }

        return k;
    };

    int maxPLen = 0;
    std::vector<uint32_t> js(g->maxDu);
double st = clock();
    for(uint32_t u = 0; u < g->n1; u++) {
// printf("%u\n", u);fflush(stdout);
        if(g->deg1(u) <= 1) continue;
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }

        candLtemp.clear();
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            for(uint32_t j = g->pV[v + 1] - 1; j >= g->pV[v]; j--) {
                uint32_t w = g->e2[j];
                if(w == u) {
                    js[i - g->pU[u]] = j;
                    break;
                }

                if(cnt[w] == 0) {
                    candLtemp.push_back(w);
                }
                cnt[w]++;
            }
        }

        int pLL = 0;
        for(auto w: candLtemp) {
            if(cnt[w] > 1) {
                candL.changeTo(w, pLL++);
            }
            cnt[w] = 0;
        }
// printf("u:%u %d %d\n", u, pLL, pR);fflush(stdout);
        if(pLL == 0) continue;

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];

            int pL = 0;
            // auto st = g->e2.begin() + g->pV[v];
            // auto ed = g->e2.begin() + g->pV[v + 1];
            // uint32_t j = std::upper_bound(st, ed, u) - g->e2.begin();
            uint32_t j = js[i - g->pU[u]] + 1;
            for(; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(candL.idx(w) < pLL) candL.changeTo(w, pL++);
            }
            candR.changeTo(v, --pR);
            if(pL == 0) continue;
            if(pR == 0) break;

            for(int i = 1; i < pR; i++) {
                candR.changeToByPos(i, i - 1);
            }
            

            // if(pL < minPQ - 1 || pR < minPQ - 1) continue;
// printf("%u %u, %d %d \n", u, v, pL, pR);
            int maxPLength = computeDP(pL, pR);
            maxPLen = std::max(maxPLen, maxPLength);

            for(int j = 0; j < pU[pL]; j++) {
                for(int k = 1; k < maxPLength && k < minPQ; k++) {
                    if(dpU[k][j] < 0.5) break;
                    sumW[k] += dpU[k][j];
                }
            }
        }
    }
printf("dp first %.2fs\n", (clock() - st) / CLOCKS_PER_SEC);
    printf("maxPlen %d\n",  maxPLen);
    for(int i = 1; i < maxPLen && i < minPQ; i++) {
        printf("sumW:%.0f ", sumW[i]);
    }
    printf("\n\n");fflush(stdout);
    
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> uiDistribution(0, 1);
    std::vector<int> stackL(g->maxDv+1), stackR(g->maxDu+1);
    std::vector<double> sumCXi(g->maxDu + 5), sumCYi(g->maxDv + 5);

    std::vector<int> maxZL(g->maxDu + 1), maxZR(g->maxDv + 1);

    // if(sumW > 0)
    for(uint32_t u = 0; u < g->n1; u++) {
// if(u > 106000) {
    // printf("u:%u\n", u);fflush(stdout);
// }
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }

        candLtemp.clear();
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            for(uint32_t j = g->pV[v + 1] - 1; j >= g->pV[v]; j--) {
                uint32_t w = g->e2[j];
                if(w == u) break;

                if(cnt[w] == 0) {
                    candLtemp.push_back(w);
                }
                cnt[w]++;
            }
        }

        int pLL = 0;
        for(auto w: candLtemp) {
            if(cnt[w] > 1) {
                candL.changeTo(w, pLL++);
            }
            cnt[w] = 0;
        }

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];

            int pL = 0;
            auto st = g->e2.begin() + g->pV[v];
            auto ed = g->e2.begin() + g->pV[v + 1];
            uint32_t j = std::upper_bound(st, ed, u) - g->e2.begin();
            for(; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(candL.idx(w) < pLL) candL.changeTo(w, pL++);
            }
            candR.changeTo(v, --pR);
            for(int i = 1; i < pR; i++) {
                candR.changeToByPos(i, i - 1);
            }

            // if(pL < minPQ - 1 || pR < minPQ - 1) continue;

            int maxPLength = computeDP(pL, pR);

// printf("    Graph: %u %u:\n", u, v);
// printf("U:\n");
// for(int u = 0; u < pL; u++) {
//     printf("%d:", candL[u]);
//     for(int i = pU[u]; i < pU[u + 1]; i++) {
//         printf("%d-%d ", i, candR[e1[i]]);
//     }
//     printf("\n");
// }
// printf("V:\n");
// for(int v = 0; v < pR; v++) {
//     printf("%d:", candR[v]);
//     for(int i = pV[v]; i < pV[v + 1]; i++) {
//         printf("%d-%d ", i, candL[e2[i]]);
//     }
//     printf("\n");
// }

            for(int len = 2; len <= maxPLength && len < minPQ && sumW[len - 1] > 0; len++) {//(len-len)-bipath
                double sumWtemp = 0.0;
                for(int j = 0; j < pU[pL]; j++) {
                    sumWtemp += dpU[len - 1][j];
                }

                uint64_t sampleSize = std::ceil(T*sumWtemp / sumW[len - 1]);
// if(sampleSize > T) {
//     printf("%d %llu %.0f %.0f\n", len, sampleSize, sumWtemp, sumW[len - 1]);
//     fflush(stdout);
// }
// assert(sampleSize <= T);
                if(sampleSize == 0) continue;

                sumCXi[0] = 0.0;
                // std::fill(sumCXi.begin(), sumCXi.begin() + pR + 2, 0.0);
                // std::fill(sumCYi.begin(), sumCYi.begin() + pL + 2, 0.0);
// printf("len %d: sumW:%.0f spSize:%llu\n", len, sumWtemp, sampleSize);
                while(sampleSize--) {
                    double tmp = 0.0;
                    double r = uiDistribution(generator); 
                    int preE = 0, preU = 0, preV = 0;
                    stackL.clear();
                    stackR.clear();

                    for(int u = 0; u < pL; u++) {
                        for(int i = pU[u]; i < pU[u + 1]; i++) {
                            tmp += dpU[len - 1][mapUtoV[i]];
                            if(tmp + 1e-8 >= r * sumWtemp) {
                                stackL.push_back(u);
                                stackR.push_back(e1[i]);

                                preU = u;
                                preV = e1[i];
                                preE = i;
                                break;
                            }
                        }

                        if(stackL.size() > 0) break;
                    }
                    assert(stackL.size() > 0);

                    for(int i = 1; i < len - 1; i++) {
                        double r = uiDistribution(generator);
                        tmp = 0.0;
                    
                        for(int j = pV[preV]; j < pV[preV + 1]; j++) {
                            if(e2[j] <= preU) continue;
                            tmp += dpV[len - i - 1][mapVtoU[j]];
                            if(tmp + 1e-8 >= r * dpU[len - i][mapUtoV[preE]]) {
                                int u = e2[j];
                                stackL.push_back(u);
                                preU = u;
                                preE = j;
                                break;
                            }
                        }

                        r = uiDistribution(generator);
                        tmp = 0.0;
                        for(int j = pU[preU]; j < pU[preU + 1]; j++) {
                            if(e1[j] <= preV) continue;
                            tmp += dpU[len - i - 1][mapUtoV[j]];
                            if(tmp + 1e-8 >= r * dpV[len - i - 1][mapVtoU[preE]]) {
                                int v = e1[j];
                                stackR.push_back(v);
                                preV = v;
                                preE = j;
                                break;
                            }
                        }
                    }

                    if(stackL.size() < len - 1) {
                        continue;
                    }
                    assert(stackL.size() == len - 1);
// printf("candL:");
// for(int i = 0; i < len - 1; i++) printf("%d ", candL[stackL[i]]);
// printf("  ");
// printf("candR:");
// for(int i = 0; i < len - 1; i++) printf("%d ", candR[stackR[i]]);
// printf("\n");

                    bool connect = true;
                    for(int i = 0; i < len - 1; i++) {
                        for(int j = 0; j < len - 1; j++) {
                            if(i == j) continue;
                            if(i > 0 && i == j + 1) continue; 

                            if(!g->connectUV(candL[stackL[i]], candR[stackR[j]])) {
                                connect = false;
                                break;
                            }
                        }
                        if(!connect) break;
                    }
                    if(!connect) continue;
#ifdef PQWEDGE
 
#endif

                    int rSize = 0, kk = 0;
                    for(int j = 0; j < pR; j++) {
                        int v = candR[j];
                        if(kk < len - 1 && stackR[kk] == j) {
                            kk++;
                            continue;
                        }
                        bool f = true;
                        for(int k = 0; k < len - 1; k++) {
                            if(!g->connectUV(candL[stackL[k]], v)) {
                                f = false;
                                break;
                            }
                        }

                        if(f) rSize++;
                        // if(pR - j - 1 + rSize < maxPQ - minPQ) break;
                    }

                    // for(int x = 0; x <= rSize && x < minPQ; x++) {
                    //     sumCXi[x] += C[rSize][x];
                    // }
                    sumCXi[0] += C[rSize][0];
                    maxZR[len] = std::max(maxZR[len], rSize);

//                     int lSize = 0;
//                     kk = 0;
//                     for(int j = 0; j < pL; j++) {
// // if(u == 3 && v == 4 && candL[stackL[0]] == 4 && candR[stackR[0]] == 5) {
// //     printf("x:%d %d %u\n", j, stackL[kk], candL[stackL[kk]]);
// // }
//                         if(kk < len - 1 && stackL[kk] == j) {
//                             kk++;
//                             continue;
//                         }
// // if(u == 3 && v == 4 && candL[stackL[0]] == 4 && candR[stackR[0]] == 5) {
// //     printf("x:%d\n", j);
// // }
//                         int u = candL[j];
//                         bool f = true;
//                         for(int k = 0; k < len - 1; k++) {
//                             if(!g->connectUV(u, candR[stackR[k]])) {
//                                 f = false;
//                                 break;
//                             }
//                         }
//                         if(f) lSize++;
//                         // if(pR - j - 1 + rSize < maxPQ - minPQ) break;
//                     }

//                     for(int x = 1; x <= lSize && x < minPQ; x++) {
//                         sumCYi[x] += C[lSize][x];
//                     }
//                     maxZL[len] = std::max(maxZL[len], lSize);
// printf("lr %d %d\n", lSize, rSize);
                }

// for(int i = 0; i <= pR - (len-1); i++) {
//     printf("%.0f ", sumCXi[i]);
// }
// printf("CX\n");
// for(int i = 0; i <= pL - (len-1); i++) {
//     printf("%.0f ", sumCYi[i]);
// }
// printf("CY\n");

                sampleSize = std::ceil(T * sumWtemp / sumW[len - 1]);
                // if(ansAll[len].size() < pR) ansAll[len].resize((pR + 2)*2);
                if(len <= pR + 1 && len < minPQ && sumCXi[0] > 0.5)
                    ansAll[len] += sumCXi[0] * sumWtemp / sampleSize / C[len - 1][len - 1];

                // for(int x = 0; x + len <= pR + 1 && x + len < minPQ; x++) {
                //     if(sumCXi[x] < 0.5) break;
                //     ansAll[len][len + x] += sumCXi[x] * sumWtemp / sampleSize / C[len + x - 1][len - 1];
                // }

                // for(int x = 1; x + len <= pL + 1 && x + len < minPQ; x++) {
                //     if(sumCYi[x] < 0.5) break;
                //     ansAll[x + len][len] += sumCYi[x] * sumWtemp / sampleSize / C[len + x - 1][len - 1];
                // }
// printf("    Ans[3][2]:%.0f %.0f\n", ansAll[3][2], sumCYi[1] * sumWtemp / sampleSize / C[3][2]);
            }
        }
    }

    // for(int x = 2; x < (int)ansAll.size() && x < minPQ; x++) {
    //     for(int y = 2; y < (int)ansAll[x].size() && y < minPQ; y++) {
    //         if(ansAll[x][y] < 0.5) break;
    //         printf("%d-%d: %.0f\n", x, y, ansAll[x][y]);
    //     }
    // }
    for(int x = 2; x < (int)ansAll.size() && x < minPQ; x++) {
        // for(int y = 2; y < (int)ansAll[x].size() && y < minPQ; y++) {
            if(ansAll[x] < 0.5) break;
            printf("%d-%d: %.0f\n", x, x, ansAll[x]);
        // }
    }
    for(int i = 2; i < minPQ && i < maxPLen; i++) {
        printf("%d:maxZL %u maxZR %u\n", i, maxZL[i], maxZR[i]);
    }
    fflush(stdout);

    delete [] bufferForDpU;
    delete [] bufferForDpV;
    delete [] dpU;
    delete [] dpV;
}

