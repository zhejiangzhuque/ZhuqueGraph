#include "pivotAndPath.h"
#include <chrono>
#include <random>
#include <cassert>

void pivotAndPath::countingV5(uint64_t T) {
    std::vector<uint32_t> vs1(g->n1), vs2(g->n1);
    vs1.clear();
    vs2.clear();

    constexpr int bar = 1000;
    std::vector<uint32_t> degs(g->n2);
    for(uint32_t v = 0; v < g->n2; v++) {
        degs[v] = g->deg2(v);
    }

    // uint64_t upperBoundEdgesForSampling = 0;
    for(uint32_t u = 0; u < g->n1; u++) {
        uint32_t sumDeg = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];

            --degs[v];
            sumDeg += degs[v] * (g->pU[u + 1] - i - 1);
        }
        
        if(sumDeg >= bar) {
            vs2.push_back(u);
            // upperBoundEdgesForSampling += g->deg1(u);
        }
        else vs1.push_back(u);
    }

    printf("bar %d\n", bar);
    printf("sample size %llu\n", T);
    // printf("upperBound %llu\n", upperBoundEdgesForSampling);

    printf("exact:%u\n", (uint32_t)vs1.size());
    printf("sample:%u\n", (uint32_t)vs2.size());
    fflush(stdout);

double t = clock();
    if(vs1.size() > 0) pivot(vs1);
double d = clock();

printf("exact done, %.5fs\n", (d - t) / CLOCKS_PER_SEC);fflush(stdout);
    if(vs2.size() > 0) samplev5(vs2, T);
double e = clock();

printf("appro done, %.5fs\n", (e - d) / CLOCKS_PER_SEC);fflush(stdout);
}

void pivotAndPath::samplev5(std::vector<uint32_t> & nodes, uint64_t T) {
    uint32_t maxDuv = std::max(g->maxDu, g->maxDv);

    // const int minPQ = 5;
    printf("H:%d\n", minPQ);
    printf("T:%llu\n", T);
    // uint32_t maxE = std::min((1u<<21), g->maxDu * g->maxDv + 5);
    uint32_t maxE = 0;
    for(uint32_t u = 0; u < g->n1; u++) {
        uint32_t edges = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            auto st = g->e2.begin() + g->pV[v];
            auto ed = g->e2.begin() + g->pV[v + 1];
            edges += ed - std::lower_bound(st, ed, u) + 1;
        }
        maxE = std::max(maxE, edges);
    }
    printf("maxE %u, maxDuv %u\n", maxE, maxDuv);
    // uint32_t maxE = 105000;

    // candL.resize(g->n1);
    // candR.resize(g->n2);
    std::vector<uint32_t> candLtemp(g->maxDv);

    // double * bufferForDpU = new double[minPQ * maxE];
    // double * bufferForDpV = new double[minPQ * maxE];

    // double ** dpU = new double *[minPQ];
    // for(uint32_t k = 0; k < minPQ; k++) {
    //     dpU[k] = bufferForDpU + k * maxE;
    // }
    std::vector<std::vector<double>> dpU(minPQ);
    for(uint32_t k = 0; k < minPQ; k++) dpU[k].resize(maxE);

    for(uint32_t e = 0; e < maxE; e++) {
        dpU[0][e] = 0;
        dpU[1][e] = 1;
    }

    // double ** dpV = new double*[minPQ];
    // for(uint32_t k = 0; k < minPQ; k++) {
    //     dpV[k] = bufferForDpV + k * maxE;
    // }
    std::vector<std::vector<double>> dpV(minPQ);
    for(uint32_t k = 0; k < minPQ; k++) dpV[k].resize(maxE);
    for(uint32_t e = 0; e < maxE; e++) {
        dpV[0][e] = 1;
    }

    // ansAll.resize(minPQ + 5);
    // for(uint32_t i = 0; i <= minPQ; i++) {
    //     ansAll[i].resize(minPQ + 5);
    // }

    //(minPQ-minPQ)-bipath
    std::vector<double> sumW(minPQ + 1);
    std::vector<double> ddp(maxE);
    
    //subGraph
    std::vector<std::vector<uint32_t>> sg(g->n1);
    std::vector<int> cnt(g->n1);
    std::vector<std::vector<uint32_t>> nxtLayerEdgesL(g->n1), nxtLayerEdgesR(g->maxDu + 5);
    //subgraph
    std::vector<int> pU(g->n1 + 5), e1(maxE), pV(g->maxDu + 5), e2(maxE);
    std::vector<int> mapUtoV(maxE), mapVtoU(maxE);

    //first compute dp
    printf("start dp first v5\n");fflush(stdout);

// int aa = 0, bb = 0, cc = 0;
    auto computeDP = [&](int pL, int pR) {
        for(int u = 0; u < pL; u++) {
            for(int i = pU[u]; i < pU[u + 1]; i++) {
                dpV[1][i] = pU[u + 1] - i - 1;
            }
        }

        int minLR = std::min(pL, pR);

        int k = 2;

        for(int v = 0; v < pR; v++) {
            nxtLayerEdgesR[v].clear();
            for(int i = pV[v] + 1; i < pV[v + 1]; i++) {
                nxtLayerEdgesR[v].push_back(i);
            }
        }
        for(int u = 0; u < pL; u++) nxtLayerEdgesL[u].clear();
        
// int cntt = 0;
        for(; k <= minLR && k < minPQ; k++) {
            memset(dpU[k].data(), 0, sizeof(double)*pU[pL]);
            memset(dpV[k].data(), 0, sizeof(double)*pV[pR]);
            // std::fill(ddp.begin(), ddp.begin() + pU[pL], 0);
            memset(ddp.data(), 0, sizeof(double)* pU[pL]);
            
            for(int v = 0; v < pR; v++) {
                // for(int i = pV[v] + 1; i < pV[v + 1]; i++) {
                if(nxtLayerEdgesR[v].size() == 0) {
                    // memset(dpU[k].data() + pV[v], 0, sizeof(double)*(pV[v + 1] - pV[v]));
                    continue;
                }
// if(nxtLayerEdgesR[v].size() > pV[v + 1] - pV[v]) {
//     printf("k %d v %d %u\n", k, v, candR[v]);
// }
// assert(nxtLayerEdgesR[v].size() <= pV[v + 1] - pV[v]);

                for(auto i : nxtLayerEdgesR[v]) {
                    // int u = e2[i];
                    // if(dpV[k - 1][mapVtoU[i]] < 0.5) continue;

                    ddp[pV[v]] += dpV[k - 1][mapVtoU[i]]; 
                    ddp[i] -= dpV[k - 1][mapVtoU[i]];
// cntt += i - pV[v];
                }
                
                dpU[k][pV[v]] = ddp[pV[v]] < 0 ? 0 : ddp[pV[v]];
                if(ddp[pV[v]] >= 1) {
                    nxtLayerEdgesL[e2[pV[v]]].push_back(mapVtoU[pV[v]]);
                    // assert(e2[pV[v]] < pL);
                }
                    
                // for(int e = pV[v] + 1; e < pV[v + 1]; e++) {
                for(int e = pV[v] + 1; e < pV[v + 1]; e++) {
                    dpU[k][e] = dpU[k][e - 1] + ddp[e];
                    if(dpU[k][e] >= 1) {
// for(auto i : nxtLayerEdgesL[e2[e]]) assert(i != mapVtoU[e]);
                        nxtLayerEdgesL[e2[e]].push_back(mapVtoU[e]);
// if(e2[e] >= pL) {
//     printf("%d %d %d\n", e, e2[e], pL);fflush(stdout);
// }
                        assert(e2[e] < pL);
                    }
                        
                }
                // for(int e = lastEdge + 1; e < pV[v + 1]; e++) dpU[k][e] = 0;

                nxtLayerEdgesR[v].clear();
            }
            // for(int v = 0; v < pR; v++) {
            //     dpU[k][pV[v]] = ddp[pV[v]] < 0 ? 0 : ddp[pV[v]];
            //     for(int e = pV[v] + 1; e < pV[v + 1]; e++) {
            //         dpU[k][e] = dpU[k][e-1] + ddp[e];
            //     }
            // }

            // std::fill(ddp.begin(), ddp.begin() + pU[pL], 0);
            memset(ddp.data(), 0, sizeof(double)* pU[pL]);
            for(int u = 0; u < pL; u++) {
                // for(int i = pU[u] + 1; i < pU[u + 1]; i++) {
                if(nxtLayerEdgesL[u].size() == 0) {
                    // memset(dpV[k].data() + pU[u], 0, sizeof(double) * (pU[u + 1] - pU[u]));
                    continue;
                }
// assert(nxtLayerEdgesL[u].size() <= pU[u + 1] - pU[u]);
                for(auto i: nxtLayerEdgesL[u]) {
                    // int v = e1[i];
                    // if(dpU[k][mapUtoV[i]] < 0.5) continue;

                    ddp[pU[u]] += dpU[k][mapUtoV[i]];
                    ddp[i] -= dpU[k][mapUtoV[i]];
// cntt += i - pU[u];
                }

                dpV[k][pU[u]] = ddp[pU[u]] < 0? 0 : ddp[pU[u]];
                if(ddp[pU[u]] >= 1) {
// for(auto i : nxtLayerEdgesR[e1[pU[u]]]) assert(i != mapUtoV[pU[u]]);
                    nxtLayerEdgesR[e1[pU[u]]].push_back(mapUtoV[pU[u]]);
                }
                    
                // nxtLayerEdgesR[u] = 
                // for(int e = pU[u] + 1; e < pU[u + 1]; e++) {
                for(int e = pU[u] + 1; e < pU[u + 1]; e++) {
                    dpV[k][e] = dpV[k][e-1] + ddp[e];
                    if(dpV[k][e] >= 1) {
// for(auto i : nxtLayerEdgesR[e1[e]]) assert(i != mapUtoV[e]);
                        nxtLayerEdgesR[e1[e]].push_back(mapUtoV[e]);
// assert(nxtLayerEdgesR[e1[e]].size() <= pV[e1[e] + 1] - pV[e1[e]]);
                    }
                        
                }
                // for(int e = lastEdge + 1; e < pU[u + 1]; e++) dpV[k][e] = 0;

                nxtLayerEdgesL[u].clear();
            }
// printf("%d %d\n", cntt, pU[pL]);
// cc++;
// if(cntt <= 0.33*pU[pL]) {
//     aa++;
// }
// else if(cntt <= 0.66*pU[pL]) {
//     bb++;
// }
        }

        return k;
    };
// g->print();

   
double stT = clock();
    int maxPLen = 0;
    for(auto u : nodes) {
// if(u % 1000 == 0) {
// printf("    nodes %u\n", u);fflush(stdout);
// }
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }
        if(pR < 2) continue;

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

        std::sort(candLtemp.begin(), candLtemp.end());

        int pL = 1;
        candL.changeTo(u, 0);
        for(auto w: candLtemp) {
            if(cnt[w] > 1) {
                candL.changeTo(w, pL++);
            }
            cnt[w] = 0;
        }

        if(pL == 1) continue;

        std::fill(pV.begin(), pV.begin() + pR + 1, 0);
        // for(int i = 0; i < pL; i++) sg[i].clear();

        for(int i = 0; i < pR; i++) {
            uint32_t v = candR[i];

            for(uint32_t j = g->pV[v + 1] - 1; j >= g->pV[v]; j--) {
                uint32_t w = g->e2[j];
                if(w < u) break;

                uint32_t wId = candL.idx(w);
                if(wId >= pL) continue;

                sg[wId].push_back(i);

                pV[i + 1]++;

                if(j == 0) break;
            }
        }

        pU[0] = 0;
        for(int i = 0; i < pL; i++) {
            pU[i + 1] = pU[i]  + sg[i].size();
            memcpy(e1.data() + pU[i], sg[i].data(), sizeof(uint32_t) * sg[i].size());
        }
assert(pU[pL] < maxE);
        if(pU[pL] == 0) continue;

        for(int i = 0; i < pR; i++) pV[i + 1] += pV[i];
// assert(pU[pL] == pV[pR]);
        for(int u = 0; u < pL; u++) {
            for(int i = pU[u]; i < pU[u + 1]; i++) {
                int v = e1[i];
                e2[pV[v]] = u;
                
                mapUtoV[i] = pV[v];
                mapVtoU[pV[v]] = i;

                pV[v]++;
            }
        }
        for(int v = pR; v >= 1; v--) pV[v] = pV[v - 1];
        pV[0] = 0;
// assert(pU[pL] == pV[pR]);

        for(int i = 0; i < pL; i++) {
            sg[i].clear();
        }
// printf("U:\n");
// for(int u = 0; u < pL; u++) {
//     printf("%d:", candL[u]);
//     for(int i = pU[u]; i < pU[u + 1]; i++) {
//         printf("%d ", candR[e1[i]]);
//     }
//     printf("/");
//     for(auto v : sg[u]) printf("%d ", candR[v]);
//     printf("\n");
// }
// printf("V:\n");
// for(int v = 0; v < pR; v++) {
//     printf("%d:", candR[v]);
//     for(int i = pV[v]; i < pV[v + 1]; i++) {
//         printf("%d ", candL[e2[i]]);
//     }
//     printf("\n");
// }
        int maxPLength = computeDP(pL, pR);

        maxPLen = std::max(maxPLen, maxPLength);

        for(int j = pU[0]; j < pU[1]; j++) {
            for(int k = 2; k < maxPLength && k < minPQ; k++) {
                sumW[k] += dpU[k][j];
            }
        }
    }
double edT = clock();
printf("first dp time %f\n", (edT - stT) / CLOCKS_PER_SEC);
// printf("%d %d %d\n", aa, bb, cc);
// return ;
    printf("maxPlen %d\n\nsumW ",  maxPLen);
    for(int i = 2; i < maxPLen; i++) {
        printf("%d:%.0f ", i, sumW[i]);
    }
    printf("\n\n");fflush(stdout);
    
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> uiDistribution(0, 1);
    std::vector<int> stackL(minPQ + 5), stackR(minPQ + 5);
    std::vector<double> sumCXi(g->maxDu + 5), sumCYi(g->n1 + 1);

    std::vector<int> maxZL(g->maxDu + 1), maxZR(g->maxDv + 1);

    // if(sumW > 0)
    for(auto u : nodes) {
// if(u > 106000) {
//     printf("u:%u\n", u);fflush(stdout);
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
        std::sort(candLtemp.begin(), candLtemp.end());

        int pL = 1;
        candL.changeTo(u, 0);
        for(auto w: candLtemp) {
            if(cnt[w] > 1) {
                candL.changeTo(w, pL++);
            }
            cnt[w] = 0;
        }
        
        std::fill(pV.begin(), pV.begin() + pR + 1, 0);
        for(int i = 0; i < pR; i++) {
            uint32_t v = candR[i];

            for(uint32_t j = g->pV[v + 1] - 1; j >= g->pV[v]; j--) {
                uint32_t w = g->e2[j];
                if(w < u) break;

                uint32_t wId = candL.idx(w);
                if(wId >= pL) continue;
                sg[wId].push_back(i);

                pV[i + 1]++;

                if(j == 0) break;
            }
        }

        pU[0] = 0;
        for(int i = 0; i < pL; i++) {
            pU[i + 1] = pU[i]  + sg[i].size();
            memcpy(e1.data() + pU[i], sg[i].data(), sizeof(uint32_t) * sg[i].size());
        }

        for(int i = 0; i < pR; i++) pV[i + 1] += pV[i];
        for(int u = 0; u < pL; u++) {
            for(int i = pU[u]; i < pU[u + 1]; i++) {
                int v = e1[i];
                e2[pV[v]] = u;
                
                mapUtoV[i] = pV[v];
                mapVtoU[pV[v]] = i;

                pV[v]++;
            }
        }
        for(int v = pR; v >= 1; v--) pV[v] = pV[v - 1];
        pV[0] = 0;

        for(int i = 0; i < pL; i++) {
            sg[i].clear();
        }

// printf("U:\n");
// for(int u = 0; u < pL; u++) {
//     printf("%d:", candL[u]);
//     for(int i = pU[u]; i < pU[u + 1]; i++) {
//         printf("%d ", candR[e1[i]]);
//     }
//     printf("\n");
// }
// printf("V:\n");
// for(int v = 0; v < pR; v++) {
//     printf("%d:", candR[v]);
//     for(int i = pV[v]; i < pV[v + 1]; i++) {
//         printf("%d ", candL[e2[i]]);
//     }
//     printf("\n");
// }
// printf("%u, %d %d ", u, pL, pR);
        int maxPLength = computeDP(pL, pR);
// printf("\n dpU: \n");
// for(int u = 0; u < pL; u++)
// for(int j = pU[u]; j < pU[u + 1]; j++) {
// printf("%u-%u:", candL[u], candR[e1[j]]);
// for(int k = 1; k < maxPLength && k < minPQ; k++) {
//     printf("%.0f ", dpU[k][mapUtoV[j]]);
// }
// printf("\n");
// }
// printf("\n dpV: \n");
// for(int v = 0; v < pR; v++)
// for(int j = pV[v]; j < pV[v + 1]; j++) {
// printf("%u-%u:", candL[e2[j]], candR[v]);
// for(int k = 1; k < maxPLength && k < minPQ; k++) {
//     printf("%.0f ", dpV[k][mapVtoU[j]]);
// }
// printf("\n");
// }
// fflush(stdout);

        for(int len = 2; len < maxPLength && len <= minPQ && sumW[len] > 0.5; len++) {//(len-len)-bipath
            double sumWtemp = 0.0;
            for(int j = pU[0]; j < pU[1]; j++) {
                sumWtemp += dpU[len][mapUtoV[j]];
            }

            uint64_t sampleSize = std::ceil(sumWtemp / sumW[len] * T);

            if(sampleSize == 0) continue;

            std::fill(sumCXi.begin(), sumCXi.begin() + pR + 1, 0.0);
            std::fill(sumCYi.begin(), sumCYi.begin() + pL + 1, 0.0);
// printf("len %d: sumW:%.0f spSize:%llu\n", len, sumWtemp, sampleSize);
            while(sampleSize--) {
                double tmp = 0.0;
                double r = uiDistribution(generator); 
                int preE = 0, preU = 0, preV = 0;
                stackL.clear();
                stackR.clear();

                stackL.push_back(0);
                preU = 0;
                for(int i = pU[0]; i < pU[1]; i++) {
                    tmp += dpU[len][mapUtoV[i]];
                    if(tmp + 1e-8 >= r * sumWtemp) {        
                        stackR.push_back(e1[i]);
                        
                        preV = e1[i];
                        preE = i;
                        break;
                    }
                }
// if(stackR.size() == 0) {
//     printf("err:\n");

//     for(int i = pU[0]; i < pU[1]; i++) {
//         printf("%.0f ", dpU[len][mapUtoV[i]]);
//     }
//     printf("\n");
//     printf("r:%.2f sumW:%.2f r*sumW:%.2f\n", r, sumWtemp, r * sumWtemp);
//     fflush(stdout);
// }
                assert(stackR.size() > 0);

                for(int i = 1; i < len; i++) {
                    double r = uiDistribution(generator);
                    tmp = 0.0;

                    auto st = e2.begin() + pV[preV];
                    auto ed = e2.begin() + pV[preV + 1];
                    int j = std::upper_bound(st, ed, preU) - e2.begin();
                    for(; j < pV[preV + 1]; j++) {
                        tmp += dpV[len - i][mapVtoU[j]];
                        if(tmp + 1e-8 >= r * dpU[len - i + 1][mapUtoV[preE]]) {
                            int u = e2[j];
                            stackL.push_back(u);
                            preU = u;
                            preE = j;
                            break;
                        }
                    }

                    r = uiDistribution(generator);
                    tmp = 0.0;
                    st = e1.begin() + pU[preU];
                    ed = e1.begin() + pU[preU + 1];
                    j = std::upper_bound(st, ed, preV) - e1.begin();
                    for(; j < pU[preU + 1]; j++) {
                        // if(e1[j] <= preV) continue;
                        tmp += dpU[len - i][mapUtoV[j]];
                        if(tmp + 1e-8 >= r * dpV[len - i][mapVtoU[preE]]) {
                            int v = e1[j];
                            stackR.push_back(v);
                            preV = v;
                            preE = j;
                            break;
                        }
                    }
                }

                if(stackL.size() < len) {
                    continue;
                }
                assert(stackL.size() == len);
// printf("candL:");
// for(int i = 0; i < len; i++) printf("%d ", candL[stackL[i]]);
// printf("  ");
// printf("candR:");
// for(int i = 0; i < len; i++) printf("%d ", candR[stackR[i]]);
// printf("\n");

                bool connect = true;
                for(int i = 0; i < len; i++) {
                    for(int j = 0; j < len; j++) {
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

                int rSize = 0, kk = 0;
                for(int j = 0; j < pR; j++) {
                    int v = candR[j];
                    if(kk < len && stackR[kk] == j) {
                        kk++;
                        continue;
                    }
                    bool f = true;
                    for(int k = 0; k < len; k++) {
                        if(!g->connectUV(candL[stackL[k]], v)) {
                            f = false;
                            break;
                        }
                    }

                    if(f) rSize++;
                    // if(pR - j - 1 + rSize < maxPQ - minPQ) break;
                }

                for(int x = 0; x <= rSize && x < minPQ; x++) {
                    sumCXi[x] += C[rSize][x];
                }
                maxZR[len] = std::max(maxZR[len], rSize);

                int lSize = 0;
                kk = 0;
                for(int j = 0; j < pL; j++) {
                    if(kk < len && stackL[kk] == j) {
                        kk++;
                        continue;
                    }

                    int u = candL[j];
                    bool f = true;
                    for(int k = 0; k < len; k++) {
                        if(!g->connectUV(u, candR[stackR[k]])) {
                            f = false;
                            break;
                        }
                    }
                    if(f) lSize++;
                    // if(pR - j - 1 + rSize < maxPQ - minPQ) break;
                }

                for(int x = 1; x <= lSize && x < minPQ; x++) {
                    sumCYi[x] += C[lSize][x];
                }
                maxZL[len] = std::max(maxZL[len], lSize);
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

            sampleSize = std::ceil(sumWtemp / sumW[len] * T);
            // if(ansAll[len].size() < pR) ansAll[len].resize((pR + 2)*2);

            for(int x = 0; x + len <= pR && x + len < minPQ; x++) {
                if(sumCXi[x] < 0.5) break;
                ansAll[len][len + x] += sumCXi[x] * sumWtemp / sampleSize / C[len + x][len];
            }

            for(int x = 1; x + len <= pL && x + len < minPQ; x++) {
                if(sumCYi[x] < 0.5) break;
                ansAll[x + len][len] += sumCYi[x] * sumWtemp / sampleSize / C[len + x - 1][len - 1];
            }
// printf("    ansAll[3][2] %.0f\n", ansAll[3][2]);
        }
    }
edT = clock();
printf("second dp time %f\n", (edT - stT) / CLOCKS_PER_SEC);

    for(int x = 2; x < (int)ansAll.size() && x < minPQ; x++) {
        for(int y = 2; y < (int)ansAll[x].size() && y < minPQ; y++) {
            if(ansAll[x][y] < 0.5) break;
            printf("%d-%d: %.0f\n", x, y, ansAll[x][y]);
        }
    }
    for(int i = 2; i < minPQ && i < maxPLen; i++) {
        printf("%d:maxZL %u maxZR %u\n", i, maxZL[i], maxZR[i]);
    }
    printf("ed\n");
    fflush(stdout);

}