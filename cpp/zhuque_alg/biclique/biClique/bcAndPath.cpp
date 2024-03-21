#include "bcAndPath.h"
#include <cmath>
#include <random>

void bcAndPath::counting(uint64_t T, double realV, uint32_t bar) {
    std::vector<uint32_t> vs1(g->n1), vs2(g->n1);
    vs1.clear();
    vs2.clear();

    // uint32_t bar = 0;
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

    printf("bcAndPivot bar %d\n", bar);
    printf("sample size %llu\n", T);
    // printf("upperBound %llu\n", upperBoundEdgesForSampling);

    printf("exact:%u\n", (uint32_t)vs1.size());
    printf("sample:%u\n", (uint32_t)vs2.size());
    fflush(stdout);

double t = clock();
    if(vs1.size() > 0) exactCount(vs1);
double d = clock();

printf("exact done, %.5fs, ans %.0f\n", (d - t) / CLOCKS_PER_SEC, ans);
fflush(stdout);
    if(vs2.size() > 0) sample(vs2, T);
double e = clock();

printf("appro done, %.5fs, ans %.0f\n", (e - d) / CLOCKS_PER_SEC, ans);
printf("realV %f\n", realV);
printf("error: %f\n", std::abs(ans - realV) / realV);
fflush(stdout);
}

void bcAndPath::sample(std::vector<uint32_t> & nodes, uint64_t T) {
    printf("v2 maxDu:%u maxDv:%u\n", g->maxDu, g->maxDv);
    // printf("H")
    // fflush(stdout);
    uint32_t maxDuv = std::max(g->maxDu, g->maxDv);

    // const int minPQ = 11;
    printf("H=min(p,q):%d\n", minPQ);
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

    // ansAll.resize(minPQ + 1);
    // for(uint32_t i = 0; i <= minPQ; i++) {
    //     ansAll[i].resize(minPQ + 1);
    // }

    //subgraph
    std::vector<int> pU(g->maxDv+5), e1(maxE), pV(g->maxDu+5), e2(maxE);
    std::vector<int> mapUtoV(maxE), mapVtoU(maxE);
    std::vector<int> cnt(g->n1);
    std::vector<uint32_t> candLtemp(g->maxDv);
    //(minPQ-minPQ)-bipath
    double sumW = 0.0;
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
        if(pU[pL] == 0) return 1;

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
    for(auto u:nodes) {
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

            if(maxPLength < minPQ) continue;

            for(int j = 0; j < pU[pL]; j++) {
                // for(int k = 1; k < maxPLength && k < minPQ; k++) {
                    // if(dpU[minPQ-1][j] < 0.5) continue;

                    // if(k == minPQ-1)
                        sumW += dpU[minPQ - 1][j];
                // }
            }
        }
    }
printf("dp fitsrt %.2fs\n", (clock() - st) / CLOCKS_PER_SEC);
    printf("maxPlen %d\n",  maxPLen);
    if(maxPLen < minPQ) {
        printf("count 0\n");
        return ;
    }
    // for(int i = 1; i < maxPLen && i < minPQ; i++) {
    //     printf("sumW:%.0f ", sumW[i]);
    // }
    printf("sumW:%.0f ", sumW);
    printf("\n\n");fflush(stdout);
    
    std::default_random_engine generator;
    std::uniform_real_distribution<double> uiDistribution(0, 1);
    std::vector<int> stackL(g->maxDv+1), stackR(g->maxDu+1);
    std::vector<double> sumCXi(g->maxDu + 5), sumCYi(g->maxDv + 5);
    double ansTmp = 0.0;
    
    if(sumW > 0)
    for(auto u:nodes) {
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
            if(maxPLength < minPQ) continue;
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

            for(int len = minPQ; len <= minPQ && len <= maxPLength; len++) {//(len-1,len-1)-bipath
                double sumWtemp = 0.0;
                for(int j = 0; j < pU[pL]; j++) {
                    sumWtemp += dpU[len - 1][j];
                }

                uint64_t sampleSize = std::ceil(T*sumWtemp / sumW);
// if(sampleSize > T) {
//     printf("%d %llu %.0f %.0f\n", len, sampleSize, sumWtemp, sumW[len - 1]);
//     fflush(stdout);
// }
// assert(sampleSize <= T);
                if(sampleSize == 0) continue;

                std::fill(sumCXi.begin(), sumCXi.begin() + pR + 2, 0.0);
                std::fill(sumCYi.begin(), sumCYi.begin() + pL + 2, 0.0);
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

                    if(q >= p) {
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
                        int x = std::max(p,q)-len;

                        if(rSize >= x)
                            sumCXi[x] += C[rSize][x];
                    }
                    else {
                        int lSize = 0;
                        int kk = 0;
                        for(int j = 0; j < pL; j++) {
    // if(u == 3 && v == 4 && candL[stackL[0]] == 4 && candR[stackR[0]] == 5) {
    //     printf("x:%d %d %u\n", j, stackL[kk], candL[stackL[kk]]);
    // }
                            if(kk < len - 1 && stackL[kk] == j) {
                                kk++;
                                continue;
                            }
    // if(u == 3 && v == 4 && candL[stackL[0]] == 4 && candR[stackR[0]] == 5) {
    //     printf("x:%d\n", j);
    // }
                            int u = candL[j];
                            bool f = true;
                            for(int k = 0; k < len - 1; k++) {
                                if(!g->connectUV(u, candR[stackR[k]])) {
                                    f = false;
                                    break;
                                }
                            }
                            if(f) lSize++;
                            // if(pR - j - 1 + rSize < maxPQ - minPQ) break;
                        }

                        int x = std::max(p,q)-len;

                        if(lSize >= x)
                            sumCYi[x] += C[lSize][x];
                        // for(int x = 1; x <= lSize && x < minPQ; x++) {
                        //     sumCYi[x] += C[lSize][x];
                        // }
                    }

                    // for(int x = 0; x <= rSize && x < minPQ; x++) {
                    //     sumCXi[x] += C[rSize][x];
                    // }
                    

                    
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

                sampleSize = std::ceil(T * sumWtemp / sumW);
                // if(ansAll[len].size() < pR) ansAll[len].resize((pR + 2)*2);
                int x = std::max(p,q)-len;

                if(q >= p) {
                    if(sumCXi[x] < 0.5) continue;
                    ansTmp += sumCXi[x] * sumWtemp / sampleSize / C[std::max(p,q) - 1][len - 1];
                }
                else {
                    if(sumCYi[x] < 0.5) continue;
                    ansTmp += sumCYi[x] * sumWtemp / sampleSize / C[std::max(p,q) - 1][len - 1];
                }
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
    // ans += ansTmp / C[std::max(p,q) - 1][minPQ - 1];
    ans += ansTmp;
    printf("%d-%d: %.0f\n", p, q, ans);
    fflush(stdout);

    delete [] bufferForDpU;
    delete [] bufferForDpV;
    delete [] dpU;
    delete [] dpV;
}

bool bcAndPath::costEstimate() {
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

void bcAndPath::exactCount(std::vector<uint32_t> & nodes) {
    
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

    // layerBasedListing(0, g->n1, g->n2);
    int pH = g->n1, l = 0;
    for(auto u: nodes) {
        if(H.lists[u].size() < uint32_t(p - l - 1)) {
            continue;
        }

        int pSS = 0;
        if(l == 0) {
            for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
                S.changeTo(g->e1[i], pSS++);
            }
        }

        if(pSS < q) continue;
        
        int pHH = 0;

        for(int i = 0; i < H.d[0][u]; i++) {
            auto v = H.lists[u][i];
            if(H.nodes.idx(v) < (uint32_t)pH) {
                H.nodes.changeTo(v, pHH++);
            }
        }

        // if(1 < p)
        for(int i = 0; i < pHH; i++) {
            // uint32_t u = H.nodes[l + 1][i];
            uint32_t u = H.nodes[i];
            int d = H.d[0][u];
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
            H.d[1][u] = d;
        }

        layerBasedListing(1, pHH, pSS);
    }

    printf("ans %.2f\n", ans);
    fflush(stdout);
}

void bcAndPath::layerBasedListing(int l, int pH, int pS) {
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

// auto t1 = std::chrono::steady_clock::now();

    for(int j = 0; j < pH; j++) {
// if(l == 0 && j % 200 == 0) {
//     auto t2 = std::chrono::steady_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
//     printf("%u, %lld\n", j/10000, (long long)duration.count());fflush(stdout);
//     t1 = t2;
// }

        uint32_t u = tmpNodes[l][j];
        // uint32_t u = H.nodes[l][j];
        // if(l == 0)
        if(H.lists[u].size() < uint32_t(p - l - 1)) {
            continue;
        }

        int pSS = 0;

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            if(S.idx(g->e1[i]) < (uint32_t)pS) {
                S.changeTo(g->e1[i], pSS++);
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

void bcAndPath::collect2HopNeighbors() {
    H.lists.resize(g->n1);

    std::vector<uint32_t> sum(g->n1);
    std::vector<uint32_t> tmp(g->n1);
    uint32_t l = 0;

// double twotwo = 0;
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
// twotwo += C[sum[w]][q];
                H.lists[u].push_back(w);
            }
            sum[w] = 0;
        }
        l = 0;
    }
// printf("2-2 clique %.0f\n", twotwo);
}


void bcAndPath::collect2HopNeighborsContainedCDAG() {
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
