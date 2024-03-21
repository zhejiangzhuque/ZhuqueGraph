#include "turan.h"
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>

void turan::sample(uint64_t T) {
    std::vector<uint32_t> L(g->m), R(g->m);
    double sum = 0.0;
    double ans = 0.0;
// g->print();
// printf("pq: %d %d\n", p, q);

    for(uint32_t u = 0; u < g->n1; u++) {
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            auto st = g->e2.begin() + g->pV[v];
            auto ed = g->e2.begin() + g->pV[v + 1];
            uint32_t pL = ed - std::upper_bound(st, ed, u);
            uint32_t pR = g->pU[u + 1] - i - 1;

            if(pL < q - 1) continue;
            if(pR < p - 1) continue;

            L[i] = pL;
            R[i] = pR;
            sum += C[pL][p - 1] * C[pR][q - 1];
// printf("i %u, u:%u, v:%u, %u %u, %.0f*%.0f=%.0f\n", i, u, v, L[i], R[i],
//     C[pL][p - 1], C[pR][q - 1], C[pL][p - 1] * C[pR][q - 1]);
        }
    }
    printf("sum %.0f\n", sum);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> uiDistribution(0, 1);
    std::vector<uint32_t> stackL(g->maxDv), stackR(g->maxDu);
    std::vector<bool> visL(g->maxDv), visR(g->maxDu);
    uint64_t totalCnt = 0, totalT = 0;

    for(uint32_t u = 0; u < g->n1; u++) {
// if(u % 1000 == 0) {
//     printf("%u\n", u);fflush(stdout);
// }

        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }
        if(pR < q) continue;

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, --pR);
            if(pR < q - 1) continue;

            int pL = 0;

            auto st = g->e2.begin() + g->pV[v];
            auto ed = g->e2.begin() + g->pV[v + 1];
            auto mid = std::upper_bound(st, ed, u);

            if(ed - mid + 1 < p) continue;

            uint32_t j = mid - g->e2.begin();
            for(; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                candL.changeTo(w, pL++);
            }
            assert(pL >= p - 1);

            uint64_t sampleSize = std::ceil(T*(C[pL][p - 1] * C[pR][q - 1]) / sum);
            totalT += sampleSize;
            uint64_t cnt = 0;
// if(u > 1811000) {
    // printf("%llu\n", sampleSize);fflush(stdout);
// }
            for(uint64_t j = 0; j < sampleSize; j++) {
                for(int k = 0; k < p - 1; k++) {
                    int x = (pL - k) * uiDistribution(generator);
                    while(visL[x]) {
                        x = (x + 1) % pL;
// printf("there\n"); fflush(stdout);
                    }
                    visL[x] = true;

                    stackL.push_back(x);
                }

                bool f = true;
                for(int k = 0; k < q - 1; k++) {
                    int x = (pR - k) * uiDistribution(generator);
                    while(visR[x]) {
                        x = (x + 1) % pR;
                    }
                    visR[x] = true;
                    stackR.push_back(x);
                    
                    for(int h = 0; h < p - 1; h++) {
                        if(!g->connectVU(candR[x], candL[stackL[h]])) {
                            f = false;
                            break;
                        }
                    }
                    if(!f) {
                        break;
                    }
                }

                if(f) {
                    cnt++;
                }

                for(int k = 0; k < stackL.size(); k++) visL[stackL[k]] = false;
                stackL.clear();
                for(int k = 0; k < stackR.size(); k++) visR[stackR[k]] = false;
                stackR.clear();
            }

            ans += (1.0*cnt / sampleSize) * C[pL][p - 1] * C[pR][q - 1];
// printf("i %u, u:%u, v:%u, %.0f*%.0f=%.0f, %llu/%llu=%.2f, ans %.2f\n", i, u, v, 
//     C[pL][p - 1], C[pR][q - 1], C[pL][p - 1] * C[pR][q - 1], cnt, sampleSize, (1.0*cnt / sampleSize), ans);
// for(int i = 0; i < pL; i++) {
//     printf("%u ", candL[i]);
// }printf("\n");
// for(int i = 0; i < pR; i++) {
//     printf("%u ", candR[i]);
// }printf("\n");

            totalCnt += cnt;
        }
    }

    printf("%.0f\n", 1.0 * totalCnt / totalT * sum);
    printf("%.0f\n", ans);
}