#include "colorPath.h"
#include <chrono>
#include <random>
#include <cassert>

void colorPath::testSubgraphSize() {
    candL.resize(g->n1);
    candR.resize(g->n2);

    uint32_t maxPe = 0;
    for(uint32_t u = 0; u < g->n1; u++) {
        int pR = 0;
        uint32_t edges = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];

            int pL = 0;
            auto st = g->e2.begin() + g->pV[v];
            auto ed = g->e2.begin() + g->pV[v + 1];
            uint32_t j = ed - std::upper_bound(st, ed, u);
            edges += j;
        }

        maxPe = std::max(maxPe, edges);

        printf("%u:%u\n", u, edges);
    }

    printf("%u\n", maxPe);
}



//approximate p<q
void colorPath::approximateCounting3(uint64_t T) {
    if(p > q) {
        g->swapUV();
        std::swap(p, q);
    }

    auto t1 = std::chrono::steady_clock::now();
    g->coreReduction(p, q);
    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "core reduction time:" << duration.count() << "ms" << std::endl;
    printf("%u %u %u\n", g->n1, g->n2, g->m);

    
    //p <= q
// g->print();

    printf("maxDu:%u maxDv:%u\n", g->maxDu, g->maxDv);
    uint32_t maxDuv = std::max(g->maxDu, g->maxDv);
    uint32_t maxE = g->maxDu * g->maxDv;

    double * bufferForDpU = new double[maxE * p];
    double * bufferForDpV = new double[maxE * p];

    double ** dpU = new double *[maxE];
    for(uint32_t e = 0; e < maxE; e++) {
        dpU[e] = bufferForDpU + e * p;
    }
    for(uint32_t e = 0; e < maxE; e++) {
        dpU[e][0] = 0;
        dpU[e][1] = 1;
    }

    double ** dpV = new double*[maxE];
    for(uint32_t e = 0; e < maxE; e++) {
        dpV[e] = bufferForDpV + e * p;
    }
    for(uint32_t e = 0; e < maxE; e++) {
        dpV[e][0] = 1;
    }

    candL.resize(g->n1);
    candR.resize(g->n2);
    // uint32_t subGraphEcount = std::min(g->m, g->maxDu * g->maxDv);
    std::vector<int> pU(g->maxDv+1), e1(maxE), pV(g->maxDu+1), e2(maxE);
    std::vector<int> mapUtoV(maxE), mapVtoU(maxE);
    std::vector<double> P(g->m);
    double sumW = 0.0;

    //first compute dp
// int xxx = 0;
    auto computeDP = [&](int pL, int pR) {
// printf("there %d\n", xxx);fflush(stdout);
        std::fill(pV.begin(), pV.begin() + pR + 1, 0);
        for(int j = 0; j < pL; j++) {
            uint32_t x = candL[j];
            pU[j + 1] = pU[j];
            
            for(int k = 0; k < pR; k++) {
                uint32_t y = candR[k];
                if(g->connectUV(x, y)) {
                    e1[pU[j + 1]++] = k;
                    pV[k + 1]++;
                }
            }
        }

        for(int v = 0; v < pR; v++) {
            pV[v + 1] += pV[v];
        }
        assert(pU[pL] == pV[pR]);
// printf("there %d\n", xxx);fflush(stdout);
        for(int u = 0; u < pL; u++) {
            for(int i = pU[u]; i < pU[u + 1]; i++) {
                int v = e1[i];

                e2[pV[v]] = u;

                mapUtoV[i] = pV[v];
                mapVtoU[pV[v]] = i;

                pV[v]++;
            }
        }

// printf("there %d\n", xxx++);fflush(stdout);
        for(int v = pR; v >= 1; v--) {
            pV[v] = pV[v - 1];
        }
        pV[0] = 0;

// for(int u = 0; u < pL; u++) {
//     printf("u %d:", candL[u]);
//     for(int i = pU[u]; i < pU[u + 1]; i++) {
//         printf("%d ", candR[e1[i]]);
//     }
//     printf("\n");
// }
// for(int v = 0; v < pR; v++) {
//     printf("v %d:", candR[v]);
//     for(int i = pV[v]; i < pV[v + 1]; i++) {
//         printf("%d ", candL[e2[i]]);
//     }
//     printf("\n");
// }

       
// printf("there %d\n", xxx++);fflush(stdout);
// if(903 == xxx) {
// printf("pU:%d pV:%d\n", pU[pL], pV[pR]);
// fflush(stdout);
// }    
        for(int v = 0; v < pR; v++) {
            for(int i = pV[v]; i < pV[v + 1]; i++) {
                int u = e2[i];
                int h = std::upper_bound(e1.begin() + pU[u], 
                                        e1.begin() + pU[u + 1], v) - e1.begin();
                dpV[i][1] = pU[u + 1] - h;
            }
        }
        for(int k = 2; k < p; k++) {
            for(int u = 0; u < pL; u++) {
                for(int i = pU[u]; i < pU[u + 1]; i++) {
                    int v = e1[i];
                    int j = std::upper_bound(e2.begin() + pV[v], e2.begin() + pV[v + 1], u) - e2.begin();
                    if(e2[j] == u) j++;
                    dpU[i][k] = 0;
                    
                    for(; j < pV[v + 1]; j++) {
                        dpU[i][k] += dpV[j][k - 1];
                    }
                }
            }

            for(int v = 0; v < pR; v++) {
                for(int i = pV[v]; i < pV[v + 1]; i++) {
// if(903 == xxx) {
// printf("k:%d v:%d, i:%d\n", k, v, i);
// fflush(stdout);
// }
                    int u = e2[i];
                    int j = std::upper_bound(e1.begin() + pU[u], e1.begin() + pU[u + 1], v) - e1.begin();
          
                    if(e1[j] == v) j++;
                    dpV[i][k] = 0;

                    for(; j < pU[u + 1]; j++) {
                        dpV[i][k] += dpU[j][k];
                    }
                }
            }

            
        }
    };

    for(uint32_t u = 0; u < g->n1; u++) {
// if(u % 1000 == 0) {
//     printf("u%1000:%u\n", u);fflush(stdout);
// }
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];

            int pL = 0;
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) candL.changeTo(w, pL++);
            }
            candR.changeTo(v, --pR);
            for(int i = 1; i < pR; i++) {
                candR.changeToByPos(i, i - 1);
            }
            // candR.popFirst();
            // --pR;

            if(pL < p-1 || pR < q-1) continue;

            computeDP(pL, pR);
// printf("u:%u v:%u\n", u, v);
// printf("candL Nv:");
// for(int i = 0; i < pL; i++) {
//     printf("%d ", candL[i]);
// }
// printf("\ncandR Nu:");
// for(int i = 0; i < pR; i++) {
//     printf("%d ", candR[i]);
// }
// printf("\n dpU:\n");
// for(int i = 0; i < pL; i++) for(int k = pU[i]; k < pU[i + 1]; k++) {
//     printf("%d,%d:", candL[i], candR[e1[k]]);
//     for(int j = 0; j < p; j++) {
//         printf("%.0f ", dpU[k][j]);
//     }
//     printf("\n");
// }
// printf("\n dpV:\n");
// for(int i = 0; i < pR; i++) for(int k = pV[i]; k < pV[i + 1]; k++) {
//     printf("%d,%d:", candR[i], candL[e2[k]]);
//     for(int j = 0; j < p; j++) {
//         printf("%.0f ", dpV[k][j]);
//     }
//     printf("\n");
// }

            for(int j = 0; j < pU[pL]; j++) {
                P[i] += dpU[j][p - 1];
            }
            // P[i] = dpU[][p - 1];
            sumW += P[i];
// printf("P[i], %u %u, v, %.0f\n", u, v, P[i]);
        }

        // sumW += P[u];
    }
printf("first %.2f\n", sumW);fflush(stdout);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> uiDistribution(0, 1);
    std::vector<int> stackL(p), stackR(p);
    double sumCXi = 0.0;

    double realT = 0, ans2 = 0.0;

    if(sumW > 0)
    for(uint32_t u = 0; u < g->n1; u++) {
// if(u % 1000 == 0) {
//     printf("%u\n", u);fflush(stdout);
// }

        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];

            int pL = 0;
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) candL.changeTo(w, pL++);
            }
            candR.changeTo(v, --pR);
            for(int i = 1; i < pR; i++) {
                candR.changeToByPos(i, i - 1);
            }

            if(pL < p-1 || pR < q-1) continue;

            uint64_t sampleSize = std::ceil(P[i] / sumW * T);

            if(sampleSize == 0) continue;

            computeDP(pL, pR);
// printf("dp, L %d\n", pU[pL + 1]);fflush(stdout);

// printf("\n dpU:\n");
// for(int i = 0; i < pL; i++) for(int k = pU[i]; k < pU[i + 1]; k++) {
//     printf("%d,%d:", candL[i], candR[e1[k]]);
//     for(int j = 0; j < p; j++) {
//         printf("%.0f ", dpU[k][j]);
//     }
//     printf("\n");
// }
            double sumWtemp = P[i];
// printf("spsize %llu, sumTmp %.0f\n", sampleSize, sumWtemp);
// fflush(stdout);
        
            realT += sampleSize;
            double preSumCXi = sumCXi;
// printf("u:%u v:%u\n", u, v);
// printf("candL Nv:");
// for(int i = 0; i < pL; i++) {
//     printf("%d ", candL[i]);
// }
// printf("\ncandR Nu:");
// for(int i = 0; i < pR; i++) {
//     printf("%d ", candR[i]);
// }
            
            while(sampleSize--) {
// printf("sp %d\n", sampleSize + 1);fflush(stdout);
                double tmp = 0.0;
                double r = uiDistribution(generator); 
                int preE = 0, preU = 0, preV = 0;
                stackL.clear();
                stackR.clear();

                for(int u = 0; u < pL; u++) {
                    for(int i = pU[u]; i < pU[u + 1]; i++) {
                        tmp += dpU[i][p - 1];
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

                for(int i = 1; i < p - 1; i++) {
                    double r = uiDistribution(generator);
                    tmp = 0.0;
                
                    for(int j = pV[preV]; j < pV[preV + 1]; j++) {
                        if(e2[j] <= preU) continue;
                        tmp += dpV[j][p - i - 1];
                        if(tmp + 1e-8 >= r * dpU[preE][p - i]) {
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
                        tmp += dpU[j][p - i - 1];
                        if(tmp + 1e-8 >= r * dpV[preE][p - i - 1]) {
                            int v = e1[j];
                            stackR.push_back(v);
                            preV = v;
                            preE = j;
                            break;
                        }
                    }
                }

                // int x = (pU[preU+1] - pU[preU])*(uiDistribution(generator)+1e-8);
                // stackR.push_back(e1[x + pU[preU]]);

                assert(stackL.size() == p - 1);
                // if(stackL.size() < p - 1) {
                //     realT--;

                //     continue;
                // }

                bool connect = true;
                for(int i = 0; i < p - 1; i++) {
                    for(int j = 0; j < p - 1; j++) {
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
// printf("L:");
// for(int i = 0; i < p - 1; i++) printf("%d ", candL[stackL[i]]);
// printf("R:");
// for(int i = 0; i < p - 1; i++) printf("%d ", candR[stackR[i]]);
// printf("\n");

                int rSize = 0, kk = 0;
                for(int j = 0; j < pR; j++) {
                    int v = candR[j];
                    if(stackR[kk] == j) {
                        kk++;
                        continue;
                    }
                    bool f = true;
                    for(int k = 0; k < p - 1; k++) {
                        if(!g->connectUV(candL[stackL[k]], v)) {
                            f = false;
                            break;
                        }
                    }

                    if(f) rSize++;
                    if(pR - j - 1 + rSize < q - p) break;
                }

                if(rSize >= q - p) {
                    sumCXi += C[rSize][q - p];
// printf("add %.0f \n", C[rSize][q - p]);
                }
            }
            
            ans2 += (sumCXi - preSumCXi) / std::ceil(P[i] / sumW * T) * P[i] /C[q-1][p-1];
// printf("est %d, %u %u, sumCXi %.0f, %.0f\n", i, u, v, (sumCXi - preSumCXi), 
//     (sumCXi - preSumCXi) / std::ceil(P[i] / sumW * T)/ C[q-1][p-1]*P[i]);
        }
// printf("ed\n");fflush(stdout);
    }
printf("ed, sumCXi %.0f, T %.0f\n", sumCXi, realT);fflush(stdout);

    double ans = sumCXi / realT / C[q-1][p-1] * sumW;
    // ans2 /= C[q-1][p-1];

    printf("extimate answer1:%.0f, ans2:%.0f\n", ans, ans2);
    fflush(stdout);

    delete [] bufferForDpU;
    delete [] bufferForDpV;
    delete [] dpU;
    delete [] dpV;
}

//approximate input p-x, x>=q, input p only
void colorPath::approximateCounting4(uint64_t T) {
    auto t1 = std::chrono::steady_clock::now();
    g->coreReduction(p, q);
    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "core reduction time:" << duration.count() << "ms" << std::endl;
    printf("%u %u %u\n", g->n1, g->n2, g->m);
// g->print();
    printf("maxDu:%u maxDv:%u\n", g->maxDu, g->maxDv);
    uint32_t maxDuv = std::max(g->maxDu, g->maxDv);
    uint32_t maxE = std::min((1u<<23), maxDuv * 200);
    int minPQ = std::min(p, q);
    int maxPQ = std::max(p, q);

    double * bufferForDpU = new double[maxE * minPQ];
    double * bufferForDpV = new double[maxE * minPQ];

    double ** dpU = new double *[maxE];
    for(uint32_t e = 0; e < maxE; e++) {
        dpU[e] = bufferForDpU + e * minPQ;
    }
    for(uint32_t e = 0; e < maxE; e++) {
        dpU[e][0] = 0;
        dpU[e][1] = 1;
    }

    double ** dpV = new double*[maxE];
    for(uint32_t e = 0; e < maxE; e++) {
        dpV[e] = bufferForDpV + e * minPQ;
    }
    for(uint32_t e = 0; e < maxE; e++) {
        dpV[e][0] = 1;
    }

    candL.resize(g->n1);
    candR.resize(g->n2);

    ansAll.resize(maxDuv);

    //subgraph
    std::vector<int> pU(g->maxDv+1), e1(maxE), pV(g->maxDu+1), e2(maxE);
    std::vector<int> mapUtoV(maxE), mapVtoU(maxE);
    //(minPQ-minPQ)-bipath
    std::vector<double> P(g->m);
    double sumW = 0.0;

    //first compute dp
// int xxx = 0;
    auto computeDP = [&](int pL, int pR) {
// printf("there %d\n", xxx);fflush(stdout);
        std::fill(pV.begin(), pV.begin() + pR + 1, 0);
        for(int j = 0; j < pL; j++) {
            uint32_t x = candL[j];
            pU[j + 1] = pU[j];
            
            for(int k = 0; k < pR; k++) {
                uint32_t y = candR[k];
                if(g->connectUV(x, y)) {
                    e1[pU[j + 1]++] = k;
                    pV[k + 1]++;
                }
            }
        }
        assert(pU[pL] < maxE);

        for(int v = 0; v < pR; v++) {
            pV[v + 1] += pV[v];
        }
        assert(pU[pL] == pV[pR]);
// printf("there %d\n", xxx);fflush(stdout);
        for(int u = 0; u < pL; u++) {
            for(int i = pU[u]; i < pU[u + 1]; i++) {
                int v = e1[i];

                e2[pV[v]] = u;

                mapUtoV[i] = pV[v];
                mapVtoU[pV[v]] = i;

                pV[v]++;
            }
        }

// printf("there %d\n", xxx++);fflush(stdout);
        for(int v = pR; v >= 1; v--) {
            pV[v] = pV[v - 1];
        }
        pV[0] = 0;

        for(int v = 0; v < pR; v++) {
            for(int i = pV[v]; i < pV[v + 1]; i++) {
                int u = e2[i];
                int h = std::upper_bound(e1.begin() + pU[u], 
                                        e1.begin() + pU[u + 1], v) - e1.begin();
                dpV[i][1] = pU[u + 1] - h;
            }
        }
        for(int k = 2; k < minPQ; k++) {
            for(int u = 0; u < pL; u++) {
                for(int i = pU[u]; i < pU[u + 1]; i++) {
                    int v = e1[i];
                    int j = std::upper_bound(e2.begin() + pV[v], e2.begin() + pV[v + 1], u) - e2.begin();
                    if(e2[j] == u) j++;
                    dpU[i][k] = 0;
                    
                    for(; j < pV[v + 1]; j++) {
                        dpU[i][k] += dpV[j][k - 1];
                    }
                }
            }

            for(int v = 0; v < pR; v++) {
                for(int i = pV[v]; i < pV[v + 1]; i++) {
// if(903 == xxx) {
// printf("k:%d v:%d, i:%d\n", k, v, i);
// fflush(stdout);
// }
                    int u = e2[i];
                    int j = std::upper_bound(e1.begin() + pU[u], e1.begin() + pU[u + 1], v) - e1.begin();
          
                    if(e1[j] == v) j++;
                    dpV[i][k] = 0;

                    for(; j < pU[u + 1]; j++) {
                        dpV[i][k] += dpU[j][k];
                    }
                }
            }

            
        }
    };

    for(uint32_t u = 0; u < g->n1; u++) {
// if(u % 1000 == 0) {
//     printf("u%1000:%u\n", u % 1000);fflush(stdout);
// }
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];

            int pL = 0;
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) candL.changeTo(w, pL++);
            }
            candR.changeTo(v, --pR);
            for(int i = 1; i < pR; i++) {
                candR.changeToByPos(i, i - 1);
            }
            // candR.popFirst();
            // --pR;

            if(pL < minPQ - 1 || pR < minPQ - 1) continue;

            computeDP(pL, pR);
// printf("u:%u v:%u\n", u, v);
// printf("candL Nv:");
// for(int i = 0; i < pL; i++) {
//     printf("%d ", candL[i]);
// }
// printf("\ncandR Nu:");
// for(int i = 0; i < pR; i++) {
//     printf("%d ", candR[i]);
// }
// printf("\n dpU:\n");
// for(int i = 0; i < pL; i++) for(int k = pU[i]; k < pU[i + 1]; k++) {
//     printf("%d,%d:", candL[i], candR[e1[k]]);
//     for(int j = 0; j < p; j++) {
//         printf("%.0f ", dpU[k][j]);
//     }
//     printf("\n");
// }
// printf("\n dpV:\n");
// for(int i = 0; i < pR; i++) for(int k = pV[i]; k < pV[i + 1]; k++) {
//     printf("%d,%d:", candR[i], candL[e2[k]]);
//     for(int j = 0; j < p; j++) {
//         printf("%.0f ", dpV[k][j]);
//     }
//     printf("\n");
// }

            for(int j = 0; j < pU[pL]; j++) {
                P[i] += dpU[j][minPQ - 1];
            }

            sumW += P[i];
// printf("P[i], %u %u, v, %.0f\n", u, v, P[i]);
        }
    }
    
    printf("first %.2f\n", sumW);fflush(stdout);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> uiDistribution(0, 1);
    std::vector<int> stackL(minPQ+1), stackR(minPQ+1);
    std::vector<double> sumCXi(g->maxDu);

    double realT = 0;
    std::vector<double> ans2(maxDuv);

    if(sumW > 0)
    for(uint32_t u = 0; u < g->n1; u++) {
// if(u % 1000 == 0) {
//     printf("%u\n", u);fflush(stdout);
// }

        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];

            int pL = 0;
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) candL.changeTo(w, pL++);
            }
            candR.changeTo(v, --pR);
            for(int i = 1; i < pR; i++) {
                candR.changeToByPos(i, i - 1);
            }

            if(pL < minPQ - 1 || pR < minPQ - 1) continue;

            uint64_t sampleSize = std::ceil(P[i] / sumW * T);

            if(sampleSize == 0) continue;

            computeDP(pL, pR);
// printf("dp, L %d\n", pU[pL + 1]);fflush(stdout);

// printf("\n dpU:\n");
// for(int i = 0; i < pL; i++) for(int k = pU[i]; k < pU[i + 1]; k++) {
//     printf("%d,%d:", candL[i], candR[e1[k]]);
//     for(int j = 0; j < p; j++) {
//         printf("%.0f ", dpU[k][j]);
//     }
//     printf("\n");
// }
            double sumWtemp = P[i];
// printf("spsize %llu, sumTmp %.0f\n", sampleSize, sumWtemp);
// fflush(stdout);
        
            realT += sampleSize;
            // double preSumCXi = sumCXi;
            std::fill(sumCXi.begin(), sumCXi.begin() + pR + 1, 0.0);
// printf("u:%u v:%u\n", u, v);
// printf("candL Nv:");
// for(int i = 0; i < pL; i++) {
//     printf("%d ", candL[i]);
// }
// printf("\ncandR Nu:");
// for(int i = 0; i < pR; i++) {
//     printf("%d ", candR[i]);
// }
            
            while(sampleSize--) {
// printf("sp %d\n", sampleSize + 1);fflush(stdout);
                double tmp = 0.0;
                double r = uiDistribution(generator); 
                int preE = 0, preU = 0, preV = 0;
                stackL.clear();
                stackR.clear();

                for(int u = 0; u < pL; u++) {
                    for(int i = pU[u]; i < pU[u + 1]; i++) {
                        tmp += dpU[i][minPQ - 1];
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

                for(int i = 1; i < minPQ - 1; i++) {
                    double r = uiDistribution(generator);
                    tmp = 0.0;
                
                    for(int j = pV[preV]; j < pV[preV + 1]; j++) {
                        if(e2[j] <= preU) continue;
                        tmp += dpV[j][minPQ - i - 1];
                        if(tmp + 1e-8 >= r * dpU[preE][minPQ - i]) {
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
                        tmp += dpU[j][minPQ - i - 1];
                        if(tmp + 1e-8 >= r * dpV[preE][minPQ - i - 1]) {
                            int v = e1[j];
                            stackR.push_back(v);
                            preV = v;
                            preE = j;
                            break;
                        }
                    }
                }

                // int x = (pU[preU+1] - pU[preU])*(uiDistribution(generator)+1e-8);
                // stackR.push_back(e1[x + pU[preU]]);

                assert(stackL.size() == minPQ - 1);
                // if(stackL.size() < p - 1) {
                //     realT--;
                //     continue;
                // }

                bool connect = true;
                for(int i = 0; i < minPQ - 1; i++) {
                    for(int j = 0; j < minPQ - 1; j++) {
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
// printf("L:");
// for(int i = 0; i < p - 1; i++) printf("%d ", candL[stackL[i]]);
// printf("R:");
// for(int i = 0; i < p - 1; i++) printf("%d ", candR[stackR[i]]);
// printf("\n");

                int rSize = 0, kk = 0;
                for(int j = 0; j < pR; j++) {
                    int v = candR[j];
                    if(stackR[kk] == j) {
                        kk++;
                        continue;
                    }
                    bool f = true;
                    for(int k = 0; k < minPQ - 1; k++) {
                        if(!g->connectUV(candL[stackL[k]], v)) {
                            f = false;
                            break;
                        }
                    }

                    if(f) rSize++;
                    // if(pR - j - 1 + rSize < maxPQ - minPQ) break;
                }

                for(int x = 0; x <= rSize; x++) {
                    sumCXi[x] += C[rSize][x];
                }
                // if(rSize >= maxPQ - minPQ) {
                //     sumCXi[0] += C[rSize][maxPQ - minPQ];
                // }
            }

            sampleSize = std::ceil(P[i] / sumW * T);
            for(int x = 0; minPQ + x <= pR; x++) {
                ans2[x] += sumCXi[x] * P[i] / sampleSize / C[minPQ + x - 1][minPQ - 1];
            }
            // ans2 += (sumCXi - preSumCXi) / std::ceil(P[i] / sumW * T) * P[i] /C[maxPQ-1][minPQ-1];

// printf("est %d, %u %u, sumCXi %.0f, %.0f\n", i, u, v, (sumCXi - preSumCXi), 
//     (sumCXi - preSumCXi) / std::ceil(P[i] / sumW * T)/ C[q-1][p-1]*P[i]);
        }
    }
    // printf("ed, sumCXi %.0f, T %.0f\n", sumCXi, realT);fflush(stdout);

    // double ans = sumCXi / realT / C[maxPQ-1][minPQ-1] * sumW;
    // ans2 /= C[q-1][p-1];

    // printf("extimate answer1:%.0f, ans2:%.0f\n", ans, ans2);
    for(int x = 0; x < ans2.size(); x++) {
        if(ans2[x] < 0.5) break;
        printf("%d-%d: %.1f\n", minPQ, minPQ + x, ans2[x]);
    }
    fflush(stdout);

    delete [] bufferForDpU;
    delete [] bufferForDpV;
    delete [] dpU;
    delete [] dpV;
}


void colorPath::approximateCountingAll(uint64_t T) {
// g->print();
    printf("maxDu:%u maxDv:%u\n", g->maxDu, g->maxDv);
    fflush(stdout);
    uint32_t maxDuv = std::max(g->maxDu, g->maxDv);

    const int minPQ = 100;
    uint32_t maxE = std::min((1u<<20), g->maxDu * g->maxDv);
    // uint32_t maxE = 105000;

    candL.resize(g->n1);
    candR.resize(g->n2);

// double st = clock();
//     for(uint32_t u = 0; u < g->n1; u++) {
// printf("%u\n", u);fflush(stdout);
//         int pR = 0;
//         for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//             uint32_t v = g->e1[i];
//             candR.changeTo(v, pR++);
//         }

//         for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//             uint32_t v = g->e1[i];

//             int pL = 0;
//             for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
//                 uint32_t w = g->e2[j];
//                 if(w > u) candL.changeTo(w, pL++);
//             }
//             candR.changeTo(v, --pR);

//             uint32_t sumD = 0;
//             for(int j = 0; j < pL; j++) {
//                 uint32_t x = candL[j];
//                 uint32_t k = std::upper_bound(g->e1.begin() + g->pU[x], 
//                     g->e1.begin() + g->pU[x + 1], v) - g->e1.begin();
                
//                 for( ; k < g->pU[x + 1]; k++) {
//                     if(candR.idx(g->e1[k]) < pR) {
//                         sumD++;
//                     }
//                 }
//                 // for(int k = 0; k < pR; k++) {
//                 //     uint32_t y = candR[k];
//                 //     if(g->connectUV(x, y)) {
//                 //         sumD++;
//                 //     }
//                 // }
//             }
//             maxE = std::max(maxE, sumD);
//         }
//     }
//     printf("maxE:%u\n", maxE);fflush(stdout);
// printf("get maxE:%f\n", (clock() - st) / CLOCKS_PER_SEC);

    double * bufferForDpU = new double[maxE * minPQ];
    double * bufferForDpV = new double[maxE * minPQ];

    double ** dpU = new double *[maxE];
    for(uint32_t e = 0; e < maxE; e++) {
        dpU[e] = bufferForDpU + e * minPQ;
    }
    for(uint32_t e = 0; e < maxE; e++) {
        dpU[e][0] = 0;
        dpU[e][1] = 1;
    }

    double ** dpV = new double*[maxE];
    for(uint32_t e = 0; e < maxE; e++) {
        dpV[e] = bufferForDpV + e * minPQ;
    }
    for(uint32_t e = 0; e < maxE; e++) {
        dpV[e][0] = 1;
    }

    ansAll.resize(minPQ + 1);
    for(uint32_t i = 0; i <= minPQ; i++) {
        ansAll[i].resize(minPQ + 1);
    }

    //subgraph
    std::vector<int> pU(g->maxDv+1), e1(maxE), pV(g->maxDu+1), e2(maxE);
    std::vector<int> mapUtoV(maxE), mapVtoU(maxE);
    //(minPQ-minPQ)-bipath
    std::vector<double> sumW(minPQ);

    //first compute dp
printf("start dp first\n");fflush(stdout);
// uint32_t realMaxE = 0;

    auto computeDP = [&](int pL, int pR) {
        std::fill(pV.begin(), pV.begin() + pR + 1, 0);
        for(int j = 0; j < pL; j++) {
            uint32_t x = candL[j];
            pU[j + 1] = pU[j];
            
            for(int k = 0; k < pR; k++) {
                uint32_t y = candR[k];
                if(g->connectUV(x, y)) {
                    e1[pU[j + 1]++] = k;
                    pV[k + 1]++;
                }
            }
        }
        assert(pU[pL] < maxE);

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
// printf("Graph:\n");
// for(int u = 0; u < pL; u++) {
//     printf("u %d:", candL[u]);
//     for(int i = pU[u]; i < pU[u + 1]; i++) {
//         printf("%d ", candR[e1[i]]);
//     }
//     printf("\n");
// }
// for(int v = 0; v < pR; v++) {
//     printf("v %d:", candR[v]);
//     for(int i = pV[v]; i < pV[v + 1]; i++) {
//         printf("%d ", candL[e2[i]]);
//     }
//     printf("\n");
// }


        //init dpV
        // for(int v = 0; v < pR; v++) {
        //     for(int i = pV[v]; i < pV[v + 1]; i++) {
        //         int u = e2[i];
        //         int h = std::upper_bound(e1.begin() + pU[u], 
        //                                 e1.begin() + pU[u + 1], v) - e1.begin();
        //         dpV[i][1] = pU[u + 1] - h;
        //     }
        // }
        for(int u = 0; u < pL; u++) {
            for(int i = pU[u]; i < pU[u + 1]; i++) {
                // int v = e1[i];
                dpV[mapUtoV[i]][1] = pU[u + 1] - i - 1;
            }
        }

        int minLR = std::min(pL, pR);
        int k = 2;
        for(; k <= minLR && k < minPQ; k++) {
            bool f = false;

            for(int u = 0; u < pL; u++) {
                for(int i = pU[u]; i < pU[u + 1]; i++) {
                    dpU[i][k] = 0;

                    int v = e1[i];
                    int j = std::upper_bound(e2.begin() + pV[v], e2.begin() + pV[v + 1], u) - e2.begin();
                    // if(e2[j] == u) j++;
                    
                    for(; j < pV[v + 1]; j++) {
                        dpU[i][k] += dpV[j][k - 1];
                    }

                    if(dpU[i][k] > 0) f = true;
                }
            }

            for(int v = 0; v < pR; v++) {
                for(int i = pV[v]; i < pV[v + 1]; i++) {
                    int u = e2[i];
                    int j = std::upper_bound(e1.begin() + pU[u], e1.begin() + pU[u + 1], v) - e1.begin();
          
                    // if(e1[j] == v) j++;
                    dpV[i][k] = 0;

                    for(; j < pU[u + 1]; j++) {
                        dpV[i][k] += dpU[j][k];
                    }

                    if(dpV[i][k] > 0) f = true;
                }
            }

            if(!f) break;
        }

        return k;
    };

    int maxPLen = 0;
    for(uint32_t u = 0; u < g->n1; u++) {
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];

            int pL = 0;
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) candL.changeTo(w, pL++);
            }
            candR.changeTo(v, --pR);
            for(int i = 1; i < pR; i++) {
                candR.changeToByPos(i, i - 1);
            }

            // if(pL < minPQ - 1 || pR < minPQ - 1) continue;
// printf("Graph: %u %u:\n", u, v);
            int maxPLength = computeDP(pL, pR);
            maxPLen = std::max(maxPLen, maxPLength);

            for(int j = 0; j < pU[pL]; j++) {
                for(int k = 1; k < maxPLength && k < minPQ; k++) {
                    sumW[k] += dpU[j][k];
                }
            }
// printf("sumW2:%.0f\n", sumW[1]);
        }
    }

    for(int i = 1; i < maxPLen; i++) {
        printf("sumW:%.0f ", sumW[i]);
    }
    printf("\n\n");fflush(stdout);
    
    std::default_random_engine generator;
    std::uniform_real_distribution<double> uiDistribution(0, 1);
    std::vector<int> stackL(g->maxDv+1), stackR(g->maxDu+1);
    std::vector<double> sumCXi(g->maxDu + 5), sumCYi(g->maxDv + 5);

    // if(sumW > 0)
    for(uint32_t u = 0; u < g->n1; u++) {
// if(u > 106000) {
//     printf("u:%u\n", u);fflush(stdout);
// }
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];

            int pL = 0;
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) candL.changeTo(w, pL++);
            }
            candR.changeTo(v, --pR);
            for(int i = 1; i < pR; i++) {
                candR.changeToByPos(i, i - 1);
            }


            // if(pL < minPQ - 1 || pR < minPQ - 1) continue;
// printf("Graph: %u %u:\n", u, v);
            int maxPLength = computeDP(pL, pR);

            for(int len = 2; len <= maxPLength && len < minPQ && sumW[len - 1] > 0.5; len++) {//(len-len)-bipath
                double sumWtemp = 0.0;
                for(int j = 0; j < pU[pL]; j++) {
                    sumWtemp += dpU[j][len - 1];
                }

                uint64_t sampleSize = std::ceil(sumWtemp / sumW[len - 1] * T);

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
                            tmp += dpU[i][len - 1];
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
                            tmp += dpV[j][len - i - 1];
                            if(tmp + 1e-8 >= r * dpU[preE][len - i]) {
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
                            tmp += dpU[j][len - i - 1];
                            if(tmp + 1e-8 >= r * dpV[preE][len - i - 1]) {
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
// printf("\n");
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

                    for(int x = 0; x <= rSize && x < minPQ; x++) {
                        sumCXi[x] += C[rSize][x];
                    }

                    int lSize = 0;
                    kk = 0;
                    for(int j = 0; j < pL; j++) {
                        if(kk < len  - 1 && stackL[kk] == j) {
                            kk++;
                            continue;
                        }

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

                    for(int x = 1; x <= lSize && x < minPQ; x++) {
                        sumCYi[x] += C[lSize][x];
                    }
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

                sampleSize = std::ceil(sumWtemp / sumW[len - 1] * T);
                // if(ansAll[len].size() < pR) ansAll[len].resize((pR + 2)*2);

                for(int x = 0; x + len <= pR+1 && x + len < minPQ; x++) {
                    if(sumCXi[x] < 0.5) break;
// if(len >= ansAll.size()) {
//     printf("len11 %d, %d\n", len, (int)ansAll.size());
//     fflush(stdout);
// }
// if(len + x >= ansAll[len].size()) {
//     printf("len1 %d, %d, %d\n", len, len + x, (int)ansAll[len].size());
//     fflush(stdout);
// }
// assert(len < ansAll.size());
// assert(len + x < ansAll[len].size());
                    ansAll[len][len + x] += sumCXi[x] * sumWtemp / sampleSize / C[len + x - 1][len - 1];
                }

                for(int x = 1; x + len <= pL +1 && x + len < minPQ; x++) {
                    if(sumCYi[x] < 0.5) break;
                    // if(ansAll[x + len].size() <= len) ansAll[x + len].resize((len + 2)*2);
// if(len >= ansAll[x + len].size()) {
//     printf("len %d, %u\n", len, (int)ansAll[x + len].size());
//     fflush(stdout);
// }
// assert(len < ansAll[x + len].size());
                    ansAll[x + len][len] += sumCYi[x] * sumWtemp / sampleSize / C[len + x - 1][len - 1];
                }
            }
        }
    }

    for(int x = 2; x < (int)ansAll.size() && x < minPQ; x++) {
        for(int y = 2; y < (int)ansAll[x].size() && y < minPQ; y++) {
            if(ansAll[x][y] < 0.5) break;
            printf("%d-%d: %.0f\n", x, y, ansAll[x][y]);
        }
    }
    fflush(stdout);

    delete [] bufferForDpU;
    delete [] bufferForDpV;
    delete [] dpU;
    delete [] dpV;
}


void colorPath::approximateCountingAllVersion2(uint64_t T) {
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
    for(uint32_t i = 0; i <= minPQ; i++) {
        ansAll[i].resize(minPQ + 1);
    }

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

                    for(int x = 0; x <= rSize && x < minPQ; x++) {
                        sumCXi[x] += C[rSize][x];
                    }
                    maxZR[len] = std::max(maxZR[len], rSize);

                    int lSize = 0;
                    kk = 0;
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

                sampleSize = std::ceil(T * sumWtemp / sumW[len - 1]);
                // if(ansAll[len].size() < pR) ansAll[len].resize((pR + 2)*2);

                for(int x = 0; x + len <= pR + 1 && x + len < minPQ; x++) {
                    if(sumCXi[x] < 0.5) break;
                    ansAll[len][len + x] += sumCXi[x] * sumWtemp / sampleSize / C[len + x - 1][len - 1];
                }

                for(int x = 1; x + len <= pL + 1 && x + len < minPQ; x++) {
                    if(sumCYi[x] < 0.5) break;
                    ansAll[x + len][len] += sumCYi[x] * sumWtemp / sampleSize / C[len + x - 1][len - 1];
                }
// printf("    Ans[3][2]:%.0f %.0f\n", ansAll[3][2], sumCYi[1] * sumWtemp / sampleSize / C[3][2]);
            }
        }
    }

    for(int x = 2; x < (int)ansAll.size() && x < minPQ; x++) {
        for(int y = 2; y < (int)ansAll[x].size() && y < minPQ; y++) {
            if(ansAll[x][y] < 0.5) break;
            printf("%d-%d: %.0f\n", x, y, ansAll[x][y]);
        }
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





void colorPath::approximateCountingAllVersion3(uint64_t T) {
// g->print();
    uint32_t maxDuv = std::max(g->maxDu, g->maxDv);

    const int minPQ = 100;
    uint32_t maxE = std::min((1u<<21), g->maxDu * g->maxDv);
    // uint32_t maxE = 105000;

    candL.resize(g->n1);
    candR.resize(g->n2);
    std::vector<uint32_t> candLtemp(g->maxDv);

    double * bufferForDpU = new double[minPQ * maxE];
    double * bufferForDpV = new double[minPQ * maxE];

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
    for(uint32_t i = 0; i <= minPQ; i++) {
        ansAll[i].resize(minPQ + 1);
    }

    //subgraph
    std::vector<int> pU(g->n1), e1(maxE), pV(g->maxDu+5), e2(maxE);
    std::vector<int> mapUtoV(maxE), mapVtoU(maxE);
    //(minPQ-minPQ)-bipath
    std::vector<double> sumW(minPQ);
    std::vector<double> ddp(maxE);

    //first compute dp
printf("start dp first v3\n");fflush(stdout);

    auto computeDP = [&](int pL, int pR) {
        for(int u = 0; u < pL; u++) {
            for(int i = pU[u]; i < pU[u + 1]; i++) {
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
                    if(dpV[k - 1][mapVtoU[i]] == 0) continue;

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

            std::fill(ddp.begin(), ddp.begin() + pU[pL], 0);
            // memset(dpV[k], 0, sizeof(double) * pU[pL]);
            for(int u = 0; u < pL; u++) {
                for(int i = pU[u] + 1; i < pU[u + 1]; i++) {
                    int v = e1[i];

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
    for(uint32_t u = 0; u < g->n1; u++) {
// printf("    nodes %u\n", u);fflush(stdout);
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }
// if(u == 38288) {
//     printf("there\n");fflush(stdout);
// }
        candLtemp.clear();
        int pL = 0;
        for(int i = 0; i < pR; i++) {
            uint32_t v = candR[i];

            auto st = g->e2.begin() + g->pV[v];
            auto ed = g->e2.begin() + g->pV[v + 1];
            uint32_t j = std::lower_bound(st, ed, u) - g->e2.begin();
            for(; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(candL.idx(w) >= pL) {
// assert(w >= u);
                    candL.changeTo(w, pL++);
// assert(candL.idx(w) < pL);
                    candLtemp.push_back(w);
                }
            }
        }
// if(u == 38288) {
//     printf("there pL %d\n", pL);fflush(stdout);
// }
        std::sort(candLtemp.begin(), candLtemp.end());
        pL = 0;
        for(auto w : candLtemp) candL.changeTo(w, pL++);
// assert(candL[0] == u);
// if(u == 38288) {
//     printf("there pL %d pR %d\n", pL, pR);
    
//     for(int b = 0; b < pR; b++) {
//         int edges = 0;
//         for(int a = 0; a < pL; a++) {
//             if(g->connectUV(candL[a], candR[b])) edges++;
//         }
//         printf("edges %d\n", edges);
//     }
    

//     fflush(stdout);
// }


        pV[0] = 0;
        std::fill(pU.begin(), pU.begin() + pL + 1, 0);
        for(int i = 0; i < pR; i++) {
            uint32_t v = candR[i];
            pV[i + 1] = pV[i];

            auto st = g->e2.begin() + g->pV[v];
            auto ed = g->e2.begin() + g->pV[v + 1];
            uint32_t j = std::lower_bound(st, ed, u) - g->e2.begin();
// assert(j >= g->pV[v]);
// if(u == 38288) {
// //     edges 1
// // edges 587
// // edges 824
// // edges 3085
// assert(0 == pV[0]);
// if(i == 0) assert(g->pV[v + 1] - j == 1);
// else if(i == 1) assert(g->pV[v + 1] - j == 587);
// else if(i == 2) assert(g->pV[v + 1] - j == 824);
// else if(i == 3) assert(g->pV[v + 1] - j == 3085);
// printf("right here\n");
// fflush(stdout);
// }
            for(; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                int h = candL.idx(w);
                e2[pV[i + 1]++] = h;
                pU[h + 1]++;
            }
// if(u == 38288) {
// //     edges 1
// // edges 587
// // edges 824
// // edges 3085
// assert(pV[0] == 0);
// if(i == 0) assert(pV[i + 1] - pV[i] == 1);
// else if(i == 1) assert(pV[i + 1] - pV[i]== 587);
// else if(i == 2) assert(pV[i + 1] - pV[i] == 824);
// else if(i == 3) assert(pV[i + 1] - pV[i] == 3085);
// printf("right here\n");
// fflush(stdout);
// }
        }


        for(int w = 0; w < pL; w++) {
            pU[w + 1] += pU[w];
// if(u == 38288) {
//     printf("u%d : %d\n", w, pU[w]);
// }
        }
// if(pU[pL] != pV[pR]) {
//     printf("pU[pL] %d, pV[pR] %d\n", pU[pL], pV[pR]);
// }
// // assert(pU[pL] == pV[pR]);
// if(u == 38288) {
//     printf("there pV[pR] %d\n", pV[pR]);
//     for(int v = 0; v < pR; v++) {
//         printf("%d: %d %d\n", v, pV[v], pV[v + 1]);
//     }
//     printf("pU %d\n", pU[pL]);
//     fflush(stdout);
// }

        for(int v = 0; v < pR; v++) {
// if(u == 38288) {
// printf("v:%d\n", v);fflush(stdout);
// }
            for(int i = pV[v]; i < pV[v + 1]; i++) {
// printf("i:%d\n", i);fflush(stdout);
// if(u == 38288 && 7008 == i) {
//     printf("%d, v%d, pv %d %d, e2sz %d, u %d, pL %d\n", 
//         i, v, pV[v], pV[v + 1], (int)e2.size(), e2[i], pL);fflush(stdout);
// }
                int w = e2[i];
// if(u == 38288 && 7008 == i) {
//     printf("%d %d w %d, pU[w]:%d\n", i, u, w, pU[w]);fflush(stdout);
// }
                e1[pU[w]] = v;
// if(u == 38288 && 7008 == i) {
//     printf("%d \n", i);fflush(stdout);
// }   
                mapUtoV[pU[w]] = i;
                mapVtoU[i] = pU[w];
// if(u == 38288 && 7008 == i) {
//     printf("%d \n", i);fflush(stdout);
// }
                pU[w]++;
            }
        }

        for(int u = pL; u >= 1; u--) pU[u] = pU[u - 1];
        pU[0] = 0;
// if(u == 38288) {
//     printf("there\n");fflush(stdout);
// }
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
            // if(pL < minPQ - 1 || pR < minPQ - 1) continue;
// printf("Graph: %u %u:\n", u, v);
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


        maxPLen = std::max(maxPLen, maxPLength);

        for(int j = pU[0]; j < pU[1]; j++) {
            for(int k = 2; k < maxPLength && k < minPQ; k++) {
                sumW[k] += dpU[k][j];
            }
        }
    }

    printf("maxPlen %d\n\nsumW ",  maxPLen);
    for(int i = 2; i < maxPLen; i++) {
        printf("%d:%.0f ", i, sumW[i]);
    }
    printf("\n\n");fflush(stdout);
    
    std::default_random_engine generator(1000);
    std::uniform_real_distribution<double> uiDistribution(0, 1);
    std::vector<int> stackL(g->maxDv+1), stackR(g->maxDu+1);
    std::vector<double> sumCXi(g->maxDu), sumCYi(g->maxDv);

    // if(sumW > 0)
    for(uint32_t u = 0; u < g->n1; u++) {
// if(u > 106000) {
//     printf("u:%u\n", u);fflush(stdout);
// }
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }

        candLtemp.clear();
        int pL = 0;
        for(int i = 0; i < pR; i++) {
            uint32_t v = candR[i];

            auto st = g->e2.begin() + g->pV[v];
            auto ed = g->e2.begin() + g->pV[v + 1];
            uint32_t j = std::lower_bound(st, ed, u) - g->e2.begin();
            for(; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(candL.idx(w) >= pL) {
                    candL.changeTo(w, pL++);
                    candLtemp.push_back(w);
                }
            }
        }
        std::sort(candLtemp.begin(), candLtemp.end());
        pL = 0;
        for(auto w : candLtemp) candL.changeTo(w, pL++);

        pV[0] = 0;
        std::fill(pU.begin(), pU.begin() + pL + 1, 0);
        for(int i = 0; i < pR; i++) {
            uint32_t v = candR[i];
            pV[i + 1] = pV[i];

            auto st = g->e2.begin() + g->pV[v];
            auto ed = g->e2.begin() + g->pV[v + 1];
            uint32_t j = std::lower_bound(st, ed, u) - g->e2.begin();
            for(; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                int u = candL.idx(w);
                e2[pV[i + 1]++] = u;
                pU[u + 1]++;
            }
        }
        
        for(int u = 1; u <= pL; u++) pU[u] += pU[u - 1];
        for(int v = 0; v < pR; v++) {
            for(int i = pV[v]; i < pV[v + 1]; i++) {
                int u = e2[i];
                e1[pU[u]] = v;
                
                mapUtoV[pU[u]] = i;
                mapVtoU[i] = pU[u];

                pU[u]++;
            }
        }

        for(int u = pL; u >= 1; u--) pU[u] = pU[u - 1];
        pU[0] = 0;

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

        for(int len = 2; len < maxPLength && len < minPQ && sumW[len] > 0.5; len++) {//(len-len)-bipath
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

    for(int x = 2; x < (int)ansAll.size() && x < minPQ; x++) {
        for(int y = 2; y < (int)ansAll[x].size() && y < minPQ; y++) {
            if(ansAll[x][y] < 0.5) break;
            printf("%d-%d: %.8f\n", x, y, ansAll[x][y]);
        }
    }
    fflush(stdout);

    delete [] bufferForDpU;
    delete [] bufferForDpV;
    delete [] dpU;
    delete [] dpV;
}



void colorPath::approximateCountingAllVersion4(uint64_t T) {
// g->print();
    uint32_t maxDuv = std::max(g->maxDu, g->maxDv);

    const int minPQ = 100;

    uint32_t maxE = g->m;
    candL.resize(g->n1);
    candR.resize(g->n2);

    double * bufferForDpU = new double[minPQ * maxE];
    double * bufferForDpV = new double[minPQ * maxE];

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
    for(uint32_t i = 0; i <= minPQ; i++) {
        ansAll[i].resize(minPQ + 1);
    }

    //subgraph
    std::vector<int> mapUtoV(maxE), mapVtoU(maxE);
    for(uint32_t u = 0; u < g->n1; u++) {
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];

            auto st = g->e2.begin() + g->pV[v];
            auto ed = g->e2.begin() + g->pV[v + 1];
            uint32_t j = std::lower_bound(st, ed, u) - g->e2.begin();

            mapUtoV[i] = j;
            mapVtoU[j] = i;
        }
    }
    //(minPQ-minPQ)-bipath
    std::vector<double> sumW(minPQ);
    std::vector<double> ddp(maxE);

    //first compute dp
printf("start dp first\n");fflush(stdout);
    uint32_t pL = g->n1;
    uint32_t pR = g->n2;

    for(uint32_t u = 0; u < pL; u++) {
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            dpV[1][i] = g->pU[u + 1] - i - 1;
        }
    }

    uint32_t minLR = std::min(pL, pR);
    int k = 2;
    for(; k <= minLR && k < minPQ; k++) {
        bool f = false;
        std::fill(ddp.begin(), ddp.begin() + g->pU[pL], 0);
        for(uint32_t v = 0; v < pR; v++) {
            for(uint32_t i = g->pV[v] + 1; i < g->pV[v + 1]; i++) {
                uint32_t u = g->e2[i];
                if(dpV[k - 1][mapVtoU[i]] < 0.5) continue;

                ddp[g->pV[v]] += dpV[k - 1][mapVtoU[i]]; 
                ddp[i] -= dpV[k - 1][mapVtoU[i]];
// if(k == 2)
// printf("%u-%u: %u %u, %u-%u, %.1f\n", 
//     u, v, g->pV[v], i, g->e2[g->pV[v]], v, dpV[k - 1][mapVtoU[i]]);
            }
        }
        for(int v = 0; v < pR; v++) {
            dpU[k][g->pV[v]] = ddp[g->pV[v]] < 0 ? 0 : ddp[g->pV[v]];
            if(dpU[k][g->pV[v]] > 0.5) f = true;
            for(uint32_t e = g->pV[v] + 1; e < g->pV[v + 1]; e++) {
                dpU[k][e] = dpU[k][e-1] + ddp[e];
                if(dpU[k][e] > 0.5) f = true;
            }
        }

// printf("%f %f\n", dpU[2][mapUtoV[0]], ddp[g->pV[1]]);

        std::fill(ddp.begin(), ddp.begin() + g->pU[pL], 0);
        // memset(dpV[k], 0, sizeof(double) * pU[pL]);
        for(uint32_t u = 0; u < pL; u++) {
            for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
                uint32_t v = g->e1[i];
                if(dpU[k][mapUtoV[i]] < 0.5) continue;

                ddp[g->pU[u]] += dpU[k][mapUtoV[i]];
                ddp[i] -= dpU[k][mapUtoV[i]];
            }
        }
        for(uint32_t u = 0; u < pL; u++) {
            dpV[k][g->pU[u]] = ddp[g->pU[u]] < 0? 0 : ddp[g->pU[u]];
            if(dpV[k][g->pU[u]] > 0.5) f = true;
            for(uint32_t e = g->pU[u] + 1; e < g->pU[u + 1]; e++) {
                dpV[k][e] = dpV[k][e-1] + ddp[e];
                if(dpV[k][e] > 0.5) f = true;
            }
        }

        if(!f) break;
    }
    int maxPLen = k;
    

    for(int j = 0; j < maxE; j++) {
        for(int k = 2; k < maxPLen && k < minPQ; k++) {
            if(dpU[k][j] < 0.5) break;
            sumW[k] += dpU[k][j];
        }
    }

// g->print();

// printf("dpU:\n");
// for(uint32_t u = 0; u < g->n1; u++) {
//     for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//         printf("%u-%u:", u, g->e1[i]);
//         for(int j = 1; j < maxPLen; j++) {
//             if(dpU[j][mapUtoV[i]] < 0.5) break;
//             printf("%.0f ", dpU[j][mapUtoV[i]]);
//         }
//         printf("\n");
//     }
// }

// printf("dpV:\n");
// for(uint32_t u = 0; u < g->n1; u++) {
//     for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//         printf("%u-%u:", u, g->e1[i]);
//         for(int j = 1; j < maxPLen; j++) {
//             if(dpV[j][i] < 0.5) break;
//             printf("%.0f ", dpV[j][i]);
//         }
//         printf("\n");
//     }
// }

    printf("maxPlen %d\n\nsumW ",  maxPLen);
    for(int i = 2; i < maxPLen; i++) {
        printf("%d:%.0f ", i, sumW[i]);
    }
    printf("\n\n");fflush(stdout);
    
    std::default_random_engine generator;
    std::uniform_real_distribution<double> uiDistribution(0, 1);
    std::vector<int> stackL(minPQ + 1), stackR(minPQ + 1);
    std::vector<double> sumCXi(minPQ + 1), sumCYi(minPQ + 1);

    for(uint32_t u = 0; u < g->n1; u++) {
// printf("    u: %u\n", u);fflush(stdout);
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
// printf("e: %u-%u\n", u, v);fflush(stdout);
            int pL = 0;
            auto st = g->e2.begin() + g->pV[v];
            auto ed = g->e2.begin() + g->pV[v + 1];
            uint32_t j = std::upper_bound(st, ed, u) - g->e2.begin();
            for(; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) candL.changeTo(w, pL++);
            }
            candR.changeTo(v, --pR);
            for(int i = 1; i < pR; i++) {
                candR.changeToByPos(i, i - 1);
            }
// printf("candL:");
// for(int i = 0; i < pL; i++) printf("%u ", candL[i]);
// printf("  ");
// printf("candR:");
// for(int i = 0; i < pR; i++) printf("%u ", candR[i]);
// printf("\n");
            for(int len = 2; len < maxPLen && len < minPQ && sumW[len] > 0.5; len++) {//(len-len)-bipath
                if(dpU[len][mapUtoV[i]] < 0.5) break;
                double sumWtemp = dpU[len][mapUtoV[i]];

                uint64_t sampleSize = std::ceil(sumWtemp / sumW[len] * T);
// printf("len %d: sumW:%.0f spSize:%llu\n", len, sumWtemp, sampleSize);
                if(sampleSize == 0) continue;

                std::fill(sumCXi.begin(), sumCXi.begin() + pR + 1, 0.0);
                std::fill(sumCYi.begin(), sumCYi.begin() + pL + 1, 0.0);

// int tt = 0;
                while(sampleSize--) {
                    double tmp = 0.0;
                    double r = uiDistribution(generator); 
                    int preE = i, preU = u, preV = v;
                    stackL.clear();
                    stackR.clear();
                    stackL.push_back(u);
                    stackR.push_back(v);

                    bool connect = true;
                    for(int i = 1; i < len; i++) {
                        double r = uiDistribution(generator);
                        tmp = 0.0;
                        
                        auto st = g->e2.begin() + g->pV[preV];
                        auto ed = g->e2.begin() + g->pV[preV + 1];
                        int j = std::upper_bound(st, ed, preU) - g->e2.begin();
                        for(; j < g->pV[preV + 1]; j++) {
                            // if(e2[j] <= preU) continue;
                            tmp += dpV[len - i][mapVtoU[j]];
                            if(tmp + 1e-8 >= r * dpU[len - i + 1][mapUtoV[preE]]) {
                                uint32_t w = g->e2[j];
                                stackL.push_back(w);
                                preU = w;
                                preE = j;
                                break;
                            }
                        }
                        for(int j = 0; j < i - 1; j++) {
                            if(!g->connectUV(preU, stackR[j])) {
                                connect = false;
// printf("err1stackL:");
// for(int x = 0; x < stackL.size(); x++) printf("%u ", stackL[x]);
// printf("  ");
// printf("stackR:");
// for(int x = 0; x < stackR.size(); x++) printf("%u ", stackR[x]);
// printf("\n");
                                break;
                            }
                        }
                        if(!connect) break;

                        r = uiDistribution(generator);
                        tmp = 0.0;
                        st = g->e1.begin() + g->pU[preU];
                        ed = g->e1.begin() + g->pU[preU + 1];
                        j = std::upper_bound(st, ed, preV) - g->e1.begin();
                        for(; j < g->pU[preU + 1]; j++) {
                            // if(e1[j] <= preV) continue;
                            tmp += dpU[len - i][mapUtoV[j]];
                            if(tmp + 1e-8 >= r * dpV[len - i][mapVtoU[preE]]) {
                                uint32_t v = g->e1[j];
                                stackR.push_back(v);
                                preV = v;
                                preE = j;
                                break;
                            }
                        }
                        for(int j = 0; j < i; j++) {
                            if(!g->connectUV(stackL[j], preV)) {
                                connect = false;
// printf("err2stackL:");
// for(int x = 0; x < stackL.size(); x++) printf("%u ", stackL[x]);
// printf("  ");
// printf("stackR:");
// for(int x = 0; x < stackR.size(); x++) printf("%u ", stackR[x]);
// printf("\n");
                                break;
                            }
                        }
                        if(!connect) break;
                    }

                    if(!connect) continue;

                    if(stackL.size() < len) {
                        continue;
                    }
                    assert(stackL.size() == len);
// printf("stackL:");
// for(int i = 0; i < len; i++) printf("%u ", stackL[i]);
// printf("  ");
// printf("stackR:");
// for(int i = 0; i < len; i++) printf("%u ", stackR[i]);
// printf("\n");

                    int rSize = 0, kk = 1;
                    for(int j = 0; j < pR; j++) {
                        if(kk < len && stackR[kk] == candR[j]) {
                            kk++;
                            continue;
                        }
                        bool f = true;
                        for(int k = 1; k < len; k++) {
                            if(!g->connectUV(candR[j], stackL[k])) {
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

                    int lSize = 0;
                    kk = 1;
                    for(int j = 0; j < pL; j++) {
                        if(kk < len && stackL[kk] == candL[j]) {
                            kk++;
                            continue;
                        }

                        bool f = true;
                        for(int k = 1; k < len; k++) {
                            if(!g->connectUV(candL[j], stackR[k])) {
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
// printf("lr %d %d\n", lSize, rSize);
// if(rSize > 0) tt++;
                }
// printf("%d\n", tt);

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

                for(int x = 0; x + len <= pR + 1 && x + len < maxPLen; x++) {
                    if(sumCXi[x] < 0.5) break;
                    ansAll[len][len + x] += sumCXi[x] * sumWtemp / sampleSize / C[len + x - 1][len - 1];
                }

                for(int x = 1; x + len <= pL + 1 && x + len < maxPLen; x++) {
                    if(sumCYi[x] < 0.5) break;
                    ansAll[x + len][len] += sumCYi[x] * sumWtemp / sampleSize / C[len + x - 1][len - 1];
                }
// printf("    ansAll[2][3] %.0f\n", ansAll[2][3]);
            }
        }
    }

    for(int x = 2; x < (int)ansAll.size() && x < minPQ; x++) {
        for(int y = 2; y < (int)ansAll[x].size() && y < minPQ; y++) {
            if(ansAll[x][y] < 0.5) break;
            printf("%d-%d: %.0f\n", x, y, ansAll[x][y]);
        }
    }
    fflush(stdout);

    delete [] bufferForDpU;
    delete [] bufferForDpV;
    delete [] dpU;
    delete [] dpV;
}


// void colorPath::approximateCounting(uint64_t T) {
//     candL.resize(g->n1);
//     candR.resize(g->n2);

//     std::vector<uint32_t> edgesFirstLevel(g->n1);
//     std::vector<uint32_t> sampleSize(g->n1);
//     std::vector<double> delta(g->n1);
//     double sumDelta;

//     std::vector<uint32_t> col(g->m);
//     std::vector<uint32_t> colMap(g->m);
//     /**
//     *** Pr(clique) = path/pathv * C(ev,pq)/C(E,pq) 
//     *** ans = \sum_v sum(Xv)/[C(ev,pq)/C(E,pq)]
//     **/

//     //first level, edge count, sample size
//     for(uint32_t u = 0; u < g->n1; u++) {
//         int pR = 0;
//         for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//             uint32_t v = g->e1[i];
//             candR.changeTo(v, pR++);
//         }
//         uint32_t m = 0;

//         for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//             uint32_t v = g->e1[i];

//             int pL = 0;
//             for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
//                 uint32_t w = g->e2[j];
//                 if(w > u) candL.changeTo(w, pL++);
//             }
//             candR.changeTo(v, --pR);//ordered still

//             //build sg, for edge count
//             // buildSg();//linear
//             for(int l = 0; l < pL; l++) {
//                 uint32_t x = candL[l];
//                 for(uint32_t j = g->pU[x]; j < g->pU[x + 1]; j++) {
//                     uint32_t y = g->e1[j];
//                     if(candR.idx(y) < pR) {
//                         m++;
//                     }
//                 }
//             }
//             edgesFirstLevel[u] += m;
//             delta[u] += m;
//         }

//         delta[u] /= g->n1;
//         sumDelta += delta[u];
//     }

//     for(uint32_t u = 0; u < g->n1; u++) {
//         sampleSize[u] = T * delta[u] / sumDelta;
//     }

//     //color
//     assert(g->n1 > 0);
//     uint32_t maxColor = 0;
//     for(uint32_t u = g->n1 - 1; u >= 0; u--) {
//         int pR = 0;
//         for(uint32_t i = g->pU[u + 1] - 1; i >= g->pU[u]; i--) {
//             uint32_t v = g->e1[i];
//             candR.changeTo(v, pR++);
//         }

//         for(uint32_t i = g->pU[u + 1] - 1; i >= g->pU[u]; i--) {
//             uint32_t v = g->e1[i];
            
//             int pL = 0;
//             for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
//                 uint32_t w = g->e2[j];
//                 if(w > u) candL.changeTo(w, pL++);
//             }
//             candR.changeTo(v, --pR);//ordered still

//             for(int l = 0; l < pL; l++) {
//                 uint32_t x = candL[l];
//                 for(uint32_t j = g->pU[x]; j < g->pU[x + 1]; j++) {
//                     uint32_t y = g->e1[j];
//                     if(candR.idx(y) < pR) {
//                         colMap[col[j]] = i;
//                     }
//                 }
//             }

//             uint32_t cc = 0;
//             while(cc < g->m && colMap[cc] == i) cc++;
//             col[i] = cc;
//             maxColor = std::max(cc, maxColor);
//         }
//     }

//     //dp 1
//     std::vector<uint32_t> pU(g->m), pV(g->m), e(g->m), tmpe(g->m);

//     for(uint32_t u = 0; u < g->n1; u++) {
//         int pR = 0;
//         for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//             uint32_t v = g->e1[i];
//             candR.changeTo(v, pR++);
//         }

//         for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//             uint32_t v = g->e1[i];

//             int pL = 0;
//             for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
//                 uint32_t w = g->e2[j];
//                 if(w > u) candL.changeTo(w, pL++);
//             }
//             candR.changeTo(v, --pR);//ordered still

//             //build sg
//             std::fill(pU.begin(), pU.begin() + maxColor + 1, 0);
//             tmpe.clear();
//             uint32_t pe = 0;
//             for(int j = 0; j < pL; j++) {
//                 uint32_t x = candL[j];
//                 for(uint32_t l = g->pU[x]; l < g->pU[x + 1]; l++) {
//                     uint32_t y = g->e1[l];
//                     if(candR.idx(y) < pR) {
//                         pU[col[l] + 1]++;
//                         tmpe.push_back(l);
//                     }
//                 }
//             }
// //count of colors of the edges
//             for(uint32_t c = 0; c < maxColor; c++) {
//                 pU[c + 1] += pU[c];
//             }
//             if(e.capacity() < tmpe.size()) e.resize(tmpe.size());
//             for(int j = 0; j < tmpe.size(); j++) {
//                 e[pU[col[tmpe[j]]]++] = j;
//             }

            

//             // for()

//             for(int k = 1; k < p + q; k++) {
//                 for(int j = 0; j < tmpe.size(); j++) {
//                     for(int l = pU[j]; l < pU[j + 1]; l++) {

//                     }
//                 }
//             }

//             // dynamicProgramming();//mm'
//         }
//     }

// }


// void colorPath::approximateCounting2(uint64_t T) {
//     auto t1 = std::chrono::steady_clock::now();
//     g->coreReduction(p, q);
//     auto t2 = std::chrono::steady_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
//     std::cout << "core reduction time:" << duration.count() << "ms" << std::endl;

//     if(p > q) {
//         g->swapUV();
//         std::swap(p, q);
//     }
//     //p <= q

//     double * bufferForDpU = new double[g->maxDv * (p + q)];
//     double * bufferForDpV = new double[g->maxDu * (p + q)];

//     double ** dpU = new double*[g->maxDv];
//     for(uint32_t u = 0; u < g->maxDv; u++) {
//         dpU[u] = bufferForDpU + u * (p + q);
//     }
//     for(uint32_t u = 0; u < g->maxDv; u++) {
//         dpU[u][0] = 1;
//     }

//     double ** dpV = new double*[g->maxDu];
//     for(uint32_t v = 0; v < g->maxDu; v++) {
//         dpV[v] = bufferForDpV + v * (p + q);
//     }
//     for(uint32_t v = 0; v < g->maxDu; v++) {
//         dpV[v][0] = 1;
//     }

//     candL.resize(g->n1);
//     candR.resize(g->n2);
//     // uint32_t subGraphEcount = std::min(g->m, g->maxDu * g->maxDv);
//     std::vector<uint32_t> pU(g->maxDv), e1(g->m);
//     std::vector<uint32_t> e2[g->maxDu];
//     std::vector<std::vector<double>> w;
//     std::vector<double> sumW;

//     w.resize(p + q + 1);
//     sumW.resize(p + q + 1);
//     for(int k = 4; k <= p + q; k++) {
//         w[k].resize(g->n1);
//     }

//     for(uint32_t u = 0; u < g->n1; u++) {
//         int pR = 0;
//         for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//             uint32_t v = g->e1[i];
//             candR.changeTo(v, pR++);
//         }
//         uint32_t m = 0;

//         for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//             uint32_t v = g->e1[i];

//             int pL = 0;
//             for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
//                 uint32_t w = g->e2[j];
//                 if(w > u) candL.changeTo(w, pL++);
//             }
//             candR.changeTo(v, --pR);

//             if(pL < p || pR < q) continue;

//             std::fill(pU.begin(), pU.begin() + pL + 1, 0);
//             for(int k = 0; k < pR; k++) e2[k].clear();
//             // std::fill(pV.begin(), pV.begin() + pR + 1, 0);
//             for(int j = 0; j < pL; j++) {
//                 uint32_t x = candL[j];
//                 for(int k = 0; k < pR; k++) {
//                     uint32_t y = candR[k];
//                     if(g->connectUV(x, y)) {
//                         e1[pU[j]++] = k;
//                         e2[k].push_back(j);
//                     }
//                 }
//             }

//             for(int k = 1; k < p + q; k++) {
//                 for(int u = 0; u < pL; u++) {
//                     for(int j = pU[u]; j < pU[u + 1]; j++) {
//                         int v = e1[j];
//                         dpU[u][k] += dpV[v][k - 1];
//                     }
//                 }

//                 for(int v = 0; v < pR; v++) {
//                     for(int j = 0; j < e2[v].size(); j++) {
//                         int u = e2[v][j];
//                         dpV[v][k] += dpU[u][k - 1];
//                     }
//                 }
//             }
//         }
    
//         for(int pp = 2; pp <= p; pp++) {
//             for(int qq = 2; qq <= q; qq++) {
//                 // w[pp + qq][u] = pp >= qq ? dpU[u][pp+qq] : dpV[v][pp+qq];
//                 // sumW[pp + qq] 
//             }
//         }
//     }

//     std::vector<std::vector<uint32_t>> sampleSize;

//     sampleSize.resize(p + 1 + q);

//     // for(int k = 4; k <= p + q; k++) {
//     //     for(uint32_t u = 0; u < g->n1; u++) {
//     //         sampleSize[u] = T * sumW[k][u] / ;
//     //     }
//     // }

//     delete [] bufferForDpU;
//     delete [] bufferForDpV;
//     delete [] dpU;
//     delete [] dpV;
// }