#include "rawEdgePivot.h"

#include <chrono>
#include <random>



void rawEdgePivot::approximateCountMaximalPivot(double rate) {
    candL.resize(g->n1);
    candR.resize(g->n2);
    tmpNodesL.resize(std::max(g->maxDu, g->maxDv) + 5);
    tmpNodesR.resize(std::max(g->maxDu, g->maxDv) + 5);

    ansAll.resize(std::max(g->maxDu, g->maxDv) + 5);
    printf("maxD:%u\n", std::max(g->maxDu, g->maxDv));
    fflush(stdout);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> uiDistribution(0, 1);

    for(uint32_t u = 0; u < g->n1; u++) {
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }
        if(pR == 0) continue;

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, --pR);

            if(uiDistribution(generator) > rate) continue;

            int pL = 0;

            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) candL.changeTo(w, pL++);
            }

#ifdef DEBUG
printf("st choose p %u-%u\n", u, v);
#endif

            pivotCount(1, pL, pR, treePath{0, 1, 0, 1});
        }

    }
    for(int i = 2; i < (int)ansAll.size(); i++) {
        for(int j = 2; j < (int)ansAll[i].size(); j++) {
            if(ansAll[i][j] > 0.0) {
                printf("%d-%d:%.0f\n", i, j, ansAll[i][j] / rate);
            }
        }
    }
    fflush(stdout);
}

void rawEdgePivot::exactCountMaximalPivotV2() {
    candL.resize(g->n1);
    candR.resize(g->n2);
    tmpNodesL.resize(g->n1);
    tmpNodesR.resize(std::max(g->maxDu, g->maxDv) + 5);

    ansAll.resize(std::max(g->maxDu, g->maxDv) + 5);
    printf("maxD:%u\n", std::max(g->maxDu, g->maxDv));
    fflush(stdout);

    std::vector<uint32_t> candLtemp(g->maxDv);
    std::vector<int> cnt(g->n1);

    for(uint32_t u = 0; u < g->n1; u++) {
// if(u % 1000 == 0) {
// printf("    nodes %u\n", u);fflush(stdout);
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

        int pL = 0;
        for(auto w: candLtemp) {
            candL.changeTo(w, pL++);
            cnt[w] = 0;
        }
        
        pivotCount(1, pL, pR, treePath{0, 1, 0, 0});
    }

    for(int i = 2; i < (int)ansAll.size(); i++) {
        for(int j = 2; j < (int)ansAll[i].size(); j++) {
            if(ansAll[i][j] > 0.0) {
                printf("%d-%d:%.0f\n", i, j, ansAll[i][j]);
            }
        }
    }
    fflush(stdout);
}

void rawEdgePivot::exactCountMaximalPivot() {
#ifdef DEBUG
g->print();
#endif

// g->coreReductionFast22();

    candL.resize(g->n1);
    candR.resize(g->n2);
    tmpNodesL.resize(std::max(g->maxDu, g->maxDv) + 5);
    tmpNodesR.resize(std::max(g->maxDu, g->maxDv) + 5);

    ansAll.resize(std::max(g->maxDu, g->maxDv) + 5);

#ifdef PQWEDGE
ansWedge.resize(std::max(g->maxDu, g->maxDv) + 5);
// Pl.resize(g->maxDv);
// Hl.resize(g->maxDv);
// Pr.resize(g->maxDu);
// Hr.resize(g->maxDu);
// Pl.clear();
// Pr.clear();
// Hl.clear();
// Hr.clear();
// sumDegPl = sumDegPr = sumDegHl = sumDegHr = 0;
#endif

    printf("maxD:%u\n", std::max(g->maxDu, g->maxDv));
    fflush(stdout);

// std::vector<int> tmpDegsL(g->maxDv), tmpDegsR(g->maxDu);
// auto t1 = std::chrono::steady_clock::now();
    for(uint32_t u = 0; u < g->n1; u++) {
// double t1 = clock();
// auto t1 = std::chrono::steady_clock::now();
// if(u % 1000 == 0) {
//     auto t2 = std::chrono::steady_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
//     printf("%u, %u, %lld\n", u, g->deg1(u), (long long)duration.count());fflush(stdout);
//     t1 = t2;
// }
// printf("%u ", u);fflush(stdout);

        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }
        if(pR == 0) continue;

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            int pL = 0;
// printf("%u ", i);fflush(stdout);
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) candL.changeTo(w, pL++);
            }

            candR.changeTo(v, --pR);
            // for(int i = 1; i < pR; i++) {
            //     candR.changeToByPos(i, i - 1);
            // }
#ifdef DEBUG
printf("st choose p %u-%u\n", u, v);
#endif
// auto t1 = std::chrono::steady_clock::now();
// printf("%u ", i);fflush(stdout);
#ifdef PQWEDGE
sumDegHl = g->deg1(u);
sumDegHr = g->deg2(v);
sumDegPl = 0.0;
sumDegPr = 0.0;
#endif
            pivotCount(1, pL, pR, treePath{0, 1, 0, 1});

// if(u > 1950000) {

// auto t2 = std::chrono::steady_clock::now();
// auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
// long long ts = (long long)duration.count();

// if(ts < 10) continue; 

// int edges = 0, edgesOfedges = 0;

// for(int j = 0; j < pR; j++) tmpDegsR[j] = 0;
// for(int i = 0; i < pL; i++) {
//     tmpDegsL[i] = 0;
//     for(int j = 0; j < pR; j++) {
//         if(g->connectUV(candL[i], candR[j])) {
//             edges++;
//             tmpDegsL[i]++;
//             tmpDegsR[j]++;
//         }
//     }
// }
// for(int i = 0; i < pL; i++) {
//     for(int j = 0; j < pR; j++) {
//         if(g->connectUV(candL[i], candR[j])) {
//             edgesOfedges += tmpDegsR[j];
//         }
//     }
// }

// printf("%u %lld e/l:%.2f e/r:%.2f e/(l+r):%.2f e/()l*r:%.2f ee/e:%.2f\n", 
//     u, ts, 1.0*edges/pL, 1.0*edges/pR, 
//     1.0*edges/(pL + pR), 1.0*edges/(pL + pR) ,
//     1.0*edgesOfedges/edges);
// fflush(stdout);
// }

        }

// double t2 = clock();
// auto t2 = std::chrono::steady_clock::now();
// auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
// uint32_t sumDeg = 0;
// for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//     uint32_t v = g->e1[i];
//     sumDeg += g->e2.begin()+g->pV[v + 1]-std::upper_bound(g->e2.begin()+g->pV[v], 
//         g->e2.begin()+g->pV[v + 1], u);
//     // sumDeg += g->deg2(v);
// }
// printf("%u, %u, %u, %.2f, %.2f\n", u, g->deg1(u), sumDeg, (1.0*sumDeg)/g->deg1(u),
// // (long long)duration.count());
// (t2 - t1));
// t1 = t2;

    }

    // printf("\t");
    // for(int i = 1; i <= q; i++) {
    //     printf("%d\t", i);
    // }
    // printf("\n");
    // for(int i = 1; i <= p; i++) {
    //     printf("%d\t", i);
    //     for(int j = 1; j <= q; j++) {
    //         printf("%.0f ", ansAll[i][j]);
    //     }
    //     printf("\n");
    // }
    for(int i = 2; i < (int)ansAll.size(); i++) {
        for(int j = 2; j < (int)ansAll[i].size(); j++) {
            if(ansAll[i][j] > 0.0) {
#ifdef PQWEDGE
                printf("%d-%d:%.8f\n", i, j, ansAll[i][j]*i*j / ansWedge[i][j]);
#else
                printf("%d-%d:%.0f\n", i, j, ansAll[i][j]);
#endif
            }
        }
    }

    fflush(stdout);

#ifdef COEI_BASELINE
g->print();
g->n[0] = g->n1;
g->n[1] = g->n2;

constexpr int ls = 2, rs = 2;
int clique = 0;
int wedge = 0;
auto checkSize = [&](uint32_t i, uint32_t j, int ls, int rs)->bool {
    uint32_t sz = 0;
    for(uint32_t u = 0; u < g->n[0]; u++) {
        if((i & (1<<u)) == 0) continue;
        sz++;
    }
    if(sz != ls) return false;
    sz = 0;
    for(uint32_t v = 0; v < g->n[1]; v++) {
        if((j & (1<<v)) == 0) continue;
        sz++;
    }
    if(sz != rs) return false;

    return true;
};

auto checkClique = [&](uint32_t i, uint32_t j) {
    for(uint32_t u = 0; u < g->n[0]; u++) {
        if((i & (1<<u)) == 0) continue;
        for(uint32_t v = 0; v < g->n[1]; v++) {
            if((j & (1<<v)) == 0) continue;
            if(!g->connect(u, v, 0)) return false; 
        }
    }

    for(uint32_t v = 0; v < g->n[1]; v++) {
        if((j & (1<<v)) == 0) continue;
        for(uint32_t u = 0; u < g->n[0]; u++) {
            if((i & (1<<u)) == 0) continue;
            if(!g->connect(u, v, 0)) return false; 
        }
    }

    return true;
};

for(uint32_t i = (1<<g->n[0])-1; i > 0; i--) {
    for(uint32_t j = (1<<g->n[1])-1; j > 0; j--) {
        if(!checkSize(i, j, ls, rs)) continue;
        if(checkClique(i, j)) clique++;
    }
}
printf("baseline_clique_%d_%d:%d\n", ls, rs, clique);
printf("alg_clique_%d_%d:%.0f\n", ls, rs, ansAll[ls][rs]);
auto updateWedgel = [&](uint32_t i, uint32_t j) { 
    for(uint32_t v = 0; v < g->n[1]; v++) {
        if((j & (1<<v)) == 0) continue;
        wedge += g->deg2(v) - ls+1;
    }
};
auto updateWedger = [&](uint32_t i, uint32_t j) { 
    for(uint32_t u = 0; u < g->n[0]; u++) {
        if((i & (1<<u)) == 0) continue;
        wedge += g->deg1(u) - rs+1;
    }
};

for(uint32_t i = (1<<g->n[0])-1; i > 0; i--) {
    for(uint32_t j = (1<<g->n[1])-1; j > 0; j--) {
        if(checkSize(i, j, ls - 1, rs) && checkClique(i, j)) {
            updateWedgel(i, j);
        }
        if(checkSize(i, j, ls, rs - 1) && checkClique(i, j)) {
            updateWedger(i, j);
        }
    }
}
printf("baseline_wedge_%d_%d:%d\n", ls, rs, wedge);
printf("alg_wedge_%d_%d:%.0f\n", ls, rs, ansWedge[ls][rs]);
printf("baseline_cluster_%d_%d:%f\n", ls, rs, 1.0*clique/wedge);
printf("alg_cluster_%d_%d:%f\n", ls, rs, ansAll[ls][rs]/ansWedge[ls][rs]);

#endif

}



void rawEdgePivot::pivotCount(int l, int pL, int pR, treePath t) {
// printf("l:%d, pL:%d, pR:%d, %d :%d :%d :%d\n",
//     l, pL, pR, t.p1, t.h1, t.p2, t.h2);
// fflush(stdout);
#ifdef DEBUG
printf("l:%d, pL:%d, pR:%d, %d :%d :%d :%d\n",
    l, pL, pR, t.p1, t.h1, t.p2, t.h2);
for(int i = 0; i < pL; i++) {
printf("%u ", candL[i]);
}

printf("\n");
for(int i = 0; i < pR; i++) {
    printf("%u ", candR[i]);
}
printf("\n");
#endif
#ifdef DEBUG
int x = 3, y = 1;
#endif
// if(l > 10) return;

    if(pL == 0 && pR == 0) {
        // t.p1 += pL;
        // t.p2 += pR;
#ifdef DEBUG
printf("adding %d %d %d %d\n", t.p1, t.h1, t.p2, t.h2);
#endif
        for(int ll = 0; ll <= t.p1 && ll < maxD2; ll++) {
            for(int r = 0; r <= t.p2 && r < maxD2; r++) {
                if(r + t.h2 >= (int)ansAll[ll + t.h1].size()) {
                    ansAll[ll + t.h1].resize(r + t.h2 + 4);
                }
                ansAll[ll + t.h1][r + t.h2] += C[t.p1][ll] * C[t.p2][r];
#ifdef PQWEDGE
if(r + t.h2 + 1 >= (int)ansWedge[ll + t.h1].size())
    ansWedge[ll + t.h1].resize(r + t.h2 + 6);

//for each v in h1, it must be selected. \sum_v(dv - (p1+h2))*C()*C()=\sum_vdv-|t.h1|*(p1+h2)
ansWedge[ll + t.h1][r + t.h2 + 1] += C[t.p1][ll] * C[t.p2][r] 
    * (sumDegHl - t.h1*(r + t.h2));

if(r + t.h2 >= (int)ansWedge[ll + t.h1 + 1].size())
    ansWedge[ll + t.h1 + 1].resize(r + t.h2 + 5);

ansWedge[ll + t.h1 + 1][r + t.h2] += C[t.p1][ll] * C[t.p2][r] 
    * (sumDegHr - t.h2*(ll + t.h1));
// printf("type1:%u %u %u %u, %u %u\n", t.p1, t.p2, t.h1, t.h2, ll, r);
#endif
            }
        }

#ifdef PQWEDGE
//for each v in p1, when it is selected, dv-(r+t.h2)
for(int ll = 0; ll < t.p1 && ll < maxD2; ll++) {
    for(int r = 0; r <= t.p2 && r < maxD2; r++) {
if(r + t.h2 + 1 >= (int)ansWedge[ll + t.h1 + 1].size())
    ansWedge[ll + t.h1 + 1].resize(r + t.h2 + 5);

        ansWedge[ll + t.h1 + 1][r + t.h2 + 1] += C[t.p1-1][ll] * C[t.p2][r] 
            * (sumDegPl - t.p1*(r + t.h2));
    }
}
//for each v in p2, when it is selected, dv-(l+t.h1)
for(int ll = 0; ll <= t.p1 && ll < maxD2; ll++) {
    for(int r = 0; r < t.p2 && r < maxD2; r++) {
if(r + t.h2 + 1 >= (int)ansWedge[ll + t.h1 + 1].size())
    ansWedge[ll + t.h1 + 1].resize(r + t.h2 + 5);
        ansWedge[ll + t.h1 + 1][r + t.h2 + 1] += C[t.p1][ll] * C[t.p2-1][r] 
            * (sumDegPr - t.p2*(ll + t.h1));
    }
}
// printf("type2:%u %u %u %u\n", t.p1, t.p2, t.h1);
// printf("wedge22 %u:%.f\n", l, ansWedge[2][2]);
#endif
          

        
#ifdef DEBUG
if((int)ansAll[x].size() > y)
printf("ansAll[%d][%d] %.2f\n", x, y, ansAll[x][y]);
else printf("ansAll[%d][%d] %.2f\n", x, y, 0.0);
#endif
        return ;
    }

    int pivotL = 0;
    int maxDeg = 0;
    for(int i = 0; i < pL; i++) {
        int deg = 0;

        uint32_t u = candL[i];

        if((uint32_t)pR*3 < g->deg1(u)) {
            for(int j = 0; j < pR; j++) {
                if(g->connectUV(u, candR[j])) {
                    deg++;
                }
            }
        }
        else {
            for(uint32_t j = g->pU[u]; j < g->pU[u + 1]; j++) {
                if(candR.idx(g->e1[j]) < (uint32_t)pR) {
                    deg++;
                }
            }
        }

        if(deg > maxDeg) {
            maxDeg = deg;
            pivotL = u;
        }
    }
#ifdef DEBUG
printf("pivotL:%d, %d\n", pivotL, maxDeg);
#endif

    if(maxDeg == 0) {
#ifdef PQWEDGE
double tmpSumDegPl = 0.0;
for(int j = 0; j < pL; j++) tmpSumDegPl += g->deg1(candL[j]);
sumDegPl += tmpSumDegPl;
#endif
        pivotCount(l + 1, 0, 0, {t.p1 + pL, t.h1, t.p2, t.h2});
#ifdef PQWEDGE
sumDegPl -= tmpSumDegPl;
#endif
        
        for(int i = 1; i <= pR; i++) {
            // pivotCount(l + 1, 0, 0, {t.p1, t.h1, t.p2 + pR, t.h2 + i})*C[pR][i];
            t.h2 += i;
#ifdef DEBUG
printf("xadding %d %d %d %d, %d/%d\n", t.p1, t.h1, t.p2, t.h2, i, pR);
#endif

            for(int l = 0; l <= t.p1 && l < maxD2; l++) {
                for(int r = 0; r <= t.p2 && r < maxD2; r++) {
                    if(r + t.h2 >= (int)ansAll[l + t.h1].size()) {
                        ansAll[l + t.h1].resize(r + t.h2 + 5);
                    }
                    ansAll[l + t.h1][r + t.h2] 
                        += C[t.p1][l] * C[t.p2][r] * C[pR][i];
                }
            }

#ifdef DEBUG
if((int)ansAll[x].size() > y)
printf("ansAll[%d][%d] %.2f\n", x, y, ansAll[x][y]);
else printf("ansAll[%d][%d] %.2f\n", x, y, 0.0);
#endif
            t.h2 -= i;
        }

#ifdef PQWEDGE
//at least one vertex in Pr should be selected
for(uint32_t i = 0; i < pR; i++) sumDegPr += g->deg2(candR[i]);
t.p2 += pR;
for(uint32_t i = 0; i < pR; i++) {
    t.h2 += 1;
    sumDegHr += g->deg2(candR[i]);
    sumDegPr -= g->deg2(candR[i]);
    t.p2--;

    for(int ll = 0; ll <= t.p1 && ll < maxD2; ll++) {
        for(int r = 0; r <= t.p2 && r < maxD2; r++) {

        if(r + t.h2 + 1 >= (int)ansWedge[ll + t.h1].size())
            ansWedge[ll + t.h1].resize(r + t.h2 + 6);

        //for each v in h1, it must be selected. \sum_v(dv - (p1+h2))*C()*C()=\sum_vdv-|t.h1|*(p1+h2)
        ansWedge[ll + t.h1][r + t.h2 + 1] += C[t.p1][ll] * C[t.p2][r] 
            * (sumDegHl - t.h1*(r + t.h2));

        if(r + t.h2 >= (int)ansWedge[ll + t.h1 + 1].size())
            ansWedge[ll + t.h1 + 1].resize(r + t.h2 + 5);

        ansWedge[ll + t.h1 + 1][r + t.h2] += C[t.p1][ll] * C[t.p2][r] 
            * (sumDegHr - t.h2*(ll + t.h1));
// printf("type11:%u %u %u %u, %u %u\n", t.p1, t.p2, t.h1, t.h2, ll, r);
        }

    }
// printf("wedge221_e=0, %u:%.f\n", l, ansWedge[2][2]);

    //for each v in p1, when it is selected, dv-(r+t.h2)
    for(int ll = 0; ll < t.p1 && ll < maxD2; ll++) {
        for(int r = 0; r <= t.p2 && r < maxD2; r++) {
            ansWedge[ll + t.h1 + 1][r + t.h2 + 1] += C[t.p1-1][ll] * C[t.p2][r] 
                * (sumDegPl - t.p1*(r + t.h2));
        }
    }
    //for each v in p2, when it is selected, dv-(l+t.h1)
    for(int ll = 0; ll <= t.p1 && ll < maxD2; ll++) {
        for(int r = 0; r < t.p2 && r < maxD2; r++) {
            ansWedge[ll + t.h1 + 1][r + t.h2 + 1] += C[t.p1][ll] * C[t.p2-1][r] 
                * (sumDegPr - t.p2*(ll + t.h1));
        }
    }
// printf("type21:%u %u %u %u\n", t.p1, t.p2, t.h1);
// printf("wedge221_e=0, %u:%.f\n", l, ansWedge[2][2]);
    sumDegHr -= g->deg2(candR[i]);
    t.h2--;
}

// printf("wedge22_e=0, %u:%.f\n", l, ansWedge[2][2]);
#endif
    
        return;
    }

    int pRR = 0;
    if((uint32_t)pR*3 < g->deg1(pivotL)) {
        for(int j = 0; j < pR; j++) {
            if(g->connectUV(pivotL, candR[j])) {
                candR.changeToByPos(j, pRR++);
            }
        }
    }
    else {
        for(uint32_t j = g->pU[pivotL]; j < g->pU[pivotL + 1]; j++) {
            if(candR.idx(g->e1[j]) < (uint32_t)pR) {
                candR.changeTo(g->e1[j], pRR++);
            }
        }
    }

    int pivotR = 0;
    maxDeg = 0;
    for(int i = 0; i < pRR; i++) {
        int deg = 0;
        int v = candR[i];

        if((uint32_t)pL*3 < g->deg2(v)) {
            for(int j = 0; j < pL; j++) {
                if(g->connectVU(v, candL[j])) {
                    deg++;
                }
            }
        }
        else {
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                if(candL.idx(g->e2[j]) < (uint32_t)pL) {
                    deg++;
                }
            }
        }

        if(deg > maxDeg) {
            maxDeg = deg;
            pivotR = v;
        }
    }

    // candL.changeTo(pivotL, --pL);
    if(pRR > 0) candR.changeTo(pivotR, --pRR);
    // candR.changeTo(pivotR, --pR);

    int pLL = 0;
    if((uint32_t)pL*3 <= g->deg2(pivotR)) {
        for(int j = 0; j < pL; j++) {
            if(g->connectVU(pivotR, candL[j])) {
                candL.changeToByPos(j, pLL++);
            }
        }
    }
    else {
        for(uint32_t j = g->pV[pivotR]; j < g->pV[pivotR + 1]; j++) {
            if(candL.idx(g->e2[j]) < (uint32_t)pL) {
                candL.changeTo(g->e2[j], pLL++);
            }
        }
    }
    if(pLL > 0) candL.changeTo(pivotL, --pLL);

#ifdef DEBUG
printf("choose p %u-%u\n", pivotL, pivotR);
#endif
#ifdef PQWEDGE
sumDegPl += g->deg1(pivotL);
sumDegPr += g->deg2(pivotR);
#endif
    pivotCount(l + 1, pLL, pRR, 
        treePath{t.p1 + 1, t.h1, t.p2 + 1, t.h2});
#ifdef PQWEDGE
sumDegPl -= g->deg1(pivotL);
sumDegPr -= g->deg2(pivotR);
#endif

    if((int)tmpNodesL[l].size() < pL - pLL) tmpNodesL[l].resize((pL - pLL)*1.5);
    candL.copy(tmpNodesL[l].data(), pL, pLL + 1);
    if((int)tmpNodesR[l].size() < pR - pRR) tmpNodesR[l].resize((pR - pRR)*1.5);
    candR.copy(tmpNodesR[l].data(), pR, pRR + 1);

    int tmpLSize = pL - pLL - 1;
    int tmpRSize = pR - pRR - 1;
#ifdef DEBUG
printf("tmpLRSize:%d,%d, pLR:%d,%d\n", tmpLSize, tmpRSize, pL, pR);
#endif
    // std::sort(tmpNodesL[l].begin(), tmpNodesL[l].begin() + tmpLSize);
    // std::sort(tmpNodesR[l].begin(), tmpNodesR[l].begin() + tmpRSize);

#ifdef PQWEDGE
for(uint32_t i = 0; i < pL; i++) sumDegPl += g->deg1(candL[i]);
#endif

    for(int i = 0; i < tmpLSize; i++) {
        uint32_t u = tmpNodesL[l][i];
        candL.changeTo(u, --pL);
        pRR = 0;

        for(uint32_t j = g->pU[u]; j < g->pU[u + 1]; j++) {
            uint32_t v = g->e1[j];
            if(candR.idx(v) < (uint32_t)pR) {
                candR.changeTo(v, pRR++);
            }
        }
#ifdef DEBUG
printf("testing h %u, pRR %d\n", u, pRR);
#endif
        t.p1 += pL;
        t.h1 += 1;
        // for(int j = 1; j <= tmpRSize - i; j++) {
            // t.h2 += j;
#ifdef DEBUG
printf("xadding %d %d %d %d, %d/%d\n", t.p1, t.h1, t.p2, t.h2, i, pR);
#endif
#ifdef PQWEDGE
sumDegPl -= g->deg1(u);
sumDegHl += g->deg1(u);
#endif
        for(int l = 0; l <= t.p1 && l < maxD2; l++) {
            for(int r = 0; r <= t.p2 && r < maxD2; r++) {
                if(r + t.h2 >= (int)ansAll[l + t.h1].size()) {
                    ansAll[l + t.h1].resize(r + t.h2 + 5);
                }
                ansAll[l + t.h1][r + t.h2] 
                    += C[t.p1][l] * C[t.p2][r];
#ifdef PQWEDGE
if(r + t.h2 + 1 >= (int)ansWedge[l + t.h1].size())
    ansWedge[l + t.h1].resize(r + t.h2 + 6);

//for each v in h1, it must be selected. \sum_v(dv - (p1+h2))*C()*C()=\sum_vdv-|t.h1|*(p1+h2)
ansWedge[l + t.h1][r + t.h2 + 1] += C[t.p1][l] * C[t.p2][r] 
    * (sumDegHl - t.h1*(r + t.h2));

if(r + t.h2 >= (int)ansWedge[l + t.h1 + 1].size())
    ansWedge[l + t.h1 + 1].resize(r + t.h2 + 5);

ansWedge[l + t.h1 + 1][r + t.h2] += C[t.p1][l] * C[t.p2][r] 
    * (sumDegHr - t.h2*(l + t.h1));

//for each v in p1, when it is selected
// printf("%u %u %u %u\n", t.p1, t.p2, t.h1, t.h2);
#endif
            }
        }
#ifdef PQWEDGE
//for each v in p1, when it is selected, dv-(r+t.h2)
for(int ll = 0; ll < t.p1 && ll < maxD2; ll++) {
    for(int r = 0; r <= t.p2 && r < maxD2; r++) {
        if(r + t.h2 + 1 >= (int)ansWedge[ll + t.h1 + 1].size())
            ansWedge[ll + t.h1 + 1].resize(r + t.h2 + 5);
        ansWedge[ll + t.h1 + 1][r + t.h2 + 1] += C[t.p1-1][ll] * C[t.p2][r] 
            * (sumDegPl - t.p1*(r + t.h2));
    }
}
//for each v in p2, when it is selected, dv-(l+t.h1)
for(int ll = 0; ll <= t.p1 && ll < maxD2; ll++) {
    for(int r = 0; r < t.p2 && r < maxD2; r++) {
        if(r + t.h2 + 1 >= (int)ansWedge[ll + t.h1 + 1].size())
            ansWedge[ll + t.h1 + 1].resize(r + t.h2 + 5);
        ansWedge[ll + t.h1 + 1][r + t.h2 + 1] += C[t.p1][ll] * C[t.p2-1][r] 
            * (sumDegPr - t.p2*(ll + t.h1));
    }
}
#endif
            // t.h2 -= j;
        // }
        t.p1 -= pL;
        t.h1 -= 1;
#ifdef PQWEDGE
sumDegHl -= g->deg1(u);
#endif

        for(uint32_t j = g->pU[u]; j < g->pU[u + 1]; j++) {
            uint32_t v = g->e1[j];
            if(candR.idx(v) < (uint32_t)pR) {
                pLL = 0;
                for(int k = 0; k < pL; k++) {
                    if(g->connectVU(v, candL[k])) {
                        candL.changeToByPos(k, pLL++);
                    }
                }
                candR.changeTo(v, --pRR);
#ifdef DEBUG
printf("choose h %u-%u\n", u, v);
#endif
#ifdef PQWEDGE
sumDegHl += g->deg1(u);
sumDegHr += g->deg2(v);
#endif
                pivotCount(l + 1, pLL, pRR, 
                    treePath{t.p1, t.h1 + 1, t.p2, t.h2 + 1});
#ifdef PQWEDGE
sumDegHl -= g->deg1(u);
sumDegHr -= g->deg2(v);
#endif
            }
        }
    }

#ifdef PQWEDGE
for(uint32_t i = 0; i < pR; i++) sumDegPr += g->deg2(candR[i]);
#endif
    for(int i = 0; i < tmpRSize; i++) {
        uint32_t v = tmpNodesR[l][i];
        candR.changeTo(v, --pR);
        pLL = 0;

        for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
            uint32_t u = g->e2[j];
            if(candL.idx(u) < (uint32_t)pL) {
                candL.changeTo(u, pLL++);
            }
        }

        t.p2 += pR;//pivotR + [i+1,tmpRSize-1]
        t.h2 += 1;
        // for(int j = 1; j <= tmpRSize - i; j++) {
            // t.h2 += j;
#ifdef DEBUG
printf("xadding %d %d %d %d, %d/%d\n", t.p1, t.h1, t.p2, t.h2, i, pR);
#endif
#ifdef PQWEDGE
sumDegPr -= g->deg2(v);
sumDegHr += g->deg2(v);
#endif
        for(int l = 0; l <= t.p1 && l < maxD2; l++) {
            for(int r = 0; r <= t.p2 && r < maxD2; r++) {
                if(r + t.h2 >= (int)ansAll[l + t.h1].size()) {
                    ansAll[l + t.h1].resize(r + t.h2 + 5);
                }
                ansAll[l + t.h1][r + t.h2] 
                    += C[t.p1][l] * C[t.p2][r];
#ifdef PQWEDGE
if(r + t.h2 + 1 >= (int)ansWedge[l + t.h1].size())
    ansWedge[l + t.h1].resize(r + t.h2 + 6);

//for each v in h1, it must be selected. \sum_v(dv - (p1+h2))*C()*C()=\sum_vdv-|t.h1|*(p1+h2)
ansWedge[l + t.h1][r + t.h2 + 1] += C[t.p1][l] * C[t.p2][r] 
    * (sumDegHl - t.h1*(r + t.h2));

if(r + t.h2 >= (int)ansWedge[l + t.h1 + 1].size())
    ansWedge[l + t.h1 + 1].resize(r + t.h2 + 5);

ansWedge[l + t.h1 + 1][r + t.h2] += C[t.p1][l] * C[t.p2][r] 
    * (sumDegHr - t.h2*(l + t.h1));

//for each v in p1, when it is selected

#endif
            }
        }
#ifdef PQWEDGE
//for each v in p1, when it is selected, dv-(r+t.h2)
for(int ll = 0; ll < t.p1 && ll < maxD2; ll++) {
    for(int r = 0; r <= t.p2 && r < maxD2; r++) {
        ansWedge[ll + t.h1 + 1][r + t.h2 + 1] += C[t.p1-1][ll] * C[t.p2][r] 
            * (sumDegPl - t.p1*(r + t.h2));
    }
}
//for each v in p2, when it is selected, dv-(l+t.h1)
for(int ll = 0; ll <= t.p1 && ll < maxD2; ll++) {
    for(int r = 0; r < t.p2 && r < maxD2; r++) {
        ansWedge[ll + t.h1 + 1][r + t.h2 + 1] += C[t.p1][ll] * C[t.p2-1][r] 
            * (sumDegPr - t.p2*(ll + t.h1));
    }
}
#endif
            // t.h2 -= j;
        // }
        t.p2 -= pR;
        t.h2 -= 1;
#ifdef PQWEDGE
sumDegHr -= g->deg2(v);
#endif

        for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
            uint32_t u = g->e2[j];
            if(candL.idx(u) < (uint32_t)pL && g->connectUV(u, pivotR)) {
                pRR = 0;
                for(int k = 0; k < pR; k++) {
                    if(g->connectUV(u, candR[k])) {
                        candR.changeToByPos(k, pRR++);
                    }
                }
                candL.changeTo(u, --pLL);
#ifdef DEBUG
printf("choose h %u-%u\n", u, v);
#endif
#ifdef PQWEDGE
sumDegHl += g->deg1(u);
sumDegHr += g->deg2(v);
#endif
                pivotCount(l + 1, pLL, pRR, 
                    treePath{t.p1, t.h1 + 1, t.p2, t.h2 + 1});
#ifdef PQWEDGE
sumDegHl -= g->deg1(u);
sumDegHr -= g->deg2(v);
#endif
            }
        }
    }
}
