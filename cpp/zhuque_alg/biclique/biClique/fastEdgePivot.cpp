#include "fastEdgePivot.h"

void fastEdgePivot::exactCountMaximalPivotFast() {
#ifdef DEBUG
g->print();
#endif

    candL.resize(g->n1);
    candR.resize(g->n2);
    uint32_t maxD = std::max(g->maxDu, g->maxDv);
    tmpNodesL.resize(maxD);
    tmpNodesR.resize(maxD);

    ansAll.resize(g->maxDv + 1);
    ansAll[1].resize(g->maxDu + 1);
    // for(uint32_t i = 1; i <= 10 && i <= g->maxDv; i++) {
    //     ansAll[i].resize(g->maxDu + 1);
    // }
    // for(uint32_t i = 11; i <= g->maxDv; i++) {
    //     ansAll[i].resize(10);
    // }

    printf("maxD:%u\n", std::max(g->maxDu, g->maxDv));
    fflush(stdout);

    for(uint32_t u = 0; u < g->n1; u++) {
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }
        if(pR == 0) continue;

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            int pL = 0;

            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) candL.changeTo(w, pL++);
            }

            candR.changeTo(v, --pR);
#ifdef DEBUG
printf("st choose h %u-%u\n", u, v);
#endif
            pivotCountFast(1, pL, pR, treePath{0, 1, 0, 1});
        }
    }

    for(int i = 1; i < (int)ansAll.size(); i++) {
        for(int j = 1; j < (int)ansAll[i].size(); j++) {
            if(ansAll[i][j] > 0.0) {
                printf("%d-%d:%.0f\n", i, j, ansAll[i][j]);
            }
        }
    }
}

#ifdef DEBUG
int x = 3, y = 1;
#endif

void fastEdgePivot::pivotCountFast(int l, int pL, int pR, treePath t) {
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

    auto update = [&]() {
        for(int ll = 0; ll <= t.p1; ll++) {
            for(int r = 0; r <= t.p2; r++) {
                if(r + t.h2 >= (int)ansAll[ll + t.h1].size()) {
                    ansAll[ll + t.h1].resize(r + t.h2 + 1);
                }
                // ansAll[ll + t.h1][r + t.h2] += Cnm(t.p1, ll) * Cnm(t.p2, r);
                ansAll[ll + t.h1][r + t.h2] += C[t.p1][ll] * C[t.p2][r];
            }
        }
    };

    if(pL == 0 && pR == 0) {
#ifdef DEBUG
printf("adding %d %d %d %d\n", t.p1, t.h1, t.p2, t.h2);
#endif
        update();
#ifdef DEBUG
if((int)ansAll[x].size() > y)
printf("ansAll[%d][%d] %.2f\n", x, y, ansAll[x][y]);
else printf("ansAll[%d][%d] %.2f\n", x, y, 0.0);
#endif
        return ;
    }

    // if(pL == 1 && pR == 1) {
    //     if(g->connectUVFast(candL[0], candR[0])) {
    //         t.p1++;
    //         t.p2++;
    //         update();
    //     }
    //     else {
    //         t.p1++;
    //         update();
    //         t.p1--;
    //         t.h2++;
    //         update();
    //     }

    //     return ;
    // }

    int pivotL = 0;
    int maxDeg = 0;
    for(int i = 0; i < pL; i++) {
        int deg = 0;

        uint32_t u = candL[i];

        if((uint32_t)pR < g->deg1(u)) {
            for(int j = 0; j < pR; j++) {
                if(g->connectUVFast(u, candR[j])) {
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
        t.p1 += pL;
        update();
        t.p1 -= pL;

        // pivotCountFast(l + 1, 0, 0, {t.p1 + pL, t.h1, t.p2, t.h2});
        for(int i = 1; i <= pR; i++) {
            // pivotCount(l + 1, 0, 0, {t.p1, t.h1, t.p2 + pR, t.h2 + i})*C[pR][i];
            t.h2 += i;
#ifdef DEBUG
printf("xadding %d %d %d %d, %d/%d\n", t.p1, t.h1, t.p2, t.h2, i, pR);
#endif
            for(int l = 0; l <= t.p1; l++) {
                for(int r = 0; r <= t.p2; r++) {
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
    
        return;
    }

    int pRR = 0;
    if((uint32_t)(pR+pR) < g->deg1(pivotL)) {
        for(int j = 0; j < pR; j++) {
            if(g->connectUVFast(pivotL, candR[j])) {
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

        if((uint32_t)(pL+pL) < g->deg2(v)) {
            for(int j = 0; j < pL; j++) {
                if(g->connectVUFast(v, candL[j])) {
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
    if((uint32_t)(pL+pL) <= g->deg2(pivotR)) {
        for(int j = 0; j < pL; j++) {
            if(g->connectVUFast(pivotR, candL[j])) {
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
    pivotCountFast(l + 1, pLL, pRR, 
        treePath{t.p1 + 1, t.h1, t.p2 + 1, t.h2});


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

    for(int i = 0; i < tmpLSize; i++) {
        uint32_t u = tmpNodesL[l][i];
        candL.changeTo(u, --pL);
        pRR = 0;

        // if((uint32_t)pR < g->deg1(u)) {
            // for(int i = 0; i < pR; i++) {
            //     if(g->connectUVFast(u, candR[i])) {
            //         candR.changeTo(candR[i], pRR++);
            //     }
            // }
        // }
        // else {
            for(uint32_t j = g->pU[u]; j < g->pU[u + 1]; j++) {
                uint32_t v = g->e1[j];
                if(candR.idx(v) < (uint32_t)pR) {
                    candR.changeTo(v, pRR++);
                }
            }
        // }
        
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
        for(int l = 0; l <= t.p1; l++) {
            for(int r = 0; r <= t.p2; r++) {
                if(r + t.h2 >= (int)ansAll[l + t.h1].size()) {
                    ansAll[l + t.h1].resize(r + t.h2 + 5);
                }
                // ansAll[l + t.h1][r + t.h2] += Cnm(t.p1, l) * Cnm(t.p2, r);
                ansAll[l + t.h1][r + t.h2] += C[t.p1][l] * C[t.p2][r];
            }
        }
            // t.h2 -= j;
        // }
        t.p1 -= pL;
        t.h1 -= 1;

        for(uint32_t j = g->pU[u]; j < g->pU[u + 1]; j++) {
            uint32_t v = g->e1[j];
            if(candR.idx(v) < (uint32_t)pR) {
                pLL = 0;
                for(int k = 0; k < pL; k++) {
                    if(g->connectVUFast(v, candL[k])) {
                        candL.changeToByPos(k, pLL++);
                    }
                }
                candR.changeTo(v, --pRR);
#ifdef DEBUG
printf("choose h %u-%u\n", u, v);
#endif
                pivotCountFast(l + 1, pLL, pRR, 
                    treePath{t.p1, t.h1 + 1, t.p2, t.h2 + 1});
            }
        }
    }

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
        for(int l = 0; l <= t.p1; l++) {
            for(int r = 0; r <= t.p2; r++) {
                if(r + t.h2 >= (int)ansAll[l + t.h1].size()) {
                    ansAll[l + t.h1].resize(r + t.h2 + 5);
                }
                // ansAll[l + t.h1][r + t.h2] += Cnm(t.p1, l) * Cnm(t.p2, r);
                ansAll[l + t.h1][r + t.h2] += C[t.p1][l] * C[t.p2][r];
            }
        }
            // t.h2 -= j;
        // }
        t.p2 -= pR;
        t.h2 -= 1;

        for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
            uint32_t u = g->e2[j];
            if(candL.idx(u) < (uint32_t)pL && g->connectUVFast(u, pivotR)) {
                pRR = 0;
                for(int k = 0; k < pR; k++) {
                    if(g->connectUVFast(u, candR[k])) {
                        candR.changeToByPos(k, pRR++);
                    }
                }
                candL.changeTo(u, --pLL);
#ifdef DEBUG
printf("choose h %u-%u\n", u, v);
#endif
                pivotCountFast(l + 1, pLL, pRR, 
                    treePath{t.p1, t.h1 + 1, t.p2, t.h2 + 1});
            }
        }
    }
}
