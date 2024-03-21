#include "edgePivotSpecificPQ.h"

void edgePivotSpecificPQ::sparseCounting(uint32_t bar) {
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
    // printf("sample size %llu\n", T);
    // printf("upperBound %llu\n", upperBoundEdgesForSampling);

    printf("exact:%u\n", (uint32_t)vs1.size());
    printf("sample:%u\n", (uint32_t)vs2.size());
    fflush(stdout);

double t = clock();
if(vs1.size() > 0) {

    candL.resize(g->n1);
    candR.resize(g->n2);
    uint32_t maxD = std::max(g->maxDu, g->maxDv);
    tmpNodesL.resize(maxD);
    tmpNodesR.resize(maxD);

    printf("maxDDD:%u\n", std::max(g->maxDu, g->maxDv));
    // fflush(stdout);

    for(auto u : vs1) { 
        if(g->deg1(u) == 0) continue;
// printf("%u\n", u);fflush(stdout);
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }
// printf("%u\n", u);fflush(stdout);
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            int pL = 0;

            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) candL.changeTo(w, pL++);
            }
// printf("%u\n", u);fflush(stdout);
            candR.changeTo(v, --pR);
#ifdef DEBUG
printf("st choose h %u-%u\n", u, v);fflush(stdout);
#endif
            pivotCountFast(1, pL, pR, treePath{0, 1, 0, 1});
        }
    }

    printf("done\n");

}
double d = clock();

printf("exact done, %.5fs, ans %.0f\n", (d - t) / CLOCKS_PER_SEC, ansAll[p][q]);
}

void edgePivotSpecificPQ::exactCountMaximalPivotFast() {
#ifdef DEBUG
g->print();
#endif

    candL.resize(g->n1);
    candR.resize(g->n2);
    uint32_t maxD = std::max(g->maxDu, g->maxDv);
    tmpNodesL.resize(maxD);
    tmpNodesR.resize(maxD);

    printf("maxDDD:%u\n", std::max(g->maxDu, g->maxDv));
    fflush(stdout);

    for(uint32_t u = 0; u < g->n1; u++) { 
        if(g->deg1(u) == 0) continue;
// printf("%u\n", u);fflush(stdout);
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }
// printf("%u\n", u);fflush(stdout);
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            int pL = 0;

            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) candL.changeTo(w, pL++);
            }
// printf("%u\n", u);fflush(stdout);
            candR.changeTo(v, --pR);
#ifdef DEBUG
printf("st choose h %u-%u\n", u, v);fflush(stdout);
#endif
            pivotCountFast(1, pL, pR, treePath{0, 1, 0, 1});
        }
    }

    for(int i = 1; i <= p; i++) {
        for(int j = 1; j <= q; j++) {
            if(ansAll[i][j] > 0.0) {
                printf("%d-%d:%.0f\n", i, j, ansAll[i][j]);
            }
        }
    }
}



void edgePivotSpecificPQ::pivotCountFast(int l, int pL, int pR, treePath t) {
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
    if(t.h1 > p || t.h2 > q) return;

    auto update = [&]() {
        for(int ll = 0; ll + t.h1 <= p; ll++) {
            for(int r = 0; r + t.h2 <= q; r++) {
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
        // pivotCountFast(l + 1, 0, 0, {t.p1 + pL, t.h1, t.p2, t.h2});
        t.p1 += pL;
        update();
        t.p1 -= pL;

        // t.p2 += pR;
        // for(int l = 0; l + t.h1 <= p; l++) {
        //     for(int r = 1; r + t.h2 <= q; r++) {
        //         ansAll[l + t.h1][r + t.h2] += C[t.p1][l] * C[t.p2][r];
        //     }
        // }
        // t.p2 -= pR;
        for(int i = 1; i <= pR; i++) {
            // pivotCount(l + 1, 0, 0, {t.p1, t.h1, t.p2 + pR, t.h2 + i})*C[pR][i];
            t.h2 += i;
            if(t.h2 > q) break;
#ifdef DEBUG
printf("xadding %d %d %d %d, %d/%d\n", t.p1, t.h1, t.p2, t.h2, i, pR);
#endif
            for(int l = 0; l + t.h1 <= p; l++) {
                for(int r = 0; r + t.h2 <= q; r++) {
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
        for(int l = 0; l + t.h1 <= p; l++) {
            for(int r = 0; r + t.h2 <= q; r++) {
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
        for(int l = 0; l + t.h1 <= p; l++) {
            for(int r = 0; r + t.h2 <= q; r++) {
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
