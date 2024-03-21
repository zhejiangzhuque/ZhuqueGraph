#include "rawEdgePivot2.h"

#include <chrono>
#include <random>

// #define DENSE_DEBUG

void rawEdgePivot2::addEdgesType2and3(Dinic * dinic) {
    this->dinic = dinic;
    PLtmp.resize(p);
    PRtmp.resize(q);
    PLtmp.clear();
    PRtmp.clear();
    comNeiL.resize(g->n1);
    comNeiR.resize(g->n2);
    // idxComNeiL.resize(p + q + 20);
    // idxComNeiR.resize(p + q + 20);

    for(uint32_t u = 0; u < g->n1; u++) {
        int pR = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, pR++);
        }
        if(pR == 0) continue;
        
HL.push_back(u);
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            int pL = 0;

            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) candL.changeTo(w, pL++);
            }

            candR.changeTo(v, --pR);

            HR.push_back(v);
            pivotForDensestSubgraph(1, pL, pR, treePath{0, 1, 0, 1});
            HR.pop_back();
        }
HL.pop_back();
    }

    printf("cliqueId:%u\n", cliqueId);
}


void rawEdgePivot2::listAllPminus1QBicliquesL(int i, int l, treePath t, uint32_t idxl) {
    if(l + t.h1 == p - 1) {
        listAllPminus1QBicliquesR(0, 0, t, idxl);
        return;
    }

    if(t.p1 - i + l + t.h1 < p - 1) return;

    for(int j = i; j < t.p1; j++) {
        PLtmp.push_back(PL[j]);
        listAllPminus1QBicliquesL(j + 1, l + 1, t, idxl);
        PLtmp.pop_back();
    }
}
void rawEdgePivot2::listAllPminus1QBicliquesR(int i, int r, treePath t, uint32_t idxl) {
    if(r + t.h2 == q) {
        uint32_t i = 0, j = 0;
        for(auto u : PLtmp) comNeiL.changeTo(u, i++);
        for(auto u : HL) comNeiL.changeTo(u, i++);

        if(i == idxl) return;

        j = i;
        while(j < idxl) {
            dinic->addEdge(comNeiL[j++], cliqueId + g->n1 + g->n2, 1);
        }

        for(auto u : PLtmp) dinic->addEdge(cliqueId + g->n1 + g->n2, u, INF);
        for(auto u : HL) dinic->addEdge(cliqueId + g->n1 + g->n2, u, INF);
        for(auto v : PRtmp) dinic->addEdge(cliqueId + g->n1 + g->n2, v + g->n1, INF);
        for(auto v : HR) dinic->addEdge(cliqueId + g->n1 + g->n2, v + g->n1, INF);
       
#ifdef DENSE_DEBUG
printf("cliqueId_%u:", cliqueId+ g->n1 + g->n2);
printf("L:");
for(auto u : PLtmp) printf(" %u", u);
for(auto u : HL) printf(" %u", u);
printf(" R:");
for(auto v : PRtmp) printf(" %u", v);
for(auto v : HR) printf(" %u", v);
printf("\n");
#endif
        cliqueId++;

        return;
    }

    if(t.p2 - i + r + t.h2 < q) return;

    for(int j = i; j < t.p2; j++) {
        PRtmp.push_back(PR[j]);
        uint32_t newIdxL = updateComNeiL(idxl, PR[j]);
        listAllPminus1QBicliquesR(j + 1, r + 1, t, newIdxL);
        PRtmp.pop_back();
    }
}
void rawEdgePivot2::listAllPQminus1BicliquesL(int i, int l, treePath t, uint32_t idxr) {
    if(l + t.h1 == p) {
        listAllPQminus1BicliquesR(0, 0, t, idxr);
        return;
    }

    if(t.p1 - i + l + t.h1 < p) return;

    for(int j = i; j < t.p1; j++) {
        PLtmp.push_back(PL[j]);
        uint32_t newIdxR = updateComNeiR(idxr, PL[j]);
        listAllPQminus1BicliquesL(j + 1, l + 1, t, newIdxR);
        PLtmp.pop_back();
    }
}
void rawEdgePivot2::listAllPQminus1BicliquesR(int i, int r, treePath t, uint32_t idxr) {
    if(r + t.h2 == q - 1) {
        uint32_t i = 0, j = 0;
        for(auto v : PRtmp) comNeiR.changeTo(v, i++);
        for(auto v : HR) comNeiR.changeTo(v, i++);

        if(i == idxr) return;
        j = i;

        while(j < idxr) {
            dinic->addEdge(comNeiR[j++] + g->n1, cliqueId + g->n1 + g->n2, 1);
        }

        for(auto u : PLtmp) dinic->addEdge(cliqueId + g->n1 + g->n2, u, INF);
        for(auto u : HL) dinic->addEdge(cliqueId + g->n1 + g->n2, u, INF);
        for(auto v : PRtmp) dinic->addEdge(cliqueId + g->n1 + g->n2, v + g->n1, INF);
        for(auto v : HR) dinic->addEdge(cliqueId + g->n1 + g->n2, v + g->n1, INF);
#ifdef DENSE_DEBUG
printf("cliqueId_%u:", cliqueId+ g->n1 + g->n2);
printf("L:");
for(auto u : PLtmp) printf(" %u", u);
for(auto u : HL) printf(" %u", u);
printf(" R:");
for(auto v : PRtmp) printf(" %u", v);
for(auto v : HR) printf(" %u", v);
printf("\n");
#endif
        cliqueId++;

        return;
    }

    if(t.p2 - i + r + t.h2 < q - 1) return;

    for(int j = i; j < t.p2; j++) {
        PRtmp.push_back(PR[j]);
        listAllPQminus1BicliquesR(j + 1, r + 1, t, idxr);
        PRtmp.pop_back();
    }
}


void rawEdgePivot2::pivotForDensestSubgraph(int l, int pL, int pR, treePath t) {
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

//prune the branches can not (p-1,q) or (p, q-1)
    if(t.h1 > p || t.h2 > q) return;
    if(t.h1 == p && t.h2 == q) return;
    if(pL + t.h1 + t.p1 < p-1) return;
    if(pR + t.h2 + t.p2 < q-1) return;

    if(pL == 0 && pR == 0) {
#ifdef DENSE_DEBUG
printf("adding %d %d %d %d\n", t.p1, t.h1, t.p2, t.h2);
// printf("adding %d %d %d %d, pq %u %u, pqcnt %.0f\n", t.p1, t.h1, t.p2, t.h2, p, q, 
//     C[t.p1][p - t.h1] * C[t.p2][q - t.h2]);
// printf("PL.sz %u, PR.sz %u, %u %u\n", PL.size(), PR.size(), HL.size(), HR.size());
assert(PL.size() == t.p1);
assert(PR.size() == t.p2);
assert(HL.size() == t.h1);
assert(HR.size() == t.h2);

fflush(stdout);
printf("pri cliqueId:%u, edges %u\n", cliqueId, dinic->getEdges());
#endif
        //list all (p-1,q)-bicliques that contain all h nodes
        uint32_t idxl = g->n1, idxr = g->n2;
        if(t.h1 <= p - 1) {
            for(uint32_t i = 0; i < t.h2; i++) idxl = updateComNeiL(idxl, HR[i]);
#ifdef DENSE_DEBUG
printf("HcomNeiL:%u,,", idxl);
for(uint32_t i = 0; i < idxl; i++) printf("%u ", comNeiL[i]);printf("\n");
#endif
            listAllPminus1QBicliquesL(0, 0, t, idxl);
// printf("aftL cliqueId:%u, should add %.0f, edges %u\n", 
//     cliqueId, C[t.p1][p-1-t.h1]*C[t.p2][q-t.h2], dinic->getEdges());
        }
        //list all (p,q-1)-bicliques that contain all h nodes
        if(t.h2 <= q - 1) {
            for(uint32_t i = 0; i < t.h1; i++) idxr = updateComNeiR(idxr, HL[i]);
#ifdef DENSE_DEBUG
printf("HcomNeiR:%u,,", idxr);
for(uint32_t i = 0; i < idxr; i++) printf("%u ", comNeiR[i]);printf("\n");
#endif
            listAllPQminus1BicliquesL(0, 0, t, idxr);
// printf("aftR cliqueId:%u, should add %.0f, edges %u\n", 
//     cliqueId, C[t.p1][p-t.h1]*C[t.p2][q-1-t.h2], dinic->getEdges());

        }    

// printf("cliqueId:%u, %u\n", cliqueId, cliqueId + g->n1 +g->n2);


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
        for(uint32_t i = 0; i < pL; i++) PL.push_back(candL[i]);
        pivotForDensestSubgraph(l + 1, 0, 0, {t.p1 + pL, t.h1, t.p2, t.h2});
        for(uint32_t i = 0; i < pL; i++) PL.pop_back();

        if(pR == 0) return;
        for(int i = pR - 1; i >= 0; i--) {
            PR.push_back(candR[i]);
        }
        
        for(int i = 0; i < pR; i++) {
            HR.push_back(candR[i]);
            PR.pop_back();
// printf("%u %u\n", PR.size(), t.p2 + pR - i - 1);fflush(stdout);
            pivotForDensestSubgraph(l + 1, 0, 0, {t.p1, t.h1, t.p2 + pR - i - 1, t.h2 + 1});
            HR.pop_back();
        }

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

PL.push_back(pivotL);
PR.push_back(pivotR);
    pivotForDensestSubgraph(l + 1, pLL, pRR, 
        treePath{t.p1 + 1, t.h1, t.p2 + 1, t.h2});
PL.pop_back();
PR.pop_back();


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
    for(int i = 0; i < pLL + 1; i++) {
        uint32_t u = candL[i];
        PL.push_back(u);
    }
    for(int i = tmpLSize - 1; i >= 0; i--) {
        uint32_t u = tmpNodesL[l][i];
        PL.push_back(u);
    }
    for(int i = 0; i < tmpLSize; i++) {
        uint32_t u = tmpNodesL[l][i];
        PL.pop_back();
        HL.push_back(u);
        pivotForDensestSubgraph(l + 1, 0, 0, treePath{t.p1 + pL - i - 1, t.h1 + 1, t.p2, t.h2});
        HL.pop_back();
    }
    for(int i = 0; i < pLL + 1; i++) PL.pop_back();

    for(int i = 0; i < pRR + 1; i++) {
        uint32_t u = candR[i];
        PR.push_back(u);
    }
    for(int i = tmpRSize - 1; i >= 0; i--) {
        uint32_t u = tmpNodesR[l][i];
        PR.push_back(u);
    }
    for(int i = 0; i < tmpRSize; i++) {
        uint32_t v = tmpNodesR[l][i];
        PR.pop_back();
        HR.push_back(v);
        pivotForDensestSubgraph(l + 1, 0, 0, treePath{t.p1, t.h1, t.p2 + pR - i - 1, t.h2 + 1});
        HR.pop_back();
    }
    for(int i = 0; i < pRR + 1; i++) PR.pop_back();



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
        // t.p1 += pL;
        // t.h1 += 1;
        // for(int l = 0; l <= t.p1 && l < maxD2; l++) {
        //     for(int r = 0; r <= t.p2 && r < maxD2; r++) {
        //         if(r + t.h2 >= (int)ansAll[l + t.h1].size()) {
        //             ansAll[l + t.h1].resize(r + t.h2 + 5);
        //         }
        //         ansAll[l + t.h1][r + t.h2] 
        //             += C[t.p1][l] * C[t.p2][r];
        //     }
        // }
        // t.p1 -= pL;
        // t.h1 -= 1;


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
HL.push_back(u);
HR.push_back(v);
                pivotForDensestSubgraph(l + 1, pLL, pRR, 
                    treePath{t.p1, t.h1 + 1, t.p2, t.h2 + 1});
HL.pop_back();
HR.pop_back();
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

//         t.p2 += pR;//pivotR + [i+1,tmpRSize-1]
//         t.h2 += 1;
//         // for(int j = 1; j <= tmpRSize - i; j++) {
//             // t.h2 += j;
// #ifdef DEBUG
// printf("xadding %d %d %d %d, %d/%d\n", t.p1, t.h1, t.p2, t.h2, i, pR);
// #endif
//         for(int l = 0; l <= t.p1 && l < maxD2; l++) {
//             for(int r = 0; r <= t.p2 && r < maxD2; r++) {
//                 if(r + t.h2 >= (int)ansAll[l + t.h1].size()) {
//                     ansAll[l + t.h1].resize(r + t.h2 + 5);
//                 }
//                 ansAll[l + t.h1][r + t.h2] 
//                     += C[t.p1][l] * C[t.p2][r];
//             }
//         }
//             // t.h2 -= j;
//         // }
//         t.p2 -= pR;
//         t.h2 -= 1;

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
HL.push_back(u);
HR.push_back(v);
                pivotForDensestSubgraph(l + 1, pLL, pRR, 
                    treePath{t.p1, t.h1 + 1, t.p2, t.h2 + 1});
HL.pop_back();
HR.pop_back();
            }
        }
    }
}


void rawEdgePivot2::localCount() {
//     if(p == 1) {
//         for(uint32_t u = 0; u < g->n1; u++) {
//             if(g->deg1(u) >= q) (*localAns)[u] = C[g->deg1(u)][q];
//             for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
//                 uint32_t v = g->e1[i] + g->n1;
//                 (*localAns)[v] += C[g->deg1(u) - 1][q - 1];
//             }
//         }

//         return;
//     }
//     else if(q == 1) {
// g->print();
// printf("pq %u %u\n", p, q);
//         for(uint32_t vv = 0; vv < g->n2; vv++) {
//             uint32_t v = vv + g->n1;
//             if(g->deg1(vv) >= p) (*localAns)[vv] = C[g->deg1(vv)][p];
//             for(uint32_t i = g->pV[vv]; i < g->pV[vv + 1]; i++) {
//                 uint32_t u = g->e2[i];
//                 (*localAns)[u] += C[g->deg2(vv) - 1][p - 1];
//             }
//         }
//         return;
//     }
// g->print();
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

            HL.push_back(u); HR.push_back(v);
            pivotCount(1, pL, pR, treePath{0, 1, 0, 1});
            HL.pop_back(); HR.pop_back();
        }
    }

    // printf("totalCount:%.0f\n", ans);
}

void rawEdgePivot2::pivotCount(int l, int pL, int pR, treePath t) {
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
    if(t.h1 > p || t.h2 > q) return;
    if(pL + t.h1 + t.p1 < p) return;
    if(pR + t.h2 + t.p2 < q) return;

    if(pL == 0 && pR == 0) {
        // t.p1 += pL;
        // t.p2 += pR;
#ifdef DEBUG
printf("adding %d %d %d %d\n", t.p1, t.h1, t.p2, t.h2);
#endif
        // for(int ll = 0; ll <= t.p1 && ll < maxD2; ll++) {
        //     for(int r = 0; r <= t.p2 && r < maxD2; r++) {
        //         if(r + t.h2 >= (int)ansAll[ll + t.h1].size()) {
        //             ansAll[ll + t.h1].resize(r + t.h2 + 4);
        //         }
        //         ansAll[ll + t.h1][r + t.h2] += C[t.p1][ll] * C[t.p2][r];
        //     }
        // }
        ans += C[t.p1][p - t.h1] * C[t.p2][q - t.h2];
#ifdef DEBUGDENSE
printf("adding %d %d %d %d, pq %u %u\n", t.p1, t.h1, t.p2, t.h2, p, q);
printf("PL.sz %u, PR.sz %u, %u %u\n", PL.size(), PR.size(), HL.size(), HR.size());
assert(PL.size() == t.p1);
assert(PR.size() == t.p2);
assert(HL.size() == t.h1);
assert(HR.size() == t.h2);
#endif
        if(t.h1 < p)
        for(auto u : PL) {
            (*localAns)[u] += C[t.p1-1][p-1 - t.h1] * C[t.p2][q - t.h2];
#ifdef DEBUGDENSE
printf("addpl %u:%.0f\n", u, C[t.p1-1][p-1 - t.h1] * C[t.p2][q - t.h2]);
#endif
        }
// printf("adding %d %d %d %d, %.0f, %.0f\n",
//      t.p1, t.h1, t.p2, t.h2, C[t.p1][p - t.h1] * C[t.p2][q - t.h2], ans);
// for(int i = 0; i < PR.size(); i++) printf("%u ", PR[i] + g->n1);printf("\n");fflush(stdout);
// for(auto v : PR) printf("%u ", v + g->n1);printf("\n");fflush(stdout);
// printf("%d %d %d %d", t.p1, p - t.h1, t.p2-1, q-1 - t.h2);printf("\n");fflush(stdout);
        if(t.h2 < q)
        for(auto v : PR) {
            (*localAns)[v + g->n1] += C[t.p1][p - t.h1] * C[t.p2-1][q-1 - t.h2];
#ifdef DEBUGDENSE
printf("addpr %u:%.0f\n", v, C[t.p1][p - t.h1] * C[t.p2-1][q-1 - t.h2]);
#endif
        }
        for(auto u : HL) {
            (*localAns)[u] += C[t.p1][p - t.h1] * C[t.p2][q - t.h2];
#ifdef DEBUGDENSE
printf("addhl %u:%.0f\n", u, C[t.p1][p - t.h1] * C[t.p2][q - t.h2]);
#endif
        }
        for(auto v : HR) {
            (*localAns)[v + g->n1] += C[t.p1][p - t.h1] * C[t.p2][q - t.h2];
#ifdef DEBUGDENSE
printf("addhr %u:%.0f\n", v, C[t.p1][p - t.h1] * C[t.p2][q - t.h2]);
#endif
        }

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
        for(uint32_t i = 0; i < pL; i++) PL.push_back(candL[i]);
        pivotCount(l + 1, 0, 0, {t.p1 + pL, t.h1, t.p2, t.h2});
        for(uint32_t i = 0; i < pL; i++) PL.pop_back();

        
//         for(int i = 1; i <= pR; i++) {
//             // pivotCount(l + 1, 0, 0, {t.p1, t.h1, t.p2 + pR, t.h2 + i})*C[pR][i];
//             t.h2 += i;
// #ifdef DEBUG
// printf("xadding %d %d %d %d, %d/%d\n", t.p1, t.h1, t.p2, t.h2, i, pR);
// #endif

//             // for(int l = 0; l <= t.p1 && l < maxD2; l++) {
//             //     for(int r = 0; r <= t.p2 && r < maxD2; r++) {
//             //         if(r + t.h2 >= (int)ansAll[l + t.h1].size()) {
//             //             ansAll[l + t.h1].resize(r + t.h2 + 5);
//             //         }
//             //         ansAll[l + t.h1][r + t.h2] 
//             //             += C[t.p1][l] * C[t.p2][r] * C[pR][i];
//             //     }
//             // }
            
// #ifdef DEBUG
// if((int)ansAll[x].size() > y)
// printf("ansAll[%d][%d] %.2f\n", x, y, ansAll[x][y]);
// else printf("ansAll[%d][%d] %.2f\n", x, y, 0.0);
// #endif
//             t.h2 -= i;
//         }
        if(pR == 0) return;
        for(int i = pR - 1; i >= 0; i--) PR.push_back(candR[i]);
        for(int i = 0; i < pR; i++) {
            HR.push_back(candR[i]);
            PR.pop_back();
// printf("%u %u\n", PR.size(), t.p2 + pR - i - 1);fflush(stdout);
            pivotCount(l + 1, 0, 0, {t.p1, t.h1, t.p2 + pR - i - 1, t.h2 + 1});
            HR.pop_back();
        }

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

PL.push_back(pivotL);
PR.push_back(pivotR);
    pivotCount(l + 1, pLL, pRR, 
        treePath{t.p1 + 1, t.h1, t.p2 + 1, t.h2});
PL.pop_back();
PR.pop_back();


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
    for(int i = 0; i < pLL + 1; i++) {
        uint32_t u = candL[i];
        PL.push_back(u);
    }
    for(int i = tmpLSize - 1; i >= 0; i--) {
        uint32_t u = tmpNodesL[l][i];
        PL.push_back(u);
    }
    for(int i = 0; i < tmpLSize; i++) {
        uint32_t u = tmpNodesL[l][i];
        PL.pop_back();
        HL.push_back(u);
        pivotCount(l + 1, 0, 0, treePath{t.p1 + pL - i - 1, t.h1 + 1, t.p2, t.h2});
        HL.pop_back();
    }
    for(int i = 0; i < pLL + 1; i++) PL.pop_back();

    for(int i = 0; i < pRR + 1; i++) {
        uint32_t u = candR[i];
        PR.push_back(u);
    }
    for(int i = tmpRSize - 1; i >= 0; i--) {
        uint32_t u = tmpNodesR[l][i];
        PR.push_back(u);
    }
    for(int i = 0; i < tmpRSize; i++) {
        uint32_t v = tmpNodesR[l][i];
        PR.pop_back();
        HR.push_back(v);
        pivotCount(l + 1, 0, 0, treePath{t.p1, t.h1, t.p2 + pR - i - 1, t.h2 + 1});
        HR.pop_back();
    }
    for(int i = 0; i < pRR + 1; i++) PR.pop_back();



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
        // t.p1 += pL;
        // t.h1 += 1;
        // for(int l = 0; l <= t.p1 && l < maxD2; l++) {
        //     for(int r = 0; r <= t.p2 && r < maxD2; r++) {
        //         if(r + t.h2 >= (int)ansAll[l + t.h1].size()) {
        //             ansAll[l + t.h1].resize(r + t.h2 + 5);
        //         }
        //         ansAll[l + t.h1][r + t.h2] 
        //             += C[t.p1][l] * C[t.p2][r];
        //     }
        // }
        // t.p1 -= pL;
        // t.h1 -= 1;


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
HL.push_back(u);
HR.push_back(v);
                pivotCount(l + 1, pLL, pRR, 
                    treePath{t.p1, t.h1 + 1, t.p2, t.h2 + 1});
HL.pop_back();
HR.pop_back();
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

//         t.p2 += pR;//pivotR + [i+1,tmpRSize-1]
//         t.h2 += 1;
//         // for(int j = 1; j <= tmpRSize - i; j++) {
//             // t.h2 += j;
// #ifdef DEBUG
// printf("xadding %d %d %d %d, %d/%d\n", t.p1, t.h1, t.p2, t.h2, i, pR);
// #endif
//         for(int l = 0; l <= t.p1 && l < maxD2; l++) {
//             for(int r = 0; r <= t.p2 && r < maxD2; r++) {
//                 if(r + t.h2 >= (int)ansAll[l + t.h1].size()) {
//                     ansAll[l + t.h1].resize(r + t.h2 + 5);
//                 }
//                 ansAll[l + t.h1][r + t.h2] 
//                     += C[t.p1][l] * C[t.p2][r];
//             }
//         }
//             // t.h2 -= j;
//         // }
//         t.p2 -= pR;
//         t.h2 -= 1;

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
HL.push_back(u);
HR.push_back(v);
                pivotCount(l + 1, pLL, pRR, 
                    treePath{t.p1, t.h1 + 1, t.p2, t.h2 + 1});
HL.pop_back();
HR.pop_back();
            }
        }
    }
}
