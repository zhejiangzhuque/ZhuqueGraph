#include "BK.h"

void BK::collect2HopNeighborsContained(int p, int q) {
    H.contained.resize(g->n1);
    H.lists.resize(g->n1);
    H.nodes.resize(g->n1);

    std::vector<uint32_t> sum(g->n1);
    std::vector<uint32_t> tmp(g->n1);
    uint32_t l = 0;

    std::vector<uint32_t> sum2(g->n1);
    std::vector<uint32_t> tmp2(g->n1);
    uint32_t l2 = 0;

// int twotwo = 0;
    for(uint32_t u = 0; u < g->n1; u++) {
        H.contained[u] = false;

        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) {
                    if(sum[w] == 0) tmp[l++] = w;
                    sum[w]++;
                }
                else if(w < u) {
                    if(sum2[w] == 0) tmp2[l2++] = w;
                    sum2[w]++;
                }
            }
        }

        for(uint32_t i = 0; i < l; i++) {
            uint32_t w = tmp[i];
            if(sum[w] >= (uint32_t)q) {
// twotwo += C[sum[w]][q];
                H.lists[u].push_back(w);
            }
            sum[w] = 0;
        }
        l = 0;

        for(uint32_t i = 0; i < l2; i++) {
            uint32_t w = tmp2[i];
            if(sum2[w] == g->pU[u + 1] - g->pU[u]) {
                H.contained[u] = true;
// printf("%u contained by %u\n", u, w);
            }
            sum2[w] = 0;
        }
        l2 = 0;
    }

    for(uint32_t u = 0; u < g->n1; u++) {
        std::sort(H.lists[u].begin(), H.lists[u].end());
    }
}

void BK::exactCountMaximal(int p, int q) {
    // g->coreReduction(p, q);
    // if(costEstimate(p, q)) {
    //     g->swapUV();
    //     std::swap(p, q);
    //     printf("swap\n");
    // }
#ifdef DEBUG
    g->print();
#endif

    collect2HopNeighborsContained(p, q);
    printf("collect 2/contain\n"); fflush(stdout);
#ifdef DEBUG
    H.print();
#endif
    S.resize(g->n2);
    Q.resize(g->n1);

    uint32_t maxDegree = 0;
    for(uint32_t u = 0; u < g->n1; u++) {
        maxDegree = std::max(maxDegree, (uint32_t)H.lists[u].size());
    }
    printf("maxDegree of H %u\n", maxDegree);

    tmpNodes.resize(maxDegree + 1);
    tmpNodes[0].resize(g->n1);
    for(uint32_t i = 1; i < maxDegree + 1; i++) {
        tmpNodes[i].resize(maxDegree);
    }

    ans = 0;

    for(uint32_t u = 0; u + p - 1 < g->n1; u++) {
        if(H.contained[u]) continue;

        int pH = 0;
        for(auto v: H.lists[u]) H.nodes.changeTo(v, pH++);

        int pS = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            S.changeTo(g->e1[i], pS++);
        }
#ifdef DEBUG
printf("choose %u, cands:", u);
for(int i = 0; i < pH; i++) printf("%d ", H.nodes[i]);
printf("\n");
#endif

        countFormMaximalClique(1, 1, pH, pS, 0, p, q);
#ifdef DEBUG
printf("      ans %.2f\n", ans);
#endif
    }

    printf("ans %.2f\n", ans);
    fflush(stdout);
}

void BK::countFormMaximalClique(int l, int L, int pH, int pS, int pQ, int p, int q) {
#ifdef DEBUG
printf("l:%d L:%d pH:%d pS:%d pQ:%d\n", l, L, pH, pS, pQ);
#endif
    H.nodes.copy(tmpNodes[l].data(), pH);
    int ppH = pH;//deleted from H

    for(int j = 0; j < pH; j++) {
        uint32_t u = tmpNodes[l][j];

#ifdef DEBUG
printf("l%u,testing %u, idx %u, ppH %d, ", l,u, H.nodes.idx(u), ppH);
#endif

        if(H.nodes.idx(u) >= (uint32_t)ppH || H.lists[u].size() + 1 < uint32_t(p - L)) {
#ifdef DEBUG
printf(" test end\n");
#endif
            continue;
        }
        H.nodes.changeTo(u, --ppH);

        int pSS = 0;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            if(S.idx(g->e1[i]) < (uint32_t)pS) {
                S.changeTo(g->e1[i], pSS++);
            }
        }

        if(pSS < q) {
#ifdef DEBUG
printf(" test end\n");
#endif
            continue;
        }
#ifdef DEBUG
printf("pSS:%u ", pSS);
#endif
        int pQQ = 0;
        bool isFirst = true;
        for(int i = 0; i < pQ; i++) {
            int nQ = 0;

            for(int j = 0; j < pSS; j++) {
                if(g->connectUV(Q[i], S[j])) {
                    nQ++;   
                }
            }

            if(nQ == pSS) {
                isFirst = false;
                break;
            }
            else if(nQ >= q) {
                Q.changeToByPos(i, pQQ++);
            }
        }

        if(!isFirst) {
#ifdef DEBUG
printf("is not first\n");
#endif
            continue;
        }
#ifdef DEBUG
printf("choose %u, cands:\n", u);
#endif

        int pHH = 0;
        int addL = 0;
        for(auto v: H.lists[u]) {
#ifdef DEBUG
printf("consider %u,idx %u,ppH %d,", v, H.nodes.idx(v), ppH);
#endif
            if(H.nodes.idx(v) >= (uint32_t)ppH) {
#ifdef DEBUG
printf("not choose\n");
#endif
                continue;
            }

            int nV = 0;
            for(int j = 0; j < pSS; j++) {
                if(g->connectUV(v, S[j])) {
                    nV++;        
                }
            }
#ifdef DEBUG
printf("nV:%u ", nV);
#endif
            if(nV == pSS) {
                addL++;
#ifdef DEBUG
printf("addL\n");
#endif
                nV = 0;
                for(int j = pSS; j < pS; j++) {
                    if(g->connectUV(v, S[j])) {
                        nV++;
                    }
                }
                if(nV == 0) {
                    //del v from H
                    H.nodes.changeTo(v, --ppH);
                    // Q.changeTo(v, pQ++);
                }
            }
            else if(nV >= q) {
#ifdef DEBUG
printf("pHH\n");
#endif
                H.nodes.changeTo(v, pHH++);
            }
        }

        if(L + 1 + addL >= p && pSS >= q) {
#ifdef DEBUG
printf("adding:(%.2f - %.2f)*%.2f\n", C[L + 1 + addL][p], C[L][p], C[pSS][q]);
#endif
            // ans += (C[L + 1 + addL][p] - C[L][p]) * C[pSS][q];
            ans++;
        }

// printf("\n");
        
        countFormMaximalClique(l + 1, L + 1 + addL, pHH, pSS, pQQ, p, q);
        Q.changeTo(u, pQ++);
    }
#ifdef DEBUG
printf("return l:%d\n", l);
#endif
}


