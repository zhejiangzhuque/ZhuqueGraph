#include "biplex.h"
#define DEBUG

void biplex::run() {
g->print();printf("\n");

    plexEnum(0, {0, 0}, {g->n[0], g->n[1]}, {0, 0});

    printf("total cnt : %llu\n", totalMaximalBiPlex);
}

void biplex::plexEnum(uint32_t deep, parameter r, parameter c, parameter x) {
#ifdef DEBUG
printf("\n");
for(unsigned i = 0; i < deep; i++) printf(" ");
printf("deep %u\n", deep);fflush(stdout);
printf("R_l:");
for(uint32_t i = 0; i < r[0]; i++) printf("%u ", R[0][i]);printf("\n");
printf("R_r:");
for(uint32_t i = 0; i < r[1]; i++) printf("%u ", R[1][i]);printf("\n");
printf("C_l:");
for(uint32_t i = 0; i < c[0]; i++) printf("%u ", C[0][i]);printf("\n");
printf("C_r:");
for(uint32_t i = 0; i < c[1]; i++) printf("%u ", C[1][i]);printf("\n");
printf("X_l:");
for(uint32_t i = 0; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("X_r:");
for(uint32_t i = 0; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
#endif
    
    for(int t = 0; t <= 1; t++) {
        uint32_t z = t ^ 1;
        if(c[t] == 0) {
            if(x[t] == 0) {
                parameter newC = c, newR = r;
                for(uint32_t i = 0; i < newC[z]; ) {
                    uint32_t u = C[z][i];

                    if(nonNeiC.getCntNonNeiIdx(z, deep, i) == 0) {
                        R[z].changeTo(u, newR[z]++);
                        C[z].changeTo(u, --newC[z]);
                    }
                    else i++;
                }
                impPlexBranch(deep, newR, newC, x, t);
            }
            return;
        }
    }
#ifdef DEBUG
// printf("deep2 %u\n", deep);fflush(stdout);
#endif
    //find the pivot with maximal degree in cand
    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {C[0][0], C[1][0]};

    for(uint32_t t = 0; t < 2; t++) {
        uint32_t z = t ^ 1;//t + z = 1

        for(uint32_t i = 0; i < c[t]; i++) {//scan all vertices in C
            uint32_t u = C[t][i];
            uint32_t e = 0;

            if(g->deg(u, t) > c[z]) {
                for(uint32_t j = 0; j < c[z]; j++) {
                    uint32_t v = C[z][j];
                    if(g->connect(u, v, t)) {
                        e++;
                    }
                }
            }
            else {
                for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    if(C[z].pos(v) < c[z]) {//in C
                        e++;
                    }
                }
            }

            if(e > maxE[z]) {
                maxE[z] = e;
                maxI[t] = i;
                maxUV[t] = u;
            }
        }
    }
#ifdef DEBUG
// printf("deep3 %u\n", deep);fflush(stdout);
#endif
    //t is the side of pivot, and S[z].cSize() - maxE[z] is smaller than S[t].cSize() - maxE[t]
    uint32_t t = 1;
    if(c[0] - maxE[0] >= c[1] - maxE[1]) {
        t = 0;
    }
    uint32_t z = t ^ 1;

    if(maxE[0] == 0) {
#ifdef DEBUG
printf("no edge\n");
#endif
        plexEnum(deep + 1, r, {0, c[1]}, x);
#ifdef DEBUG
printf("    return1 deep %u\n", deep);
#endif
        plexEnum(deep + 1, r, {c[0], 0}, x);
#ifdef DEBUG
printf("    return2 deep %u\n", deep);
#endif
        return;
    }

    //ws[deep] is the nonNeighbors of pivot in cand
#ifdef DEBUG
printf("pivot %u, t %u, degree %u\n", maxUV[t], t, maxE[z]);fflush(stdout);
#endif
    //run pivot
    R[t].changeTo(maxUV[t], r[t]++);
    C[t].changeTo(maxUV[t], --c[t]);

    nonNeiC.resizeNextDeep(deep, c);
    nonNeiR.resizeNextDeep(deep, r);
    nonNeiX.resizeNextDeep(deep, x);
    //scan c x r
    // for(uint32_t i = 0; i < c[z]; i++) if(!g->connect(maxUV[t], C[z][i], t)) nonNei[z].addNonNei(C[z][i], maxUV[t]);
    // for(uint32_t i = 0; i < x[z]; i++) if(!g->connect(maxUV[t], X[z][i], t)) nonNei[z].addNonNei(X[z][i], maxUV[t]);
    for(uint32_t i = 0; i < r[z]; i++) {
        nonNeiR.cntNonNei[0][deep + 1][i] = nonNeiR.cntNonNei[0][deep][i];
        nonNeiR.cntNonNei[1][deep + 1][i] = nonNeiR.cntNonNei[1][deep][i];
        if(!g->connect(maxUV[t], R[z][i], t))
            nonNeiR.addNonNeiIdx(t, R[z][i], maxUV[t], deep + 1, i);
    }

    //make sure newC[0] points to U
    parameter newC = plexUpdate(t, deep, r, maxUV[t], c, C, nonNeiC);
    parameter newX = plexUpdate(t, deep, r, maxUV[t], x, X, nonNeiX);
#ifdef DEBUG
// printf("deep4.1 %u\n", deep);fflush(stdout);
#endif
    uint32_t wsSize = c[z] - newC[z];
#ifdef DEBUG
// printf("deep4.1 %u, wsSize %u\n", deep, wsSize);fflush(stdout);
#endif
    if(ws[deep].size() < wsSize) {
#ifdef DEBUG
// printf("deep4.1 %u, wsSize %u\n", deep, wsSize);fflush(stdout);
#endif
        ws[deep].resize(wsSize * 2);
    }
#ifdef DEBUG
// printf("deep4.2 %u\n", deep);fflush(stdout);
#endif
    memcpy(ws[deep].data(), C[z].begin() + newC[z], sizeof(uint32_t) * wsSize);
#ifdef DEBUG
// printf("deep5 %u\n", deep);fflush(stdout);
#endif
    
    plexEnum(deep + 1, r, newC, newX);
#ifdef DEBUG
printf("    return3 deep %u\n", deep);
#endif

    // for(uint32_t i = 0; i < c[z]; i++) if(!g->connect(maxUV[t], C[z][i], t)) nonNei[z].pop(C[z][i]);
    // for(uint32_t i = 0; i < x[z]; i++) if(!g->connect(maxUV[t], X[z][i], t)) nonNei[z].pop(X[z][i]);
    // for(uint32_t i = 0; i < r[z]; i++) if(!g->connect(maxUV[t], R[z][i], t)) nonNei[z].pop(R[z][i]);

    //rm pivot
    X[t].changeTo(maxUV[t], x[t]++);
    // R[t].changeTo(maxUV[t], --r[t]);
    --r[t];

    for(uint32_t i = 0; i < wsSize; i++) {
        uint32_t w = ws[deep][i];

        R[z].changeTo(w, r[z]++);
        C[z].changeTo(w, --c[z]);

        parameter newC = plexUpdate(z, deep, r, w, c, C, nonNeiC);
        parameter newX = plexUpdate(z, deep, r, w, x, X, nonNeiX);

        //scan c x r
        // for(uint32_t i = 0; i < newC[t]; i++) if(!g->connect(w, C[t][i], z)) nonNei[t].addNonNei(C[t][i], w);
        // for(uint32_t i = 0; i < newX[t]; i++) if(!g->connect(w, X[t][i], z)) nonNei[t].addNonNei(X[t][i], w);
        for(uint32_t i = 0; i < r[t]; i++)
            if(!g->connect(w, R[t][i], z))
                nonNeiR.addNonNeiIdx(t, R[t][i], w, deep, i);

        plexEnum(deep + 1, r, newC, newX);
#ifdef DEBUG
printf("    return deep4 %u\n", deep);
#endif
        // for(uint32_t i = 0; i < newC[t]; i++) if(!g->connect(w, C[t][i], z)) nonNei[t].pop(C[t][i]);
        // for(uint32_t i = 0; i < newX[t]; i++) if(!g->connect(w, X[t][i], z)) nonNei[t].pop(X[t][i]);
        // for(uint32_t i = 0; i < r[t]; i++) if(!g->connect(w, R[t][i], z)) nonNei[t].pop(R[t][i]);

        X[z].changeTo(w, x[z]++);
        // R[z].changeTo(w, --r[z]);
        --r[z];
    }
}

void biplex::impPlexBranch(uint32_t deep, parameter r, parameter c, parameter x, uint32_t t) {
#ifdef DEBUG
for(unsigned i = 0; i < deep; i++) printf(" ");
printf("imp_deep %u\n", deep);fflush(stdout);
printf("R_l:");
for(uint32_t i = 0; i < r[0]; i++) printf("%u ", R[0][i]);printf("\n");
printf("R_r:");
for(uint32_t i = 0; i < r[1]; i++) printf("%u ", R[1][i]);printf("\n");
printf("C_l:");
for(uint32_t i = 0; i < c[0]; i++) printf("%u ", C[0][i]);printf("\n");
printf("C_r:");
for(uint32_t i = 0; i < c[1]; i++) printf("%u ", C[1][i]);printf("\n");
printf("X_l:");
for(uint32_t i = 0; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("X_r:");
for(uint32_t i = 0; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
fflush(stdout);
#endif
    
    uint32_t z = t ^ 1;
    //move C[z] into R[z], let X[t]=0

    if(c[t] == 0) {
        if(x[t] + x[z] == 0 && r[t] > 0 && r[z] > 0) {
            totalMaximalBiPlex++;
#ifdef DEBUG
printf("get a answer\n");
#endif
        }

        return ;
    }
    
    uint32_t pivot = C[z][0], pivotIdx = 0;
    uint32_t nd = nonNeiC.getCntNonNeiIdx(z, deep, 0);
    for(uint32_t i = 1; i < c[z]; i++) {
        uint32_t v = C[z][i];
        if(nd > nonNeiC.getCntNonNeiIdx(z, deep, i)) {
            nd = nonNeiC.getCntNonNeiIdx(z, deep, i);
            pivot = v;
            pivotIdx = i;
        }
    }
    
    //build Cv/Pv
    uint32_t wsSize = 0;

    uint32_t costkd = 0, cntNonNeiPivot = nonNeiC.getCntNonNeiIdx(z, deep, pivotIdx);
    for(uint32_t i = 0; i < cntNonNeiPivot; i++) {
        uint32_t v = nonNeiC.getBuffer(z, pivot)[i];
        costkd += g->deg(v, t);
    }
    if(costkd < c[z] * cntNonNeiPivot) {
        wsSize = 0;
        for(uint32_t i = 0; i < cntNonNeiPivot; i++) {
            uint32_t v = nonNeiC.getBuffer(z, pivot)[i];
            uint32_t tmp = c[z];
            
            for(uint32_t j = g->p[t][v]; j < g->p[t][v + 1]; j++) {
                uint32_t w = g->e[t][j];

                if(wsSize <= C[z].pos(w) && C[z].pos(w) < c[z]) {
                    C[z].changeTo(w, --tmp);
                }
            }

            wsSize = tmp;
        }
    }
    else {
        wsSize = c[z];
        for(uint32_t j = 0; j < c[z]; j++) {
            uint32_t u = C[z][j];
            bool connectAll = true;

            for(uint32_t i = 0; i < cntNonNeiPivot; i++) {
                uint32_t v = nonNeiC.getBuffer(z, pivot)[i];
                if(!g->connect(u, v, z)) {
                    connectAll = false;
                    break;
                }
            }

            if(connectAll) {
                C[z].changeToByPos(j, --wsSize);
            }
        }
    }

    if(ws[deep].size() < wsSize) {
        ws[deep].resize(wsSize * 2);
    }
    memcpy(ws[deep].data(), C[z].begin(), sizeof(uint32_t) * wsSize);

    for(uint32_t i = 0; i < wsSize; i++) {
        uint32_t v = ws[deep][i];

        R[z].changeTo(v, r[z]++);
        C[z].changeTo(v, --c[z]);

        parameter newC = plexUpdate(z, deep, r, v, c, C, nonNeiC);
        parameter newX = plexUpdate(z, deep, r, v, x, X, nonNeiX);

        for(uint32_t i = 0; i < r[t]; i++) 
            if(!g->connect(v, R[t][i], z)) 
                nonNeiR.addNonNeiIdx(t, R[t][i], v, deep, i);

        // plexEnum(deep + 1, r, newC, newX);
        impPlexBranch(deep + 1, r, newC, newX, t);

        X[z].changeTo(v, x[z]++);
        // R[z].changeTo(v, --r[z]);
        --r[z];
    }
}

biplex::parameter biplex::plexUpdate(uint32_t t, uint32_t deep, parameter r, uint32_t u, 
    parameter c, LinearSet * C, nonNeiMatainer & nonNei) 
{    
    uint32_t z = t ^ 1;
    parameter newC;
    newC[t] = c[t];
    newC[z] = 0;

    for(uint32_t i = 0; i < c[z]; i++) {
        uint32_t v = C[z][i];
        if(g->connect(u, v, t)) {
            C[z].changeTo(v, newC[z]++);
        }
        // else {
        //     if(nonNei.getCntNonNeiIdx(z, deep, i) + 1 < k) {
        //         nonNei.addNonNeiIdx(z, v, u, deep + 1, i);
        //         C[z].changeTo(v, newC[z]++);
        //     }
        // }
    }

    for(uint32_t i = 0; i < nonNei.getCntNonNei(t, deep, u, C); i++) {
        uint32_t v = nonNei.getBuffer(t, u)[i];
        
        if(nonNei.getCntNonNei(z, deep, v, C) == k) {
            if(g->deg(v, z) < c[t]) {
                uint32_t tmp = 0;
                for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    if(C[t].pos(w) < newC[t]) {
                        C[t].changeTo(w, tmp++);
                    }
                }
                newC[t] = tmp;
            }
            else {
                for(uint32_t j = 0; j < c[t]; j++) {
                    uint32_t w = C[t][j];
                    
                    if(!g->connect(v, w, z)) {
                        C[t].changeTo(w, --newC[t]);
                    }
                }
            }
        }
    }

    return newC;
}
