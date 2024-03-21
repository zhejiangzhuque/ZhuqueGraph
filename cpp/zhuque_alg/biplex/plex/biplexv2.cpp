#include "biplexv2.h"
// #define DEBUG
// #define DEBUG2
// #define BASELINE
#define LIMITEDNUM
// #define ANSPRINT
#define DEEPTHRESHOLD 40

#ifdef BASELINE
void print(uint32_t i) {
    if(i == 0) return;
    print(i>>1);
    printf("%d", i&1);
}
#endif

void biplexv2::run() {
#ifdef DEBUG2
g->print();printf("\n");

fflush(stdout);
#endif
#ifdef DEBUG
g->print();printf("\n");

fflush(stdout);
#endif
#ifdef DEEPTHRESHOLD
// printf("DEEPTHRESHOLD %d\n", DEEPTHRESHOLD);
#endif

    uint32_t t = 0, z = 1;
    if(od != "two") {
        if(g->maxDu < g->maxDv) {
            t = 1; z = 0;
        }
    }
    else if(g->core[0] < g->core[1]) {t=1; z=0;}
    printf("t=%d,z=%d\n", t, z);
    
    if(std::min(ls, rs) < 2*k + 1) {
    // if (false) {
        // printf("all\n");
        plexEnum(0, {0, 0}, {g->n[0], g->n[1]}, {0, 0}, {0, 0});
    }
    else {
        parameter r{0, 0}, c{0, 0}, x{0, 0}, sx{0, 0};
        r[t] = 1;
        std::vector<uint32_t> deg(g->n[t]), degZ(g->n[z]);
        std::vector<uint32_t> stk(g->n[t]), stkZ(g->n[z]);
        uint32_t l = 0, lz = 0;
        uint32_t testlen = 0;
        for(uint32_t u = 0; u < g->n[t]; u++) {
            c[t] = 0;
            x[t] = 0;
            c[z] = 0;
            if (g->cores[t][u]+k < lrs[z]) continue;
            R[t].changeTo(u, 0);
            l = lz = 0;
            testlen++;

            for(uint32_t i = g->p[t][u]; i < g->p[t][u + 1]; i++) {
                uint32_t v = g->e[t][i];
                if (g->cores[z][v]+k < lrs[t]) continue;
                C[z].changeTo(v, c[z]++);

                if(g->p[z][v + 1] > 0)
                for(uint32_t j = g->p[z][v + 1] - 1; j >= g->p[z][v]; j--) {
                    uint32_t w = g->e[z][j];

                    if(w > u) {
                        deg[w]++;
                        
                        if(deg[w] + 2*k == lrs[z] && g->cores[t][w]+k >= lrs[z]) {
                            // stk[l++] = w;
                            C[t].changeTo(w, c[t]++);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }



            for(uint32_t i = g->p[t][u]; i < g->p[t][u + 1]; i++) {
                uint32_t v = g->e[t][i];
                if (g->cores[z][v]+k < lrs[t]) continue;
                if(g->p[z][v + 1] > 0)
                for(uint32_t j = g->p[z][v + 1] - 1; j >= g->p[z][v]; j--) {
                    uint32_t w = g->e[z][j];

                    if(w > u) {
                        if(C[t].pos(w) < c[t]) degZ[v]++;
                    }
                    else break;

                    if(j == 0) break;
                }
            }

            for(uint32_t i = g->p[t][u]; i < g->p[t][u + 1]; i++) {
                uint32_t v = g->e[t][i];
                if(degZ[v] + k + 1 < lrs[t] &&  (g->cores[z][v]+k >= lrs[t]) ) {
                    stkZ[lz++] = v;
                }
            }

            while(lz > 0 || l > 0) {
                //remove z
                for(uint32_t ii = 0; ii < lz; ii++) {
                    uint32_t v = stkZ[ii];
                    C[z].changeTo(v, --c[z]);
                    degZ[v] = 0;

                    if(g->deg(v, z) < c[t]) {
                        if(g->p[z][v + 1] > 0)
                        for(uint32_t j = g->p[z][v + 1] - 1; j >= g->p[z][v]; j--) {
                            uint32_t w = g->e[z][j];
                            if(deg[w] > 0) {
                                --deg[w];
                                if(deg[w] + 2*k + 1 == lrs[z] && g->cores[t][w]+k >= lrs[z]) {
                                    stk[l++] = w;
                                }
                            }

                            if(j == 0) break;
                        }
                    }
                    else {
                        for(uint32_t i = 0; i < c[t]; i++) {
                            uint32_t w = C[t][i];

                            if(g->connect(v, w, z)) {
                                if(deg[w] > 0) {
                                    --deg[w];
                                    if(deg[w] + 2*k + 1 == lrs[z] ) {
                                        stk[l++] = w;
                                    }
                                }
                            }
                        }
                    }
                }

                //remove t
                lz = 0;
                for(uint32_t ii = 0; ii < l; ii++) {
                    uint32_t u = stk[ii];
                    C[t].changeTo(u, --c[t]);
                    deg[u] = 0;
                    
                    for(uint32_t i = 0; i < c[z]; i++) {
                        uint32_t v = C[z][i];

                        if(g->connect(u, v, t)) {
                            if(degZ[v] > 0) {
                                --degZ[v];
                                if(degZ[v] + k + 2 == lrs[t]) {
                                    stkZ[lz++] = v;
                                }
                            }
                        }
                    }
                }

                l = 0;
            }

#ifdef DEBUG2
printf("C[t]: ");
for(uint32_t ii = 0; ii < c[t]; ii++) {
    uint32_t u = C[t][ii];
    printf("%u ", u);
}

printf("C[z]: ");
for(uint32_t ii = 0; ii < c[z]; ii++) {
    uint32_t u = C[z][ii];
    printf("%u-%u ", u, degZ[u]);
}

printf("\n");
fflush(stdout);
#endif

            for(uint32_t ii = 0; ii < c[t]; ii++) {
                uint32_t uu = C[t][ii];
                for(uint32_t i = g->p[t][uu]; i < g->p[t][uu + 1]; i++) {
                    uint32_t v = g->e[t][i];
                    if (g->cores[z][v]+k < lrs[t]) continue;
                    degZ[v]++;
                    if(C[z].pos(v) >= c[z] && degZ[v] + k == lrs[t]) {
                        C[z].changeTo(v, c[z]++);
                        
                        nonNei.addNonNei(z, v, u);
                    }
                }
            }
#ifdef DEBUG2
printf("C[t]: ");
for(uint32_t ii = 0; ii < c[t]; ii++) {
    uint32_t u = C[t][ii];
    printf("%u ", u);
}

printf("C[z]: ");
for(uint32_t ii = 0; ii < c[z]; ii++) {
    uint32_t u = C[z][ii];
    printf("%u-%u ", u, degZ[u]);
}
printf("%u-%u ", 2, degZ[2]);

printf("\n");
fflush(stdout);
#endif

            for(uint32_t ii = 0; ii < c[t]; ii++) {
                uint32_t u = C[t][ii];
                for(uint32_t i = g->p[t][u]; i < g->p[t][u + 1]; i++) {
                    uint32_t v = g->e[t][i];
                    // if(degZ[v] + k <= lrs[t]) 
                        degZ[v] = 0;
                }
            }

            for(uint32_t ii = 0; ii < c[t]; ii++) deg[C[t][ii]] = 0;

            for(uint32_t ii = 0; ii < c[z]; ii++) {
                uint32_t v = C[z][ii];

                for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];

                    if(w < u) {
                        deg[w]++;
                        
                        if(deg[w] + 2*k == lrs[z]) {
                            X[t].changeTo(w, x[t]++);
                        }
                    }
                    else break;
                }
            }

            for(uint32_t ii = 0; ii < c[z]; ii++) {
                uint32_t v = C[z][ii];

                for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];

                    if(w < u) deg[w] = 0;
                    else break;
                }
            }

            plexEnum(0, r, c, x, sx);
#ifdef LIMITEDNUM
            if(totalMaximalbiplexv2 >= outPutT) break;
#endif
            for(uint32_t ii = 0; ii < c[z]; ii++) nonNei.clear(z, C[z][ii]);
            for(uint32_t ii = 0; ii < c[z]; ii++) degZ[C[z][ii]] = 0;
        }
        printf("testlen : %u\n", testlen);
    }

    printf("total cnt : %llu\n", totalMaximalbiplexv2);

//print the right answer
#ifdef BASELINE
auto check = [&](uint32_t i, uint32_t j)->bool {
    uint32_t sz = 0;
    for(uint32_t u = 0; u < g->n[0]; u++) {
        if((i & (1<<u)) == 0) continue;
        sz++;
    }
    if(sz < ls) return false;
    sz = 0;
    for(uint32_t v = 0; v < g->n[1]; v++) {
        if((j & (1<<v)) == 0) continue;
        sz++;
    }
    if(sz < rs) return false;

    for(uint32_t u = 0; u < g->n[0]; u++) {
        if((i & (1<<u)) == 0) continue;
        int nonNei = 0;
        for(uint32_t v = 0; v < g->n[1]; v++) {
            if((j & (1<<v)) == 0) continue;
            if(!g->connect(u, v, 0)) nonNei++; 
        }

        if(nonNei > k) {
            return false;
        }
    }

    for(uint32_t v = 0; v < g->n[1]; v++) {
        if((j & (1<<v)) == 0) continue;
        int nonNei = 0;
        for(uint32_t u = 0; u < g->n[0]; u++) {
            if((i & (1<<u)) == 0) continue;
            if(!g->connect(u, v, 0)) nonNei++; 
        }

        if(nonNei > k) {
            return false;
        }
    }
    return true;
};

int realCnt = 0;
for(uint32_t i = (1<<g->n[0])-1; i > 0; i--) {
    for(uint32_t j = (1<<g->n[1])-1; j > 0; j--) {
        if(!check(i, j)) continue;

        bool isMaximal = true;
        for(uint32_t u = 0; u < g->n[0]; u++) {
            if((i & (1<<u)) > 0) continue;
            if(check(i|(1<<u), j)) {
                isMaximal = false;

                break;
            }
        }

        if(isMaximal)
        for(uint32_t v = 0; v < g->n[1]; v++) {
            if((j & (1<<v)) > 0) continue;
            if(check(i, j|(1<<v))) {
                isMaximal = false;

                break;
            }
        }

        if(isMaximal) {
            realCnt++;

// print(i);printf("\n");
// print(j);printf("\n");
// printf("L:");
// for(uint32_t u = 0; u < g->n[0]; u++) if((i & (1<<u)) > 0) printf("%u ", u);printf("\n");
// printf("R:");
// for(uint32_t v = 0; v < g->n[1]; v++) if((j & (1<<v)) > 0) printf("%u ", v);printf("\n");

        }
    }
}
printf("realCnt:%d\n", realCnt);
#endif
}

void biplexv2::plexEnum(uint32_t deep, parameter r, parameter c, parameter x, parameter sx) {

// printf("deep:%u cnt:%u\n", deep, totalMaximalbiplexv2);fflush(stdout);

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
printf("NoNeiC_l:");
for(uint32_t i = 0; i < c[0]; i++) printf("%u ", nonNei.getCntNonNei(0, C[0][i]));printf("\n");
printf("NoNeiC_r:");
for(uint32_t i = 0; i < c[1]; i++) printf("%u ", nonNei.getCntNonNei(1, C[1][i]));printf("\n");
printf("X_l:");
for(uint32_t i = sx[0]; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("X_r:");
for(uint32_t i = sx[1]; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
printf("aX_l:");
for(uint32_t i = 0; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("aX_r:");
for(uint32_t i = 0; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
fflush(stdout);
#endif
// if(deep == 0)
// printf("c[0]: %u %u %u %u\n", c[0], c[1], r[0], r[1]);
    if(c[0] + r[0] < ls) return;
    if(c[1] + r[1] < rs) return;
    
    // if(c[0] + r[0] == 0 || c[1] + r[1] == 0) return;
    if(c[0] + c[1] == 0) {
        if(x[0] - sx[0] + x[1] - sx[1] == 0) {
            totalMaximalbiplexv2++;
        if (totalMaximalbiplexv2 == outPutT) {printf("done\n"); fflush(stdout);} 

#ifdef ANSPRINT
printf("R_l:");
for(uint32_t i = 0; i < r[0]; i++) printf("%u ", R[0][i]);printf("\n");
printf("R_r:");
for(uint32_t i = 0; i < r[1]; i++) printf("%u ", R[1][i]);printf("\n");

fflush(stdout);
#endif

#ifdef DEBUG
printf("answer\n");
#endif
        }
        return;
    }

    //find the pivot with maximal degree in cand
    uint32_t maxE[2] = {g->n[1] + 1, g->n[0] + 1}, maxUV[2] = {g->n[0], g->n[1]};

#ifdef DEEPTHRESHOLD
    bool hasCommonNeiInC[2] = {false, false};
    bool hasCommonNeiInX[2] = {false, false};
    if(deep <= DEEPTHRESHOLD)
    for(uint32_t t = 0; t < 2; t++) {
        uint32_t z = t ^ 1;//t + z = 1
        for(uint32_t i = 0; i < c[t]; i++) {
            if(nonNei.getCntNonNei(t, C[t][i]) == 0) {
                hasCommonNeiInC[t] = true;
                break;
            }
        }
        for(uint32_t i = sx[t]; i < x[t]; i++) {
            if(nonNei.getCntNonNei(t, X[t][i]) == 0) {
                hasCommonNeiInX[t] = true;
                break;
            }
        }
    }
#endif

    auto findPivotC = [&](uint32_t t, uint32_t z, uint32_t u, uint32_t i) {
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

        if(maxE[t] == g->n[z] + 1 || e > maxE[t]) {
            maxE[t] = e;
            maxUV[t] = u;
        }
        return e;
    };

    auto findPivotX = [&](uint32_t t, uint32_t z, uint32_t u, uint32_t i) {
        uint32_t e = 0;

        if(g->deg(u, t) > x[z]) {
            for(uint32_t j = 0; j < x[z]; j++) {
                uint32_t v = X[z][j];
                if(g->connect(u, v, t)) {
                    e++;
                }
            }
        }
        else {
            for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                if(X[z].pos(v) < x[z]) {//in X
                    e++;
                }
            }
        }

        if(maxE[t] == g->n[z] + 1 || e > maxE[t]) {
            maxE[t] = e;
            maxUV[t] = u;
        }
        return e;
    };


    for(uint32_t t = 0; t < 2; t++) {
        uint32_t z = t ^ 1;//t + z = 1
#ifdef DEEPTHRESHOLD
        if(hasCommonNeiInC[t]) {
            for(uint32_t i = 0; i < c[t]; i++) if(nonNei.getCntNonNei(t, C[t][i]) == 0) {
                uint32_t ui = C[t][i];
                uint32_t e = findPivotC(t, z, ui, i);
                if (e + k + r[z] < lrs[z]){
                    C[t].changeTo(ui, --c[t]); i--;
                }
            }
        }
        else
#endif
        for(uint32_t i = 0; i < c[t]; i++) {
            uint32_t ui = C[t][i];
            uint32_t e = findPivotC(t, z, ui, i);
            if (e + k + r[z] < lrs[z]+ nonNei.getCntNonNei(t, ui)){
                C[t].changeTo(ui, --c[t]); i--;
            }
        }

#ifdef DEEPTHRESHOLD
        if(hasCommonNeiInX[t]) {
            for(uint32_t i = sx[t]; i < x[t]; i++) if(nonNei.getCntNonNei(t, X[t][i]) == 0) {
                findPivotX(t, z, X[t][i], i);
            }
        }
        else
#endif
        for(uint32_t i = sx[t]; i < x[t]; i++) {
            findPivotX(t, z, X[t][i], i);
        }
    }

    parameter newSC = {0, 0}, newC = c, newSX = sx, newX = x;

    auto computePuv = [&](uint32_t t, uint32_t z, uint32_t u) {
        newSC[t] = 0;
        if(nonNei.getCntNonNei(t, u) == 0) {
            if(C[t].pos(u) < c[t]) {//pivot in C,同一边不用枚举
                newSC[t] = 1;
                C[t].changeTo(u, 0);
            }
            // else {//pivot in X 都不用枚举

            // }
            return;
        }

        for(uint32_t i = 0; i < nonNei.getCntNonNei(t, u); i++) {
            uint32_t v = nonNei.getBuffer(t, u)[i];
            if(g->deg(v, z) < c[t]) {
                uint32_t tmp = c[t];
                for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    uint32_t pw = C[t].pos(w);
                    if(newSC[t] <= pw && pw < c[t]) {
                        C[t].changeTo(w, --tmp);
                    }
                }
                newSC[t] = tmp;
            }
            else {
                uint32_t tmp = c[t];
                for(uint32_t j = newSC[t]; j < tmp; ) {
                    uint32_t w = C[t][j];
                    
                    if(g->connect(v, w, z)) {
                        C[t].changeTo(w, --tmp);
                    }
                    else j++;
                }
                newSC[t] = tmp;
            }
        }
    };
    
    if(maxUV[0] != g->n[0]) computePuv(0, 1, maxUV[0]);
    if(maxUV[1] != g->n[1]) computePuv(1, 0, maxUV[1]);

    uint32_t Pu[2] = {c[0] - newSC[0], c[1] - newSC[1]};
    
    uint32_t t = 0;
    if(maxUV[0] == g->n[0] ) t = 1;
    else if(maxUV[1] == g->n[1]) t = 0;
    else if(Pu[0] + maxE[0] < Pu[1] + maxE[1]) t = 1;
    uint32_t z = t ^ 1;

    if(Puv[t][deep].size() < c[t] - Pu[t]) {
        Puv[t][deep].resize((c[t] - Pu[t])*2);
    }
    memcpy(Puv[t][deep].data(), C[t].begin(), sizeof(uint32_t) * newSC[t]);

    if(g->deg(maxUV[t], t) > c[z]) {
        uint32_t tmp = c[z];
        for(uint32_t j = 0; j < tmp; ) {
            uint32_t v = C[z][j];
            if(g->connect(maxUV[t], v, t)) {
                C[z].changeTo(v, --tmp);
            }
            else j++;
        }
        newC[z] = tmp;
    }
    else {
        uint32_t u = maxUV[t];
        uint32_t tmp = c[z];
        for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            if(C[z].pos(v) < c[z]) {//in C
                C[z].changeTo(v, --tmp);
            }
        }
        newC[z] = tmp;
    }
    if(Puv[z][deep].size() < newC[z]) {
        Puv[z][deep].resize(newC[z] * 2);
    }
    memcpy(Puv[z][deep].data(), C[z].begin(), sizeof(uint32_t) * newC[z]);

    uint32_t PuvSize[2];
    PuvSize[t] = newSC[t];
    PuvSize[z] = newC[z];

#ifdef DEBUG
printf("t %u, pivot %u\n", t, maxUV[t]);
for(uint32_t ii = 0; ii < 2; ii++) {
    uint32_t tt = ii == 0? t : z;
    uint32_t zz = tt ^ 1;

    for(uint32_t i = 0; i < PuvSize[tt]; i++) {
        uint32_t w = Puv[tt][deep][i];
        printf("%u ", w);
    }
    printf("\n");
    fflush(stdout);
}
#endif

    for(uint32_t ii = 0; ii < 2; ii++) {
        uint32_t tt = ii == 0? z : t;
        uint32_t zz = tt ^ 1;
        for(uint32_t i = 0; i < PuvSize[tt]; i++) {
            uint32_t w = Puv[tt][deep][i];

#ifdef DEBUG
printf("deep %u tt %u c[tt] %u, w %u\n", deep, tt, c[tt], w);
printf(" R_l:");
for(uint32_t i = 0; i < r[0]; i++) printf("%u ", R[0][i]);printf("\n");
printf(" R_r:");
for(uint32_t i = 0; i < r[1]; i++) printf("%u ", R[1][i]);printf("\n");
printf("C_l:");
for(uint32_t i = 0; i < c[0]; i++) printf("%u ", C[0][i]);printf("\n");
printf("C_r:");
for(uint32_t i = 0; i < c[1]; i++) printf("%u ", C[1][i]);printf("\n");
printf("NoNeiC_l:");
for(uint32_t i = 0; i < c[0]; i++) printf("%u ", nonNei.getCntNonNei(0, C[0][i]));printf("\n");
printf("NoNeiC_r:");
for(uint32_t i = 0; i < c[1]; i++) printf("%u ", nonNei.getCntNonNei(1, C[1][i]));printf("\n");
printf("X_l:");
for(uint32_t i = sx[0]; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("X_r:");
for(uint32_t i = sx[1]; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
printf("aX_l:");
for(uint32_t i = 0; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("aX_r:");
for(uint32_t i = 0; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
printf("\n");
fflush(stdout);
#endif

            C[tt].changeTo(w, --c[tt]);
            R[tt].changeTo(w, r[tt]++);

            parameter newC = plexUpdateC(tt, deep, r, w, c);
            parameter newSX = plexUpdateX(tt, deep, r, w, x, sx);

            for(uint32_t i = 0; i < newC[zz]; i++) {
                if(!g->connect(w, C[zz][i], tt)) {
                    nonNei.addNonNei(zz, C[zz][i], w);
                }
            }
            for(uint32_t i = newSX[zz]; i < x[zz]; i++) {
                if(!g->connect(w, X[zz][i], tt)) {
                    nonNei.addNonNei(zz, X[zz][i], w);
                }
            }
            for(uint32_t i = 0; i < r[zz]; i++) {
                if(!g->connect(w, R[zz][i], tt)) {
                    nonNei.addNonNei(zz, R[zz][i], w);
                }
            }

            plexEnum(deep + 1, r, newC, x, newSX);
#ifdef LIMITEDNUM
            if(totalMaximalbiplexv2 >= outPutT) return;
#endif
            for(uint32_t i = 0; i < newC[zz]; i++) {
                if(!g->connect(w, C[zz][i], tt)) {
                    nonNei.pop(zz, C[zz][i]);
                }
            }
            for(uint32_t i = newSX[zz]; i < x[zz]; i++) {
                if(!g->connect(w, X[zz][i], tt)) {
                    nonNei.pop(zz, X[zz][i]);
                }
            }
            for(uint32_t i = 0; i < r[zz]; i++) {
                if(!g->connect(w, R[zz][i], tt)) {
                    nonNei.pop(zz, R[zz][i]);
                }
            }

            R[tt].changeTo(w, --r[tt]);
            X[tt].changeTo(w, x[tt]++);
        }
    }
    
    for(uint32_t ii = 0; ii < 2; ii++) {
        uint32_t tt = ii == 0? z : t;
        uint32_t zz = tt ^ 1;
        for(uint32_t i = 0; i < PuvSize[tt]; i++) {
            uint32_t w = Puv[tt][deep][i];
            X[tt].changeTo(w, --x[tt]);
        }
    }
}

biplexv2::parameter biplexv2::plexUpdateC(uint32_t t, uint32_t deep, parameter r, uint32_t u, 
    parameter c) 
{    
    uint32_t z = t ^ 1;
    parameter newC;
    newC[t] = c[t];
    newC[z] = 0;

    for(uint32_t i = 0; i < c[z]; i++) {
        uint32_t v = C[z][i];
        if(g->connect(u, v, t)) {
            C[z].changeTo(v, newC[z]++);
// printf("inUpdateC:%u C %u, newC %u\n", u, v, newC[z]);
        }
        else if(nonNei.getCntNonNei(t, u) < k) {
            if(nonNei.getCntNonNei(z, v) + 1 <= k) {
                // nonNei.addNonNeiIdx(z, v, u, deep + 1, i);
                C[z].changeTo(v, newC[z]++);
// printf("inUpdateCnonNei:%u C %u, newC %u\n", u, v, newC[z]);
// printf("nonNei.getCntNonNei(z, v) %u\n", nonNei.getCntNonNei(z, v));
            }
        }
    }

    if(newC[t]> 0)
    for(uint32_t i = 0; i < nonNei.getCntNonNei(t, u); i++) {
        uint32_t v = nonNei.getBuffer(t, u)[i];
        
        if(nonNei.getCntNonNei(z, v) + 1 == k) {
            if(g->deg(v, z) < c[t]) {
                uint32_t tmp = 0;
                for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    if(C[t].pos(w) < newC[t]) {
// printf("change %u %u\n", w, tmp);
                        C[t].changeTo(w, tmp++);
                    }
                }
                newC[t] = tmp;
            }
            else {
                for(uint32_t j = 0; j < newC[t]; ) {
                    uint32_t w = C[t][j];
                    
                    if(!g->connect(v, w, z)) {
// printf("change2 %u %u\n", w, newC[t]);
                        C[t].changeTo(w, --newC[t]);
                    }
                    else j++;
                }
            }
        }
    }

    return newC;
}

biplexv2::parameter biplexv2::plexUpdateX(uint32_t t, uint32_t deep, parameter r, uint32_t u, 
                        parameter x, parameter sx)
{
    uint32_t z = t ^ 1;
    parameter newSX;
    newSX[t] = sx[t];
    newSX[z] = x[z];

    for(uint32_t i = sx[z]; i < newSX[z]; ) {
        uint32_t v = X[z][i];
        if(g->connect(u, v, t)) {
            X[z].changeTo(v, --newSX[z]);
        }
        else if(nonNei.getCntNonNei(t, u) < k && nonNei.getCntNonNei(z, v) + 1 <= k) {
            X[z].changeTo(v, --newSX[z]);
        }
        else i++;
    }

    if(x[t] - sx[t] > 0)
    for(uint32_t i = 0; i < nonNei.getCntNonNei(t, u); i++) {
        uint32_t v = nonNei.getBuffer(t, u)[i];  
        if(nonNei.getCntNonNei(z, v) + 1 == k) {
            if(g->deg(v, z) < x[t] - newSX[t]) {
                uint32_t tmp = x[t];
                for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    if(newSX[t] <= X[t].pos(w) && X[t].pos(w) < tmp) {
                        X[t].changeTo(w, --tmp);
                    }
                }
                newSX[t] = tmp;
            }
            else {
                uint32_t tmp = x[t];
                for(uint32_t j = newSX[t]; j < tmp; ) {
                    uint32_t w = X[t][j];
                    
                    if(g->connect(v, w, z)) {
                        X[t].changeTo(w, --tmp);
                    }
                    else j++;
                }
                newSX[t] = tmp;
            }
        }
    }

    return newSX;
}


void biplexv2::plexEnumAll(uint32_t deep, parameter r, parameter c, parameter x, parameter sx) {

// printf("deep:%u cnt:%u\n", deep, totalMaximalbiplexv2);fflush(stdout);

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
printf("NoNeiC_l:");
for(uint32_t i = 0; i < c[0]; i++) printf("%u ", nonNei.getCntNonNei(0, C[0][i]));printf("\n");
printf("NoNeiC_r:");
for(uint32_t i = 0; i < c[1]; i++) printf("%u ", nonNei.getCntNonNei(1, C[1][i]));printf("\n");
printf("X_l:");
for(uint32_t i = sx[0]; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("X_r:");
for(uint32_t i = sx[1]; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
printf("aX_l:");
for(uint32_t i = 0; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("aX_r:");
for(uint32_t i = 0; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
fflush(stdout);
#endif
// if(deep == 0)
// printf("c[0]: %u %u %u %u\n", c[0], c[1], r[0], r[1]);
    if(c[0] + r[0] < ls) return;
    if(c[1] + r[1] < rs) return;
    
    // if(c[0] + r[0] == 0 || c[1] + r[1] == 0) return;
    if(c[0] + c[1] == 0) {
        if(x[0] - sx[0] + x[1] - sx[1] == 0) {
            totalMaximalbiplexv2++;
        if (totalMaximalbiplexv2 == outPutT) {printf("done\n"); fflush(stdout);} 

#ifdef ANSPRINT
printf("R_l:");
for(uint32_t i = 0; i < r[0]; i++) printf("%u ", R[0][i]);printf("\n");
printf("R_r:");
for(uint32_t i = 0; i < r[1]; i++) printf("%u ", R[1][i]);printf("\n");

fflush(stdout);
#endif

#ifdef DEBUG
printf("answer\n");
#endif
        }
        return;
    }

    //find the pivot with maximal degree in cand
    uint32_t maxE[2] = {g->n[1] + 1, g->n[0] + 1}, maxUV[2] = {g->n[0], g->n[1]};

#ifdef DEEPTHRESHOLD
    bool hasCommonNeiInC[2] = {false, false};
    bool hasCommonNeiInX[2] = {false, false};
    if(deep <= DEEPTHRESHOLD)
    for(uint32_t t = 0; t < 2; t++) {
        uint32_t z = t ^ 1;//t + z = 1
        for(uint32_t i = 0; i < c[t]; i++) {
            if(nonNei.getCntNonNei(t, C[t][i]) == 0) {
                hasCommonNeiInC[t] = true;
                break;
            }
        }
        for(uint32_t i = sx[t]; i < x[t]; i++) {
            if(nonNei.getCntNonNei(t, X[t][i]) == 0) {
                hasCommonNeiInX[t] = true;
                break;
            }
        }
    }
#endif

    auto findPivotC = [&](uint32_t t, uint32_t z, uint32_t u, uint32_t i) {
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

        if(maxE[t] == g->n[z] + 1 || e > maxE[t]) {
            maxE[t] = e;
            maxUV[t] = u;
        }
    };

    auto findPivotX = [&](uint32_t t, uint32_t z, uint32_t u, uint32_t i) {
        uint32_t e = 0;

        if(g->deg(u, t) > x[z]) {
            for(uint32_t j = 0; j < x[z]; j++) {
                uint32_t v = X[z][j];
                if(g->connect(u, v, t)) {
                    e++;
                }
            }
        }
        else {
            for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                if(X[z].pos(v) < x[z]) {//in X
                    e++;
                }
            }
        }

        if(maxE[t] == g->n[z] + 1 || e > maxE[t]) {
            maxE[t] = e;
            maxUV[t] = u;
        }
    };


    for(uint32_t t = 0; t < 2; t++) {
        uint32_t z = t ^ 1;//t + z = 1
#ifdef DEEPTHRESHOLD
        if(hasCommonNeiInC[t]) {
            for(uint32_t i = 0; i < c[t]; i++) if(nonNei.getCntNonNei(t, C[t][i]) == 0) {
                findPivotC(t, z, C[t][i], i);
            }
        }
        else
#endif
        for(uint32_t i = 0; i < c[t]; i++) {
            findPivotC(t, z, C[t][i], i);
        }

#ifdef DEEPTHRESHOLD
        if(hasCommonNeiInX[t]) {
            for(uint32_t i = sx[t]; i < x[t]; i++) if(nonNei.getCntNonNei(t, X[t][i]) == 0) {
                findPivotX(t, z, X[t][i], i);
            }
        }
        else
#endif
        for(uint32_t i = sx[t]; i < x[t]; i++) {
            findPivotX(t, z, X[t][i], i);
        }
    }

    parameter newSC = {0, 0}, newC = c, newSX = sx, newX = x;

    auto computePuv = [&](uint32_t t, uint32_t z, uint32_t u) {
        newSC[t] = 0;
        if(nonNei.getCntNonNei(t, u) == 0) {
            if(C[t].pos(u) < c[t]) {//pivot in C,同一边不用枚举
                newSC[t] = 1;
                C[t].changeTo(u, 0);
            }
            // else {//pivot in X 都不用枚举

            // }
            return;
        }

        for(uint32_t i = 0; i < nonNei.getCntNonNei(t, u); i++) {
            uint32_t v = nonNei.getBuffer(t, u)[i];
            if(g->deg(v, z) < c[t]) {
                uint32_t tmp = c[t];
                for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    uint32_t pw = C[t].pos(w);
                    if(newSC[t] <= pw && pw < c[t]) {
                        C[t].changeTo(w, --tmp);
                    }
                }
                newSC[t] = tmp;
            }
            else {
                uint32_t tmp = c[t];
                for(uint32_t j = newSC[t]; j < tmp; ) {
                    uint32_t w = C[t][j];
                    
                    if(g->connect(v, w, z)) {
                        C[t].changeTo(w, --tmp);
                    }
                    else j++;
                }
                newSC[t] = tmp;
            }
        }
    };
    
    if(maxUV[0] != g->n[0]) computePuv(0, 1, maxUV[0]);
    if(maxUV[1] != g->n[1]) computePuv(1, 0, maxUV[1]);

    uint32_t Pu[2] = {c[0] - newSC[0], c[1] - newSC[1]};
    
    uint32_t t = 0;
    if(maxUV[0] == g->n[0] ) t = 1;
    else if(maxUV[1] == g->n[1]) t = 0;
    else if(Pu[0] + maxE[0] < Pu[1] + maxE[1]) t = 1;
    uint32_t z = t ^ 1;

    if(Puv[t][deep].size() < c[t] - Pu[t]) {
        Puv[t][deep].resize((c[t] - Pu[t])*2);
    }
    memcpy(Puv[t][deep].data(), C[t].begin(), sizeof(uint32_t) * newSC[t]);

    if(g->deg(maxUV[t], t) > c[z]) {
        uint32_t tmp = c[z];
        for(uint32_t j = 0; j < tmp; ) {
            uint32_t v = C[z][j];
            if(g->connect(maxUV[t], v, t)) {
                C[z].changeTo(v, --tmp);
            }
            else j++;
        }
        newC[z] = tmp;
    }
    else {
        uint32_t u = maxUV[t];
        uint32_t tmp = c[z];
        for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            if(C[z].pos(v) < c[z]) {//in C
                C[z].changeTo(v, --tmp);
            }
        }
        newC[z] = tmp;
    }
    if(Puv[z][deep].size() < newC[z]) {
        Puv[z][deep].resize(newC[z] * 2);
    }
    memcpy(Puv[z][deep].data(), C[z].begin(), sizeof(uint32_t) * newC[z]);

    uint32_t PuvSize[2];
    PuvSize[t] = newSC[t];
    PuvSize[z] = newC[z];

#ifdef DEBUG
printf("t %u, pivot %u\n", t, maxUV[t]);
for(uint32_t ii = 0; ii < 2; ii++) {
    uint32_t tt = ii == 0? t : z;
    uint32_t zz = tt ^ 1;

    for(uint32_t i = 0; i < PuvSize[tt]; i++) {
        uint32_t w = Puv[tt][deep][i];
        printf("%u ", w);
    }
    printf("\n");
    fflush(stdout);
}
#endif

    for(uint32_t ii = 0; ii < 2; ii++) {
        uint32_t tt = ii == 0? z : t;
        uint32_t zz = tt ^ 1;
        for(uint32_t i = 0; i < PuvSize[tt]; i++) {
            uint32_t w = Puv[tt][deep][i];

#ifdef DEBUG
printf("deep %u tt %u c[tt] %u, w %u\n", deep, tt, c[tt], w);
printf(" R_l:");
for(uint32_t i = 0; i < r[0]; i++) printf("%u ", R[0][i]);printf("\n");
printf(" R_r:");
for(uint32_t i = 0; i < r[1]; i++) printf("%u ", R[1][i]);printf("\n");
printf("C_l:");
for(uint32_t i = 0; i < c[0]; i++) printf("%u ", C[0][i]);printf("\n");
printf("C_r:");
for(uint32_t i = 0; i < c[1]; i++) printf("%u ", C[1][i]);printf("\n");
printf("NoNeiC_l:");
for(uint32_t i = 0; i < c[0]; i++) printf("%u ", nonNei.getCntNonNei(0, C[0][i]));printf("\n");
printf("NoNeiC_r:");
for(uint32_t i = 0; i < c[1]; i++) printf("%u ", nonNei.getCntNonNei(1, C[1][i]));printf("\n");
printf("X_l:");
for(uint32_t i = sx[0]; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("X_r:");
for(uint32_t i = sx[1]; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
printf("aX_l:");
for(uint32_t i = 0; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("aX_r:");
for(uint32_t i = 0; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
printf("\n");
fflush(stdout);
#endif

            C[tt].changeTo(w, --c[tt]);
            R[tt].changeTo(w, r[tt]++);

            parameter newC = plexUpdateC(tt, deep, r, w, c);
            parameter newSX = plexUpdateX(tt, deep, r, w, x, sx);

            for(uint32_t i = 0; i < newC[zz]; i++) {
                if(!g->connect(w, C[zz][i], tt)) {
                    nonNei.addNonNei(zz, C[zz][i], w);
                }
            }
            for(uint32_t i = newSX[zz]; i < x[zz]; i++) {
                if(!g->connect(w, X[zz][i], tt)) {
                    nonNei.addNonNei(zz, X[zz][i], w);
                }
            }
            for(uint32_t i = 0; i < r[zz]; i++) {
                if(!g->connect(w, R[zz][i], tt)) {
                    nonNei.addNonNei(zz, R[zz][i], w);
                }
            }

            plexEnum(deep + 1, r, newC, x, newSX);
#ifdef LIMITEDNUM
            if(totalMaximalbiplexv2 >= outPutT) return;
#endif
            for(uint32_t i = 0; i < newC[zz]; i++) {
                if(!g->connect(w, C[zz][i], tt)) {
                    nonNei.pop(zz, C[zz][i]);
                }
            }
            for(uint32_t i = newSX[zz]; i < x[zz]; i++) {
                if(!g->connect(w, X[zz][i], tt)) {
                    nonNei.pop(zz, X[zz][i]);
                }
            }
            for(uint32_t i = 0; i < r[zz]; i++) {
                if(!g->connect(w, R[zz][i], tt)) {
                    nonNei.pop(zz, R[zz][i]);
                }
            }

            R[tt].changeTo(w, --r[tt]);
            X[tt].changeTo(w, x[tt]++);
        }
    }
    
    for(uint32_t ii = 0; ii < 2; ii++) {
        uint32_t tt = ii == 0? z : t;
        uint32_t zz = tt ^ 1;
        for(uint32_t i = 0; i < PuvSize[tt]; i++) {
            uint32_t w = Puv[tt][deep][i];
            X[tt].changeTo(w, --x[tt]);
        }
    }
}
