#include "BCE.h"

// void BCE::run() {
//     // for(uint32_t i = 0; i < g->n[0]; i++) {
//     //     for(uint32_t j = g->p[0][i]; j < g->p[0][i + 1]; j++) {
//     //         printf("%u %u\n", i, g->e[0][j]);
//     //     }
//     // }

// // g->print();

//     bbranch(0); 
//     // bbranch2(0, {0, 0}, {0, 0}, {g->n[0], g->n[1]})

//     printf("maxBiCliqueCount: %llu\n", totalMaximalCount);
// }
// #define BASELINE
constexpr bool DEBUG = false;
constexpr bool BACKTRACK_DEBUG = false;


void BCE::run() {
    // for(uint32_t i = 0; i < g->n[0]; i++) {
    //     for(uint32_t j = g->p[0][i]; j < g->p[0][i + 1]; j++) {
    //         printf("%u %u\n", i, g->e[0][j]);
    //     }
    // }

// g->print();

#ifdef BASELINE
auto print = [&](uint32_t x, uint32_t y) {
    printf("L:");
    for(uint32_t i = 0; i < g->n[0]; i++) {
        if(x & (1 << i)) printf("%u ", i);
    }
    printf("\n");
    printf("R:");
    for(uint32_t i = 0; i < g->n[1]; i++) {
        if(y & (1 << i)) printf("%u ", i);
    }
    printf("\n");printf("\n");
};
uint32_t ls = 1, rs = 1;
auto check = [&](uint32_t x, uint32_t y)->bool {
    uint32_t cnt = 0;
    for(uint32_t i = 0; i < g->n[0]; i++) {
        if(x & (1 << i)) cnt++;
    }
    if(cnt < ls) return false;
    cnt = 0;
    for(uint32_t i = 0; i < g->n[1]; i++) {
        if(y & (1 << i)) cnt++;
    }
    if(cnt < rs) return false;
    for(uint32_t i = 0; i < g->n[0]; i++) if(x & (1 << i)){
        for(uint32_t j = 0; j < g->n[1]; j++) if(y & (1 << j)) {
            if(!g->connect(i, j, 0)) return false;
        }
    }
    return true;
};
uint32_t ans = 0;
for(uint32_t x = (1<<g->n[0])-1; x >= 1; x--) {
    for(uint32_t y = (1<<g->n[1])-1; y >= 1; y--) {
        if(!check(x, y)) continue;
        
        bool isMaximal = true;
        for(uint32_t i = 0; i < g->n[0]; i++) {
            if(x & (1 << i)) continue;
            if(check(x | (1 << i), y)) {
                isMaximal = false;
                break;
            }
        }
        if(isMaximal)
        for(uint32_t j = 0; j < g->n[1]; j++) {
            if(y & (1 << j)) continue;
            if(check(x, y | (1 << j))) {
                isMaximal = false;
                break;
            }
        }
        if(isMaximal) {
            ans++;
            print(x, y);
        }
    }
}
printf("ans %u\n", ans);
g->print();
#endif

    int32_t t = 0, z = 1;
    // if (g->maxCu2 > g->maxCv2) {t=1;z=0;}
    if (g->maxDu < g->maxDv) {t=1;z=0;}
    printf("t=%d,z=%d\n",t,z);
    for(uint32_t u = 0; u < g->n[t]; u++) {
        S[t].c = S[t].r = g->n[t];
        S[z].c = S[z].r = g->n[z];
        realRSize[t] = 1;
        realRSize[z] = 0;
        realR[t].clear(); realR[z].clear();
realR[t].push_back(u);
        for(uint32_t i = g->p[t][u]; i < g->p[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            S[z].swapByPos(--S[z].c, S[z].pos(v));
        }
        S[z].x = S[z].c;

        for(uint32_t i = g->p[t][u]; i < g->p[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->p[z][v + 1] > 0)
            for(uint32_t j = g->p[z][v + 1] - 1; j >= g->p[z][v]; j--) {
                uint32_t w = g->e[z][j];
                
                if(w > u) {
                    uint32_t pw = S[t].pos(w);
                    if(pw < S[t].c) {
                        S[t].swapByPos(--S[t].c, pw);
                    }
                }
                else break;

                if(j == 0) break;
            }
        }
        S[t].x = S[t].c;

        for(uint32_t i = g->p[t][u]; i < g->p[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->p[z][v + 1] > 0)
            for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                uint32_t w = g->e[z][j];
                
                if(w < u) {
                    uint32_t pw = S[t].pos(w);
                    if(pw < S[t].x) {
                        S[t].swapByPos(--S[t].x, pw);
                    }
                }
                else break;

                if(j == 0) break;
            }
        }

        if(S[z].CIsEmpty()) continue;

        if(S[z].XIsEmpty()) {
            bool f = true;
            for(uint32_t i = S[t].x; i < S[t].r; i++) {
                bool connectAll = true;
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    if(!g->connect(S[t][i], S[z][j], t)) {
                        connectAll = false;
                        break;
                    }
                }
                if(connectAll) {
                    f = false;
                    break;
                }
            }
            if(f) {
if(DEBUG) {
printf("t=%u add0\n", t);
printf("Output:");
for(uint32_t i = 0; i <= 1; i++) {
    printf(" / ");
    for(auto u: realR[i]) {
        printf("%u ", u);
    }
    if(i==z) {
        for(uint32_t i = S[z].c; i<S[z].r;i++) {
            printf("%u ", S[z][i]);
        }
    }
}
printf("\n");
}
                totalMaximalCount++;
            }
        }

        if(S[t].cSize() > 0)
           bbranch(0); 
    }

    
    // bbranch2(0, {0, 0}, {0, 0}, {g->n[0], g->n[1]})

    printf("maxBiCliqueCount: %llu\n", totalMaximalCount);
}


/**
 *  deep : search depth
 *  S : X and C
 *  realRSize[0],realRSize[1] : the size of R
 *  totalMaximalCount : the answer
 **/
void BCE::bbranch(uint32_t deep) {


if(BACKTRACK_DEBUG){
for(uint32_t i = 0; i < deep; i++) printf(" ");
printf("st deep %u\n", deep);
S[0].print();
printf("\n");
S[1].print();
printf("\n");
printf("realR:");
for(uint32_t i = 0; i <= 1; i++) {
    printf(" / ");
    for(auto u: realR[i]) {
        printf("%u ", u);
    }
}
printf("\n");

fflush(stdout);
}
    if(S[0].CIsEmpty() && S[1].CIsEmpty()) {
if(DEBUG) printf("all empty\n");
//         if(S[0].XIsEmpty() && S[1].XIsEmpty()) {
//             if(realRSize[0] > 0 && realRSize[1] > 0)
//                 totalMaximalCount++;
// if(DEBUG) {

// printf("        one ans\n");
// for(uint32_t i = 0; i <= 1; i++) {
//     printf(" / ");
//     for(auto u: realR[i]) {
//         printf("%u ", u);
//     }
// }
// printf("\n");

// }

//         }
        return;
    }

    //find the pivot with maximal degree in cand
    //maxI: the ith vertex in cand is the pivot
    //maxE: the degree of the pivot
    //maxUV: the vertex id of the pivot
    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    maxI[0] = S[0].c;
    maxI[1] = S[1].c;
    maxUV[0] = S[0][S[0].c];
    maxUV[1] = S[1][S[1].c];
    for(uint32_t t = 0; t < 2; t++) {
        uint32_t z = t ^ 1;//t + z = 1

        for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
        // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
            uint32_t u = S[t][i];
            uint32_t e = 0;

            if(g->deg(u, t) > S[z].r - S[z].c) {
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    uint32_t v = S[z][j];
                    if(g->connect(u, v, t)) {
                        e++;
                        deg[v]++;
                    }
                }
            }
            else {
                for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    uint32_t pv = S[z].pos(v);
                    if(S[z].c <= pv && pv < S[z].r) {//in C
                        e++;
                        deg[v]++;
                    }
                }
            }

            if(e > maxE[z]) {
                maxE[z] = e;
                maxI[t] = i;
                maxUV[t] = u;
            }
        }

        for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
        // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
            uint32_t u = S[1][i];
            uint32_t e = deg[u];
            deg[u] = 0;
            if(e > maxE[0]) {
                maxE[0] = e;
                maxI[1] = i;
                maxUV[1] = u;
            }
        }
        break;
    }

        //when there is no edge
//     if(maxE[0] == 0) {
// if(DEBUG) {
// printf("maxE== 0:\n realR[0]:");
// for(auto u: realR[0]) printf("%u ", u); printf("\n");
// printf("realR[1]:");
// for(auto u: realR[1]) printf("%u ", u); printf("\n");

// printf("S[0]:\n");
// S[0].print();
// printf("\nS[1]:\n");
// S[1].print();
// printf("\n");
// }
//         auto tryAddZSide = [&](uint32_t t) {
//             uint32_t z = t ^ 1;
//             bool canZAdd = true;

//             //add the vertices in z side, the X in z side must be empty
//             if(!S[z].XIsEmpty() || realRSize[t] == 0) return;

//             //if a vertex in X and t side is connected to all vertices in z side C
//             for(uint32_t i = S[t].x; i < S[t].c; i++) {
//                 uint32_t u = S[t][i];
//                 bool all = true;

//                 for(uint32_t j = S[z].c; j < S[z].r; j++) {
//                     uint32_t v = S[z][j];

//                     if(!g->connect(u, v, t)) {
//                         all = false;
//                         break;
//                     }
//                 }

//                 if(all) {
//                     canZAdd = false;
//                 }
//             }

//             if(canZAdd && realRSize[t] + S[z].cSize() > 0) {
// if(DEBUG) {
// printf("t=%u add1\n", t);
// printf("Output:");
// for(uint32_t i = 0; i <= 1; i++) {
//     printf(" / ");
//     for(auto u: realR[i]) {
//         printf("%u ", u);
//     }
//     if(i==z) {
//         for(uint32_t i = S[z].c; i<S[z].r;i++) {
//             printf("%u ", S[z][i]);
//         }
//     }
// }
// printf("\n");
// }
//                 totalMaximalCount++;
//             }
//             // else {
//                 // printf("canZAdd false\n");
//             // }
// if(DEBUG) {
// printf("\n");
// }
//         };

        
//         if(!S[0].CIsEmpty() && !S[1].CIsEmpty()) {
//             tryAddZSide(0);
//             tryAddZSide(1);
//         }
//         else if(S[0].CIsEmpty()) {
//             tryAddZSide(0);
//         }
//         else if(S[1].CIsEmpty()) {
//             tryAddZSide(1);
//         }

        
//         return;
//     }

    

    bool isPivotInX[2] = {false, false};
    for(uint32_t t = 0; t < 2; t++) {
        uint32_t z = t ^ 1;//t + z = 1

        for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
        // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
            uint32_t u = S[t][i];
            uint32_t e = 0;

            if(g->deg(u, t) > S[z].r - S[z].c) {
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    uint32_t v = S[z][j];
                    if(g->connect(u, v, t)) {
                        e++;
                    }
                }
            }
            else {
                for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    uint32_t pv = S[z].pos(v);
                    if(S[z].c <= pv && pv < S[z].r) {//in C
                        e++;
                    }
                }
            }

            if(e > maxE[z]) {
                isPivotInX[t] = true;
                maxE[z] = e;
                maxI[t] = i;
                maxUV[t] = u;
            }
        }
    }
    
    //t is the side of pivot, and S[z].cSize() - maxE[z] is smaller than S[t].cSize() - maxE[t]
    uint32_t t = 1;
    if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
        t = 0;
    }
    uint32_t z = t ^ 1;


if(BACKTRACK_DEBUG){
for(uint32_t i = 0; i < deep; i++) printf(" ");
printf("t %u z %u\n", t, z);
printf("t: maxI %u maxE %u maxUV %u\n", maxI[t], maxE[t], maxUV[t]);
printf("z: maxI %u maxE %u maxUV %u\n", maxI[z], maxE[z], maxUV[z]);
}
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;
    
    //put the pivot at the first of R

    if(isPivotInX[t] == false) {
        realRSize[t]++;
if(DEBUG) realR[t].push_back(maxUV[t]);

        
        S[t].swapByPos(maxI[t], --r[t]);
        S[t].r = r[t];
    }
    
    if(isPivotInX[t] == false) {
        if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
            for(uint32_t i = S[z].x; i < S[z].c; i++) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(S[z].x++, i);
                }
            }
        }
        else {
            uint32_t u = maxUV[t];
            uint32_t tmpX = S[z].c;
            for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].x <= pv && pv < S[z].c) {
                    S[z].swapByPos(--tmpX, pv);
                }
            }
            S[z].x = tmpX;
        }
    }
    
    if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
        if(S[z].r > 0)
        for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
            if(!g->connect(maxUV[t], S[z][i], t)) {
                S[z].swapByPos(--S[z].r, i);
            }

            if(i == 0) break;
        }
    }
    else {
        uint32_t u = maxUV[t];
        uint32_t tmpR = S[z].c;
        for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            uint32_t pv = S[z].pos(v);
            if(S[z].c <= pv && pv < S[z].r) {
                S[z].swapByPos(tmpR++, pv);
            }
        }
        S[z].r = tmpR;
    }

    uint32_t wsSize = r[z] - S[z].r;
    if(ws[deep].size() < wsSize) {
        ws[deep].resize(wsSize * 2);
    }
    memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);
    

if(BACKTRACK_DEBUG){
for(uint32_t i = 0; i < deep; i++) printf(" ");
printf("p-node\n");fflush(stdout);
}
    if(isPivotInX[t] == false && S[z].cSize() > 0) {
        if(S[z].XIsEmpty()) {
            bool f = true;
            for(uint32_t i = S[t].x; i < S[t].r; i++) {
                bool connectAll = true;
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    if(!g->connect(S[t][i], S[z][j], t)) {
                        connectAll = false;
                        break;
                    }
                }
                if(connectAll) {
                    f = false;
                    break;
                }
            }
            if(f) {
if(DEBUG) {
printf("t=%u add2\n", t);
printf("Output:");
for(uint32_t i = 0; i <= 1; i++) {
    printf(" / ");
    for(auto u: realR[i]) {
        printf("%u ", u);
    }
    if(i==z) {
        for(uint32_t i = S[z].c; i<S[z].r;i++) {
            printf("%u ", S[z][i]);
        }
    }
}
printf("\n");
}
                totalMaximalCount++;
            }
        }

        uint32_t newR = S[z].r;
        
        for(uint32_t i = S[z].c; i < newR; ) {
            bool degZero = true;
            for(uint32_t j = S[t].c; j < S[t].r; j++) {
                if(g->connect(S[z][i], S[t][j], z)) {
                    degZero = false;
                    break;
                }
            }
            if(degZero) {
                S[z].swapByPos(i, --newR);
            }
            else i++;
        }
        S[z].r = newR;

        if(S[t].cSize() > 0)
            bbranch(deep + 1);

        
    }
    //put pivot into X, child branches do not change R
    // S[t].swapByPos(S[t].c++, S[t].r++);

    if(isPivotInX[t] == false) {
        S[t].swapByPos(c[t]++, r[t]++);
if(DEBUG) realR[t].pop_back();
        realRSize[t]--;
    }

    // x[t] = S[t].x;
    // r[t] = S[t].r;
// if(x[t] != S[t].x) printf("xnot same\n");
// if(r[t] != S[t].r) printf("rnot same\n");
    // S[z].x = x[z];
    // S[z].r = r[z];

if(BACKTRACK_DEBUG) {
for(uint32_t i = 0; i < deep; i++) printf(" ");
printf("Deep %u After p-node: wsSize %u, pivot t%u-%u\n", deep, wsSize, t, maxUV[t]);
//initial next level parameter
        S[t].x = x[t]; S[z].x = x[z];
        S[t].c = c[t]; S[z].c = c[z];
        S[t].r = r[t]; S[z].r = r[z];
S[0].print();
printf("\n");
S[1].print();
printf("\n");
printf("realR:");
for(uint32_t i = 0; i <= 1; i++) {
    printf(" / ");
    for(auto u: realR[i]) {
        printf("%u ", u);
    }
}
printf("\n");
}

    for(uint32_t j = 0; j < wsSize; j++) {
// if(deep == 0) {
//     if(j % 100 == 0) {
//         printf("%u of %u\n", j, wsSize);
//         fflush(stdout);
//     }
// }
        uint32_t w = ws[deep][j];
if(DEBUG) realR[z].push_back(w);
        realRSize[z]++;
        S[z].swapByPos(S[z].pos(w), --r[z]);

        //get C_V
if(BACKTRACK_DEBUG) {
for(uint32_t i = 0; i < deep; i++) printf(" ");
printf("h-node deep %u, z %u, w %u\n", deep, z, w);

} 
  
        //initial next level parameter
        S[t].x = x[t]; S[z].x = x[z];
        S[t].c = c[t]; S[z].c = c[z];
        S[t].r = r[t]; S[z].r = r[z];
        
if(BACKTRACK_DEBUG) {
S[0].print();
printf("\n");
S[1].print();
printf("\n");
} 
        if(S[t].cSize() == 0) continue;   

        if(g->deg(w, z) > S[t].c - S[t].x) {
            for(uint32_t i = S[t].x; i < S[t].c; i++) {
                if(!g->connect(w, S[t][i], z)) {
                    S[t].swapByPos(S[t].x++, i);
                }
            }
        }
        else {
            uint32_t tmpX = S[t].c;
            for(uint32_t j = g->p[z][w]; j < g->p[z][w + 1]; j++) {
                uint32_t u = g->e[z][j];
                uint32_t pu = S[t].pos(u);
                if(S[t].x <= pu && pu < S[t].c) {
                    S[t].swapByPos(--tmpX, pu);
                }
            }
            S[t].x = tmpX;
        }
        
        if(g->deg(w, z) > S[t].c - S[t].r) {
            if(S[t].r > 0)
            for(uint32_t i = S[t].r - 1; i >= S[t].c; i--) {
                if(!g->connect(w, S[t][i], z)) {
                    S[t].swapByPos(--S[t].r, i);
                }

                if(i == 0) break;
            }
        }
        else {
            uint32_t tmpR = S[t].c;
            for(uint32_t j = g->p[z][w]; j < g->p[z][w + 1]; j++) {
                uint32_t u = g->e[z][j];
                uint32_t pu = S[t].pos(u);
                if(S[t].c <= pu && pu < S[t].r) {
                    S[t].swapByPos(tmpR++, pu);
                }
            }
            S[t].r = tmpR;
        }

        if(S[t].XIsEmpty()) {
            bool f = true;
            for(uint32_t i = S[z].x; i < S[z].r; i++) {
                bool connectAll = true;
                for(uint32_t j = S[t].c; j < S[t].r; j++) {
                    if(!g->connect(S[z][i], S[t][j], z)) {
                        connectAll = false;
                        break;
                    }
                }
                if(connectAll) {
                    f = false;
                    break;
                }
            }
            if(f) {
if(DEBUG) {
printf("t=%u add3\n", t);
printf("Output:");
for(uint32_t i = 0; i <= 1; i++) {
    printf(" / ");
    for(auto u: realR[i]) {
        printf("%u ", u);
    }
    if(i==t) {
        for(uint32_t i = S[t].c; i<S[t].r;i++) {
            printf("%u ", S[t][i]);
        }
    }
}
printf("\n");
}
                totalMaximalCount++;
            }
        }

        uint32_t newR = S[t].r;
        for(uint32_t i = S[t].c; i < newR; ) {
            bool degZero = true;
            for(uint32_t j = S[z].c; j < S[z].r; j++) {
                if(g->connect(S[t][i], S[z][j], t)) {
                    degZero = false;
                    break;
                }
            }
            if(degZero) {
                S[t].swapByPos(i, --newR);
            }
            else i++;
        }
        S[t].r = newR;

        if(S[z].cSize() > 0)
            bbranch(deep + 1);

        // S[t].x = x[t];
        // S[t].r = r[t];

        //put h from R into X
        // S[z].swapByPos(S[z].c++, S[z].r++);
        S[z].swapByPos(c[z]++, r[z]++);
if(DEBUG) realR[z].pop_back();

if(BACKTRACK_DEBUG) {
for(uint32_t i = 0; i < deep; i++) printf(" ");
printf("After h-node: Deep %u\n", deep);
S[0].print();
printf("\n");
S[1].print();
printf("\n");
}
        realRSize[z]--;
    }

    for(uint32_t j = 0; j < wsSize; j++) {
        uint32_t w = ws[deep][j];
if(BACKTRACK_DEBUG) {
printf("put back %u from %u to %u\n", w, S[z].pos(w), c[z]-1);
}
        S[z].swapByPos(S[z].pos(w), --c[z]);

    }

    //put pivot back into C
    S[t].swapByPos(S[t].pos(maxUV[t]), --c[t]);

if(BACKTRACK_DEBUG) printf("back deep %u\n\n", deep);
    // S[t].x = x[t]; S[z].x = x[z];
    // S[t].c = c[t]; S[z].c = c[z];
    // S[t].r = r[t]; S[z].r = r[z];
    
}   
