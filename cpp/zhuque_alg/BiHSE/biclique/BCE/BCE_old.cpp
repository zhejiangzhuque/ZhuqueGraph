#include "BCE.h"
#include "../tools/listLinearHeap.hpp"
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

void BCE::run() {
#ifdef BASELINE
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
        if(isMaximal) ans++;
    }
}
printf("ans %u\n", ans);
g->print();
#endif


    uint32_t t = 0, z = 1;
    if(od == "core") {
        if(g->maxDu < g->maxDv) {
            t = 1; z = 0;
        }
    }
    else if(g->core[0] < g->core[1]) {t=1; z=0;}
    printf("t=%d,z=%d\n", t, z);

    std::vector<uint32_t> d[2];
    d[0].resize(g->n[0]);
    d[1].resize(g->n[1]);
    std::vector<uint32_t> id(std::max(g->n[0], g->n[1]));
    std::vector<uint32_t> keys(std::max(g->n[0], g->n[1]));
    for(uint32_t u = 0; u < std::max(g->n[0], g->n[1]); u++) id[u] = u;
    std::vector<uint32_t> buffer[2];
    buffer[0].resize(g->n[0]);
    buffer[1].resize(g->n[1]);

    uint32_t maxDTmp[2];
    // uint32_t lrs[2] = {ls, rs};

uint32_t sumCandidateSets = 0, maxCandidateSets = 0;

    for(uint32_t u = 0; u < g->n[t]; u++) {
        S[t].c = S[t].r = g->n[t];
        S[z].c = S[z].r = g->n[z];
        realRSize[t] = 1;
        realRSize[z] = 0;
        // realR[t].clear(); realR[z].clear();
        buffer[t].clear(); buffer[z].clear();

        for(uint32_t i = g->p[t][u]; i < g->p[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            S[z].swapByPos(--S[z].c, S[z].pos(v));
        }
        
        maxDTmp[0] = maxDTmp[1] = 0;
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

                    d[t][w]++;
                    d[z][v]++;
                    maxDTmp[t] = std::max(maxDTmp[t], d[t][w]);
                    maxDTmp[z] = std::max(maxDTmp[z], d[z][v]);
                }
                else break;

                if(j == 0) break;
            }
        }
        
    
        //core reduction
    //if (lrs[0] > 1 && lrs[1] > 1) {

        ListLinearHeap zheap(S[z].r - S[z].c, maxDTmp[z] + 1), theap(S[t].r - S[t].c, maxDTmp[t] + 1);
        for(uint32_t i = S[z].c; i < S[z].r; i++) keys[i - S[z].c] = d[z][S[z][i]] + 1;
        for(uint32_t i = S[z].c; i < S[z].r; i++) buffer[z].push_back(S[z][i]);
// for(uint32_t i = S[z].c; i < S[z].r; i++) printf("z%u %u\n", S[z][i], keys[i - S[z].c]);
        zheap.init(S[z].r - S[z].c, maxDTmp[z] + 1, id.data(), keys.data());
        
        for(uint32_t i = S[t].c; i < S[t].r; i++) keys[i - S[t].c] = d[t][S[t][i]] + 1;
        for(uint32_t i = S[t].c; i < S[t].r; i++) buffer[t].push_back(S[t][i]);
        theap.init(S[t].r - S[t].c, maxDTmp[t] + 1, id.data(), keys.data());
        
        for(uint32_t i = S[z].c; i < S[z].r; i++) d[z][S[z][i]] = 0;
        for(uint32_t i = S[t].c; i < S[t].r; i++) d[t][S[t][i]] = 0;

// for(uint32_t i = S[t].c; i < S[t].r; i++) printf("t%u %u\n", S[t][i], keys[i - S[t].c]);
// printf("csize %u %u, u %u\n", S[t].r - S[t].c, S[z].r - S[z].c, u);
        while(true) {
            uint32_t u, deg = 0;
            bool updated = false;

            while(zheap.get_min(u, deg) && deg < lrs[t]) {
// printf("popz %u key %u\n", S[z][u + S[z].c], deg);fflush(stdout);
                zheap.pop_min(u, deg);

                for(uint32_t i = S[t].c; i < S[t].r; i++) {
                    if(g->connect(S[z][u + S[z].c], S[t][i], z)) {
                        theap.decrement(i - S[t].c);
// printf("deT %u\n", S[t][i]);fflush(stdout);
                        updated = true;
                    }
                }
            }
            while(theap.get_min(u, deg) && deg < lrs[z] + 1) {
// printf("popt %u key %u\n", S[t][u + S[t].c], deg);fflush(stdout);
                theap.pop_min(u, deg);

                for(uint32_t i = S[z].c; i < S[z].r; i++) {
                    if(g->connect(S[t][u + S[t].c], S[z][i], t)) {
                        zheap.decrement(i - S[z].c);
// printf("deZ %u\n", S[z][i]);fflush(stdout);
                        updated = true;
                    }
                }
            }

            if(!updated) {
                S[t].c = S[t].r;
                S[z].c = S[z].r;
                while(!zheap.empty()) {
                    zheap.pop_min(u, deg);
                    S[z].swapByPos(--S[z].c, S[z].pos(buffer[z][u]));
                }
                while(!theap.empty()) {
                    theap.pop_min(u, deg);
                    S[t].swapByPos(--S[t].c, S[t].pos(buffer[t][u]));
                }
                break;
            }
        }
// printf("csize %u %u, u %u\n", S[t].r - S[t].c, S[z].r - S[z].c, u);
        S[z].x = S[z].c;
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
    //}

sumCandidateSets += S[t].r - S[t].c + S[z].r - S[z].c;
maxCandidateSets = std::max(maxCandidateSets, S[t].r - S[t].c + S[z].r - S[z].c);
        bbranch(0); 
    }
printf("sumCandidateSets: %u\n", sumCandidateSets);
printf("maxCandidateSets: %u\n", maxCandidateSets);
    // bbranch2(0, {0, 0}, {0, 0}, {g->n[0], g->n[1]})

    printf("maxBiCliqueCount: %llu\n", totalMaximalCount);
}



constexpr bool DEBUG = false;
constexpr bool BACKTRACK_DEBUG = false;

/**
 *  deep : search depth
 *  S : X and C
 *  realRSize[0],realRSize[1] : the size of R
 *  totalMaximalCount : the answer
 **/
void BCE::bbranch(uint32_t deep) {


// if(BACKTRACK_DEBUG){
// for(uint32_t i = 0; i < deep; i++) printf(" ");
// printf("st deep %u\n", deep);
// S[0].print();
// printf("\n");
// S[1].print();
// printf("\n");
// printf("realR:");
// for(uint32_t i = 0; i <= 1; i++) {
//     printf(" / ");
//     for(auto u: realR[i]) {
//         printf("%u ", u);
//     }
// }
// printf("\n");

// fflush(stdout);
// }

    if(S[0].cSize() + realRSize[0] < ls) return;
    if(S[1].cSize() + realRSize[1] < rs) return;

    if(S[0].CIsEmpty() && S[1].CIsEmpty()) {
        if(S[0].XIsEmpty() && S[1].XIsEmpty()) {
            if(realRSize[0] >= lrs[0] && realRSize[1] >= lrs[1])
                totalMaximalCount++;
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

        }
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
    for(uint32_t t = 0; t < 1; t++) {
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
    }
    for(uint32_t i = S[1].c; i < S[1].r; i++) {
        uint32_t u = S[1][i];
        uint32_t e = deg[u];
        deg[u] = 0;

        if(e > maxE[0]) {
            maxE[0] = e;
            maxI[1] = i;
            maxUV[1] = u;
        }
    }
    

        //when there is no edge
    if(maxE[0] == 0) {
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
        auto tryAddZSide = [&](uint32_t t) {
            uint32_t z = t ^ 1;
            bool canZAdd = true;

            //add the vertices in z side, the X in z side must be empty
            if(!S[z].XIsEmpty() || realRSize[t] == 0 || realRSize[t] < lrs[t]) return;

            //if a vertex in X and t side is connected to all vertices in z side C
            for(uint32_t i = S[t].x; i < S[t].c; i++) {
                uint32_t u = S[t][i];
                bool all = true;

                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    uint32_t v = S[z][j];

                    if(!g->connect(u, v, t)) {
                        all = false;
                        break;
                    }
                }

                if(all) {
                    canZAdd = false;
                }
            }

            if(canZAdd && realRSize[z] + S[z].cSize() >= lrs[z]) {
// if(DEBUG) {
// printf("t=%u add\n", t);
// }
                totalMaximalCount++;
            }
            // else {
                // printf("canZAdd false\n");
            // }
// if(DEBUG) {
// printf("\n");
// }
        };

        
        if(!S[0].CIsEmpty() && !S[1].CIsEmpty()) {
            tryAddZSide(0);
            tryAddZSide(1);
        }
        else if(S[0].CIsEmpty()) {
            tryAddZSide(0);
        }
        else if(S[1].CIsEmpty()) {
            tryAddZSide(1);
        }

        
        return;
    }

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
    


// if(BACKTRACK_DEBUG){
// for(uint32_t i = 0; i < deep; i++) printf(" ");
// printf("t %u z %u\n", t, z);
// printf("t: maxI %u maxE %u maxUV %u\n", maxI[t], maxE[t], maxUV[t]);
// printf("z: maxI %u maxE %u maxUV %u\n", maxI[z], maxE[z], maxUV[z]);
// }
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;
    
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
    

// if(BACKTRACK_DEBUG){
// for(uint32_t i = 0; i < deep; i++) printf(" ");
// printf("p-node\n");fflush(stdout);
// }
        //put the pivot at the first of R

    if(isPivotInX[t] == false) {
        realRSize[t]++;
// if(DEBUG) realR[t].push_back(maxUV[t]);
  
        S[t].swapByPos(maxI[t], --r[t]);
        S[t].r = r[t];

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
   
        okDegreeSide = z;
        bbranch(deep + 1);
        
    //put pivot into X, child branches do not change R
    // S[t].swapByPos(S[t].c++, S[t].r++);

        S[t].swapByPos(c[t]++, r[t]++);
// if(DEBUG) realR[t].pop_back();
        realRSize[t]--;
    }

    // x[t] = S[t].x;
    // r[t] = S[t].r;
// if(x[t] != S[t].x) printf("xnot same\n");
// if(r[t] != S[t].r) printf("rnot same\n");
    // S[z].x = x[z];
    // S[z].r = r[z];

// if(BACKTRACK_DEBUG) {
// for(uint32_t i = 0; i < deep; i++) printf(" ");
// printf("Deep %u After p-node: wsSize %u, pivot %u-%u\n", deep, wsSize, t, maxUV[t]);
// S[0].print();
// printf("\n");
// S[1].print();
// printf("\n");
// }

    for(uint32_t j = 0; j < wsSize; j++) {
// if(deep == 0) {
//     if(j % 100 == 0) {
//         printf("%u of %u\n", j, wsSize);
//         fflush(stdout);
//     }
// }
        uint32_t w = ws[deep][j];
// if(DEBUG) realR[z].push_back(w);
        realRSize[z]++;
        S[z].swapByPos(S[z].pos(w), --r[z]);
        
        //initial next level parameter
        S[t].x = x[t]; S[z].x = x[z];
        S[t].c = c[t]; S[z].c = c[z];
        S[t].r = r[t]; S[z].r = r[z];
        

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
        
// if(BACKTRACK_DEBUG) {
// for(uint32_t i = 0; i < deep; i++) printf(" ");
// printf("h-node deep %u, %u-%u\n", deep, z, w);
// }
        bbranch(deep + 1);

        // S[t].x = x[t];
        // S[t].r = r[t];

        //put h from R into X
        // S[z].swapByPos(S[z].c++, S[z].r++);
        S[z].swapByPos(c[z]++, r[z]++);
// if(DEBUG) realR[z].pop_back();

// if(BACKTRACK_DEBUG) {
// for(uint32_t i = 0; i < deep; i++) printf(" ");
// printf("After h-node: Deep %u\n", deep);
// S[0].print();
// printf("\n");
// S[1].print();
// printf("\n");
// }
        realRSize[z]--;
    }

    for(uint32_t j = 0; j < wsSize; j++) {
        uint32_t w = ws[deep][j];
        S[z].swapByPos(S[z].pos(w), --c[z]);
    }

    //put pivot back into C
    if(isPivotInX[t] == false) {
        S[t].swapByPos(S[t].pos(maxUV[t]), --c[t]);
    }

// if(BACKTRACK_DEBUG) printf("back deep %u\n\n", deep);
    // S[t].x = x[t]; S[z].x = x[z];
    // S[t].c = c[t]; S[z].c = c[z];
    // S[t].r = r[t]; S[z].r = r[z];
    
}   

// void BCE::run() {
//     for(uint32_t u = 0; u < g->n[0]; u++) {
//         S[0].c = S[0].r = g->n[0];
//         S[1].c = S[1].r = g->n[1];
//         realRSize[0] = 1;
//         realRSize[1] = 0;
//         realR[0].clear(); realR[1].clear();

//         for(uint32_t i = g->p[0][u]; i < g->p[0][u + 1]; i++) {
//             uint32_t v = g->e[0][i];
//             S[1].swapByPos(--S[1].c, S[1].pos(v));
//         }
//         S[1].x = S[1].c;

//         for(uint32_t i = g->p[0][u]; i < g->p[0][u + 1]; i++) {
//             uint32_t v = g->e[0][i];
//             if(g->p[1][v + 1] > 0)
//             for(uint32_t j = g->p[1][v + 1] - 1; j >= g->p[1][v]; j--) {
//                 uint32_t w = g->e[1][j];
                
//                 if(w > u) {
//                     uint32_t pw = S[0].pos(w);

//                     if(pw < S[0].c) {
//                         S[0].swapByPos(--S[0].c, pw);    
//                     }
//                 }
//                 else break;

//                 if(j == 0) break;
//             }
//         }
//         S[0].x = S[0].c;

//         for(uint32_t i = g->p[0][u]; i < g->p[0][u + 1]; i++) {
//             uint32_t v = g->e[0][i];
//             if(g->p[1][v + 1] > 0)
//             for(uint32_t j = g->p[1][v]; j < g->p[1][v + 1]; j++) {
//                 uint32_t w = g->e[1][j];
                
//                 if(w < u) {
//                     uint32_t pw = S[0].pos(w);

//                     if(pw < S[0].x) {
//                         S[0].swapByPos(--S[0].x, pw);
//                     }
//                 }
//                 else break;

//                 if(j == 0) break;
//             }
//         }

//         uint32_t x[2] = {S[0].x, S[1].x};

//         bbranch(0);

        
//     }

    
//     // bbranch2(0, {0, 0}, {0, 0}, {g->n[0], g->n[1]})

//     printf("maxBiCliqueCount: %llu\n", totalMaximalCount);
// }
