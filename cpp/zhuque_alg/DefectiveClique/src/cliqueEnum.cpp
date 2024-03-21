#include "cliqueEnum.h"

CliqueEnum::CliqueEnum(/* args */)
{
    redeg = NULL;
    readj = NULL;
    recore = NULL;
}

CliqueEnum::~CliqueEnum()
{
    if (redeg != NULL) delete[] redeg; redeg=NULL;
    if (readj != NULL) delete[] readj; readj=NULL;
    if (recore != NULL) delete[] recore; recore=NULL;
}

void CliqueEnum::run() {
    clock_t tm = clock();
    printf("Run: Branch\n");
    global_time = clock();
    cur_out_size = 10;
    if (alg == 1 || alg == 2) basicEnum();
    else if (alg == 3) basicEnum2d();
    else if (alg == 4) pivotEnum();
    else if (alg == 5) pivotOrderingEnum();
    else if (alg == 6) enum2d();
    else if (alg == 7) enum2dOrdering();
    else printf("Error: alg=%d, should be in [1,7]\n", alg);
    printf("All time: %f sec\n", double(clock()-tm)/CLOCKS_PER_SEC);
}
// Basic branch
void CliqueEnum::basicEnum()
{
    vector<int> R;
    initHash();
    R.resize(md+k);
    clock_t tm = clock();
    // if (alg == 1) {
    //     vector<int> C, X;
    //     C.reserve(n);
    //     X.reserve(n);
    //     for (int i = 0; i < n; ++i)
    //         if (deg[i]>0) C.emplace_back(i);
    //     basicBranchN(R,0,0,C,C.size(),X);
    // }
    // else 
    if (alg == 1) {
        vector<upair> C, X;
        C.reserve(n);
        X.reserve(n);
        for (int i = 0; i < n; ++i)
            if (deg[i]>0) C.emplace_back(i, 0);
        basicBranchW(R,0,0,C,C.size(),X);
    }
    else if (alg == 2) {
        vector<upair> C, X;
        C.reserve(n);
        X.reserve(n);
        for (int i = 0; i < n; ++i)
            if (deg[i]>0) C.emplace_back(i, 0);
        basicCORBranch(R,0,0,C,C.size(),X);
    }
    printf("Number of def-cliques: %lld \n", cliquenums);
    printf("Maximum def-clique size: %d \n", MCliqueSize);
    printf("Running time of Enum: %.3f s \n", double(clock()-tm) / CLOCKS_PER_SEC);
}

void CliqueEnum::basicEnum2d()
{
    vector<int> R, C, X;
    vector<bool> visited; 
    initHash();
    R.resize(md+k);
    C.reserve(n);
    X.reserve(n);
    visited.resize(n,false);
    clock_t tm = clock();

    vector<upair> Cx, Xx;
    Cx.reserve(n);
    Xx.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (deg[i] > 0) {
            C.clear(); X.clear();
            Cx.clear(); Xx.clear();
            visited[i] = true;
            for (int j = 0; j < deg[i]; ++j) {
                int u = adj[i][j];
                visited[u] = true;
                if (u > i) C.emplace_back(u);
                else X.emplace_back(u);

                if (u > i) Cx.emplace_back(u, 0);
                else Xx.emplace_back(u, 0);

            }
            int csize = C.size();
            for (int j = 0; j < csize; ++j) {
                int u = C[j];
                for (int l = 0; l < deg[u]; ++l) {
                    int w = adj[u][l];
                    if (!visited[w]) {
                        if (w > i) C.emplace_back(w);
                        else X.emplace_back(w);
                        visited[w] = true;

                        if (w > i) Cx.emplace_back(w,1);
                        else Xx.emplace_back(w, 1);
                    }
                }
            }
            R[0] = i;
            visited[i] = false;
            for (auto u : C) visited[u] = false;
            for (auto u : X) visited[u] = false;
            // basicBranchN(R,1,0,C,C.size(),X);

            basicCORBranch(R,1,0,Cx,Cx.size(),Xx);
            if (cliquenums >= limited_results) break;
        }
    }
    printf("Number of def-cliques: %lld \n", cliquenums);
    printf("Maximum def-clique size: %d \n", MCliqueSize);
    printf("Running time of Enum: %.3f s \n", double(clock()-tm) / CLOCKS_PER_SEC);
}

void CliqueEnum::basicBranchN(vector<int> &R, int rsize, int nnbrs, vector<int> &C, int csize, vector<int> &X)
{
    if (rsize+csize<minsize) return;
    iterations++;
    if (csize <= 0) {
        if (X.empty() && rsize >= minsize) {
            cliquenums++;
            MCliqueSize=max(MCliqueSize,rsize);
            // if (cliquenums % cur_out_size == 0) {
            //     cur_out_size *= 10;
            //     printf("Results %lld, time: %f sec\n", cliquenums, double(clock()-global_time)/CLOCKS_PER_SEC);
            // }
        }
        return;
    }
    //if (cliquenums >= limited_results) return;
    assert(csize>0);
    assert(csize<=C.size());
    vector<int> _C, _X;
    // //_C.reserve(csize);_X.reserve(csize+X.size());
    int v = C[csize-1];
    _C.reserve(csize-1);
    _X.reserve(csize+X.size());

    int _nnbrs = nnbrs;
    for (int i = 0; i < rsize; ++i) 
        if (!is_nbr(v,R[i])) _nnbrs++;
    
    R[rsize] = v;
    for (int i = 0; i < csize; ++i) {
        int u = C[i];
        if (v != u) {
            int nu = 0;
            for (int j = 0; j <= rsize; ++j){
                if (!is_nbr(u,R[j])) {
                    nu++;
                    if (nu+_nnbrs>k) break;
                }
            }
            if (nu+_nnbrs<=k) _C.emplace_back(u);
        }
    }
    for (int i = 0; i < X.size(); ++i) {
        int u = X[i];
        if (v != u) {
            int nu = 0;
            for (int j = 0; j <= rsize; ++j){
                if (!is_nbr(u,R[j])) {
                    nu++;
                    if (nu+_nnbrs>k) break;
                }
            }
            if (nu+_nnbrs<=k) _X.emplace_back(u);
        }
    }
    basicBranchN(R,rsize+1,_nnbrs,_C,_C.size(),_X);
    X.emplace_back(v); 
    if (cliquenums >= limited_results) return;
    //vector<int> ().swap(_C); vector<int> ().swap(_X);
    basicBranchN(R,rsize,nnbrs,C,csize-1,X);
    if (cliquenums >= limited_results) return;
}

void CliqueEnum::basicBranchW(vector<int> &R, int rsize, int nnbrs, vector<upair> &C, int csize, vector<upair> &X)
{
    if (rsize+csize<minsize) return;
    iterations++;
    if (csize <= 0) {
        if (X.empty() && rsize >= minsize) {
            cliquenums++;
            MCliqueSize=max(MCliqueSize,rsize);
            // if (cliquenums % cur_out_size == 0) {
            //     cur_out_size *= 10;
            //     printf("Results %lld, time: %f sec\n", cliquenums, double(clock()-global_time)/CLOCKS_PER_SEC);
            // }
        }
        return;
    }
    //if (cliquenums >= limited_results) return;
    assert(csize>0);
    assert(csize<=C.size());
    vector<upair> _C, _X;
    // //_C.reserve(csize);_X.reserve(csize+X.size());
    auto v = C[0];
    C[0] = C[csize-1];
    C[csize-1] = v;
    _C.reserve(csize-1);
    _X.reserve(csize+X.size());

    int _nnbrs = nnbrs+v.second;
    updates(v.first,_nnbrs,C,0,csize-1,_C);
    updates(v.first,_nnbrs,X,0,X.size(),_X);
    
    R[rsize] = v.first;
    basicBranchW(R,rsize+1,_nnbrs,_C,_C.size(),_X);
    X.emplace_back(v); 
    if (cliquenums >= limited_results) return;
    basicBranchW(R,rsize,nnbrs,C,csize-1,X);
    if (cliquenums >= limited_results) return;
}

void CliqueEnum::basicCORBranch(vector<int> &R, int rsize, int nnbrs, vector<upair> &C, int csize, vector<upair> &X)
{
    if (rsize+csize<minsize) return;
    iterations++;
    if (csize <= 0) {
        if (X.empty() && rsize >= minsize) {
            cliquenums++;
            MCliqueSize=max(MCliqueSize,rsize);
        }
        return;
    }
    if (rsize+csize<minsize) return;
    assert(csize>0);
    assert(csize<=C.size());
    vector<upair> _C, _X;
    // //_C.reserve(csize);_X.reserve(csize+X.size());
    int pos = csize;
    for (int i = 0; i < pos; ++i) {
        auto u = C[i];
        if (u.second == 0) {
            pos--;
            C[i] = C[pos];
            C[pos] = u;
            i--;
        }
    }
    
    if (pos+nnbrs > k) {
        vector<upair> Co, Xo;
        Co.resize(csize);
        Xo.resize(X.size());
        for (int32 i = 0; i < csize; ++i)
            Co[i] = C[i];
        for (int32 i = 0; i < X.size(); ++i)
            Xo[i] = X[i];

        int32 p = 0;
        int32 _nnbrs = nnbrs;
        int32 colen = csize;
        int32 xolen = X.size();
        while (p + nnbrs < k && colen > 0) {
            auto v = Co[0];
            assert(v.second > 0);
            _nnbrs += v.second;
            if (v.first != C[p].first) break;
            if (_nnbrs > k || v.second == 0) break;

            _C.clear(); _X.clear();
            int32 olenn = 0;
            for (int32 i = 0; i < colen; ++i) {
                auto u = Co[i];
                if (u.first != v.first && _nnbrs+u.second <= k) {
                    if (!is_nbr(u.first, v.first)) u.second++;
                    if (_nnbrs+u.second <= k) {
                        Co[olenn++] = u;
                        _C.emplace_back(u);
                    }
                }
            }
            colen = olenn;
            olenn = 0;
            for (int32 i = 0; i < xolen; ++i) {
                auto u = Xo[i];
                if (u.first != v.first && _nnbrs+u.second <= k) {
                    if (!is_nbr(u.first, v.first)) u.second++;
                    if (_nnbrs+u.second <= k) {
                        Xo[olenn++] = u;
                        _X.emplace_back(u);
                    }
                }
            }
            xolen = olenn;
            olenn = 0;
            
            int32 _csize = _C.size();
            if (colen > 0 && Co[0].second > 0 && Co[0].first == C[p+1].first){
                _X.emplace_back(Co[0]);
                _C[0] = _C[_csize-1];
                _csize--;
            }else colen = 0;
            R[rsize + p] = v.first;
            basicCORBranch(R,rsize+1+p,_nnbrs,_C,_csize,_X);
            if (cliquenums >= limited_results) return;
            p++;
            // break;
        }

        auto v = C[0];
        C[0] = C[csize-1];
        C[csize-1] = v;

        X.emplace_back(v); 
        if (cliquenums >= limited_results) return;
        if (rsize+csize<minsize) return;
        basicCORBranch(R,rsize,nnbrs,C,csize-1,X);
        if (cliquenums >= limited_results) return;
        if (rsize+csize<minsize) return;
    }
    else {

        auto v = C[csize-1];
        int _nnbrs = nnbrs+v.second;
        bool reduction = false;
        _C.clear();
        for (int32 i = 0; i < csize-1; ++i) {
            auto u = C[i];
            if (u.second + _nnbrs <= k) {
                if (!is_nbr(v.first, u.first)) 
                    u.second++;
                if (u.second + _nnbrs <= k) {
                    _C.emplace_back(u);
                    if(u.second != 0) reduction = false;
                }
            }
        }

        updates(v.first,_nnbrs,X,0,X.size(),_X);
        
        R[rsize] = v.first;
        basicCORBranch(R,rsize+1,_nnbrs,_C,_C.size(),_X);
        X.emplace_back(v); 
        if (cliquenums >= limited_results) return;
        if (v.second == 0 && reduction) return;
        if (rsize+csize<minsize+1) return;
        basicCORBranch(R,rsize,nnbrs,C,csize-1,X);
        if (cliquenums >= limited_results) return;
    }
}


void CliqueEnum::initHash()
{
    if (cuhash.empty()) cuhash.resize(n);
    for (int i=0; i < n; ++i) {
        int d = deg[i];
        cuhash[i].reserve(d);
        for (int j = 0 ; j < d; ++j) {
            int u = adj[i][j];
            cuhash[i].insert(u);
            if (i > u && !is_nbr(u,i)) {
                printf("Can not find the nbr %d:%d\n",u,i);
                abort();
            } 
        }
    }
}

void CliqueEnum::pivotEnum()
{
    initHash();
    vector<int> R;
    vector<upair> C, X;
    R.resize(md);
    C.reserve(n);
    X.reserve(n);
    clock_t tm = clock();
    for (int i = 0; i < n; ++i) {
        if (deg[i] == 0) continue;
        C.emplace_back(i,0);
    }
    pivotBranch(R,0,0,C,C.size(),X);
    printf("Number of def-cliques: %lld \n", cliquenums);
    printf("Maximum def-clique size: %d \n", MCliqueSize);
    printf("Running time of Enum: %.3f s \n", double(clock()-tm) / CLOCKS_PER_SEC);
}

void CliqueEnum::pivotOrderingEnum()
{
    int32 *nodeset = new int32[n]();
    for (int32 i = 0; i < n; ++i) nodeset[i] = i;
    int32 maxcore = core_decompsition(nodeset, n);
    initHash();
    vector<int> R;
    vector<upair> C, X, X1;
    R.resize(md); C.reserve(n);
    X.reserve(n); X1.reserve(n);
    testC.resize(n, 0);
    clock_t tm = clock();
    for (int i = 0; i < n; ++i) {
        int v = nodeset[i];
        if (deg[v]<=0) continue;
        C.clear(); X.clear(); X1.clear();
        for (int j = 0; j < n; ++j) {
            int u = nodeset[j];
            int c = 0;
            if (deg[u] <= 0) continue;
            if (!is_nbr(v,u)) {
                c++;
                if (c > k) continue;
            }
            if (j < i) X.emplace_back(u,c);
            else if (j > i) {
                //C.emplace_back(u,c);
                if (c == 0) C.emplace_back(u,c);
                else X1.emplace_back(u,c);
            }
        }
        R[0] = v;
        combiNonbrs(R,1,0,C,X,X1,0);
        for(auto u: X1) X.emplace_back(u);
        pivotBranch(R,1,0,C,C.size(),X);
    }
    printf("Number of def-cliques: %lld \n", cliquenums);
    printf("Maximum def-clique size: %d \n", MCliqueSize);
    printf("Running time of Enum: %.3f s \n", double(clock()-tm) / CLOCKS_PER_SEC);
    delete[] nodeset;
}

void CliqueEnum::combiNonbrs(vector<int> &R, int rsize, int nonbrs, vector<upair> &C, vector<upair> &X, vector<upair> &X1, int xsp)
{
    if (nonbrs < k) {
        int xs1 = X1.size();
        int _nnbrs = 0;
        vector<upair> _C, _X, _X1;
        _C.reserve(C.size());
        _X.reserve(X.size());
        _X1.reserve(_X1.size());
        for (int i = xsp; i < xs1; ++i) {
            auto v = X1[i];
            R[rsize] = v.first;
            _nnbrs = nonbrs + v.second;
            _C.clear(); _X.clear(); _X1.clear();
            for (int j = 0; j < deg[v.first]; ++j) testC[adj[v.first][j]] = rsize;
            for (auto u : C) if (u.second+_nnbrs<=k) {
                //if (!is_nbr(u.first,v.first)) u.second++;
                if (testC[u.first] != rsize) u.second++;
                if (u.second+_nnbrs<=k) _C.emplace_back(u);
            }
            for (auto u : X) if (u.second+_nnbrs<=k) {
                //if (!is_nbr(u.first,v.first)) u.second++;
                if (testC[u.first] != rsize) u.second++;
                if (u.second+_nnbrs<=k) _X.emplace_back(u);
            }
            if (_nnbrs < k) {
                int _xsp = 0;
                for (int j = 0; j < xs1; ++j) {
                    auto u = X1[j];
                    if (_nnbrs+u.second <= k && u.first != v.first) {
                        //if (!is_nbr(u.first,v.first)) u.second++;
                        if (testC[u.first] != rsize) u.second++;
                        if (u.second+_nnbrs<=k) {
                            _X1.emplace_back(u);
                            if (j < i) _xsp++;
                        }
                    }
                }
                combiNonbrs(R, rsize+1, _nnbrs, _C, _X, _X1, _xsp);
                for(auto u : _X1) _X.emplace_back(u);
            }
            for (int j = 0; j < deg[v.first]; ++j) testC[adj[v.first][j]] = 0;
            pivotBranch(R, rsize+1, _nnbrs, _C, _C.size(), _X);
        }
    }
}

void CliqueEnum::combiNonbrs1(vector<int> &R, int rsize, int nonbrs, vector<upair> &C, vector<upair> &X, vector<upair> &X1, int xsp)
{
    if (nonbrs < k) {
        int xs1 = X1.size();
        int _nnbrs = 0;
        vector<upair> _C, _X, _X1;
        _C.reserve(C.size());
        _X.reserve(X.size());
        _X1.reserve(_X1.size());
        for (int i = xsp; i < xs1; ++i) {
            auto v = X1[i];
            R[rsize] = v.first;
            _nnbrs = nonbrs + v.second;
            _C.clear(); _X.clear(); _X1.clear();
            updates(v.first,_nnbrs,C,0,C.size(),_C);
            updates(v.first,_nnbrs,X,0,X.size(),_X);
            if (_nnbrs < k) {
                int _xsp = 0;
                for (int j = 0; j < xs1; ++j) {
                    auto u = X1[j];
                    if (_nnbrs+u.second <= k && u.first != v.first) {
                        if (!is_nbr(u.first,v.first)) u.second++;
                        if (u.second+_nnbrs<=k) {
                            _X1.emplace_back(u);
                            if (j < i) _xsp++;
                        }
                    }
                }
                combiNonbrs(R, rsize+1, _nnbrs, _C, _X, _X1, _xsp);
                for(auto u : _X1) _X.emplace_back(u);
            }
            pivotBranch(R, rsize+1, _nnbrs, _C, _C.size(), _X);
        }
    }
}

void CliqueEnum::pivotEnumTwoPhases()
{
    int32 *nodeset = new int32[n]();
    for (int32 i = 0; i < n; ++i) nodeset[i] = i;
    int32 maxcore = core_decompsition(nodeset, n);
    initIndexAndFwd(nodeset);

    vector<int> R;
    vector<upair> C, X;
    vector<bool> visited, computed;
    R.resize(md);
    C.reserve(n);
    X.reserve(n);
    visited.resize(n,false);
    computed.resize(n,false);
    executed.resize(n, false);
    clock_t tm = clock();
    int maximumCsize = 0;
    for (int i = 0; i < n && minsize>=k+2; ++i) {
        int v = nodeset[i];
        computed[v] = true;
        executed[v] = true;
        if (deg[v]<=0) continue;
        C.clear(); X.clear();
        initRoot(v, R, C, X, visited);
        R[0] = v;
        //twoHopNbrs(v,C,X,visited,computed);
        maximumCsize = max(maximumCsize,(int)C.size());
        pivotBranch(R,1,0,C,C.size(),X);
        //if (i == 0) break;
    }

    printf("Number of def-cliques: %lld \n", cliquenums);
    printf("Maximum def-clique size: %d \n", MCliqueSize);
    printf("maximumGap :\t %d \n", maximumGap);
    printf("maximumCsize :\t %d \n", maximumCsize);
    printf("Running time of Enum: %.3f s \n", double(clock()-tm) / CLOCKS_PER_SEC);
}

void CliqueEnum::pivotBranch(vector<int32> &R, int32 rsize, int32 nonbrs, vector<upair> &C, int csize, vector<upair> &X)
{
    if (rsize+csize<minsize) return;
    if (C.size() == 0) {
        if (X.size() == 0 && rsize>=minsize) {
            ++cliquenums; MCliqueSize = max(MCliqueSize, rsize);
            // if (cliquenums % cur_out_size == 0) {
            //     cur_out_size *= 10;
            //     printf("Results %lld, time: %f sec\n", cliquenums, double(clock()-global_time)/CLOCKS_PER_SEC);
            // }
        }
        return;
    }
    //if (cliquenums >= limited_results) return;
    vector<upair> _C, _X;
    _C.reserve(C.size());
    _X.reserve(C.size()+X.size());
    int32 _nonbrs = nonbrs;
    int32 pivotSize = csize;
    int32 _psize = 0;
    int32 a = 0;
    pivot(rsize, C, X, pivotSize);
    if (nonbrs < k) {
       sort(C.begin(), C.begin() + pivotSize, [&](upair a, upair b){return a.second > b.second;});
       if (C[0].second == 0) {
           int id = 0, nid = 0;
           for (int i = pivotSize; i < csize; ++i) {
               auto u = C[i];
               if (u.second > nid) {
                   id = i; nid = u.second;
               }
           }
           if (nid > 0) {
               auto u = C[id];
               C[id] = C[pivotSize];
               C[pivotSize++] = C[0];
               C[0] = u;
           }
       }
    }
    maximumGap = max(maximumGap, pivotSize);
    for (int i = 0; i < pivotSize; ++i) {
        auto v = C[i];
        R[rsize] = v.first;
        _nonbrs = nonbrs + v.second;
        updates(v.first, _nonbrs, C, i+1, C.size(), _C);
        updates(v.first, _nonbrs, X,   0, X.size(), _X);
        pivotBranch(R, rsize+1, _nonbrs, _C, _C.size(), _X);
        X.emplace_back(v);
    }
}

void CliqueEnum::pivotBranch1(vector<int32> &R, int32 rsize, int32 nonbrs, vector<upair> &C, int csize, vector<upair> &X)
{
    //if (rsize+csize<minsize) return;
    if (C.size() == 0) {
        if (X.size() == 0 && rsize>=minsize) {
            ++cliquenums; MCliqueSize = max(MCliqueSize, rsize);
            // if (R[0] == 16 && R[1] <= 6) {
            //     printf("Cli %lld:", cliquenums);
            //     for (int i = 0; i < rsize; ++i) printf(" %d", R[i]);
            //     printf("\n");
            // }
            return;
        }
        return;
    }
    vector<upair> _C, _X;
    _C.reserve(C.size());
    _X.reserve(C.size()+X.size());
    int32 _nonbrs = nonbrs;
    int32 pivotSize = csize;
    int32 _psize = 0;
    int32 a = 0;
    pivot(rsize, C, X, pivotSize);
    maximumGap = max(maximumGap, pivotSize);
    if (nonbrs < k) {
        sort(C.begin(), C.begin() + pivotSize, [&](upair a, upair b){return a.second > b.second;});
        if (C[0].second == 0) {
            int id = 0, nid = 0;
            for (int i = pivotSize; i < csize; ++i) {
                auto u = C[i];
                if (u.second > nid) {
                    id = i; nid = u.second;
                }
            }
            if (nid > 0) {
                auto u = C[id];
                C[id] = C[pivotSize];
                C[pivotSize++] = C[0];
                C[0] = u;
            }
        }
        for (int i = 0; i < pivotSize; ++i) {
            auto v = C[i];
            R[rsize] = v.first;
            _nonbrs = nonbrs + v.second;
            int nvr = 0, tc = C.size(); 
            _C.clear();
            for (int32 j = i+1; j < tc; ++j) {
                auto u = C[j];
                if (u.second + _nonbrs <= k) {
                    if (!is_nbr(v.first, u.first)) u.second++;
                    if (u.second + _nonbrs <= k) {
                        _C.emplace_back(u);
                        if (u.second == 0) nvr++;
                    }
                }
            }
            if (nvr+rsize+1+k-_nonbrs < minsize) continue;
            updates(v.first, _nonbrs, X,   0, X.size(), _X);
            pivotBranch1(R, rsize+1, _nonbrs, _C, _C.size(), _X);
            X.emplace_back(v);
        }
    }
    else {
        for (int i = 0; i < pivotSize; ++i) {
            auto v = C[i];
            R[rsize] = v.first;
            _nonbrs = nonbrs + v.second;
            updates(v.first, _nonbrs, C, i+1, C.size(), _C);
            if (_C.size()+rsize+1 < minsize) continue;
            updates(v.first, _nonbrs, X,   0, X.size(), _X);
            pivotBranch1(R, rsize+1, _nonbrs, _C, _C.size(), _X);
            X.emplace_back(v);
        }
    }
}

void CliqueEnum::twoHopNbrs(int v, vector<upair> &C,vector<upair> &X, vector<bool> &visited, vector<bool> &computed)
{
    C.clear(); X.clear();
    int d = deg[v];
    visited[v] = true;
    for (int i = 0; i < d; ++i) {
        int u = adj[v][i];
        if (computed[u]) X.emplace_back(u,0);
        else C.emplace_back(u,0);
        visited[u] = true;
    }

    for (int i = 0; i < d && k > 0; ++i) {
        int u = adj[v][i];
        int du = deg[u];
        for (int j = 0; j < du; ++j) {
            int w = adj[u][j];
            if (!visited[w]) {
                if (computed[w]) X.emplace_back(w,1);
                else C.emplace_back(w,1);
            }
            visited[w] = true;
        }
    }
    visited[v] = false;
    for (auto u: C) visited[u.first] = false;
    for (auto u: X) visited[u.first] = false;
}

void CliqueEnum::subgraph(vector<upair> &C)
{
    for (auto &a : C) {
        subdeg[a.first] = 0; //subadj[a].clear();
    }
    for (int i = 0; i < C.size(); ++i) {
        int v = C[i].first;
        for (int j = i+1; j < C.size(); ++j) {
            int u = C[j].first;
            if (is_nbr(v, u)) {
                subdeg[v]++;
                subdeg[u]++;
                //subadj[v].emplace_back(u);
                //subadj[u].emplace_back(v);
            }
        }
    }
}

void CliqueEnum::enum2d()
{
    assert(minsize>=k+2);
    int32 *nodeset = new int32[n];
    for (int32 i = 0; i < n; ++i) nodeset[i] = i;
    int32 maxcore = core_decompsition(nodeset, n);
    initIndexAndFwd(nodeset);

    clock_t tm = clock();
    int32 psize = 1;
    vector<bool> visited(n,false);
    vector<int32> R, P;
    vector<upair> C, X;
    P.resize(maxcore+k+1);
    R.resize(maxcore+k+1);
    C.reserve(maxcore+k+1);
    X.reserve(md);
    executed.resize(n, false);
    for (int32 i = 0; i < n; ++i) {
        int32 v = nodeset[i];
        executed[v] = true;
        if (core[v]+k+1 < minsize) continue;

        TwoHopNbrs(v,R,C,X,visited);
        pivotBranch(R, 1, 0, C, C.size(), X);
    }
    printf("Number of def-cliques: %lld \n", cliquenums);
    printf("Maximum def-clique size: %d \n", MCliqueSize);
    printf("Running time of Enum: %.3f s \n", double(clock()-tm) / CLOCKS_PER_SEC);
    delete[] nodeset;
}

void CliqueEnum::enum2dOrdering()
{
    assert(minsize>=k+2);
    int32 *nodeset = new int32[n]();
    for (int32 i = 0; i < n; ++i) nodeset[i] = i;
    int32 maxcore = core_decompsition(nodeset, n);
    printf("Max Core: %d\n", maxcore);
    initIndexAndFwd(nodeset);

    clock_t tm = clock();
    int32 psize = 1;
    int csize = (maxcore * md)/(minsize-k-1);
    csize = csize > n ? n:csize;
    vector<bool> visited(n,false);
    vector<int32> R;
    vector<upair> C, X, X1;
    R.resize(maxcore+k+1);
    C.reserve(maxcore+csize);
    X.reserve(maxcore+csize);
    X1.reserve(maxcore+csize);
    executed.resize(n, false);
    for (int32 i = 0; i < n; ++i) {
        int32 v = nodeset[i];
        executed[v] = true;
        if (core[v]+k+1 < minsize) continue;

        initRoot1(v, R, C, X, visited);
        if (C.size() + 1 < minsize) continue;
        R[0] = v;
        pivotBranch1(R, 1, 0, C, C.size(), X);
    }
    printf("Number of def-cliques: %lld \n", cliquenums);
    printf("Maximum def-clique size: %d \n", MCliqueSize);
    printf("Running time of Enum: %.3f s \n", double(clock()-tm) / CLOCKS_PER_SEC);
    delete[] nodeset;
}

void CliqueEnum::recursion(vector<int32> &R, int32 rsize, int32 nonbrs, vector<upair> &C, vector<upair> &X, vector<int32> &P, int32 &psize)
{
    if (C.size() + rsize < minsize ) return;
    if (C.size() == 0) {
        if (X.size() == 0 && rsize >= minsize) {
            ++cliquenums; MCliqueSize = max(MCliqueSize, rsize);
            // printf("Cli %lld:", cliquenums);
            // for (int i = 0; i < rsize; ++i) printf(" %d", R[i]);
            // printf("\n");
            return;
        }
        return;
    }
    vector<int32> _P;
    vector<upair> _C, _X;
    _C.reserve(C.size());
    _X.reserve(C.size()+X.size());
    _P.resize(C.size() + rsize);
    int32 _nonbrs = nonbrs;
    int32 pivotSize = C.size();
    int32 _psize = 0;
    int32 a = 0;
    pivot(rsize, C, X, pivotSize);
    if (nonbrs < k) {
        sort(C.begin(), C.begin() + pivotSize, [&](upair a, upair b){return a.second > b.second;});
    }
    for (int i = 0; i < pivotSize; ++i) {
        if (rsize + C.size() - a++ < minsize) break;
        int32 v = C[i].first;
        int32 nv = C[i].second;
        R[rsize] = v;
        _nonbrs = nonbrs + nv;
        updates(v, _nonbrs, C, i+1, C.size(), _C);
        updates(v, _nonbrs, X,   0, X.size(), _X);
        _P[0] = v; _psize = 1;
        recursion(R, rsize+1, _nonbrs, _C, _X, _P, _psize);
        X.emplace_back(v, nv);
    }
}

void CliqueEnum::updates(int32 v, int32 nonbrs, vector<upair> &C, int32 s, int32 t, vector<upair> &res)
{
    res.clear();
    if (s >= t) return;
    for (int32 i = s; i < t; ++i) {
        auto u = C[i];
        if (u.second + nonbrs <= k) {
            if (!is_nbr(v, u.first)) u.second++;
            if (u.second + nonbrs <= k)
                res.emplace_back(u);
        }
    }
}

void CliqueEnum::initIndexAndFwd(int32 *nodeset)
{
    cuhash.resize(n);
    fwdadj.resize(n);
    vector<bool> visited(n, false);
    for (int i = 0; i < n; ++i) {
        int v = nodeset[i];
        int d = deg[v];
        int cv = core[v];
        cuhash[v].reserve(d);
        fwdadj[v].reserve(cv);
        for (int j = 0; j < d; ++j) {
            int u = adj[v][j];
            if (!visited[u]) 
                fwdadj[v].emplace_back(u);
            cuhash[v].insert(u);
        }
        visited[v] = true;
    }
    // for (int i = 0; i < n; ++i) 
    // {
    //     int v = nodeset[i];
    //     printf("v=%d, fwdadj:",v);
    //     for (int j = 0; j < fwdadj[v].size(); ++j)
    //         printf(" %d",fwdadj[v][j]);
    //     printf("\n");
    // }
    // test
    // for (int i = 0; i < n; ++i) {
    //     int d = deg[i];
    //     for (int j = 0; j < d; ++j) {
    //         int u = adj[i][j];
    //         if (!cuhash[i].find(u)) 
    //             printf("Error: The neighbor node %d of %d is not included in the hash table.\n", u, i);
    //     }
    // }
}

inline bool CliqueEnum::is_nbr(int32 u, int32 v) {
    if (cuhash[u].getsize() <= 0) return false;
    else return cuhash[u].find(v);
}

void CliqueEnum::TwoHopNbrs(int v, vector<int32> &R, vector<upair> &C, vector<upair> &X, vector<bool> &visited)
{
    C.clear(); X.clear();
    visited[v] = true;
    bool it = true;
    for (auto u : fwdadj[v]) {
        visited[u] = true;
        C.emplace_back(u, 0);
    }
    int len = C.size();
    int32 d = deg[v];
    for (int i = 0; i < d; ++i) {
        int u = adj[v][i];
        if (!visited[u]) X.emplace_back(u, 0);
        visited[u] = true;
    }
    for (int32 i = 0; i < len && k > 0; ++i) {
        int32 u = C[i].first;
        int32 du = deg[u];
        for (int32 j = 0; j < du; ++j) {
            int32 w = adj[u][j];
            if (!visited[w]) {
                visited[w] = true;
                if (executed[w]) X.emplace_back(w, 1);
                else C.emplace_back(w,1);
            }
        }
    }
    for (auto &u : C) visited[u.first] = false;
    for (auto &u : X) visited[u.first] = false;
    //for (int i = 0; i < d; ++i) visited[adj[v][i]] = false;
    visited[v] = false;
    R[0] = v;
}

void CliqueEnum::initRoot1(int32 v, vector<int32> &R, vector<upair> &C, vector<upair> &X, vector<bool> &visited)
{
    C.clear(); X.clear();
    visited[v] = true;
    vector<vectori> temp(2);
    int p = 0, q = 1;
    bool it = true;
    temp[p].reserve(fwdadj[v].size());
    temp[q].reserve(fwdadj[v].size());
    for (auto u : fwdadj[v]) {
        visited[u] = true;
        temp[p].emplace_back(u);
    }
    while (it) {
        vector<int> &atemp = temp[p];
        vector<int> &btemp = temp[q];
        btemp.clear();
        it = false;
        for (auto u : atemp) {
            int d = 0;
            int maxc = 0, clnums = 0;
            for (auto w : atemp) if (u != w && is_nbr(u,w)) d++;
            if (d+k+2 >= minsize) btemp.emplace_back(u);
        }
        p = q; q = 1-p;
        if (btemp.size() != atemp.size()) it = true;
    }
    for (auto u : temp[p]) C.emplace_back(u,0);
    int len = C.size();
    int32 d = deg[v];
    for (int i = 0; i < d; ++i) {
        int u = adj[v][i];
        if (core[u]+k+1 >= minsize && !visited[u]){
            int ud = 0;
            for (auto w : temp[p]) if (is_nbr(u,w)) ud++;
            if (ud+k+2 >= minsize) X.emplace_back(u, 0);
        }
        visited[u] = true;
    }
    temp[q].clear();
    for (int32 i = 0; i < len && k > 0; ++i) {
        int32 u = C[i].first;
        int32 du = deg[u];
        for (int32 j = 0; j < du; ++j) {
            int32 w = adj[u][j];
            if (!visited[w] && core[w]+k+1 >= minsize) {
                visited[w] = true;
                temp[q].emplace_back(w);
                int xd = 0;
                for (auto x : temp[p]) if (is_nbr(x,w)) xd++;
                if (xd+k+1 < minsize) continue;
                if (executed[w]) X.emplace_back(w, 1);
                else C.emplace_back(w,1);
            }
        }
    }
    for (auto u : temp[q]) visited[u] = false;
    for (int i = 0; i < d; ++i) visited[adj[v][i]] = false;
    visited[v] = false;
    R[0] = v;
}

void CliqueEnum::initRoot(int32 v, vector<int32> &R, vector<upair> &C, vector<upair> &X, vector<bool> &visited)
{
    C.clear(); X.clear();
    visited[v] = true;
    vector<vectori> temp(2);
    int p = 0, q = 1;
    bool it = true;
    temp[0].reserve(fwdadj[v].size());
    temp[1].reserve(fwdadj[v].size());
    for (auto u : fwdadj[v]) {
        visited[u] = true;
        temp[0].emplace_back(u);
    }
    while (it) {
        vector<int> &atemp = temp[p];
        vector<int> &btemp = temp[q];
        btemp.clear();
        it = false;
        for (auto u : atemp) {
            int d = 0;
            int maxc = 0, clnums = 0;
            for (auto w : atemp) if (u != w && is_nbr(u,w)) d++;
            if (d+k+2 >= minsize) btemp.emplace_back(u);
        }
        p = q; q = 1-p;
        if (btemp.size() != atemp.size()) it = true;
    }
    for (auto u : temp[p]) C.emplace_back(u,0);
    int len = C.size();
    int32 d = deg[v];
    for (int i = 0; i < d; ++i) {
        int u = adj[v][i];
        if (core[u]+k+1 >= minsize && !visited[u]){
            int ud = 0;
            for (auto w : temp[p]) if (is_nbr(u,w)) ud++;
            if (ud+k+2 >= minsize) X.emplace_back(u, 0);
        }
        visited[u] = true;
    }
    temp[q].clear();
    for (int32 i = 0; i < len && k > 0; ++i) {
        int32 u = C[i].first;
        int32 du = deg[u];
        for (int32 j = 0; j < du; ++j) {
            int32 w = adj[u][j];
            if (!visited[w] && core[w]+k+1 >= minsize) {
                visited[w] = true;
                temp[q].emplace_back(w);
                int xd = 0;
                for (auto x : temp[p]) if (is_nbr(x,w)) xd++;
                if (xd+k+1 < minsize) continue;
                if (executed[w]) X.emplace_back(w, 1);
                else C.emplace_back(w,1);
            }
        }
    }
    for (auto u : temp[q]) visited[u] = false;
    for (int i = 0; i < d; ++i) visited[adj[v][i]] = false;
    visited[v] = false;
    R[0] = v;
}

bool CliqueEnum::pivot(int32 rsize, vector<upair> &C, vector<upair> &X, int &pivotSize)
{
    int32 pivot = -1;
    int32 maxd = 0;
    // for (int32 i = 0; i < X.size(); ++i) {
    //     int32 v = X[i].first;
    //     int32 nv = X[i].second;
    //     if (nv == 0) {
    //         int cnt = 0;
    //         for (int32 j = 0; j < C.size(); ++j)
    //             if (is_nbr(C[j].first, v)) cnt++;
    //         if (cnt > maxd) {pivot = v; maxd = cnt;}
    //     }
    // }
    if (rsize <= 0) {
        for (int32 i = 0; i < C.size(); ++i) {
            auto v = C[i];
            if (v.second == 0) {
                 if (deg[v.first] > maxd) {
                    pivot = v.first; maxd = deg[v.first];
                }
            }
        }
    }
    else 
    {
        for (int32 i = 0; i < C.size(); ++i) {
            int32 v = C[i].first;
            int32 nv = C[i].second;
            if (nv == 0) {
                int cnt = 0;
                for (int32 j = 0; j < C.size(); ++j)
                    if (is_nbr(C[j].first, v)) cnt++;
                if (cnt > maxd) {pivot = v; maxd = cnt;}
            }
        }
    }
    //if (pivot == -1 || rsize + maxd <= minsize) return false;
    if (pivot == -1 ) return false;
    pivotSize = C.size();
    for (int32 i = 0; i < pivotSize; ++i) {
        auto v = C[i];
        if (v.first == pivot) {
            //C[i] = C[0];
            //C[0] = v;
            continue;
        }
        else if (is_nbr(v.first, pivot)) {
            pivotSize--;
            C[i] = C[pivotSize];
            C[pivotSize] = v;
            --i;
        }
    }
    return true;
}