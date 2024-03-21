#include "polyEnum.h"

#define INR 1
#define INCR 2

PolyEnum::PolyEnum(/* args */)
{
    redeg = NULL;
    readj = NULL;
    recore = NULL;
}

PolyEnum::~PolyEnum()
{
    if (redeg != NULL) delete[] redeg; redeg=NULL;
    if (readj != NULL) delete[] readj; readj=NULL;
    if (recore != NULL) delete[] recore; recore=NULL;
}

void PolyEnum::run() {
    reids();
    clock_t tm = clock();
    printf("Run: polyEnum\n");
    polyEnum();
    printf("All time: %f sec\n", double(clock()-tm)/CLOCKS_PER_SEC);
}

void PolyEnum::updates(int32 v, int32 nonbrs, vector<upair> &C, int32 s, int32 t, vector<upair> &res)
{
    //res.clear();
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

void PolyEnum::updates1(int32 v, int32 nonbrs, vector<upair> &C, int32 s, int32 t, int32 tn, vector<upair> &Res)
{
    Res.clear();
    //if (s >= t) return;
    for (int32 i = s; i < t; ++i) {
        int32 u = C[i].first;
        int32 nu = C[i].second;
        if (nu + nonbrs <= k) {
            if (!is_nbr(v, u)) nu++;
            if (nu + nonbrs <= k){
                Res.emplace_back(u, nu);
            }
        }
    }
    int32 csize = C.size();
    for (int32 i = tn; i < csize; ++i) {
        int32 u = C[i].first;
        int32 nu = C[i].second;
        if (nu + nonbrs <= k) {
            if (!is_nbr(v, u)) nu++;
            if (nu + nonbrs <= k){
                Res.emplace_back(u, nu);
            }
        }
    }
}

inline bool PolyEnum::is_nbr(int32 u, int32 v) {
    return cuhash[u].find(v);
}

void PolyEnum::reids()
{
    int32 *nodeset = new int32[n];
    for (int32 i = 0; i < n; ++i) nodeset[i] = i;
    maxcore = core_decompsition(nodeset, n);

    vector<int> index(n);
    maps.resize(n);
    for (int i = 0, s = 0, ck = 0; i < n; ++i) {
        if (core[i] <= ck) ;
        else {
            sort(nodeset+s, nodeset+i, [&](int v, int u){return deg[v] < deg[u];});
            s = i;
            ck = core[i];
        }
        if (i+1 == n) {
            sort(nodeset+s, nodeset+n, [&](int v, int u){return deg[v] < deg[u];});
        }
    }
    // sort(nodeset, nodeset+n, [&](int v, int u){return deg[v] < deg[u];});
    for (int i = 0; i < n; ++i) {
        int v = nodeset[i];
        index[v] = i;
        maps[i] = v;
    }
    if (redeg == NULL) redeg = new int32[n]();
    if (readj == NULL) readj = new int32*[n];
    if (recore == NULL) recore = new int32[n]();

    for (int i = 0; i < n; ++i) {
        int d  = deg[i];
        for (int j = 0; j < d; ++j) {
            int u = adj[i][j];
            adj[i][j] = index[u];
            // adj[i][j] = u;
        }
        sort(adj[i], adj[i]+d);
    }

    for (int i = 0; i < n; ++i) {
        redeg[index[i]] = deg[i];
        readj[index[i]] = adj[i];
        // redeg[i] = deg[i];
        // readj[i] = adj[i];
    }
   
    cuhash.resize(n);
    for (int i = 0; i < n; ++i) {
        int d = redeg[i];
        cuhash[i].reserve(d);
        for (int j = 0; j < d; ++j) {
            int u = readj[i][j];
            cuhash[i].insert(u);
        }
    }
    delete [] nodeset;
}

void PolyEnum::polyEnum()
{
    clock_t tm = clock();
    int32 psize = 1;
    vector<bool> visited(n,false);
    vector<int32> R;
    R.resize(md+k+1);
    inR.resize(n,false);
    inRi.resize(n,0);
    int cnts = 0;
    int components = 0; 
    queue<int> Q;
    global_time = clock();
    cur_out_size = 10;
    int i = 0;
    while (i < n && redeg[i]==0) i++;
    if (i < n) {
        R[0] = i;
        upair r = polyExtendMax(R,1,0);
        polyRecursion(R, r.first, r.second, i);
        // if (cliquenums >= limited_results) break;
    }
    printf("Number of def-cliques: %lld \n", cliquenums);
    printf("Maximum def-clique size: %d \n", MCliqueSize);
    printf("Running time of PolyEnum: %.3f s \n", double(clock()-tm) / CLOCKS_PER_SEC);
}

void PolyEnum::polyRelabelVertices(int32 *nodeset)
{
    vector<int> index(n);
    maps.resize(n);
    for (int i = 0; i < n; ++i) {
        int v = nodeset[i];
        index[v] = i;
        maps[i] = v;
    }
    if (redeg == NULL) redeg = new int32[n]();
    if (readj == NULL) readj = new int32*[n];
    if (recore == NULL) recore = new int32[n]();
   
    cuhash.resize(n);
    for (int i = 0; i < n; ++i) {
        int d = deg[i];
        cuhash[i].reserve(d);
        readj[i] = adj[i];
        redeg[i] = d;
        for (int j = 0; j < d; ++j) {
            int u = readj[i][j];
            cuhash[i].insert(u);
        }
    }
}

void PolyEnum::polyRecursion(vector<int> &R, int rsize, int nonbrs, int v)
{
    ++cliquenums; MCliqueSize = max(MCliqueSize, rsize);
    // if (cliquenums % cur_out_size == 0) {
    //     cur_out_size *= 10;
    //     printf("Results %d, time: %f sec \n", cliquenums, double(clock()-global_time) / CLOCKS_PER_SEC);
    // }
    if (cliquenums >= limited_results) return;
    sort(R.begin(),R.begin()+rsize);
    //int pi =  polyGetPI(R, rsize, nonbrs);
    int pi = v;
    //printf("pi=%d\n", pi);
    int pos = 0, nextu = pi;
    for (int i = 0; i < rsize; ++i) if (R[i] == pi) {pos = i+1; break;};
    for (int i = pi+1; i < n; ++i) {
        if (i > nextu && pos < rsize) {nextu = R[pos++]; i--; continue;}
        else if (redeg[i] == 0 || nextu == i) continue;
        polyGenerateAlSat(R, rsize, nonbrs, i);
        if (cliquenums >= limited_results) return;
    }
}

int PolyEnum::polyGetPI(vector<int> &R, int rsize, int nonbrs)
{
    assert(rsize>0);
    if (rsize == 1) return R[0];
    //sort(R.begin(), R.begin()+rsize);
    int lastv = R[rsize-1];
    int _nnbr, nextu, curu = R[0];
    int cu=1, pv=1, pi=curu, cnnbr = 0;
    if (rsize>1) nextu = R[1];
    for (int i = curu; i < lastv; ++i) {
        if (redeg[i] == 0) continue;
        else if (i == curu && cu < rsize) {
            curu = R[cu++];
            if (nextu < curu) {
                checkAdd(R, pv++, cnnbr, nextu, _nnbr);
                nextu = curu;
                cnnbr = _nnbr;
            }
        }
        else if (i < nextu) {
            if (checkAdd(R, pv, cnnbr, i, _nnbr)) {
                int inn = _nnbr-cnnbr;
                while (inn+cnnbr <= k) {
                    pi = nextu;
                    checkAdd(R, pv++, cnnbr, nextu, _nnbr);
                    cnnbr = _nnbr;
                    if (!is_nbr(i,nextu))inn++;
                    if (pv >= rsize) break;
                    nextu = R[pv];
                } 
            }
        }
    }
    return pi;
}

bool PolyEnum::checkAdd(vector<int> &R, int rsize, int nonbrs, int v, int &_nonbrs)
{
    _nonbrs = nonbrs;
    for (int i = 0; i < rsize; ++i) {
        int u = R[i];
        if (!is_nbr(v, u)) {
            _nonbrs++;
            if (_nonbrs>k) return false;
        }
    }
    return true;
}

void PolyEnum::polyInitRoot(int32 v, vector<int32> &R, vector<upair> &P, vector<bool> &visited)
{
    P.clear();
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
    for (auto u : temp[p]) P.emplace_back(u,0);
    int len = P.size();
    int32 d = redeg[v];
    for (int i = 0; i < d; ++i) {
        int u = readj[v][i];
        if (recore[u]+k+1 >= minsize && !visited[u]){
            int ud = 0;
            for (auto w : temp[p]) if (is_nbr(u,w)) ud++;
            if (ud+k+2 >= minsize) P.emplace_back(u, 0);
        }
        visited[u] = true;
    }
    temp[q].clear();
    for (int32 i = 0; i < len && k > 0; ++i) {
        int32 u = P[i].first;
        int32 du = redeg[u];
        for (int32 j = 0; j < du; ++j) {
            int32 w = readj[u][j];
            if (!visited[w] && recore[w]+k+1 >= minsize) {
                visited[w] = true;
                temp[q].emplace_back(w);
                int xd = 0;
                for (auto x : temp[p]) if (is_nbr(x,w)) xd++;
                if (xd+k+1 < minsize) continue;
                P.emplace_back(w, 1);
                // else C.emplace_back(w,1);
            }
        }
    }
    for (auto u : temp[q]) visited[u] = false;
    for (int i = 0; i < d; ++i) visited[readj[v][i]] = false;
    visited[v] = false;
    R[0] = v;
    sort(P.begin(), P.end());
}

upair PolyEnum::polyExtendMax(vector<int> &R, int rsize, int nonbrs)
{
    int newrsize = rsize;
    for (int i = 0; i < rsize; ++i) inR[R[i]] = true;
    for (int i = 0; i < n; ++i) {
        if (redeg[i] == 0) continue;
        if (!inR[i]) {
            int curnnbr = 0;
            for (int j = 0; j < newrsize; ++j) {
                int u = R[j];
                if (!is_nbr(i, u)) {
                    curnnbr++;
                    if (curnnbr+nonbrs>k) break;
                }
            }
            if (nonbrs+curnnbr <= k) {
                R[newrsize++] = i;
                nonbrs += curnnbr;
                //inR[i] = true;
            }
        }
        if (nonbrs == k) {
            int v = R[0];
            int d = redeg[v];
            for (int j = 0; j < d; ++j) {
                int u = readj[v][j];
                if (u <= i || inR[u]) continue;
                int l = 0;
                for (l = 1; l < newrsize; ++l) {
                    int w = R[l];
                    if (!is_nbr(w, u)) break;
                }
                if (l >= newrsize) R[newrsize++] = u;
            }
            break;
        }
    }
    for (int i = 0; i < newrsize; ++i) inR[R[i]] = false;
    return upair(newrsize, nonbrs);
}

bool PolyEnum::pivot(int32 rsize, vector<upair> &C, vector<upair> &X, int &pivotSize)
{
    int32 pivot = -1;
    int32 maxd = 0;
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

    if (pivot == -1 ) return false;
    pivotSize = C.size();
    for (int32 i = 0; i < pivotSize; ++i) {
        auto v = C[i];
        if (v.first == pivot) {
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

void PolyEnum::recursiveReduc(vector<int> &R, int rsize, int nonbrs, vector<upair> &C, vector<upair> &X, vector<upair> &CR, int deep)
{
    int csize = C.size();
    if (C.empty()) {
        if (X.empty()) {
            int _nonbrs=0, _rsize=0;
            vector<int> _R(maxcore+k+1);
            if (validExtendMax(R,rsize,nonbrs,_R,_rsize,_nonbrs,CR)) {
                polyRecursion(_R,_rsize,_nonbrs, R[0]);
            }
        }
        return;
    }
    int xsize = X.size();
    int _nonbrs = 0, psize = csize;
    vector<upair> _C, _X;
    _C.reserve(csize);
    _X.reserve(csize+X.size());
    pivot(rsize,C,X,psize);
    for (int i = 0; i < psize; ++i) {
        auto u = C[i];
        _nonbrs = nonbrs + u.second; 
        _C.clear(); _X.clear();
        updates(u.first,_nonbrs,C,i+1,csize,_C);
        updates(u.first,_nonbrs,X,0 ,xsize,_X);
        updates(u.first,_nonbrs,C,0,i,_X);
        R[rsize] = u.first;
        recursiveReduc(R,rsize+1,_nonbrs,_C,_X, CR, 0);
    }
}

bool PolyEnum::validExtendMax(vector<int> &R, int rsize, int nonbrs, vector<int> &_R, int &_rsize, int &_nonbrs,vector<upair> &CR)
{
    int mxv = 0, mir = 0, nv = 0, v = R[0];
    _rsize = 0; _nonbrs = 0;
    for (int i = 1; i < rsize; ++i) {
        int u = R[i];
        _R[_rsize++] = u;
        inRi[u] = INR;
        assert(u < v);
        if (!is_nbr(u,v)) nv++;
        mir=min(u,mir);
    }
    if (_rsize == 0) mir = 0;
    _nonbrs = nonbrs-nv;
 

    for (auto &u : CR) {
        if (inRi[u.first] == 0) inRi[u.first] = INCR;
        mxv = max(u.first,mxv);
    }
   
    assert(_nonbrs >= 0);
    assert(_nonbrs <= k);
    //_rsize = rsize-1; 
    int newnonbrs = 0;
    bool flag = true;
    //printf("_nonbrs=%d, nv=%d\n",_nonbrs, nv);
    for (int i = 0; i <= mxv; ++i) {
        if (deg[i] == 0 || inRi[i] == INR) continue;
        if (_nonbrs == k) {
            int u = _R[0];
            int du = redeg[u];
            for (int j = 0; j < du; ++j) {
                int w = readj[u][j];
                if (w < i || inRi[w] == INR) continue;
                if (w > mxv) break;
                if (checkAdd(_R,_rsize,_nonbrs,w,newnonbrs)) {
                    if (inRi[w] <= 0) {flag = false; break; }
                    _nonbrs = newnonbrs;
                    _R[_rsize++] = w;
                }
            }
            break;
        }
        else if (checkAdd(_R,_rsize,_nonbrs,i,newnonbrs)) {
            if (inRi[i] <= 0) {flag = false; break; }
            _nonbrs = newnonbrs;
            _R[_rsize++] = i; 
        }
    }
   
    if (!flag) {
        for (int i = 0; i < rsize; ++i) inRi[R[i]] = 0;
        for (auto &u : CR) inRi[u.first] = 0;
        return flag;
    }
    assert(_rsize == CR.size());
    _rsize = rsize;
    _R[rsize-1] = v; inRi[v] = INR;
    for (int i = mir+1; i < n; ++i) {
        if (deg[i] == 0) continue;
        else if (nonbrs == k) {
            int u = _R[0];
            int du = redeg[u];
            for (int j = 0; j < du; ++j) {
                int w = readj[u][j];
                if (w < i || inRi[w] == INR) continue;
                if (checkAdd(_R,_rsize,nonbrs,w,newnonbrs)) {
                    if (w < v) {flag = false; break; }
                    nonbrs = newnonbrs;
                    _R[_rsize++] = w;
                }
            }
            break;
        }
        else if (inRi[i] == INR)continue;
        else if (checkAdd(_R,_rsize,nonbrs,i,newnonbrs)) {
            if (i < v) {flag = false; break; }
            nonbrs = newnonbrs;
            _R[_rsize++] = i;
        }

    }
 
    _nonbrs = nonbrs;
    for (int i = 0; i < rsize; ++i) inRi[R[i]] = 0;
    for (auto &u : CR) inRi[u.first] = 0;

    return flag;
}

void PolyEnum::recursiveInc(vector<int> &R, int rsize, int nonbrs, vector<upair> &C, vector<upair> &X, vector<upair> &CR)
{
    int csize = C.size();
    if (C.empty()) {
        recursiveReduc(R,rsize,nonbrs,C,X, CR, 0);
        return;
    }
    if (nonbrs < k) {
        sort(C.begin(), C.end(), [&](upair a, upair b){return a.second > b.second;});
    }
    int xsize = X.size();
    if (C[0].second == 0) {
        int ncsize = 0, nxsize = 0;
        for (int i = 0; i < csize; ++i) {
            auto u = C[i];
            for (int j = 0; j < csize; ++j) {
                auto w = C[j];
                if (!is_nbr(u.first, w.first)) {
                    u.second++; w.second++;
                }
            }
            if (u.second == 1) {
                R[rsize++] = u.first;
                for (int j = 0; j < xsize; ++j) {
                    auto w = X[j];
                    if (nonbrs+w.second<=k) {
                        if (!is_nbr(u.first,w.first)) w.second++;
                        if (nonbrs+w.second<=k)
                            X[nxsize++] = w;
                    }
                }
                xsize = nxsize; nxsize = 0;
            }
            else {
                u.second = 0;
                C[ncsize++] = u;
            }
        }
        csize = ncsize;
        C.resize(csize); X.resize(xsize);
        recursiveReduc(R,rsize,nonbrs,C,X, CR, 0);
        if (cliquenums >= limited_results) return;
    }
    else {
        int _nonbrs = 0;
        vector<upair> _C, _X;
        _C.reserve(csize);
        _X.reserve(csize+X.size()); 
        for (int i = 0; i < csize; ++i) {
            auto u = C[i];
            _nonbrs = nonbrs + u.second; 
            _C.clear(); _X.clear();
            if (u.second == 0) {
                for (int j=i; j<csize;++j)_C.emplace_back(C[j]);
                for (int j=0; j<i;++j)X.emplace_back(C[j]);
                recursiveInc(R,rsize,_nonbrs,_C,X, CR);
                break;
            }
            updates(u.first,_nonbrs,C,i+1,csize,_C);
            updates(u.first,_nonbrs,X,0 ,xsize,_X);
            updates(u.first,_nonbrs,C,0,i,_X);
            R[rsize] = u.first;
            recursiveInc(R,rsize+1,_nonbrs,_C,_X, CR);
            if (cliquenums >= limited_results) return;
        }
    }
}

void PolyEnum::polyGenerateAlSat(vector<int> &R, int rsize, int nonbrs, int v)
{
    int ssize = 0, psize = rsize-1, nnbrP = 0;
    vector<upair> C;
    vector<int> _R;
    C.resize(rsize);
    _R.resize(rsize);
    for (int i = 0; i < rsize; ++i) {
        int u = R[i];
        if (is_nbr(v,u)) C[psize--] = upair(u,0);
        else C[ssize++] = upair(u,1);
    }
    int _rsize = 1, _nonbrs = 0;
    vector<upair> _C, _X;
    _C.reserve(rsize);
    _X.reserve(rsize); 
    _R[0] = v;

    for (int i = 0; i < ssize; ++i) {
        auto u = C[i];
        if (u.first > v) continue;
        _nonbrs = u.second; 
        _C.clear();
        for (int j = i+1; j < rsize; ++j) {
            auto w = C[j];
            if (_nonbrs+w.second <= k && w.first < v) {
                if (!is_nbr(u.first,w.first))
                    w.second++;
                if (_nonbrs+w.second <= k)
                    _C.emplace_back(w);
            }
        }
        _X.clear();
        for (int j = 0; j < i; ++j) {
            auto w = C[j];
            if (_nonbrs+w.second <= k && w.first < v) {
                if (!is_nbr(u.first,w.first))
                    w.second++;
                if (_nonbrs+w.second <= k)
                    _X.emplace_back(w);
            }
        }
        _R[_rsize] = u.first;
        recursiveInc(_R,_rsize+1,_nonbrs,_C,_X, C);
        if (cliquenums >= limited_results) return;
    }
    _nonbrs = 0; _rsize = 1;
    for (int i = ssize; i < rsize; ++i) {
        auto &u = C[i];
        if (u.first > v) continue;
        _R[_rsize++] = u.first;
        for (int j=i+1; j < rsize; ++j) {
            auto &w = C[j];
            if (w.first < v && !is_nbr(u.first,w.first))
                _nonbrs++;
        }
    }
    for (int i = 0; i < ssize; ++i) {
        auto u = C[i];
        if (u.first > v) continue;
        if (_nonbrs+u.second <= k) {
            for (int j = 1; j < _rsize; ++j) {
                int w = _R[j];
                if (!is_nbr(u.first,w)) {
                    u.second++;
                    if (_nonbrs+u.second>k)
                        break;
                }
            }
        }
 
        if (_nonbrs+u.second<=k) return;
    }
   
    vector<int> _Rn(maxcore+k+1);
    if (validExtendMax(_R,_rsize,_nonbrs,_Rn,rsize,nonbrs,C)) {   
        polyRecursion(_Rn,rsize,nonbrs,v);
        if (cliquenums >= limited_results) return;
    }
}
