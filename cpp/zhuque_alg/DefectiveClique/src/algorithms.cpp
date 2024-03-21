#include <assert.h>

#include "algorithms.h"

Algorithm::Algorithm(/* args */)
{
    n = m = md = 0;
    deg = NULL;
	adj = NULL;
	datas = NULL;
}

Algorithm::~Algorithm()
{
    //printf("algorithm=%d, mincliquesize=%d\n", algorithm, mincliquesize);
    if (deg != NULL) delete[] deg; deg = NULL;
	if (adj != NULL) delete[] adj; adj = NULL;
	if (datas != NULL) delete[] datas; datas = NULL;
}

void Algorithm::testprintGraph()
{
    for (int32 i = 0; i < n; ++i) {
        printf("nbr[%d]: deg=%d\n", i, deg[i]);
        int32 d = deg[i];
        for (int32 j = 0; j < d; ++j) {
            printf("\t%d\n", adj[i][j]);
        }
    }
}

void Algorithm::read_graph(const char *str)
{
    printf("file: %s\n", str);
    bool is_bin = false;
	clock_t stm = clock();
    if (strstr(str,".bin")) is_bin = true;
    if (is_bin) {
        FILE *in = fopen(str, "rb");
        if (in == NULL) {
            printf("No such file: %s\n", str);
            exit(1);
        }

		size_t FRead = 0;
		FRead = fread(&n, sizeof(int32), 1, in);
		FRead = fread(&m, sizeof(int32), 1, in);
		deg = new int32[n]();
		adj = new int32*[n];
		datas = new int32[m]();
		FRead = fread(deg, sizeof(int32), n, in);
        FRead = fread(datas, sizeof(int32), m, in);
		fclose(in);

		for (int32 i = 0, s = 0; i < n; ++i) { // Construct offs of all vertices
			adj[i] = datas + s;
			s += deg[i];
			md = deg[i] > md ? deg[i] : md;
		}
		printf("n = %d, m = %d, maxdeg = %d\n", n, m, md);

        // int32 testlen = 0;
        // for (int32 i = 0, s = 0; i < V; ++i) {
        //     for (int32 j = 0; j < deg[i]; ++j)
        //     {
        //         if (testlen++ > 10) break;
        //         printf("%d\t%d\t%lf\n",i, adj[i][j].first, adj[i][j].second);
        //     }
        // }
    }
    else {
        FILE *in = fopen(str, "r");
        if (in == NULL) {
            printf("No such file: %s\n", str);
            exit(1);
        }
        char8 line[128];
        fgets(line, 128, in);
        if (sscanf(line, "%d %d", &n, &m) != 2) exit(1);
        //printf("n=%d, m=%d\n", n, m);
        assert(n > 0); assert(m > 0);
        vector<pair<int32,int32>> tempE; tempE.reserve(m);

        if (deg != NULL) exit(1);
        deg = new int32[n]();
        int32 u, v, cnt = 0;
        for (int32 i = 0; i < m && (!feof(in)); ++i) {
            char *r = fgets(line, 128, in);
            //if (feof(in)) break;
            sscanf(line, "%d %d", &u, &v);
            //printf("u=%d, v=%d\n", u, v);
            if (u >= v) continue;
            assert(u < n && u >= 0);
            assert(v < n && v >= 0);
            tempE.emplace_back(u,v);
            deg[u]++; deg[v]++;
        }
        fclose(in);
        m = tempE.size();
        //printf("m=%d\n", m);
        adj = new int32*[n]();
        datas = new int32[m*2]();
        for (int32 i = 0; i < n; ++i){
            int32 d = deg[i];
            md = max(md, d);
            cnt += d; deg[i] = cnt - d;
        }
        //printf("cnt=%d\n", cnt);
        for (int32 i = 0; i < n; ++i) {
            adj[i] = datas + deg[i]; deg[i] = 0;
        }

        for (int32 i = 0; i < m; ++i) {
            u = tempE[i].first;
            v = tempE[i].second;
            assert(u < v);
            adj[v][deg[v]++] = u;
            adj[u][deg[u]++] = v;
        }
        //cnt = 0;
        for (int32 i = 0; i < n; ++i) cnt += deg[i];
        //printf("cnt=%d\n", cnt);
        printf("n = %d, m = %d, maxdeg = %d\n", n, m * 2, md);
    }
	printf("Reading time: %lf s\n", double(clock()-stm)/CLOCKS_PER_SEC);

    vector<int32> tdeg(n,0);
    for (int32 i = 0; i < n; ++i) {
        for (int32 j = 0; j < deg[i]; ++j) {
            int32 v = adj[i][j];
            if (j > 0 && adj[i][j-1] >= v) {
                printf("v=%d, j=%d, adj-1=%d, adj=%d\n",i,j,adj[i][j-1], v);
                assert(adj[i][j-1] < v);
            }
            // if (i < v) {
            //     adj[i][tdeg[i]++] = adj[i][j];
            //     adj[v][tdeg[v]++] = i;
            // }
        }
        //deg[i] = tdeg[i];
    }
}

void Algorithm::scalability(bool randomv, float scal)
{
    srand(0);
    long MAXRANDOMID = RAND_MAX * scal;
    if (randomv) {
        int32 ids = 0;
        vector<int32> randids(n,-1);
        for(int32 i = 0; i < n; ++i) {
            if (rand() < MAXRANDOMID && deg[i] > 0)
                randids[i] = ids++;
        }
        int32 *tempdata = datas;
        m = 0; md = 0;
        for(int32 i = 0; i < n; ++i) {
            int32 newid = randids[i];
            //int32 d = deg[i]; deg[i] = 0;
            if (newid < 0) continue;
            int32 newdeg = 0, d = deg[i];
            adj[newid] = tempdata;
            for (int32 j = 0; j < d; ++j) {
                int32 u = adj[i][j];
                if (randids[u] >= 0) 
                    adj[newid][newdeg++] = randids[u];
            }
            deg[newid] = newdeg;
            tempdata += newdeg;
            m += newdeg;
            md = max(md, newdeg);
        }
        n = ids;
    }
    else {
        m = 0; md = 0;
        vector<int32> tempdeg(n,0);
        for(int32 i = 0; i < n; ++i) {
            int32 d = deg[i];
            deg[i] = 0;
            for (int32 j = 0; j < d; ++j) {
                int32 u = adj[i][j];
                if (i < u) {
                    if (rand() >= MAXRANDOMID) continue;
                    adj[i][tempdeg[i]++] = adj[i][j];
                    adj[u][tempdeg[u]++] = i;
                }
            }
            deg[i] = tempdeg[i];
            m += tempdeg[i];
            md = max(md, tempdeg[i]);
        }
    }
    //printf("Scale=%.1f\%, sv=%d, n=%d, m=%d, maxdeg=%d\n", scal*100, randomv, n, m, md);
}

int Algorithm::core_decompsition(int32 *nodeset, int32 nodesize)
{
    bool flag = nodeset == NULL ? true : false;
    int32 maxcore = 0;
    int32 len     = nodesize;
    vector<int32> bin, pos, curdeg, sequence;
    pos.resize(n);
    bin.resize(md+1, 0);
    curdeg.resize(n, 0);
    sequence.resize(n);
    if (core.empty()) core.resize(n,0);
    for (int32 i = 0; i < len; ++i) flag ? bin[deg[i]]++ : bin[deg[nodeset[i]]]++;
    //for (int32 i = 0; i <= md; ++i) printf("bin[%d]=%d\n", i, bin[i]);
    for (int32 i = 1; i <= md; ++i) bin[i] += bin[i-1];
    for (int32 i = md; i > 0; --i) bin[i] = bin[i-1]; bin[0] = 0;

    //for (int32 i = 0; i <= md; ++i) printf("offs[%d]=%d\n", i, bin[i]);

    for (int32 i = 0; i < len; ++i) {
        int32 v = flag ? i : nodeset[i];
        int32 dv = deg[v];
        int32 posv = bin[dv]++;
        sequence[posv] = v;
        pos[v] = posv;
        curdeg[v] = dv;
    }

    int32 k = 0;
    for (int32 i = 0; i < len; ++i) {
        int32 v = sequence[i];
        int32 d = deg[v];
        k = max(k, curdeg[v]);
        maxcore = max(k,maxcore);
        flag ? i : nodeset[i] = v;
        core[v] = k;
        for (int32 j = 0; j < d; ++j) {
            int32 w = adj[v][j];
            int32 dw = curdeg[w]--;
            if (dw > k) {
                int32 posw = pos[w];
                int32 pdws = bin[dw-1]++;
                if (posw != pdws) {
                    sequence[posw] = sequence[pdws];
                    pos[sequence[posw]] = posw;
                    sequence[pdws] = w;
                    pos[w] = pdws;
                }
            }
        }
    }
    return maxcore;
}

int32 Algorithm::coloring(int32 *nodeset, int32 nodesize)
{
    int32 color_nums = 0, len = nodesize;
    vector<int32> bin(md+1), sequence(len);

    //core_decompsition(nodeset, nodesize);
    //for (int32 i = 0; i < len; ++i) sequence[i] = nodeset[i];
    //len = topKEtaCoreDecompsition(nodeset); 
    //for (int32 i = 0; i < len; ++i) sequence[i] = nodeset[i];

    for (int32 i = 0; i < len; ++i) bin[deg[nodeset[i]]]++;
    for (int32 i = 1; i <= md; ++i) bin[i] += bin[i-1];
    for (int32 i = md; i > 0; --i) bin[i] = bin[i-1]; bin[0] = 0;
    for (int32 i = 0; i < len; ++i) {
        int32 v = nodeset[i];
        int32 posi = bin[deg[v]]++;
        sequence[posi] = v;
    }
    if (colors.empty()) {colors.resize(n, -1);}
    for (int32 i = len-1; i >= 0; --i) {
        int32 maxc = -1, curc = -1;
        int32 u = sequence[i];
        int32 d = deg[u];
        //deg[u] = 0;
        for (int32 j = 0; j < d; ++j) {
            int32 v = adj[u][j];
            int32 cv = colors[v];
            //double p = adj[u][j].first;
            //if (p + 1e-16 > eta) adj[u][deg[u]++] = adj[u][j];
            if (cv >= 0) {bin[cv]++; maxc = max(maxc, cv);}
        }
        for (int32 j = 0; j <= maxc; ++j) {
            if (bin[j] == 0 && curc == -1) curc = j;
            bin[j] = 0;
        }
        if (curc == -1) {
            color_nums = max(color_nums, ++maxc + 1);
            colors[u] = maxc;
        }
        else colors[u] = curc;
    }
    printf("Color nums: %d\n", color_nums);
    return color_nums;
}