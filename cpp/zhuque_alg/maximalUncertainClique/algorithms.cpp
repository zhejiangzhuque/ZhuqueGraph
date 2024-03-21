#include "algorithms.h"
Algorithm::Algorithm(/* args */)
{
    V = E = 0;
	maxdeg = 0;
    algorithm = 0;
    mincliquesize = 0;
    deg = NULL;
	adj = NULL;
	data = NULL;
}

Algorithm::~Algorithm()
{
    //printf("algorithm=%d, mincliquesize=%d\n", algorithm, mincliquesize);
    if (deg != NULL) delete[] deg; deg = NULL;
	if (adj != NULL) delete[] adj; adj = NULL;
	if (data != NULL) delete[] data; data = NULL;
}

void Algorithm::testprintGraph()
{
    for (int i = 0; i < V; ++i) {
        printf("nbr[%d]: deg=%d\n", i, deg[i]);
        int d = deg[i];
        for (int j = 0; j < d; ++j) {
            printf("\t%d, %lf\n", adj[i][j].first, adj[i][j].second);
        }
    }
}

void Algorithm::read_graph(const char *str)
{
    bool is_bin = false;
	clock_t stm = clock();
    if (strstr(str,".bin")) is_bin = true;
    if (is_bin) {
        FILE *in = fopen(str, "rb");
        if (in == NULL) {
            printf("No such file: %s\n", str);
            exit(1);
        }
        printf("Graph name: %s\n", str);

		size_t FRead = 0;
		FRead = fread(&V, sizeof(int), 1, in);
		FRead = fread(&E, sizeof(long long int), 1, in);
		deg = new int[V]();
		adj = new pairs*[V];
		data = new pairs[E]();
		FRead = fread(deg, sizeof(int), V, in);
		for (int i = 0; i < E; ++i) FRead = fread(&data[i].first, sizeof(int), 1, in);
		for (int i = 0; i < E; ++i) FRead = fread(&data[i].second, sizeof(double), 1, in);
		fclose(in);
		for (int i = 0, s = 0; i < V; ++i) { // Construct offs of all vertices
			adj[i] = data + s;
			s += deg[i];
			maxdeg = deg[i] > maxdeg ? deg[i] : maxdeg;
		}
		printf("n = %d, m = %d, maxdeg = %d\n", V, E, maxdeg);

        // int testlen = 0;
        // for (int i = 0, s = 0; i < V; ++i) {
        //     for (int j = 0; j < deg[i]; ++j)
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
        printf("file: %s\n", str);
        char line[128];
        fgets(line, 128, in);
        if (sscanf(line, "%d %d", &V, &E) != 2) exit(1);
        printf("n=%d, m=%d\n", V, E);
        assert(V > 0); assert(E > 0);
        pair<int,int> *tempE = new pair<int,int>[E];
        double *probs = new double[E];

        if (deg != NULL) exit(1);
        deg = new int[V]();

        double p;
        int u, v, cnt = 0;
        int maxid = 0;
        for (int i = 0; i < E; ++i) {
            char *r = fgets(line, 128, in);
            if (feof(in)) break;
            sscanf(line, "%d %d %lf", &u, &v, &p);
            //printf("u=%d, v=%d, p=%lf\n", u, v, p);
            if (u >= v) continue;
            assert(u < V && u >= 0);
            assert(v < V && v >= 0);
            tempE[cnt].first = u;
            tempE[cnt].second = v;
            probs[cnt++] = p;
            deg[u]++; deg[v]++;
        }
        fclose(in);
        E = cnt*2;
        for (int i = 0; i < V; ++i) maxdeg = max(maxdeg, deg[i]);
        for (int i = 1; i < V; ++i) deg[i] += deg[i-1];
        for (int i = V-1; i > 0; --i) deg[i] = deg[i-1]; deg[0] = 0;

        adj = new pairs*[V];
        data = new pairs[E]();
        for (int i = 0; i < cnt; ++i) {
            u = tempE[i].first;
            v = tempE[i].second;
            
            int pu = deg[u]++, pv = deg[v]++; 
            data[pu].first = v;
            data[pu].second = probs[i];
            data[pv].first = u;
            data[pv].second = probs[i];
        }
        for (int i = V-1; i > 0 ; --i) deg[i] -= deg[i-1];
        *adj = data;
        for (int i = 1; i < V; ++i) adj[i] = adj[i-1] + deg[i-1]; 

        delete[] probs;
        delete[] tempE;
        printf("n = %d, m = %d, maxdeg = %d\n", V, E, maxdeg);
    }
	printf("Reading time: %lf sec\n", double(clock()-stm)/CLOCKS_PER_SEC);

    vector<int> tdeg(V,0);
    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < deg[i]; ++j) {
            int v = adj[i][j].first;
            if (i < v) {
                adj[i][tdeg[i]++] = adj[i][j];
                adj[v][tdeg[v]].first = i;
                adj[v][tdeg[v]++].second = adj[i][j].second;
            }
        }
        deg[i] = tdeg[i];
    }
}

void Algorithm::scalability(bool randomv, float scal)
{
    srand(100);
    long MAXRANDOMID = RAND_MAX * scal;
    if (randomv) {
        int ids = 0;
        vector<int> randids(V,-1);
        for(int i = 0; i < V; ++i) {
            if (rand() < MAXRANDOMID && deg[i] > 0)
                randids[i] = ids++;
        }
        pairs *tempdata = data;
        E = 0; maxdeg = 0;
        for(int i = 0; i < V; ++i) {
            int newid = randids[i];
            //int d = deg[i]; deg[i] = 0;
            if (newid < 0) continue;
            int newdeg = 0, d = deg[i];
            adj[newid] = tempdata;
            for (int j = 0; j < d; ++j) {
                int u = adj[i][j].first;
                if (randids[u] >= 0) {
                    adj[newid][newdeg].first = randids[u];
                    adj[newid][newdeg++].second = adj[i][j].second;
                }
            }
            deg[newid] = newdeg;
            tempdata += newdeg;
            E += newdeg;
            maxdeg = max(maxdeg, newdeg);
        }
        V = ids;
    }
    else {
        E = 0; maxdeg = 0;
        vector<int> tempdeg(V,0);
        for(int i = 0; i < V; ++i) {
            int d = deg[i];
            deg[i] = 0;
            for (int j = 0; j < d; ++j) {
                int u = adj[i][j].first;
                if (i < u) {
                    if (rand() >= MAXRANDOMID) continue;
                    adj[i][tempdeg[i]++] = adj[i][j];
                    adj[u][tempdeg[u]].first = i;
                    adj[u][tempdeg[u]++].second = adj[i][j].second;
                }
            }
            deg[i] = tempdeg[i];
            E += tempdeg[i];
            maxdeg = max(maxdeg, tempdeg[i]);
        }
    }
    printf("Scale=%.1f\%, n=%d, m=%d, maxdeg=%d\n", scal*100,  V, E, maxdeg);
}

int Algorithm::core_decompsition(int *nodeset, int nodesize)
{
    bool flag = nodeset == NULL ? true : false;
    if (core.empty()) core.resize(V, 0);
    int *sequence = new int [nodesize];
    int len = nodesize;

    int maxcore = 0;
    int *bin = new int[maxdeg+1]();
    int *pos = new int[V];
    int *curdeg = new int[V];
    memcpy(curdeg, deg, sizeof(int) * V);

    for (int i = 0; i < len; ++i) flag ? bin[curdeg[i]]++ : bin[curdeg[nodeset[i]]]++;
    //for (int i = 0; i <= maxdeg; ++i) printf("bin[%d]=%d\n", i, bin[i]);
    for (int i = 1; i <= maxdeg; ++i) bin[i] += bin[i-1];
    for (int i = maxdeg; i > 0; --i) bin[i] = bin[i-1]; bin[0] = 0;

    //for (int i = 0; i <= maxdeg; ++i) printf("offs[%d]=%d\n", i, bin[i]);

    for (int i = 0; i < len; ++i) {
        int v = flag ? i : nodeset[i];
        int posv = bin[deg[v]]++;
        sequence[posv] = v;
        pos[v] = posv;
    }

    int k = 0;
    for (int i = 0; i < len; ++i) {
        int v = sequence[i];
        int d = deg[v];
        k = max(k, curdeg[v]);
        maxcore = max(k,maxcore);
        flag ? i : nodeset[i] = v;
        core[v] = k;
        for (int j = 0; j < d; ++j) {
            int w = adj[v][j].first;
            int dw = curdeg[w]--;
            if (dw > k) {
                int posw = pos[w];
                int pdws = bin[dw-1]++;
                if (posw != pdws) {
                    sequence[posw] = sequence[pdws];
                    pos[sequence[posw]] = posw;
                    sequence[pdws] = w;
                    pos[w] = pdws;
                }
            }
        }
    }

    delete[] sequence;
    delete[] bin;
    delete[] pos;
    delete[] curdeg;
    return maxcore;
}

int Algorithm::topKEtaCore(int *nodeset)
{
    int nodesetsize = 0;
    int k = mincliquesize - 1;
    double tempsize = 0;
    queue<int> Q;
    vector<bool> visited(V);
    vector<vector<pair<double,int>>> sortedAdj(V);

    for (int i = 0; i < V; ++i) {
        if (deg[i] < k) visited[i] = true;
        //printf("deg[%d]=%d\n", i, deg[i]);
    }

    for (int i = 0; i < V; ++i) {
		int v, d = deg[i];
		double p = 1.0, q = 1.0;
        if (visited[i]) continue;
        sortedAdj[i].reserve(d);
        tempsize += (sizeof(pair<double, int>)*d)/double(1024*1024);
		for (int j = 0; j < d; ++j) {
            v = adj[i][j].first;
            p = adj[i][j].second;
            if (!visited[v] && p + 1e-16 > eta) sortedAdj[i].emplace_back(p, v);
		}
        sort(sortedAdj[i].begin(), sortedAdj[i].end(), greater< pair<double,int> >());
        d = sortedAdj[i].size();
        //printf("nbr[%d]: \n", i);
        for (int j = 0; j < d; ++j) {
            p = q * sortedAdj[i][j].first;
            //printf("\t%lf, %d\n", sortedAdj[i][j].first,  sortedAdj[i][j].second);
            if (p + 1e-16 < eta) {
                d = j; break;
            }
            q = p;
		}
        if (d < k) {Q.emplace(i), visited[i] = true;}
	}
    size_t maxQ = Q.size(); 
    while (!Q.empty()) {
        int u = Q.front(); Q.pop();
        int d = sortedAdj[u].size();
        //printf("deg[%d]=%d\n", u, deg[u]);
        for (int j = 0; j < d; ++j) {
            int v = sortedAdj[u][j].second;
            if (!visited[v]) {
                double q = 1.0, pro = 1.0;
                int dv = sortedAdj[v].size();
                int core_new = 0;
                for (int l = 0; l < dv; ++l) {
                    if (!visited[sortedAdj[v][l].second]) {
                        pro *= sortedAdj[v][l].first;
                        if (pro + 1e-16 < eta) break;
                        core_new++;
                    }
                }
                if (core_new < k) {Q.emplace(v), visited[v] = true;}
            }
        }
        maxQ = max(maxQ, Q.size());
    }
    for (int i = 0; i < V; ++i) {
        if (!visited[i]) nodeset[nodesetsize++] = i;
    }
	printf("Phase 1, left nums: %d\n", nodesetsize);
    return nodesetsize;
}

int Algorithm::topKEtaCoreDecompsition(int *nodeset)
{
    if ( topcore.empty()) topcore.resize(V, 0);
	vector<bool> visited(V, 0);
    vector<vector<pair<double,int>>> sortedAdj(V);
    int mink = mincliquesize - 1;
    double tempsize = 0;
	for (int i = 0; i < V; ++i){
		int d = deg[i];
        if (d < mink) {continue;}
        sortedAdj[i].reserve(d);
        tempsize += (sizeof(pair<double, int>)*d)/double(1024*1024);
        for (int j = 0; j < d; ++j) {
            int v = adj[i][j].first;
            double p = adj[i][j].second;
            if (deg[v] >= mink && p + 1e-16 > eta) sortedAdj[i].emplace_back(p, v);
        }
        d = sortedAdj[i].size();
        if (d < mink) continue;
        sort(sortedAdj[i].begin(), sortedAdj[i].end(), greater< pair<double,int> >());
		double p = 1.0, q = 1.0;
		for (int j = 0; j < d; ++j) {
			q = p * sortedAdj[i][j].first;
			if (q + 1e-16 < eta){
				d = j; break;
			}
			p = q;
		}
		topcore[i] = d;
	}

	vector<vector<int>> vertices(V);
	vector<int> core_index(V);
	int max_init_topcore = 0, max_topcore = 0;

	for (int i = 0; i < V; ++i)
	{
		int corei = topcore[i];
        if (corei == 0) continue;
		core_index[i] = (int)vertices[corei].size();
		vertices[corei].push_back(i);
		max_init_topcore = max(max_init_topcore, corei);
	}
	//printf("Maximum init topcore number: %d\n", max_init_topcore);
	int cnt = 0;
	for (int i = 0; i <= max_init_topcore; ++i) {
		size_t size = vertices[i].size();
		for (int j = 0; j < size; ++j) {
			int v = vertices[i][j];
			int d = sortedAdj[v].size();
            if (nodeset != NULL) nodeset[cnt++] = v;
			visited[v] = true;
			for (int k = 0; k < d; ++k) {
				int u = sortedAdj[v][k].second;
				double p = sortedAdj[v][k].first;

				int core_u = topcore[u];
				if (core_u > i) {
					double q = 1.0, pro = 1.0;
					int du = sortedAdj[u].size();
					int core_new = 0;
					for (int l = 0; l < du; ++l){
						if (!visited[sortedAdj[u][l].second]) {
							pro *= sortedAdj[u][l].first;
							if (pro + 1e-16 < eta) break;
							++core_new;
						}
					}
					core_new = max(i, core_new);
					assert(core_new <= core_u);
					topcore[u] = core_new;

					if (core_new != core_u) {
						int pos = core_index[u];
						int pos_end = (int)vertices[core_u].size() - 1;
						int w = vertices[core_u][pos_end];

						core_index[u] = (int)vertices[core_new].size();
						vertices[core_new].push_back(u);
						if (w != u) {
							vertices[core_u][pos] = w;
							core_index[w] = pos;
						}
						vertices[core_u].resize(pos_end);
						size = vertices[i].size();
					}
				}
			}
		}
		if (size > 0)
		{
			max_topcore = i;
			//printf("Iterator %d, letf_vert = %d\n", i, n - cnt);
		}
	}

	//printf("Maximum max_topcore: %d\n", max_topcore);
    return cnt;
}

int Algorithm::topKEtatriangle(int *nodeset)
{
    clock_t stm = clock();
    int nodesize = topKEtaCore(nodeset);
    if (nodesize < V) resetGraph(nodeset, nodesize);
    int k = mincliquesize - 2;
    if (k <= 0) return nodesize;
    int re = nodeset[0];

    vector<int> inset(V, -1), subdeg(V, 0);
    vector<epro>comnbrs(maxdeg);
    vector<vector<int>>readj(V);
    vector<vector<epro>> topTri(E);
    vector<vector<epro>*> topadj(V);
    
    for (int i = 0, cnt = 0; i < nodesize; ++i) {
        int u = nodeset[i];
        readj[u].reserve(deg[u]);
        topadj[u] = &topTri[cnt];
        cnt += deg[u];
        re = u;
    }
    for (int i = 0; i < nodesize; ++i) {
        int u = nodeset[i];
        int d = deg[u];
        for (int j = 0; j < d; ++j) {
            int v = adj[u][j].first;
            if (u > v) continue;
            readj[v].emplace_back(subdeg[u]++);
            readj[u].emplace_back(subdeg[v]++);
        }
    }
    double tempsize = 0;
    queue< pair<int,int> > Q;
    //printf("k=%d\n", k);
    for (int i = 0; i < nodesize; ++i) {
        int u = nodeset[i];
        int d = deg[u];
        int maxu = adj[u][d-1].first;
        vector<epro> *topadju;

        for (int j = 0; j < d; ++j)
            inset[adj[u][j].first] = j;

        for (int j = 0; j < d; ++j) {
            int v = adj[u][j].first;
            if (u > v) continue;
            int dv = deg[v];
            int comsize = 0;

            if (subdeg[u] < k+1 || subdeg[v] < k+1) {
                int rpos = readj[u][j];
                adj[u][j].second = -adj[u][j].second;
                adj[v][rpos].second = -adj[v][rpos].second;
                Q.emplace(u, j);
                subdeg[u]--; subdeg[v]--;
                continue;
            }

            for (int s = 0; s < dv; ++s) {
                int w = adj[v][s].first;
                int is = inset[w];
                if (w > maxu) break;
                else if (is >= 0) {
                    double p = adj[v][s].second * adj[u][is].second;
                    if (p + 1e-16 > eta) {
                        comnbrs[comsize].p = p;
                        comnbrs[comsize++].u = w;
                    }
                }
            }
            if (comsize < k) {
                int rpos = readj[u][j];
                adj[u][j].second = -adj[u][j].second;
                adj[v][rpos].second = -adj[v][rpos].second;
                Q.emplace(u, j); 
                subdeg[u]--; subdeg[v]--;
                continue;
            }
            sort(comnbrs.begin(), comnbrs.begin() + comsize, decreasefun);
            topadj[u][j].resize(comsize+1);
            topadju = &topadj[u][j];
            tempsize += (sizeof(pair<double, int>)*(comsize+1))/double(1024*1024);
            double p = adj[u][j].second;
            int len = 0;
            for (int l = 0; l < comsize; ++l) {
                if ( l < k ) {
                    p *= comnbrs[l].p;
                    if (p + 1e-16 < eta) break;
                }
                else {
                    p /= comnbrs[l-1].p;
                    p *= comnbrs[l].p;
                    if (p + 1e-15 < eta) break;    
                }
                (*topadju)[l+1] = comnbrs[l];
                len ++;
            }
            (*topadju)[0].u = len;
            if (len < k) {
                int rpos = readj[u][j];
                adj[u][j].second = -adj[u][j].second;
                adj[v][rpos].second = -adj[v][rpos].second;
                Q.emplace(u, j);
                subdeg[u]--; subdeg[v]--;
            }
        }

        for (int j = 0; j < d; ++j)
            inset[adj[u][j].first] = -1;
    }

    //printf("Q.size=%ld\n", Q.size());
    size_t maxQ = Q.size();
    while (!Q.empty()) {
        int u = Q.front().first, j = Q.front().second; Q.pop();
        int v = adj[u][j].first;
        int ressize = 0;
        int us = 0, ut = deg[u];
        int vs = 0, vt = deg[v];
        pairs *adju = adj[u], *adjv = adj[v];
        double puv = adj[u][j].second;

        adj[u][j].second = 0;
        adj[v][readj[u][j]].second = 0;

        while (us < ut && vs < vt) {
            int a = adju[us].first;
            int b = adjv[vs].first;
            if (a < b) us++;
            else if (a > b) vs++;
            else {
                double pu = adju[us].second;
                double pv = adjv[vs].second;
                double p = pu * pv;
                if (p + 1e-16 > eta || p -1e-16 < -eta) {
                    bool is_break = false;
                    int comsize = 0;
                    int x, y, rpos;
                    if (pu + 1e-16 > eta) {
                        rpos = readj[u][us];
                        if (a > u) {x = u; y = us;}
                        else {x = a; y = rpos;}
                        vector<epro> &nbrsvu = topadj[x][y];
                        int dm =  nbrsvu[0].u;
                        if (subdeg[u] >= k+1 && subdeg[a] >= k+1)
                        for (int m = 1; m <= dm; ++m) {
                            double q = nbrsvu[m].p;
                            if (nbrsvu[m].u == v) { nbrsvu[m].p = 0; is_break = true;}
                            else if (q + 1e-16 > eta) {
                                if (comsize ++ < k)
                                    pu *= q;
                                else if (is_break) break;
                            }
                        }
                        if (pu + 1e-16 < eta || comsize < k) {
                            adj[u][us].second = -adj[u][us].second;
                            adj[a][rpos].second = -adj[a][rpos].second;
                            Q.emplace(x, y);
                            subdeg[u]--; subdeg[a]--;
                        }
                    }

                    if (pv + 1e-16 > eta) {
                        rpos = readj[v][vs];
                        if (a > v) {x = v; y = vs;}
                        else {x = a; y = rpos;}
                        vector<epro> &nbrsvu = topadj[x][y];
                        int dm =  nbrsvu[0].u;
                        comsize = 0; is_break = false;
                        if (subdeg[v] >= k+1 && subdeg[a] >= k+1)
                        for (int m = 1; m <= dm; ++m) {
                            double q = nbrsvu[m].p;
                            if (nbrsvu[m].u == u) {nbrsvu[m].p = 0; is_break = true;}
                            else if (q + 1e-16 > eta) {
                                if (comsize ++ < k) 
                                    pv *= q;
                                else if (is_break) break;
                            }
                        }

                        if (pv + 1e-16 < eta || comsize < k) {
                            adj[v][vs].second = -adj[v][vs].second;
                            adj[a][rpos].second = -adj[a][rpos].second;
                            Q.emplace(x, y);
                            subdeg[v]--; subdeg[a]--; 
                        }
                    }
                }
                us++; vs++;
            }
        }
        adj[u][j].second = 0;
        adj[v][readj[u][j]].second = 0;
        maxQ = max(maxQ,Q.size());
    }

    int resize = 0;
    for (int i = 0; i < nodesize; ++i) {
        int u = nodeset[i];
        if (subdeg[u] >= k+1) nodeset[resize++] = u;
    }
    printf("Phase 2, left nums: %d\n", resize);

    // for (int i = 0; i < V; ++i) {
    //     int d = deg[i];
    //     if (d > 0) {
    //         for (int j = 0; j < d; ++j) {
    //             int v = adj[i][j].first;
    //             double p = adj[i][j].second;
                
    //             int r = readj[i][j];
    //             int rv = adj[v][r].first;
    //             double rp = adj[v][r].second;
    //             assert(i == rv);
    //             assert(p == rp);
    //         }
    //     }
    // }
    
    return resize;
}

int Algorithm::coloring(int *nodeset, int nodesize)
{
    int color_nums = 0, len = nodesize;
    vector<int> bin(maxdeg+1), sequence(len);
    len = topKEtaCoreDecompsition(nodeset); 

    for (int i = 0; i < len; ++i) bin[deg[nodeset[i]]]++;
    for (int i = 1; i <= maxdeg; ++i) bin[i] += bin[i-1];
    for (int i = maxdeg; i > 0; --i) bin[i] = bin[i-1]; bin[0] = 0;
    for (int i = 0; i < len; ++i) {
        int v = nodeset[i];
        int posi = bin[deg[v]]++;
        sequence[posi] = v;
    }
    if (colors.empty()) {colors.resize(V, -1);}
    for (int i = len-1; i >= 0; --i) {
        int maxc = -1, curc = -1;
        int u = sequence[i];
        int d = deg[u];
        //deg[u] = 0;
        for (int j = 0; j < d; ++j) {
            int v = adj[u][j].first;
            int cv = colors[v];
            //double p = adj[u][j].first;
            //if (p + 1e-16 > eta) adj[u][deg[u]++] = adj[u][j];
            if (cv >= 0) {bin[cv]++; maxc = max(maxc, cv);}
        }
        for (int j = 0; j <= maxc; ++j) {
            if (bin[j] == 0 && curc == -1) curc = j;
            bin[j] = 0;
        }
        if (curc == -1) {
            color_nums = max(color_nums, ++maxc + 1);
            colors[u] = maxc;
        }
        else colors[u] = curc;
    }
    //printf("Color nums: %d\n", color_nums);
    return color_nums;
}

void Algorithm::resetGraph(int *nodeset, int nodesize)
{
    vector<bool> visited(V);
    maxdeg = 0; E = 0;
    for (int i = 0; i < nodesize; ++i)  visited[nodeset[i]] = true;
    for (int i = 0; i < V; ++i) {
        if (!visited[i]) { deg[i] = 0; continue; }
        int d = deg[i];
        int dsize = 0;
        for (int j = 0; j < d; ++j) {
            int v = adj[i][j].first;
            double p = adj[i][j].second;
            if (visited[v] && p + 1e-16 > eta) {
            //if (visited[v]) {
                adj[i][dsize].first = v;
                adj[i][dsize++].second = p;
            }
        }
        deg[i] = dsize; maxdeg = max(maxdeg, dsize); E += dsize;
    }
}