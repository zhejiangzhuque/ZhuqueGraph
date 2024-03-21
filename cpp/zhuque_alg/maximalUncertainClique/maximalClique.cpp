#include "maximalClique.h"

maximalClique::maximalClique()
{
	iterations = 0;
	outs = NULL;
}

maximalClique::~maximalClique()
{
	if (outs) fclose(outs);
}

void maximalClique::setSubgraph(node *Can, list &X, vector<bool> &subset)
{
	int maxsubdeg = 0;
	for (int i = Can[0].next; i > 0; i = Can[i].next) subset[Can[i].u] = true;
	for (int i = Can[0].next; i > 0; i = Can[i].next) {
		int u = Can[i].u;
		int d = deg[u];
		subgraphdeg[u] = 0;
		subgraph[u].clear();
		for (int j = 0; j < d; ++j) {
			int v = adj[u][j].first;
			if (subset[v])
				subgraph[u].emplace_back(v, adj[u][j].second);
		}
		subgraphdeg[u] = subgraph[u].size();
		maxsubdeg = max(maxsubdeg, subgraphdeg[u]);
	}
	// for (int i = X.head; i >= 0 && i < X.size; i++) {
	// 	int u = X.Nodes[i].u;
	// 	int d = deg[u];
	// 	subgraphdeg[u] = 0;
	// 	subgraph[u].clear();
	// 	for (int j = 0; j < d; ++j) {
	// 		int v = adj[u][j].first;
	// 		if (subset[v])
	// 			subgraph[u].emplace_back(v, adj[u][j].second);
	// 	}
	// 	subgraphdeg[u] = subgraph[u].size();
	// 	maxsubdeg = max(maxsubdeg, subgraphdeg[u]);
	// }

	for (int i = Can[0].next; i > 0; i = Can[i].next) subset[Can[i].u] = false;

	if (subgraphcolor.empty()) subgraphcolor.resize(V, -1);
	
	vector<int> &bin = counts;
	vector<int> &sequence = clset;

	for (int i = Can[0].next; i > 0; i = Can[i].next) bin[subgraphdeg[Can[i].u]]++;
    for (int i = 1; i <= maxsubdeg; ++i) bin[i] += bin[i-1];
    for (int i = maxsubdeg; i > 0; --i) bin[i] = bin[i-1]; bin[0] = 0;
    for (int i = Can[0].next; i > 0; i = Can[i].next) {
        int v = Can[i].u;
        int posi = bin[subgraphdeg[v]]++;
		sequence[posi] = v;
		subgraphcolor[v] = -1;
    }

	int color_nums = 0, len = bin[maxsubdeg];
	for (int i = len-1; i >= 0; --i) {
        int maxc = -1, curc = -1;
        int u = sequence[i];
        int d = subgraphdeg[u];
        for (int j = 0; j < d; ++j) {
            int v = subgraph[u][j].first;
            int cv = subgraphcolor[v];
            if (cv >= 0) {bin[cv]++; maxc = max(maxc, cv);}
        }
        for (int j = 0; j <= maxc; ++j) {
            if (bin[j] == 0 && curc == -1) curc = j;
            bin[j] = 0;
        }
        if (curc == -1) {
            color_nums = max(color_nums, ++maxc + 1);
            subgraphcolor[u] = maxc;
        }
        else subgraphcolor[u] = curc;
		bin[i] = 0;
    }
}

int maximalClique::generateCan(node *Can, int i, int csize, double qn, int u, node *CanN)
{
	//int ns = 0, nt = deg[u];
	int ns = 0, nt = subgraphdeg[u];
	int u1, u2, cs = i, cnsize = 0;
	double pn;
    //pairs *adju = adj[u];
	vector<pairs> &adju = subgraph[u];
	CanN[0].pre = CanN[0].next = -1;
	while (cs > 0 && cs <= csize && ns < nt)
	{
		u1 = Can[cs].u;
		u2 = adju[ns].first;
		if (u1 < u2) cs = Can[cs].next;
		else if (u1 > u2) ns++;
		else {
			pn = Can[cs].p * adju[ns].second;
			cs = Can[cs].next; ns++;
			if (pn * qn + 1e-15 > eta) {
				node_push_back(CanN, ++cnsize, u1, pn);
			}
		}
	}
	return cnsize; 
}

int maximalClique::generateCan(list &Can, double qn, int u, list &CanN)
{
	int u1, u2, nt, ns = 0;
	int cs = Can.head;
	double pn;
	list_clear(CanN);

	pairs *adju = adj[u];
	nt = deg[u];
	while (cs >= 0 && cs < Can.size && ns < nt)
	{
		u1 = Can.Nodes[cs].u;
		u2 = adju[ns].first;
		if (u1 < u2) cs = Can.Nodes[cs].next;
		else if (u1 > u2) ns++;
		else {
			pn = Can.Nodes[cs].p * adju[ns].second;
			cs = Can.Nodes[cs].next; ns++;
			if (pn * qn + 1e-15 > eta)
				list_push_back(CanN, u1, pn);
		}
	}
	return CanN.size;
}


bool maximalClique::CheckColors(node *Can, int csize, int curRsize, int curMaxsize)
{
	int left = mincliquesize - curRsize - 1;
	if (left <= curMaxsize) return 0;
	int v, c, sscnums = 0, numbers = 0;
	for (int j = Can[0].next; j <= csize && j > 0; j = Can[j].next) {
		v = Can[j].u;
		c = subgraphcolor[v];
		if (counts[c] == 0) {
			clset[sscnums++] = c;
			if (sscnums > left) counts[c] = 2;
			else counts[c] = 1;
		}
		if (counts[c] == 1) numbers++;
	}
	if (sscnums <= left) {
		for (int j = Can[0].next; j <= csize && j > 0; j = Can[j].next) {
			v = Can[j].u;
			c = subgraphcolor[v];
			subset[v] = true;
			counts[c] = 0;
		}
		return 1;
	}
	if (numbers >= curMaxsize && curRsize <= 4) {
		sort(clset.begin(), clset.begin()+sscnums, greater<int>());
		//sort(clset.begin(), clset.begin()+sscnums);
		for (int j = 0; j < left; ++j) counts[clset[j]] = 1;
		for (int j = left; j < sscnums; ++j) counts[clset[j]] = 2;
		for (int j = Can[0].next; j <= csize && j > 0; j = Can[j].next) {
			v = Can[j].u;
			c = subgraphcolor[v];
			if (counts[c] == 1) subset[v] = true;
			counts[c] = 0;
		}
		return 1;
	}
	else {
		for (int j = 0; j < sscnums; ++j)
			counts[clset[j]] = 0;
		return 0;
	}
}

void maximalClique::hybrid_pivot_enumerate(int *RP, int *R, double q, int rsize, node *Can, int csize, list &X, int xsize)
{
	assert(RP[0] == 1);
	iterations++;
	if (csize + rsize < mincliquesize) {
		if (RP[0] <= csize) {
			for (int j = 1; j <= csize; ++j) RP[j+1] = Can[j].u;
			RP[0] = 1 + csize;
		}
		return;
	}
	if (csize == 0 && xsize == 0) {
		maximal_cliques++; maximum_clique_size = max(maximum_clique_size, rsize);
		// 	printf("clique %d : ", maximal_cliques);
		// 	for (int i = 0; i < rsize; ++i)
		// 		printf(" %d", R[i]);
		// 	printf(" [%lf]\n",q);
		return;
	}
	if (csize == 0) return;
	else if (csize==1) {
		int u = Can[1].u;
		double qr = Can[1].p;
		if (RP[0] == 1) RP[++RP[0]] = u;
		if (xsize == 0) {
			R[rsize++] = u; 
			q *= qr;
			LB[u] = max(LB[u], rsize+1);
			maximal_cliques++; maximum_clique_size = max(maximum_clique_size, rsize);
		}
		else {
			list Xn; list_init_size(Xn, csize+xsize);
			int xnsize = generateCan(X, q * qr, u, Xn);
			if (xnsize == 0) {
				LB[u] = max(LB[u], rsize+1);
				maximal_cliques++; maximum_clique_size = max(maximum_clique_size, rsize);
			}
			list_free(Xn);
		}
		return;
	}

	int pivot = 1, maxd = 0, c = -1, maxlb = 0;
	int u, xnsize, cnsize = 0;
	
	int pivot1 = pivot, pivot2 = pivot;
	for (int j = 1; j <= csize; ++j) {
		int v = Can[j].u;
		int d = deg[v];
		int vc = colors[v];
		int lb = LB[v];
		
		//pivot selection
		if (lb > maxlb && vc >= c) {
			pivot1 = j; c = vc; maxlb = lb;
		}
		else if (vc > c && d >= maxd) {
			pivot2 = j; c = vc; maxd = d;
		}
	}
	pivot = maxlb >= mincliquesize ? pivot1 : pivot2;

	double r, qn;
	int idu = X.head;
	int *RPN = new int[csize+1];
	int qsize = 0, *Q = new int[csize];
	node *CanN = new node[csize];
	list Xn; list_init_size(Xn, csize+xsize);
	Q[qsize++] = pivot; maxlb = 0;
	for (int i = 0; i < qsize; ++i) {
		if (csize - i + rsize < mincliquesize) break;
		int id = Q[i];
		u = Can[id].u;
		r = Can[id].p;

		qn = q * r;
		R[rsize] = u;

		int pre = Can[id].pre;
		int next = Can[id].next;
		Can[pre].next = next;
		if (next > 0) Can[next].pre = pre;

		cnsize = generateCan(Can, Can[0].next, csize, qn, u, CanN);
		xnsize = generateCan(X, qn, u, Xn);

		RPN[0]=1; RPN[1] = u;
		LB[u] = max(LB[u], rsize+1);

		hybrid_pivot_enumerate(RPN, R, qn, rsize+1, CanN, cnsize, Xn, xnsize);

		if (RPN[0] >= RP[0]) {
			RP[0] = 1 + RPN[0];
			for (int j = 1; j <= RPN[0]; ++j)
				RP[j+1] = RPN[j];
		}
		
		if (i == 0) {
			if (!CheckColors(Can, csize, rsize, RPN[0]))
			for (int j = 2; j <= RPN[0]; ++j)
				subset[RPN[j]] = true;
			for (int j = Can[0].next; j > 0; j = Can[j].next) {
				int v = Can[j].u;
				if (subset[v]) subset[v] = false;
				else Q[qsize++] = j;
			}
			maxlb = csize - qsize + 1;
			if (qsize>1) {
				idu = X.head;
				while(idu >= 0 && X.Nodes[idu].u < u)
					idu = X.Nodes[idu].next;
				list_insert(X, idu, u, r);
				idu = X.head;
			}
		}
		else if (maxlb < RPN[0]) {
			maxlb = RPN[0];
			for (int j = 2; j <= RPN[0]; ++j)
				subset[RPN[j]] = true;
			qsize = i+1;
			for (int j = Can[0].next; j > 0; j = Can[j].next) {
				int v = Can[j].u;
				if (subset[v]) subset[v] = false;
				else Q[qsize++] = j;
			}
			if (qsize>i+1) {
				while(idu >= 0 && X.Nodes[idu].u < u)
					idu = X.Nodes[idu].next;
				list_insert(X, idu, u, r);
				idu = X.head;
			}
		}
		else if (next > 0) {
			while(idu >= 0 && X.Nodes[idu].u < u)
				idu = X.Nodes[idu].next;
			list_insert(X, idu, u, r);
		}
		//list_push_back(X,u,r);

		// if (next > 0) {
		// 	while(idu >= 0 && X.Nodes[idu].u < u)
		// 		idu = X.Nodes[idu].next;
		// 	list_insert(X, idu, u, r);
		// }
	}

	delete[] RPN;
	delete[] CanN;
	delete[] Q;
	list_free(Xn);    
}


void maximalClique::run_pivot_enumerate()
{
	//pruning
	clock_t stm = clock();
	int nodesize = V, *nodeset = new int[V];
	if (algorithm == 1) nodesize = topKEtaCore(nodeset);
	else if (algorithm == 2) nodesize = topKEtatriangle(nodeset);
	else {printf("No graph reduction techniuqes!\n");}
	if (nodesize < V) resetGraph(nodeset, nodesize);
	printf("Pruning time: %lf s\n", double(clock()- stm) / CLOCKS_PER_SEC);
	stm = clock();
	int maxclmums = coloring(nodeset, nodesize);

	vector<bool> visited(V, false);
	initSubgrah();

	int R[maxclmums], RP[maxclmums];
	int csize = 0, xsize = 0;
	list X; list_init_size(X, maxdeg+1);
	node *Can = new node[maxdeg+1];
	
	maximal_cliques = 0; maximum_clique_size = 0;
	for (int i = 0; i < nodesize; ++i) {
		int d, u = nodeset[i];
		if (colors[u]+1 < mincliquesize) continue;
		pairs *adji = adj[u];
		visited[u] = true;
		list_clear(X);
		csize = 0; xsize = 0;
		d = deg[u];
		Can[0].next = Can[0].pre = -1;
		for (int j = 0; j < d; ++j) {
			int v = adji[j].first;
			double p = adji[j].second;
			if (p + 1e-15 > eta) {
				if (visited[v]) {
					list_push_back(X, v, p);
					xsize++;
				}
				else node_push_back(Can,++csize,v,p);
			}
		}
		setSubgraph(Can, X, subset);
		R[0] = u; RP[0] = 1; RP[1] = u; 
		if (csize > 0 ) hybrid_pivot_enumerate(RP, R, 1.0f, 1, Can, csize, X, xsize);
	}
	delete[] Can;
	delete[] nodeset;
	list_free(X);

    printf("Maximal cliques: %lld\n", maximal_cliques);
	printf("Maximum clique size: %d\n", maximum_clique_size);
	printf("Running time: %lf s\n", double(clock()- stm) / CLOCKS_PER_SEC);
}

void maximalClique::run()
{
	LB.resize(V, false);
	subset.resize(V, false);
	counts.resize(maxdeg+1, 0);
	clset.resize(maxdeg+1, 0);
	run_pivot_enumerate();
	printf("Iterations: %lld\n", iterations);
}