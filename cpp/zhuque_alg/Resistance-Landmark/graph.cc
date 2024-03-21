#include <fstream>
#include <sstream>
#include <iostream>
#include <queue>
#include "graph.h"

int get_id(int index, map<int, int> &dic)
{
  auto it = dic.find(index);
  if (it != dic.end())
  {
    return it->second;
  }
  else
  {
    int id = dic.size();
    dic[index] = id;
    return id;
  }
}

void read_graph(string fn, graph_t &G)
{
  ifstream ifs(fn.c_str());
  string line;
  bool fst = true;
  map<int, int> dic;

  while (getline(ifs, line))
  {
    if (line.size() > 0 && line[0] == '%')
      continue;
    istringstream iss(line);
    if (fst)
    {
      // skip this line
      fst = false;
    }
    else
    {
      int u, v;
      iss >> u >> v;
      u = get_id(u, dic);
      v = get_id(v, dic);
      //      --u, --v;
      if (G.size() <= u)
        G.resize(u + 1);
      if (G.size() <= v)
        G.resize(v + 1);
      G[u].push_back(v);
      G[v].push_back(u);
    }
  }
}

int read_graph_degree(string fn, graph_t &G)
{
  ifstream ifs(fn.c_str());
  string line;
  bool fst = true;
  map<int, int> dic;

  while (getline(ifs, line))
  {
    if (line.size() > 0 && line[0] == '%')
      continue;
    istringstream iss(line);
    if (fst)
    {
      // skip this line
      fst = false;
    }
    else
    {
      int u, v;
      iss >> u >> v;
      u = get_id(u, dic);
      v = get_id(v, dic);
      //      --u, --v;
      if (G.size() <= u)
        G.resize(u + 1);
      if (G.size() <= v)
        G.resize(v + 1);
      G[u].push_back(v);
      G[v].push_back(u);
    }
  }

  int n = G.size();
  int* degree = new int[n]();
  for(int i=0; i<n; i++) {
    degree[i] = G[i].size();
  }
  int maxdegree = 0, maxid = 0;
  for(int i=0; i<n; i++) {
    if(degree[i]>maxdegree) {
      maxdegree = degree[i];
      maxid = i;
    }
  }
  cerr << "maxid: " << maxid << " maxdegree: " << maxdegree << endl;
  delete[] degree;
  return maxid;
}

size_t num_of_edges(graph_t &G)
{
  int res = 0;
  for (size_t u = 0; u < G.size(); u++)
  {
    res += G[u].size();
  }
  return res / 2;
}

void compute_BFS_tree(graph_t &G, vector<int> &bfs_parent, vector<int> &bfs_order, int landmark) {
    // compute a BFS tree from the highest degree node
    int n = G.size();
    bfs_parent.resize(n);
    bfs_order.resize(n);
    bfs_parent[landmark] = -1;
    bfs_order[0] = landmark;

    int order_count = 0;

    vector<bool> is_visited(n);
    std::fill(is_visited.begin(), is_visited.end(), false);

    queue<int> q;
    q.push(landmark);
    is_visited[landmark] = true;

    do {
        int currentNode = q.front();
        // cerr << "currentNode: " << currentNode << endl;
        q.pop();
        for(int i=0; i<G[currentNode].size(); i++) {
            int u = G[currentNode][i];
            if(!is_visited[u]) {
              is_visited[u] = true;
              q.push(u);
              bfs_parent[u] = currentNode;
              bfs_order[++order_count] = u;
            }
        }
    } while (!q.empty());

    cerr << "order count: " << order_count << endl;
    assert(order_count == n-1);
}



void contract(graph_t &G, int u, int v)
{
  if (u == v)
    return;

  // delete v from u's neighbor set
  for (size_t i = 0; i < G[u].size(); i++)
  {
    if (G[u][i] == v)
    {
      swap(G[u][i], G[u][G[u].size() - 1]);
      G[u].pop_back();
      --i;
    }
  }
  sort(G[u].begin(), G[u].end());

  for (size_t i = 0; i < G[v].size(); i++)
  {
    int w = G[v][i];
    for (size_t j = 0; j < G[w].size(); j++)
    {
      if (G[w][j] == v)
        G[w][j] = u;
    }
    sort(G[w].begin(), G[w].end());
  }
  //  swap(G[v], G[G.size() - 1]);
  //  G.pop_back();
}

SparseMatrix<double> Laplacian(graph_t &G)
{
  int n = G.size();
  vector<Eigen::Triplet<double>> triplets;
  for (int u = 0; u < n; ++u)
  {
    triplets.push_back(Eigen::Triplet<double>(u, u, G[u].size()));
    for (int i = 0; i < G[u].size(); ++i)
    {
      triplets.push_back(Eigen::Triplet<double>(u, G[u][i], -1));
    }
  }
  SparseMatrix<double> L(n, n);
  L.setFromTriplets(triplets.begin(), triplets.end());
  return L;
}

int random_walk(graph_t &G, int s, int l)
{
  int v = s;
  for (int i = 0; i < l; ++i)
  {
    if (G[v].size() == 0)
      continue;
    v = G[v][random() % G[v].size()];
  }
  return v;
}

void sample_spanning_tree(graph_t &G, vector<int> &order, vector<int> &next_edges)
{
  int n = G.size();
  vector<bool> visit(n);
  visit[order[0]] = true;

  for (int s = 1; s < n; ++s)
  {
    if (!visit[order[s]])
    {
      int u = order[s];
      while (!visit[u])
      {
        int next = rand() % G[u].size();
        next_edges[u] = next;
        u = G[u][next];
      }

      u = order[s];
      while (!visit[u])
      {
        visit[u] = true;
        u = G[u][next_edges[u]];
      }
    }
  }
}

int pagerank_landmark(graph_t& G) {
  int n = G.size();

  vector<double> pr(n, 0);
  int T = 1000000000;
  double alpha = 0.2;

  for(int i=0; i<T; i++) {
    int u = random() % n; // uniformly sample a node
    while(true) {
      double r = (double) rand() / (RAND_MAX + 1.0);
      // cerr << r << endl;
      if(r<alpha) {
        pr[u] += 1./T;
        break;
      }
      u = G[u][random() % G[u].size()];
    }
  }

  for(int i=0; i<20; i++) {
    cerr << pr[i] << endl;
  }

  double pr_max = 0.;
  int landmark = 0;
  for(int i=0; i<n; i++) {
    if(pr[i]>pr_max) {
      pr_max = pr[i];
      landmark = i;
    }
  }
  cerr << "pagerank landmark node: " << landmark << endl;
  cerr << "pagerank value: " << pr_max << endl;
  cerr << "degree: " << G[landmark].size() << endl;
  return landmark;
}

int core_landmark(graph_t &G) {
    int n = G.size();
    cout << "Core Decomposition" << endl;
		int* deg = new int[n];
		for(int i=0; i<n; i++) deg[i] = G[i].size();

		int *rid = new int[n];
		int *id = new int[n]();
		for(int i=0; i<n; i++) ++id[deg[i]];
		for(int i=1; i<n; i++) id[i] += id[i-1];

		for(int i=0; i<n; i++) rid[i] = --id[deg[i]];
		for(int i=0; i<n; i++) id[rid[i]] = i;

		int *degree_start = new int[n+1];
		for(int i=0,j=0; i<=n; i++) {
			while(j<n&&deg[id[j]]<i) ++j;
			degree_start[i] = j;
		}

		int* core = new int[n];
		int max_core = 0;
    int landmark = 0;
		for(int i=0; i<n; i++) {
			int u = id[i];
			if(deg[u]>max_core) {
        max_core = deg[u];
        landmark = u;
      }
			core[u] = max_core;

			++degree_start[deg[u]];
			if(deg[u]==0) continue;

			degree_start[deg[u]-1] = degree_start[deg[u]];
			for(uint j=0; j<deg[u]; j++) if(rid[G[u][j]]>i) {
				uint v = G[u][j];
				uint pos1 = degree_start[deg[v]], pos2 = rid[v];
				swap(id[pos1], id[pos2]);
				rid[id[pos1]] = pos1; rid[id[pos2]] = pos2;
				++degree_start[deg[v]];
				--deg[v];
			}
		}
		cerr << "Max core: " << max_core << endl;

    cerr << "core landmark node: " << landmark << endl;
    cerr << "core value: " << max_core << endl;
    cerr << "degree: " << G[landmark].size() << endl;

    return landmark;
}

int ecc_landmark(graph_t& G, int& landmark) {
  int n = G.size();
  int ecc = n;
  vector<int> eccLowerBound(n);
  vector<int> distance(n);
  for(int i=0; i<100; i++) {
    landmark = computeBFS(G, landmark, ecc, eccLowerBound, distance);
  }

  int minEccNode = 0;
  int minEcc = n;
  for(int i=0; i<n; i++) {
    if(eccLowerBound[i]<minEcc) {
      minEcc = eccLowerBound[i];
      minEccNode = i;
    }
  }

  // for(int i=0; i<20; i++) cerr << eccLowerBound[i] << endl;

  cerr << "Min Ecc: " << minEcc << endl;

  cerr << "ecc landmark node: " << landmark << endl;
  cerr << "min ecc value: " << minEcc << endl;
  cerr << "degree: " << G[landmark].size() << endl;

  return landmark;
}

int computeBFS(graph_t& G, int landmark, int& ecc, vector<int>& eccLowerBound, vector<int>& distance) {
  int n = G.size();
  int ecc_node = 0;

  queue<int> q;
  q.push(landmark);
  vector<bool> visited(n, false);
  distance[landmark] = 0;

  do {
    int u = q.front();
    q.pop();
    eccLowerBound[u] = std::max(eccLowerBound[u], distance[u]);
    ecc_node = u;

    for(int i=0; i<G[u].size(); i++) {
      int v = G[u][i];
      if(!visited[v]) {
        q.push(v);
        visited[v] = true;
        distance[v] = distance[u] + 1;
      }
    }
  } while(!q.empty());

  return ecc_node;
}