#include <iostream>
#include <queue>
#include <stack>
#include <sys/time.h>
#include <fstream>
#include "graph.h"
#include "methods.h"
#include <iomanip> 

double get_current_time_sec_method()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

double effective_resistance_akp(graph_t &G, int s, int t, double eps, double lambda)
{
  // int l = ceil(log(4 / (eps * (1 - lambda))) / log(1.0 / lambda) / 2);
  // int r = (int)ceil(40 * l * l * log(80 * l) / (eps * eps));
  // cerr << "l: " << l << " r " << r << endl;
  double l = eps;
  double r = lambda;
  cerr << "l: " << l << " r " << r << endl;
  double delta = 0;

  for (size_t i = 0; i < l; i++)
  {
    int Xis = 0, Xit = 0, Yis = 0, Yit = 0;
    for (size_t j = 0; j < r; j++)
    {
      int v = random_walk(G, s, i);
      if (v == s)
        Xis++;
      if (v == t)
        Xit++;
    }
    for (size_t j = 0; j < r; j++)
    {
      int v = random_walk(G, t, i);
      if (v == s)
        Yis++;
      if (v == t)
        Yit++;
    }
    double deltai = (double)Xis / G[s].size() - (double)Xit / G[t].size() - (double)Yis / G[s].size() + (double)Yit / G[t].size();
    deltai /= r;
    delta += deltai;
  }
  return delta;
}

double approximate_num_spanning_trees(graph_t &G, double eps, double delta)
{
  int n = G.size(), m = num_of_edges(G);
  int r = 5; //(int)ceil(90 * 90 * 90 * exp(-3));
  vector<double> S(2 * r + 1);
  for (size_t t = 1; t <= 2 * r; t++)
  {
    S[t] = S[t - 1] + (1.0 / t);
  }
  double s = S[2 * r];
  int y = 0;
  int N = 8 * log(4 / delta) * s * s / (eps * eps);
  for (size_t i = 0; i < N; i++)
  {
    int u = random() % G.size();

    auto it = lower_bound(S.begin(), S.end(), (random() % 10000) / 10000.0 * S[2 * r]);
    int t = it - S.begin();
    int v = random_walk(G, u, t);
    if (v == u)
      ++y;
  }

  int N1 = 256 * log(1 / delta) * log(n) * log(n) / (eps * eps);
  double w = 0;
  for (size_t i = 0; i < N1; i++)
  {
    int v = random() % G.size();
    w += log(2 * G[v].size());
  }
  w /= N1;

  return -log(4 * m) / n + w - s * y / N + s / n;
}

double effective_resistance_st(graph_t &G, int s, int t, double eps, double delta)
{
  int n = G.size();
  auto G_contracted = G;
  contract(G_contracted, s, t);
  double a = approximate_num_spanning_trees(G_contracted, eps / 2, delta / 2);
  double b = approximate_num_spanning_trees(G, eps / 2, delta / 2);
  return exp(a * (n - 1) - b * n);
}

double effective_resistance_mc(graph_t &G, int s, int t, double eps, double delta)
{
  if (G[s].size() > G[t].size())
    swap(s, t);

  int X = 0, N0 = 3 * log(6) * 1 * G[s].size() / (eps * eps);
  cerr << "N0: " << N0 << endl;
  // long long cnt_sum = 0;
  for (size_t i = 0; i < N0; i++)
  {
    int v = s;
    bool visited = false;

    int cnt = 0;
    for (;;)
    {
      ++cnt;
      v = G[v][random() % G[v].size()];
      if (v == t)
        visited = true;
      if (v == s)
        break;
    }
    // cerr << cnt << endl;
    // cnt_sum += cnt;
    if (visited)
      ++X;
  }
  // cerr << cnt_sum/N0 << endl;
  return (double)N0 / (G[s].size() * X);
}

double effective_resistance_mc3(graph_t &G, int s, int t, int landmark, int N0)
{
  int X = 0;
  cerr << "N0: " << N0 << endl;
  long long cnt_sum = 0;
  int Xis = 0, Xit = 0, Yis = 0, Yit = 0;
  double delta = 0;
  for (size_t i = 0; i < N0; i++)
  {
    int v = s; // random walk from s
    for (;;)
    {
      if(v == s) Xis++;
      else if(v == t) Xit++;
      v = G[v][random() % G[v].size()];
      cnt_sum++;
      if (v == landmark)
        break;
    }
    v = t; // random walk from t
    for (;;)
    {
      if(v == s) Yis++;
      else if(v == t) Yit++;
      v = G[v][random() % G[v].size()];
      cnt_sum++;
      if (v == landmark)
        break;
    }
  }

  delta = (double)Xis / G[s].size() - (double)Xit / G[t].size() - (double)Yis / G[s].size() + (double)Yit / G[t].size();
  delta = delta/(double)N0;
  cnt_sum = cnt_sum/N0;
  cerr << "average walk length: " << cnt_sum << endl;

  return delta;
}

double effective_resistance_mc2(graph_t &G, int s, int t, double eps, double delta)
{
  if (G[s].size() > G[t].size())
    swap(s, t);

  int X = 0, M0 = log(1 / delta) * 3 / (eps * eps);
  for (size_t i = 0; i < M0; i++)
  {
    int v = s;
    for (;;)
    {
      int w = G[v][random() % G[v].size()];
      if (w == t)
      {
        if (v == s)
          ++X;
        break;
      }
      else
      {
        v = w;
      }
    }
  }
  return (double)X / M0;
}

double dot(vector<double> x, vector<double> y)
{
  assert(x.size() == y.size());
  double res = 0;
  for (int i = 0; i < x.size(); ++i)
    res += x[i] * y[i];
  return res;
}

double effective_resistance_ex2(graph_t &G, int s, int t, double eps, double lambda)
{
  int l = ceil(log(4 / (eps * (1 - lambda))) / log(1 / lambda) / 2);
  int r = ceil(l * log(80 * l) / eps);
  int n = G.size();
  vector<double> Xs(n), Xt(n), Ys(n), Yt(n);
  double delta = 0;

  for (int i = 0; i < l; ++i)
  {
    for (int j = 0; j < r; ++j)
    {
      int v = random_walk(G, s, (i + 1) / 2);
      Xs[v] += 1.0 / (r * sqrt(G[v].size()));
      v = random_walk(G, t, (i + 1) / 2);
      Xt[v] += 1.0 / (r * sqrt(G[v].size()));
      v = random_walk(G, s, i / 2);
      Ys[v] += 1.0 / (r * sqrt(G[v].size()));
      v = random_walk(G, t, i / 2);
      Yt[v] += 1.0 / (r * sqrt(G[v].size()));
    }
    delta += dot(Xs, Ys) - dot(Xs, Yt) - dot(Xt, Ys) + dot(Xt, Yt);
  }
  return delta;
}

void spanning_tree_centrality(graph_t &G, int num_samples, vector<vector<double>> &centralities)
{
  int n = G.size();
  // vector<int> order = VertexOrdering::sort(G);
  vector<pair<int, int>> degrees(n);
  for (int v = 0; v < n; ++v)
    degrees[v] = make_pair(G[v].size(), v);
  sort(degrees.rbegin(), degrees.rend());

  vector<int> order(n);
  for (int i = 0; i < n; ++i)
    order[i] = degrees[i].second;
  vector<int> next_edges(n, -1);

  for (int trial = 0; trial < num_samples; ++trial)
  {
    if (trial % 10 == 0)
      cerr << trial << endl;
    sample_spanning_tree(G, order, next_edges);
    for (int s = 1; s < n; ++s)
    {
      int v = order[s];
      centralities[v][next_edges[v]] += 1.0 / num_samples;
    }
  }
}

vector<double> effective_resistance_single_source_spanning_tree(graph_t &G, int s, int landmark, double eps, vector<int> &bfs_parent) {
    int n = G.size();
    vector<double> result(n);

    return result;
}

vector<double> single_source_power_method(graph_t &G, int s, int v, int L) {
    int n = G.size();
    vector<double> result(n, 0);
    vector<double> residual(n, 0);
    vector<double> new_residual(n, 0);
    residual[s] = 1.0;
    // int L = 2000;
    double r_sum = 1.0;
    for(int i=0; i<L; i++) {
      for(int id=0; id<n; id++) {
        if(id == v) continue;
        result[id] += residual[id];
        double increment = residual[id];
        // if(increment>0) cerr << id << " " << increment << endl;
        residual[id] = 0;
        if(increment>0) {
          for(int j=0; j<G[id].size(); j++) {
            if(G[id][j]==v) {
              // cerr << "v: " << increment << endl;
              r_sum -= increment/(double)G[id].size();
            } else {
              new_residual[G[id][j]] += increment/(double)G[id].size();
              // cerr << new_residual[G[id][j]] << endl;
            }
          }
        }
      }
      residual.swap(new_residual);
      cerr << "iter: " << i << " r_sum: " << r_sum << endl;
    }
    return result;
}

double effective_resistance_push(graph_t &G, int s, int t, int landmark, double eps) {
    double delta;
    int n = G.size();
    vector<double> X_s(n, 0);
    vector<double> X_t(n, 0);
    if(s == t) {
      return 0.;
    } else if(s == landmark) {
      X_t = single_source_local_push(G, t, landmark, eps);
      delta = X_t[t]/G[t].size();
      return delta;
    } else if (t == landmark) {
      X_s = single_source_local_push(G, s, landmark, eps);
      delta = X_s[s]/G[s].size();
      return delta;
    }
    X_s = single_source_local_push(G, s, landmark, eps);
    X_t = single_source_local_push(G, t, landmark, eps);
    
    delta = X_s[s]/G[s].size() - X_s[t]/G[t].size() + X_t[t]/G[t].size() - X_t[s]/G[s].size();

    return delta;
}

double effective_resistance_bipush(graph_t &G, int s, int t, int landmark, double eps, int N0) {
    double delta;
    int n = G.size();
    vector<double> X_s(n, 0);
    vector<double> X_t(n, 0);
    vector<double> R_s(n, 0);
    vector<double> R_t(n, 0);
    if(s == t) {
      return 0.;
    } else if(s == landmark) {
      X_t = single_source_local_push_res(G, t, landmark, eps, R_t);
      delta = X_t[t]/G[t].size();
      for(size_t i=0; i<N0; i++) {
        int v = t;
        for(;;) {
          delta += R_t[v]/G[v].size()/N0;
          v = G[v][random() % G[v].size()];
          if (v == landmark)
            break;
        }
      }
      return delta;
    } else if (t == landmark) {
      X_s = single_source_local_push_res(G, s, landmark, eps, R_s);
      delta = X_s[s]/G[s].size();
      for(size_t i=0; i<N0; i++) {
        int v = s;
        for(;;) {
          delta += R_s[v]/G[v].size()/N0;
          v = G[v][random() % G[v].size()];
          if (v == landmark)
            break;
        }
      }
      return delta;
    }
    X_s = single_source_local_push_res(G, s, landmark, eps, R_s);
    X_t = single_source_local_push_res(G, t, landmark, eps, R_t);

    delta = X_s[s]/G[s].size() - X_s[t]/G[t].size() + X_t[t]/G[t].size() - X_t[s]/G[s].size();

    // run random walks
    for (size_t i = 0; i < N0; i++)
    {
      int v = s; // random walk from s
      for (;;)
      {
        delta = delta + R_s[v]/G[v].size()/N0 - R_t[v]/G[v].size()/N0;
        v = G[v][random() % G[v].size()];
        if (v == landmark)
          break;
      }
      v = t; // random walk from t
      for (;;)
      {
        delta = delta + R_t[v]/G[v].size()/N0 - R_s[v]/G[v].size()/N0;
        v = G[v][random() % G[v].size()];
        if (v == landmark)
          break;
      }
    }

    return delta;
}

vector<double> single_source_local_push_res(graph_t &G, int s, int v, double eps, vector<double> &res) {
    int n = G.size();
    vector<double> result(n, 0);
    for(int i=0; i<n; i++) res[i]=0;
    vector<bool> isInQueue(n, false);
    queue<int> r_q;

    double push_threshold = eps;
    double r_sum = 1.0;

    res[s] = 1.0;
    r_q.push(s);
    isInQueue[s] = true;
    
    while(r_q.size()>0) {
      int current_node = r_q.front();
      r_q.pop();
      isInQueue[current_node] = false;

      result[current_node] += res[current_node];
      double increment = res[current_node];
      res[current_node] = 0;
      for(int i=0; i<G[current_node].size(); i++) {
        if(G[current_node][i] == v) {
          r_sum -= increment/G[current_node].size();
          // cerr << "r_sum: " << r_sum << endl;
        } else {
          int updateNode = G[current_node][i];
          // cerr << "updateNode: " << updateNode << endl;
          res[updateNode] += increment/G[current_node].size();
          if(!isInQueue[updateNode] && res[updateNode] >= push_threshold) {
            // cerr << "updateNode: " << updateNode << endl;
            r_q.push(updateNode);
            isInQueue[updateNode] = true;
          }
        }
      }
    }
    cerr << "r_sum: " << r_sum << endl;

    return result;
}

double sample_ust_local(graph_t &G, int s, int t, int landmark, vector<int> &bfs_parent) {
    double delta = 0.;
    int n = G.size();
    vector<bool> visit(n, false);
    vector<int> next_edges(n);
    vector<int> reverse(n, 1);
    visit[landmark] = true;

    if(!visit[s]) {
      int u = s;
      while(!visit[u]) {
        int next = rand() % G[u].size();
        next_edges[u] = next;
        u = G[u][next];
      }

      u = s;
      while(!visit[u]) {
        // cout << "s: " << u << endl;
        visit[u] = true;
        u = G[u][next_edges[u]];
      }
    }

    if(!visit[t]) {
      int u = t;
      while(!visit[u]) {
        int next = rand() % G[u].size();
        next_edges[u] = next;
        u = G[u][next];
      }

      u = t;
      while(!visit[u]) {
        // cout << "t: " << u << endl;
        visit[u] = true;
        reverse[u] = -1;
        u = G[u][next_edges[u]];
      }

      while(u != landmark) {
        u = G[u][next_edges[u]];
        visit[u] = false;
      }
    }

    // update delta
    int tmp = s;
    while(tmp != landmark) {
      // cerr << tmp  << " s " << endl;
      // cout << "tmp & bfs_parent: " << tmp << " " << bfs_parent[tmp] << " next " << next_edges[tmp] << endl; // G[tmp][next_edges[.]]
      if((G[tmp][next_edges[tmp]] == bfs_parent[tmp]) && visit[tmp] && visit[bfs_parent[tmp]]) { // edge (tmp, bfs_parent[tmp]) is sampled
        delta += reverse[tmp]*1;
      } else if((G[tmp][next_edges[bfs_parent[tmp]]] == tmp) && visit[tmp] && visit[bfs_parent[tmp]]) { // edge (bfs_parent[tmp], tmp) is sampled
        delta -= reverse[bfs_parent[tmp]]*1;
      }
      tmp = bfs_parent[tmp];
    }
    tmp = t; // s->v current - t->v current
    while(tmp != landmark) {
      // cerr << tmp << endl;
      if(G[tmp][next_edges[tmp]] == bfs_parent[tmp] && visit[tmp] && visit[bfs_parent[tmp]]) { // edge (tmp, bfs_parent[tmp]) is sampled
        delta -= reverse[tmp]*1;
      } else if(G[tmp][next_edges[bfs_parent[tmp]]] == tmp && visit[tmp] && visit[bfs_parent[tmp]]) { // edge (bfs_parent[tmp], tmp) is sampled
        delta += reverse[bfs_parent[tmp]]*1;
      }
      tmp = bfs_parent[tmp];
    }
    return delta;
}

double effective_resistance_local_tree(graph_t &G, int s, int t, int landmark, int N, vector<int> &bfs_parent) {
  double delta = 0.;
  for(int i=0; i<N; i++) {
    delta += sample_ust_local(G, s, t, landmark, bfs_parent);
  }
  delta = delta/(double)N;
  return delta;
}

vector<double> single_source_local_push(graph_t &G, int s, int v, double eps) {
    int n = G.size();
    vector<double> result(n, 0);
    vector<double> res(n, 0);
    vector<bool> isInQueue(n, false);
    queue<int> r_q;

    double push_threshold = eps;
    double r_sum = 1.0;

    res[s] = 1.0;
    r_q.push(s);
    isInQueue[s] = true;
    
    while(r_q.size()>0) {
      int current_node = r_q.front();
      r_q.pop();
      isInQueue[current_node] = false;

      result[current_node] += res[current_node];
      double increment = res[current_node];
      res[current_node] = 0;
      for(int i=0; i<G[current_node].size(); i++) {
        if(G[current_node][i] == v) {
          r_sum -= increment/G[current_node].size();
          // cerr << "r_sum: " << r_sum << endl;
        } else {
          int updateNode = G[current_node][i];
          // cerr << "updateNode: " << updateNode << endl;
          res[updateNode] += increment/G[current_node].size();
          if(!isInQueue[updateNode] && res[updateNode] >= push_threshold) {
            // cerr << "updateNode: " << updateNode << endl;
            r_q.push(updateNode);
            isInQueue[updateNode] = true;
          }
        }
      }
    }
    cerr << "r_sum: " << r_sum << endl;

    return result;
}

vector<double> single_source_random_walk(graph_t &G, int s, int v, int N) {
    cerr << v << " degree: " << G[v].size() << endl; 
    int n = G.size();
    vector<double> result(n, 0);
    long long cnt_sum = 0;

    for(size_t i = 0; i < N; i++) {
      int u = s;
      result[u] += 1./N;
      int cnt = 0;
      for (;;) {
        u = G[u][random() % G[u].size()];
        ++cnt;
        if(u == v)
          break;
        else {
          result[u] += 1./N;
        }
      }
      cnt_sum += cnt;
      // cerr << cnt << endl;
    }
    cerr << "average random walk length: " << cnt_sum/N << endl;

    return result;
}

vector<double> single_source_baseline(graph_t &G, int s, int v, int N, string& filename) {
    double t_start = get_current_time_sec_method();
    cerr << v << " degree: " << G[v].size() << endl; 
    int n = G.size();
    vector<double> result(n, 0);
    // return result;
    for(int i=0; i<n; i++) {
      result[i] = effective_resistance_bipush(G, i, s, v, 1e-4, 1000);
      cerr << i << " " << result[i] << endl;
    }
    double t_end = get_current_time_sec_method();

    ofstream fout("result/" + filename + "/baseline/" + to_string(s) + ".txt");
    for(int i=0; i<n; i++) {
      fout << setprecision(16) << result[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    fout.close();

    return result;
}

vector<double> single_source_loop_erased_walk_baseline(graph_t &G, int s, int v, int N, string& filename) {
    double t_start = get_current_time_sec_method();
    int n = G.size();
    vector<double> result(n, 0);
    long long cnt_sum = 0;

    for(size_t i = 0; i < N; i++) {
      if(i % 1000 == 0) cerr << i << endl;
      int cnt = 0;
      vector<bool> intree(n, false);
      vector<int> next(n, -1);
      intree[s] = true;
      for (int j=0; j<n; j++) {
        if(intree[j]) continue;
        int u = j;
        ++cnt;
        result[u] += 1./N;
        while(!intree[u]) {
          next[u] = G[u][random() % G[u].size()];
          u = next[u];
          result[u] += 1./N;
          ++cnt;
        }
        result[u] -= 1./N;
        u = j;
        while(!intree[u]) {
          intree[u] = true;
          u = next[u];
        }
      }
      cnt_sum += cnt;
    }
    cerr << "average random walk length: " << cnt_sum/N << endl;
    double t_end = get_current_time_sec_method();

    ofstream fout("result/" + filename + "/baseline/" + to_string(s) + ".txt");
    for(int i=0; i<n; i++) {
      result[i] = result[i] / G[i].size();
      fout << setprecision(16) << result[i] << "\n";
    }
    fout << "time : " << t_end - t_start << "\n";
    fout.close();

    return result;
}

vector<double> single_source_loop_erased_walk(graph_t &G, int s, int v, int N, string& filename) {
    double t_start = get_current_time_sec_method();
    int n = G.size();
    vector<double> result(n, 0);
    long long cnt_sum = 0;

    for(size_t i = 0; i < N; i++) {
      if(i % 1000 == 0) cerr << i << endl;
      int cnt = 0;
      vector<bool> intree(n, false);
      vector<int> next(n, -1);
      intree[s] = true;
      for (int j=0; j<n; j++) {
        if(intree[j]) continue;
        int u = j;
        ++cnt;
        result[u] += 1./N;
        while(!intree[u]) {
          next[u] = G[u][random() % G[u].size()];
          u = next[u];
          result[u] += 1./N;
          ++cnt;
        }
        result[u] -= 1./N;
        u = j;
        while(!intree[u]) {
          intree[u] = true;
          u = next[u];
        }
      }
      cnt_sum += cnt;
    }
    cerr << "average random walk length: " << cnt_sum/N << endl;
    double t_end = get_current_time_sec_method();

    ofstream fout("result/" + filename + "/looperasedwalk/" + to_string(s) + ".txt");
    for(int i=0; i<n; i++) {
      result[i] = result[i] / G[i].size();
      fout << setprecision(16) << result[i] << "\n";
    }
    fout << "time : " << t_end - t_start << "\n";
    fout.close();

    return result;
}

vector<double> single_source_monte_carlo(graph_t &G, int s, int v, int N, string& filename) {
    double t_start = get_current_time_sec_method();
    cerr << v << " degree: " << G[v].size() << endl; 
    int n = G.size();
    vector<double> result(n, 0);
    return result;
    for(int i=0; i<n; i++) {
      result[i] = effective_resistance_mc3(G, i, s, v, 10000);
      cerr << i << " " << result[i] << endl;
    }
    double t_end = get_current_time_sec_method();

    ofstream fout("result/" + filename + "/montecarlo/" + to_string(s) + ".txt");
    for(int i=0; i<n; i++) {
      fout << setprecision(16) << result[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    fout.close();

    return result;
}

vector<double> single_source_monte_carlo_new(graph_t &G, int s, int v, int N, string& filename) {
    double t_start = get_current_time_sec_method();
    cerr << v << " degree: " << G[v].size() << endl; 
    int n = G.size();
    vector<double> result(n, 0);
    // return result;
    for(int i=0; i<n; i++) {
      result[i] = effective_resistance_mc3(G, i, s, v, 10000);
      cerr << i << " " << result[i] << endl;
    }
    double t_end = get_current_time_sec_method();

    ofstream fout("result/" + filename + "/montecarlo/" + to_string(s) + ".txt");
    for(int i=0; i<n; i++) {
      fout << setprecision(16) << result[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    fout.close();

    return result;
}

vector<double> single_source_spanning_tree(graph_t &G, int s, int v, int N, string& filename, vector<int> &bfs_parent) {
    double t_start = get_current_time_sec_method();
    int n = G.size();
    vector<double> result(n, 0);
    long long cnt_sum = 0;
    for(size_t i = 0; i < N; i++) {
      if(i % 1000 == 0) cerr << i << endl;
      int cnt = 0;
      vector<bool> intree(n, false);
      vector<int> next(n, -1);
      intree[s] = true;
      for (int j=0; j<n; j++) {
        if(intree[j]) continue;
        int u = j;
        ++cnt;
        // result[u] += 1./N;
        while(!intree[u]) {
          next[u] = G[u][random() % G[u].size()];
          u = next[u];
          // result[u] += 1./N;
          ++cnt;
        }
        // result[u] -= 1./N;
        u = j;
        while(!intree[u]) {
          intree[u] = true;
          u = next[u];
        }
      }
      cnt_sum += cnt;

      // build childPtr and siblingPtr
      vector<int> childPtr(n, -1);
      vector<int> siblingPtr(n,-1);

      int visitNodes = 0;
      vector<bool> visited(n, false);
      for(int k=0; k<n; k++) {
          int u = k;
          while(!visited[u]) {
            visited[u] = true;
            ++visitNodes;
            int parentU = next[u];
            if(next[u] != -1) {
              assert(siblingPtr[u] == -1);
              if(childPtr[parentU] != -1) {
                siblingPtr[u] = childPtr[parentU];
              }
              childPtr[parentU] = u;
              u = parentU;
            } else {
              break;
            }
          }
          if(visitNodes == n) break; 
      }

      // do DFS on the sampled tree
      vector<int> tVisit(n, 0);
      vector<int> tFinish(n, 0);

      std::stack<std::pair<int, int>> stack;
      stack.push({s, childPtr[s]});

      int timestamp = 0;
      do {
        int nodeu = stack.top().first;
        int nodev = stack.top().second;

        if(nodev == -1) {
          stack.pop();
          tFinish[nodeu] = ++timestamp;
        } else {
          stack.top().second = siblingPtr[nodev];
          tVisit[nodev] = ++timestamp;
          stack.push({nodev, childPtr[nodev]});
        }
      } while (!stack.empty());

      // compare BFS tree with sampled tree and aggragate
      for(int k=0; k<n; k++) {
        int u = k;

        int p = bfs_parent[u];
        int c = u;

        while (p != -1) {
          int e1 = p, e2 = c;
          bool reverse = false;
          if(e1 != next[e2]) {
            if(e2 != next[e1]) {
              c = p;
              p = bfs_parent[p];
              continue;
            }
            std::swap(e1, e2);
            reverse = true;
          }

          if(tVisit[u] >= tVisit[e2] && tFinish[u] <= tFinish[e2]) {
            result[u] += reverse ? -1./N : 1./N;
          }
          c = p;
          p = bfs_parent[p];
        }
      }
    }
    cerr << "average random walk length: " << cnt_sum/N << endl;
    double t_end = get_current_time_sec_method();

    ofstream fout("result/" + filename + "/spanningtree/" + to_string(s) + ".txt");
    for(int i=0; i<n; i++) {
      result[i] = result[i];
      fout << setprecision(16) << result[i] << "\n";
    }
    fout << "time : " << t_end - t_start << "\n";
    fout.close();

    return result;
}

vector<double> single_source_push_index(graph_t &G, int s, int v, double eps, string& filename, vector<double> &landmark_index) {
  int n = G.size();
  vector<double> result(n, 0);

  vector<double> push_result(n, 0);
  vector<double> res(n, 0);
  vector<bool> isInQueue(n, false);
  queue<int> r_q;

  double push_threshold = eps;
  double r_sum = 1.0;

  res[s] = 1.0;
  r_q.push(s);
  isInQueue[s] = true;
    
  while(r_q.size()>0) {
    int current_node = r_q.front();
    r_q.pop();
    isInQueue[current_node] = false;

    push_result[current_node] += res[current_node];
    double increment = res[current_node];
    res[current_node] = 0;
    for(int i=0; i<G[current_node].size(); i++) {
      if(G[current_node][i] == v) {
        r_sum -= increment/G[current_node].size();
      } else {
        int updateNode = G[current_node][i];
        res[updateNode] += increment/G[current_node].size();
        if(!isInQueue[updateNode] && res[updateNode] >= push_threshold) {
          r_q.push(updateNode);
          isInQueue[updateNode] = true;
        }
      }
    }
  }
  cerr << "r_sum: " << r_sum << endl;

  for(int i=0; i<n; i++) {
    result[i] = landmark_index[s] + landmark_index[i] - 2*push_result[i]/G[i].size();
  }

  ofstream fout("result/" + filename + "/pushindex/" + to_string(s) + ".txt");
  for(int i=0; i<n; i++) {
    result[i] = result[i];
    fout << setprecision(16) << result[i] << "\n";
  }
  fout.close();

  return result;
}

vector<double> single_source_abwalk_index(graph_t &G, int s, int v, int N, string& filename, vector<double> &landmark_index) {
  int n = G.size();
  vector<double> result(n, 0);

  vector<double> walk_result(n, 0);

  long long cnt_sum = 0;
  for(size_t i = 0; i < N; i++) {
    int u = s;
    walk_result[u] += 1./N;
    int cnt = 0;
    for (;;) {
      u = G[u][random() % G[u].size()];
      ++cnt;
      if(u == v)
        break;
      else {
        walk_result[u] += 1./N;
      }
    }
    cnt_sum += cnt;
  }
  cerr << "average random walk length: " << cnt_sum/N << endl;

  for(int i=0; i<n; i++) {
    result[i] = landmark_index[s] + landmark_index[i] - 2*walk_result[i]/G[i].size();
  }

  ofstream fout("result/" + filename + "/abwalkindex/" + to_string(s) + ".txt");
  for(int i=0; i<n; i++) {
    result[i] = result[i];
    fout << setprecision(16) << result[i] << "\n";
  }
  fout.close();

  return result;
}