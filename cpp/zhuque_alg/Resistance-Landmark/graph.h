#pragma once

#include <vector>
#include <string>
#include <Eigen/Sparse>
using namespace std;
using namespace Eigen;

typedef vector<vector<int>> graph_t;
void read_graph(string fn, graph_t &G);
int read_graph_degree(string fn, graph_t &G);
int core_landmark(graph_t& G);
int pagerank_landmark(graph_t& G);
int ecc_landmark(graph_t& G, int& landmark);
void compute_BFS_tree(graph_t &G, vector<int> &bfs_parent, vector<int> &bfs_order, int landmark);
int computeBFS(graph_t& G, int landmark, int& ecc, vector<int>& eccLowerBound, vector<int>& distance);
size_t num_of_edges(graph_t &G);
void contract(graph_t &G, int u, int v);

SparseMatrix<double> Laplacian(graph_t &G);
int random_walk(graph_t &G, int s, int l);

template <typename T>
void laplacian_solver(graph_t &G, T &solver)
{
  SparseMatrix<double> L = Laplacian(G);
  solver.compute(L);
  cerr << "hoge" << endl;
}

template <typename T>
double effective_resistance_exact(graph_t &G, int s, int t, T &solver)
{
  VectorXd chi = VectorXd::Zero(G.size());
  chi[s] = 1;
  chi[t] = -1;
  VectorXd x = solver.solve(chi);
  VectorXd v1 = VectorXd::Ones(G.size());
  v1 /= v1.norm();
  x -= x.dot(v1) * v1;
  return x[s] - x[t];
}

void sample_spanning_tree(graph_t &G, vector<int> &order, vector<int> &next_edges);
