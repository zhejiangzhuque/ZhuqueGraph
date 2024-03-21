#include "graph.h"
using namespace std;

double effective_resistance_akp(graph_t &G, int s, int t, double eps, double lambda);
double effective_resistance_st(graph_t &G, int s, int t, double eps, double delta);
double effective_resistance_mc(graph_t &G, int s, int t, double eps, double delta);
double effective_resistance_mc2(graph_t &G, int s, int t, double eps, double delta);
double effective_resistance_mc3(graph_t &G, int s, int t, int N0, int landmark);
double effective_resistance_ex2(graph_t &G, int s, int t, double eps, double delta);
double effective_resistance_push(graph_t &G, int s, int t, int landmark, double eps);
double effective_resistance_bipush(graph_t &G, int s, int t, int landmark, double eps, int N0);
double effective_resistance_local_tree(graph_t &G, int s, int t, int landmark, int N, vector<int> &bfs_parent);
vector<double> effective_resistance_single_source_spanning_tree(graph_t &G, int s, int landmark, double eps, vector<int> &bfs_parent);
double sample_ust_local(graph_t &G, int s, int t, int landmark, vector<int> &bfs_parent);
vector<double> single_source_power_method(graph_t &G, int s, int v, int L);
vector<double> single_source_local_push(graph_t &G, int s, int v, double eps);
vector<double> single_source_local_push_res(graph_t &G, int s, int v, double eps, vector<double> &res);
vector<double> single_source_random_walk(graph_t &G, int s, int v, int N);
vector<double> single_source_loop_erased_walk_baseline(graph_t &G, int s, int v, int N, string& filename);
vector<double> single_source_loop_erased_walk(graph_t &G, int s, int v, int N, string& filename);
vector<double> single_source_spanning_tree(graph_t &G, int s, int v, int N, string& filename, vector<int> &bfs_parent);
vector<double> single_source_baseline(graph_t &G, int s, int v, int N, string& filename);
vector<double> single_source_monte_carlo(graph_t &G, int s, int v, int N, string& filename);
vector<double> single_source_monte_carlo_new(graph_t &G, int s, int v, int N, string& filename);
vector<double> single_source_push_index(graph_t &G, int s, int v, double eps, string& filename, vector<double> &landmark_index);
vector<double> single_source_abwalk_index(graph_t &G, int s, int v, int N, string& filename, vector<double> &landmark_index);
void spanning_tree_centrality(graph_t &G, int num_samples, vector<vector<double>> &centralities);
