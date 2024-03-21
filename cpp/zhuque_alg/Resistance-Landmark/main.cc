#include <iostream>
#include <vector>
#include <cmath>
#include <sys/time.h>
#include <fstream>
#include <iomanip> 
#include <pybind11/pybind11.h>
#include "graph.h"
#include "methods.h"
using namespace std;
using namespace Eigen;

double get_current_time_sec(){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

vector<int> generate_node(int n, int num){
  srandom(0);
  vector<int> node;
  for(int i=0; i<num; i++) {
    int v = random() % n;
    node.push_back(v);
  }
  return node;
}

vector<pair<int, int>> generate_vertex_pairs(int n, int num){
  srandom(0);
  vector<pair<int, int>> res;
  for (int i = 0; i < num; ++i)
  {
    int s = random() % n;
    int t = random() % n;
    res.push_back(make_pair(s, t));
  }
  return res;
}

vector<pair<int, int>> generate_edges(graph_t &G, int num){
  srandom(0);
  vector<pair<int, int>> res;
  int n = G.size();
  vector<int> sum(n + 1);
  for (int i = 1; i <= n; ++i)
  {
    sum[i] = sum[i - 1] + G[i - 1].size();
  }

  for (int i = 0; i < num; ++i)
  {
    int j = random() % sum[n];
    int v = upper_bound(sum.begin(), sum.end(), j) - sum.begin() - 1;
    res.push_back(make_pair(v, G[v][j - sum[v]]));
  }
  return res;
}

int algorithm(string file_name, int num_trials, int N, double rmax, string method, string in_pairs, string filename, string in_landmark){
    graph_t G;
    // read_graph(file_name, G);
    int landmark;
    landmark = read_graph_degree(file_name, G);
    if(in_landmark == "degree") {
        cout << "landmark: " << in_landmark << endl;
    } else if(in_landmark == "core") {
        landmark = core_landmark(G);
        cout << "landmark: " << in_landmark << endl;
        // return 0;
    } else if(in_landmark == "pagerank") {
        landmark = pagerank_landmark(G);
        cout << "landmark: " << in_landmark << endl;
        // return 0;
    } else if(in_landmark == "ecc") {
        landmark = ecc_landmark(G, landmark);
        cout << "landmark: " << in_landmark << endl;
    } else if(in_landmark == "random") {
        cout << "landmark: " << in_landmark << endl;
        landmark = 0;
    }
    int n = G.size();
    cout << "cout # of vertices: " << n << endl;
    cerr << "# of vertices: " << n << endl;

    vector<int> bfs_parent;
    vector<int> bfs_order;
    compute_BFS_tree(G, bfs_parent, bfs_order, landmark);

    vector<pair<int, int>> pairs;
    vector<int> node;
    if (in_pairs == "vertex-pairs")
    {
        pairs = generate_vertex_pairs(n, num_trials);
    }
    else if (in_pairs == "edges")
    {
        pairs = generate_edges(G, num_trials);
    }
    else if (in_pairs == "node")
    {
        node = generate_node(n, num_trials);
        ifstream infile("query/" + filename + "-node.query");
        if(infile.good()) {
            cerr << "query file already generated!" << endl;
        } else {
            ofstream fout("query/" + filename + "-node.query");
            for(int i=0; i<num_trials; i++) {
                fout << node[i] << "\n";
            }
            fout.close();
        }
        infile.close();
    }

    if(in_pairs == "vertex-pairs") {
        ifstream infile("query/" + filename + "-pair.query");
        if(infile.good()) {
            cerr << "pair query file already generated!" << endl;
        } else {
            ofstream fout("query/" + filename + "-pair.query");
            for(int i=0; i<pairs.size(); i++) {
                int s = pairs[i].first, t = pairs[i].second;
                fout << s << " " << t << "\n";
            }
            fout.close();
        }
        infile.close();
    }

    if(method == "source_landmark") {
        int s = landmark;
        cout << "landmark node: " << landmark << endl;
        vector<double> gt;
        vector<double> lewalk;
        vector<double> spanningtree;
        ifstream infile("result/" + filename + "/baseline/" + to_string(landmark) + ".txt");
        if(infile.good()) {
            for(int i=0; i<n; i++) {
                double result;
                infile >> result;
                gt.push_back(result);
            }
        } else {
            gt = single_source_loop_erased_walk_baseline(G, s, landmark, 100000, filename);
        }
        infile.close();

        double max_error = 0.;
        double l1_error = 0.;
        double time = 0.;
        double start_time = get_current_time_sec();
        lewalk = single_source_loop_erased_walk(G, s, landmark, N, filename);
        double end_time = get_current_time_sec();
        time = end_time - start_time;
        for(int i=0; i<n; i++) {
            l1_error += abs(lewalk[i] - gt[i]);
            if(abs(lewalk[i] - gt[i]) > max_error) max_error = abs(lewalk[i] - gt[i]);
        }
        cout << "loop erased walk " << "N: " << N << endl;
        cout << "max error: " << max_error << endl;
        cout << "l1 error: " << l1_error << endl;
        cout << "time: " << time << endl;
        return 0;
    } else if(method == "source") {
        ofstream fout_lewalk("result/" + filename + "-lewalk-single-source-result.txt");
        ofstream fout_spanning_tree("result/" + filename + "-spanning-tree-single-source-result.txt");
        ofstream fout_pushp("result/" + filename + "-pushp-single-source-result.txt");
        ofstream fout_abwalkp("result/" + filename + "-abwalkp-single-source-result.txt");
        ifstream infile("query/" + filename + "-node.query");
        double avg_loop_erased_walk_time, avg_loop_erased_walk_max, avg_loop_erased_walk_l1 = 0.;
        double avg_spanning_tree_time, avg_spanning_tree_max, avg_spanning_tree_l1 = 0.;
        double avg_push_index_time, avg_push_index_max, avg_push_index_l1 = 0.;
        double avg_abwalk_index_time, avg_abwalk_index_max, avg_abwalk_index_l1 = 0.;
        double avg_monte_carlo_time, avg_monte_carlo_max, avg_monte_carlo_l1 = 0.;
        for(int i=0; i<num_trials; i++) {
            int s;
            for(int j=0; j<6; j++) infile >> s;

            // ground-truth
            vector<double> landmark_index_base;
            // load index
            ifstream infile2("result/" + filename + "/baseline/" + to_string(landmark) + ".txt");
            for(int i=0; i<n; i++) {
                double result;
                infile2 >> result;
                landmark_index_base.push_back(result);
            }
            vector<double> gt;
            ifstream infile1("result/" + filename + "/baseline/" + to_string(s) + ".txt");
            if(infile1.good()) {
                cerr << "ground-truth has been generated!" << endl;
                for(int i=0; i<n; i++) {
                    double result;
                    infile1 >> result;
                    gt.push_back(result);
                }
            } else {
                cout << "ground truth: " << to_string(s) << endl;
                gt = single_source_push_index(G, s, landmark, 1e-6, filename, landmark_index_base);
            }
            infile1.close();

            // online method
            // query n Monte Carlo
            vector<double> monte_carlo;
            double max_error = 0.;
            double l1_error = 0.;
            double time = 0.;
            double start_time = get_current_time_sec();
            monte_carlo = single_source_monte_carlo(G, s, landmark, 10000, filename);
            double end_time = get_current_time_sec();
            time = end_time - start_time;
            for(int i=0; i<n; i++) {
                l1_error += abs(monte_carlo[i] - gt[i]);
                if(abs(monte_carlo[i] - gt[i]) > max_error) max_error = abs(monte_carlo[i] - gt[i]);
            }
            cout << "n monte carlo" << " N: " << endl;
            cout << "max error: " << max_error << endl;
            cout << "l1 error: " << l1_error << endl;
            cout << "time: " << time << endl;
            avg_monte_carlo_time += time;
            avg_monte_carlo_max += max_error;
            avg_monte_carlo_l1 += l1_error;

            // loop erased walk
            vector<double> lewalk;
            vector<double> spanningtree;
            l1_error = 0., max_error = 0.;
            start_time = get_current_time_sec();
            lewalk = single_source_loop_erased_walk(G, s, landmark, N, filename);
            end_time = get_current_time_sec();
            time = end_time - start_time;
            for(int i=0; i<n; i++) {
                l1_error += abs(lewalk[i] - gt[i]);
                if(abs(lewalk[i] - gt[i]) > max_error) max_error = abs(lewalk[i] - gt[i]);
            }
            cout << "loop erased walk" << " N: " << endl;
            cout << "max error: " << max_error << endl;
            cout << "l1 error: " << l1_error << endl;
            cout << "time: " << time << endl;
            avg_loop_erased_walk_time += time;
            avg_loop_erased_walk_max += max_error;
            avg_loop_erased_walk_l1 += l1_error;
            fout_lewalk << l1_error << " " << time << endl;

            // spanning tree
            vector<int> bfs_parent_s;
            vector<int> bfs_order_s;
            compute_BFS_tree(G, bfs_parent_s, bfs_order_s, s);

            l1_error = 0., max_error = 0.;
            start_time = get_current_time_sec();
            spanningtree = single_source_spanning_tree(G, s, landmark, N, filename, bfs_parent_s);
            end_time = get_current_time_sec();
            time = end_time - start_time;
            for(int i=0; i<n; i++) {
                l1_error += abs(spanningtree[i] - gt[i]);
                if(abs(spanningtree[i] - gt[i]) > max_error) max_error = abs(spanningtree[i] - gt[i]);
            }
            cout << "spanning tree " << "N: " << N << endl;
            cout << "max error: " << max_error << endl;
            cout << "l1 error: " << l1_error << endl;
            cout << "time: " << time << endl;
            avg_spanning_tree_time += time;
            avg_spanning_tree_max += max_error;
            avg_spanning_tree_l1 += l1_error;
            fout_spanning_tree << l1_error << " " << time << endl;

            // index based method

            vector<double> push_index;
            vector<double> abwalk_index;

            vector<double> landmark_index;
            // load index
            ifstream infile3("result/" + filename + "/looperasedwalk/" + to_string(landmark) + ".txt");
            for(int i=0; i<n; i++) {
                double result;
                infile3 >> result;
                landmark_index.push_back(result);
            }

            l1_error = 0., max_error = 0.;
            start_time = get_current_time_sec();
            push_index = single_source_push_index(G, s, landmark, rmax, filename, landmark_index);
            end_time = get_current_time_sec();
            time = end_time - start_time;
            for(int i=0; i<n; i++) {
                l1_error += abs(push_index[i] - gt[i]);
                if(abs(push_index[i] - gt[i]) > max_error) max_error = abs(push_index[i] - gt[i]);
            }
            cout << "push index" << " rmax: " << rmax << endl;
            cout << "max error: " << max_error << endl;
            cout << "l1 error: " << l1_error << endl;
            cout << "time: " << time << endl;
            avg_push_index_time += time;
            avg_push_index_max += max_error;
            avg_push_index_l1 += l1_error;
            fout_pushp << l1_error << " " << time << endl;

            l1_error = 0., max_error = 0.;
            start_time = get_current_time_sec();
            abwalk_index = single_source_abwalk_index(G, s, landmark, N, filename, landmark_index);
            end_time = get_current_time_sec();
            time = end_time - start_time;
            for(int i=0; i<n; i++) {
                l1_error += abs(abwalk_index[i] - gt[i]);
                if(abs(abwalk_index[i] - gt[i]) > max_error) max_error = abs(abwalk_index[i] - gt[i]);
            }
            cout << "abwalk index" << " N: " << N << endl;
            cout << "max error: " << max_error << endl;
            cout << "l1 error: " << l1_error << endl;
            cout << "time: " << time << endl;
            avg_abwalk_index_time += time;
            avg_abwalk_index_max += max_error;
            avg_abwalk_index_l1 += l1_error;
            fout_abwalkp << l1_error << " " << time << endl;
        }

        cout << "------------------avg_performance----------------------" << endl;
        cout << "N: " << N << endl;
        cout << "rmax: " << rmax << endl;
        cout << "Monte Carlo" << endl;
        cout << "avg time: " << avg_monte_carlo_time / num_trials << endl;
        cout << "avg max error: " << avg_monte_carlo_max / num_trials << endl;
        cout << "avg l1 error: " << avg_monte_carlo_l1 / num_trials << endl;

        cout << "loop erased walk" << endl;
        cout << "avg time: " << avg_loop_erased_walk_time / num_trials << endl;
        cout << "avg max error: " << avg_loop_erased_walk_max / num_trials << endl;
        cout << "avg l1 error: " << avg_loop_erased_walk_l1 / num_trials << endl;

        cout << "spanning tree" << endl;
        cout << "avg time: " << avg_spanning_tree_time / num_trials << endl;
        cout << "avg max error: " << avg_spanning_tree_max / num_trials << endl;
        cout << "avg l1 error: " << avg_spanning_tree_l1 / num_trials << endl;

        cout << "push index" << endl;
        cout << "avg time: " << avg_push_index_time / num_trials << endl;
        cout << "avg max error: " << avg_push_index_max / num_trials << endl;
        cout << "avg l1 error: " << avg_push_index_l1 / num_trials << endl;

        cout << "abwalk index" << endl;
        cout << "avg time: " << avg_abwalk_index_time / num_trials << endl;
        cout << "avg max error: " << avg_abwalk_index_max / num_trials << endl;
        cout << "avg l1 error: " << avg_abwalk_index_l1 / num_trials << endl;

        fout_lewalk.close();
        fout_abwalkp.close();
        fout_spanning_tree.close();
        fout_pushp.close();

        return 0;
    } else if(method == "monte-carlo") {
        ofstream fout_mc("result/" + filename + "-mc-single-source-result.txt");
        ifstream infile("query/" + filename + "-node.query");
        double avg_monte_carlo_time, avg_monte_carlo_max, avg_monte_carlo_l1 = 0.;
        for(int i=0; i<num_trials; i++) {
            int s;
            infile >> s;

            vector<double> landmark_index_base;
            // load index
            ifstream infile2("result/" + filename + "/baseline/" + to_string(landmark) + ".txt");
            for(int i=0; i<n; i++) {
                double result;
                infile2 >> result;
                landmark_index_base.push_back(result);
            }

            vector<double> gt;
            ifstream infile1("result/" + filename + "/baseline/" + to_string(s) + ".txt");
            if(infile1.good()) {
                cerr << "ground-truth has been generated!" << endl;
                for(int i=0; i<n; i++) {
                    double result;
                    infile1 >> result;
                    gt.push_back(result);
                }
            } else {
                cout << "ground truth: " << to_string(s) << endl;
                gt = single_source_push_index(G, s, landmark, 1e-7, filename, landmark_index_base);
            }
            infile1.close();

            vector<double> monte_carlo;
            double max_error = 0.;
            double l1_error = 0.;
            double time = 0.;
            l1_error = 0., max_error = 0.;
            double start_time = get_current_time_sec();
            monte_carlo = single_source_monte_carlo_new(G, s, landmark, 10000, filename);
            double end_time = get_current_time_sec();
            time = end_time - start_time;
            for(int i=0; i<n; i++) {
                l1_error += abs(monte_carlo[i] - gt[i]);
                if(abs(monte_carlo[i] - gt[i]) > max_error) max_error = abs(monte_carlo[i] - gt[i]);
            }
            cout << "monte carlo" << " rmax: " << rmax << endl;
            cout << "max error: " << max_error << endl;
            cout << "l1 error: " << l1_error << endl;
            cout << "time: " << time << endl;
            avg_monte_carlo_time += time;
            avg_monte_carlo_max += max_error;
            avg_monte_carlo_l1 += l1_error;
            fout_mc << l1_error << " " << time << endl;
        }
        cout << "Monte Carlo" << endl;
        cout << "avg time: " << avg_monte_carlo_time / num_trials << endl;
        cout << "avg max error: " << avg_monte_carlo_max / num_trials << endl;
        cout << "avg l1 error: " << avg_monte_carlo_l1 / num_trials << endl;
        fout_mc.close();
        return 0;
    }

    double* gt = new double[num_trials]();
    ifstream gt_in("result/" + filename + "/gt/single-pair-result.txt");
    if(gt_in.good()) {
        cerr << "ground-truth has been generated!" << endl;
        for(int i=0; i<num_trials; i++) {
            gt_in >> gt[i];
        }
    } else {
        ifstream infile1("query/" + filename + "-pair.query");
        ofstream fout1("result/" + filename + "/gt/single-pair-result.txt");
        for(int i=0; i<pairs.size(); i++) {
            int s, t;
            infile1 >> s >> t;
            cerr << s << " " << t << endl;
            double start_time = get_current_time_sec();
            double result = effective_resistance_bipush(G, s, t, landmark, 1e-5, 10000);
            double end_time = get_current_time_sec();
            fout1 << setprecision(16) << result << endl;
        }
        infile1.close();
        fout1.close();
    }
    gt_in.close();

    double max_error = 0.;
    double avg_error = 0.;

    double avg_time = 0.;

    max_error = 0.;
    avg_error = 0.;
    avg_time = 0.;

    if(method == "akp") {
        ifstream infile3("query/" + filename + "-pair.query");
        ofstream fout3("result/" + filename + "-akp-single-pair-result.txt");
        for(int i=0; i<pairs.size(); i++) {
            int s, t;
            infile3 >> s >> t;
            double start_time = get_current_time_sec();
            double result = effective_resistance_akp(G, s, t, 100, 10000);
            double end_time = get_current_time_sec();
            fout3 << setprecision(16) << abs(result-gt[i]) << " " << end_time-start_time << endl;
            if(abs(result-gt[i])>max_error) max_error = abs(result-gt[i]);
            avg_error += abs(result-gt[i]);
            avg_time += end_time-start_time;
        }
        infile3.close();
        cout << "akp" << endl;
        cout << "max error: " << max_error << endl;
        cout << "avg error: " << avg_error / num_trials << endl;
        cout << "avg time: " << avg_time / num_trials << endl;
    } else if(method == "commute") {
        ifstream infile3("query/" + filename + "-pair.query");
        ofstream fout3("result/" + filename + "-commute-single-pair-result.txt");
        for(int i=0; i<pairs.size(); i++) {
            int s, t;
            infile3 >> s >> t;
            double start_time = get_current_time_sec();
            double result = effective_resistance_mc(G, s, t, 0.1, 0.1);
            double end_time = get_current_time_sec();
            fout3 << setprecision(16) << abs(result-gt[i]) << " " << end_time-start_time << endl;
            if(abs(result-gt[i])>max_error) max_error = abs(result-gt[i]);
            avg_error += abs(result-gt[i]);
            avg_time += end_time-start_time;
        }
        infile3.close();
        cout << "commute" << endl;
        cout << "max error: " << max_error << endl;
        cout << "avg error: " << avg_error / num_trials << endl;
        cout << "avg time: " << avg_time / num_trials << endl;
    } else if(method == "abwalk") {
        ifstream infile3("query/" + filename + "-pair.query");
        ofstream fout3("result/" + filename + "-abwalk-single-pair-result.txt");
        for(int i=0; i<pairs.size(); i++) {
            int s, t;
            infile3 >> s >> t;
            double start_time = get_current_time_sec();
            double result = effective_resistance_mc3(G, s, t, landmark, N);
            double end_time = get_current_time_sec();
            fout3 << setprecision(16) << abs(result-gt[i]) << " " << end_time-start_time << endl;
            if(abs(result-gt[i])>max_error) max_error = abs(result-gt[i]);
            avg_error += abs(result-gt[i]);
            avg_time += end_time-start_time;
        }
        infile3.close();
        cout << "abwalk" << " N: " << N << endl;
        cout << "max error: " << max_error << endl;
        cout << "avg error: " << avg_error / num_trials << endl;
        cout << "avg time: " << avg_time / num_trials << endl;
    } else if(method == "localtree") {
        ifstream infile3("query/" + filename + "-pair.query");
        ofstream fout3("result/" + filename + "-localtree-single-pair-result.txt");
        for(int i=0; i<pairs.size(); i++) {
            int s, t;
            infile3 >> s >> t;
            double start_time = get_current_time_sec();
            double result = effective_resistance_local_tree(G, s, t, landmark, N, bfs_parent);
            double end_time = get_current_time_sec();
            fout3 << setprecision(16) << abs(result-gt[i]) << " " << end_time-start_time << endl;
            if(abs(result-gt[i])>max_error) max_error = abs(result-gt[i]);
            avg_error += abs(result-gt[i]);
            avg_time += end_time-start_time;
        }
        infile3.close();
        cout << "localtree" << " N: " << N << endl;
        cout << "max error: " << max_error << endl;
        cout << "avg error: " << avg_error / num_trials << endl;
        cout << "avg time: " << avg_time / num_trials << endl;
    } else if(method == "bipush") {
        ifstream infile3("query/" + filename + "-pair.query");
        ofstream fout3("result/" + filename + "-bipush-single-pair-result.txt");
        for(int i=0; i<pairs.size(); i++) {
            int s, t;
            infile3 >> s >> t;
            double start_time = get_current_time_sec();
            double result = effective_resistance_bipush(G, s, t, landmark, rmax, N);
            double end_time = get_current_time_sec();
            fout3 << setprecision(16) << abs(result-gt[i]) << " " << end_time-start_time << endl;
            if(abs(result-gt[i])>max_error) max_error = abs(result-gt[i]);
            avg_error += abs(result-gt[i]);
            avg_time += end_time-start_time;
        }
        infile3.close();
        cout << "bipush" << " N: " << N << " rmax: " << rmax << endl;
        cout << "max error: " << max_error << endl;
        cout << "avg error: " << avg_error / num_trials << endl;
        cout << "avg time: " << avg_time / num_trials << endl;
    } else if(method == "push") {
        ifstream infile3("query/" + filename + "-pair.query");
        ofstream fout3("result/" + filename + "-push-single-pair-result-4.txt");
        for(int i=0; i<pairs.size(); i++) {
            int s, t;
            infile3 >> s >> t;
            double start_time = get_current_time_sec();
            double result = effective_resistance_push(G, s, t, landmark, rmax);
            double end_time = get_current_time_sec();
            fout3 << setprecision(16) << abs(result-gt[i]) << " " << end_time-start_time << endl;
            if(abs(result-gt[i])>max_error) max_error = abs(result-gt[i]);
            avg_error += abs(result-gt[i]);
            avg_time += end_time-start_time;
        }
        infile3.close();
        cout << "push" << " rmax: " << rmax << endl;
        cout << "max error: " << max_error << endl;
        cout << "avg error: " << avg_error / num_trials << endl;
        cout << "avg time: " << avg_time / num_trials << endl;
    }

    return 0;
}

PYBIND11_MODULE(landmark, m){
m.def("algorithm", &algorithm,
    pybind11::arg("file_name")="",pybind11::arg("num_trials")=1000,pybind11::arg("N")=1e4,pybind11::arg("rmax")=1e-4,
    pybind11::arg("method")="source_landmark",pybind11::arg("in_pairs")="vertex-pairs",pybind11::arg("filename")="",
    pybind11::arg("landmark")="degree");
}