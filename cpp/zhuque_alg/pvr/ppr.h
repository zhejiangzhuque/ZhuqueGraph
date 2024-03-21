#ifndef PPR_H
#define PPR_H

#include <vector>
#include <algorithm>
#include <queue>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <future>
#include <string>
#include <sstream>
#include "Graph.h"
#include "Random.h"
#include "alias.h"
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <sys/time.h> 
#include <time.h>

bool maxScoreCmp(const pair<int, double>& a, const pair<int, double>& b){
    return a.second > b.second;
}

class pqcompare
{
  bool reverse;
public:
  pqcompare(const bool& revparam=false)
    {reverse=revparam;}
  bool operator() (const pair<int, double>& lhs, const pair<int, double>&rhs) const
  {
    if (reverse) return (lhs.second > rhs.second);
    else return (lhs.second < rhs.second);
  }
};

void RandomWalk(double walk_num, Alias &alias, Random &R, Graph& g, int* vert_count){
    for(double i = 0; i < walk_num; i++){
        int tempNode = alias.generateRandom_t(R);
        vert_count[tempNode]++;
        while(R.drand_t() > 0.2){
            int length = g.getOutSize(tempNode);
            if(length > 0){   
                int r = R.generateRandom_t() % length;
                tempNode = g.getOutVert(tempNode, r);
            }
            vert_count[tempNode]++;
        }
    }
}

class PPR
{
friend void ppr_t_PowerMethod(PPR* ppr, vector<int> nodeList, int iterations);
public:
    double avg_mc_time;
    double avg_mc_max_err;
    double avg_mc_l1_err;
    double average_walk_time;
    double walk_forest_time_ratio;
    int error_num;
    int k;
    Graph g;
    Random R;
    int vert;
    double alpha;
    string target_filename;
    string data_path;
    double* vert_count;
    double* r; 
    unsigned NUM_CORES;
    int** multiVertCount;
    double* resultList;
    double avg_L1_error;
    double avg_max_error;
    double max_max_error;
    double avg_avg_error;
    void PowerMethodMulti(int iterations, int node_count, int num_thread);
    //const static int NUMTHREAD = 20;
    Random* Rs;

    PPR(string datapath, string name, double input_alpha) {
        avg_mc_time = 0;
        avg_mc_max_err = 0;
        avg_mc_l1_err = 0;
        avg_L1_error = 0;
        avg_max_error = 0;
        max_max_error = 0;
        avg_avg_error = 0;
        average_walk_time = 0.;
        walk_forest_time_ratio = 0.;
        target_filename = name;
        data_path = datapath;
        string filename = datapath + "/" + name + ".txt";
        cout << "filename: " << filename << endl;
        g.inputGraph(filename);
        cout << "edge num: " << g.m << endl;
        vert = g.n;
        alpha = input_alpha;
        srand(unsigned(time(0)));
        R = Random(unsigned(rand()));
        resultList = new double[vert];
        r = new double[vert];
        for(int i =0 ; i < vert; i++){
            resultList[i] = 0;
            r[i] = 0;
        }
        NUM_CORES = std::thread::hardware_concurrency();
        assert(NUM_CORES >= 2);
        cout << "thread core: " << NUM_CORES << endl;
        multiVertCount = new int*[NUM_CORES];
        Rs = new Random[NUM_CORES];
        for(int i = 0; i < NUM_CORES; i++){
            Rs[i] = Random(unsigned(rand()));
            multiVertCount[i] = new int[vert];
            for(int j = 0; j < vert; j++){
                multiVertCount[i][j] = 0;
            }
        }
        cout << "init done! " << endl;

        stringstream ssin;
        ssin << "ppr-answer/" << target_filename << "wf-ratio-" << alpha << ".txt";
        string inputFile = ssin.str();
        if(!is_file_exist(inputFile)){
            cout << "please generate walk time ratio" << endl;
            for(int i=0; i<100; i++) {
                cout << "error*****************************" << endl;
            }
        }
        ifstream fin(inputFile);
        fin >> walk_forest_time_ratio;
        fin.close();
    }
    ~PPR() {
        for(int i = 0; i < NUM_CORES; i++){
            delete[] multiVertCount[i];
        }
        delete[] multiVertCount;
        delete[] r;
        delete[] Rs;
    }

    bool is_file_exist(string fileName)
    {
        ifstream infile(fileName);
        return infile.good();
    }

    //取s点的groundtruth
    vector<int> getRealTopK(int s, int k){
        stringstream ss;
        ss << "ppr-answer/" << target_filename << "/" << s << "-" << alpha << ".txt";
        string infile = ss.str();
        ifstream real(infile);
        vector<int> realList;
        vector<double> simList;
        for(int i = 0; i < vert; i++){
            int tempId; 
            double tempSim;
            real >> tempId >> tempSim;
            if(i >= k && tempSim < simList[k-1]){
                break; 
             } 
             realList.push_back(tempId);
             simList.push_back(tempSim);
         }
         real.close();
         return realList;
     }

    unordered_map<int, double> getRealTopKMap(int s, int k){
        unordered_map<int, double> answer_map;
        stringstream ss;
        ss << "ppr-answer/" << target_filename << "/" << s << "-" << alpha << ".txt";
        string infile = ss.str();
        ifstream real(infile);
        double k_Sim = 0;
        for(int i = 0; i < vert; i++){
            int tempId;
            double tempSim;
            real >> tempId >> tempSim;
            if(i == k - 1){
                k_Sim = tempSim;
            }
            if(i >= k && tempSim < k_Sim){
                break;
            }
            answer_map[tempId] = tempSim;
        }
        real.close();
        return answer_map;
    }

    //在一个大的pair<int, double>数组中取TopK大个（按照double的值大小）
    vector<pair<int, double> > getTopK(vector<pair<int, double> > target, int k){
        typedef priority_queue<pair<int, double>, vector<pair<int,double> >, pqcompare> pq;
        pq upper_pq(pqcompare(true));

        double UpperBound = 0;
        for(int i = 0; i < target.size(); i++){ 
            if(i < k){
                upper_pq.push(target[i]);
                if(i == k - 1)
                    UpperBound = upper_pq.top().second;
            }
            else{
                if(target[i].second > UpperBound){
                    upper_pq.pop();
                    upper_pq.push(target[i]);   
                    UpperBound = upper_pq.top().second;
                }
            }
        }

        vector<pair<int, double> > answer;
        for(int i = 0; i < k; i++){
            answer.push_back(upper_pq.top());
            upper_pq.pop();           
        }
        return answer;
    }

    int* loop_erased_walk() {
        bool* intree = new bool[vert];
        int* next = new int[vert];
        int* root = new int[vert];
        for(int i=0; i<vert; i++) {
            intree[i] = false;
            next[i] = -1;
            root[i] = -1;
        }
        for(int i=0; i<vert; i++) {
            int u = i;
            while(!intree[u]) {
                if(R.drand() < alpha) {
                    intree[u] = true;
                    root[u] = u;
                } else {
                    int length = g.getOutSize(u);
                    if(length == 0){                    
                        next[u] = u;                   
                    }
                    else{
                        int tempIndex = R.generateRandom() % length;
                        next[u] = g.getOutVert(u, tempIndex);
                    }
                    u = next[u];
                }
            }
            int r = root[u];
            u = i;
            while(!intree[u]) {
                root[u] = r;
                intree[u] = true;
                u = next[u];
            }
        }
        return root;
    }

    double* naive_monte_carlo(int s, double epsilon, double& time) {
        cout << "node: " << s << endl;
        double* resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        if(g.getOutSize(s) == 0){
            resultList[s] = 1;
            return resultList;
        }
        double walk_num = vert*log(vert)/epsilon/epsilon;
        clock_t t0 = clock();
        for(double i = 0; i < walk_num; i++){
            int tempNode = s;
            resultList[tempNode] += alpha;
            while(R.drand() > alpha){
                int length = g.getOutSize(tempNode);
                if(length == 0)
                    tempNode = tempNode;
                else{
                    int r = R.generateRandom() % length;
                    tempNode = g.getOutVert(tempNode, r);
                }
                resultList[tempNode] += alpha;
            }
        }
        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        for(int i = 0; i < vert; i++){
            resultList[i] /= (double) walk_num;
        }
        return resultList;
    }

    double* PPW_k(int s, double epsilon, double& time, int batch_size, double rd_ratio, int power_iterations) {
        clock_t t0 = clock();
        cout << "node: " << s << endl;
        double* resultList = new double[vert];
        double* former_resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        if(g.getOutSize(s) == 0){
            resultList[s] = 1;
            return resultList;
        }
        double walk_num = vert*log(vert)/epsilon/100;
        int K = 1000;
        int walk_num_each_batch = walk_num / batch_size;

        double random_walk_this_round;
        while(batch_size>0) {
            clock_t t_build_start = clock();
            double* res = new double[vert]();
            for(int i=0; i<vert; i++) {
                if(resultList[i] != 0) {
                    if(g.getOutSize(i) == 0) {
                        res[i] += (1-alpha)/alpha*resultList[i];
                    } else {
                        for(int j=0; j<g.getOutSize(i); j++) {
                            int tmp_node = g.getOutVert(i,j);
                            res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                res[i] -= resultList[i]/alpha;
            }
            res[s] += 1.;
            clock_t t_build_end = clock();
            cout << "compute r time: " << (t_build_end - t_build_start) / (double) CLOCKS_PER_SEC << endl;
            clock_t t_alias_start = clock();
            double r_sum = 0.;
            double max_r = 0.;
            vector<pair<int, double> > aliasP;
            
            for(int i = 0; i < vert; i++){
                if(res[i] != 0) {
                    r_sum += abs(res[i]);
                    aliasP.push_back(pair<int, double>(i, abs(res[i])));
                    if(abs(res[i])>max_r) max_r = abs(res[i]);
                }
            }
            cout << "rsum: " << r_sum << endl;
            cout << "rmax: " << max_r << endl;
            Alias alias = Alias(aliasP);
            clock_t t_alias_end = clock();
            cout << "bulid alias time: " << (t_alias_end - t_alias_start) / (double) CLOCKS_PER_SEC << endl;

            clock_t t_walk_start = clock();
            for(int i = 0; i < walk_num_each_batch; i++) {
                int tempNode = alias.generateRandom(R);
                int sourceNode = tempNode;
                if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                while(R.drand() > alpha){
                    int length = g.getOutSize(tempNode);
                    if(length == 0){                    
                        tempNode = tempNode;                   
                    }
                    else{
                        int tempIndex = R.generateRandom() % length;
                        tempNode = g.getOutVert(tempNode, tempIndex);
                    }
                    if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                    if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                }
                // if(res[sourceNode]>0) resultList[tempNode] += r_sum / walk_num_each_batch;
                // if(res[sourceNode]<0) resultList[tempNode] -= r_sum / walk_num_each_batch;
            }
            clock_t t_walk_end = clock();
            cout << "random walk time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;
            cout << "random walk number: " << walk_num_each_batch << endl;

            random_walk_this_round = t_walk_end - t_walk_start;

            clock_t t_power_start = clock();
            for(int k=0; k<power_iterations; k++) {
                for(int i=0; i<vert; i++) {
                    former_resultList[i] = resultList[i];
                    resultList[i] = 0;
                }
                for(int i=0; i<vert; i++) {
                    if(former_resultList[i]>0) {
                        int length = g.getOutSize(i);
                        if(length == 0) {
                            resultList[i] += (1-alpha)*former_resultList[i];
                        } else {
                            for(int j=0; j<length; j++) {
                                int tmp_node = g.getOutVert(i,j);
                                resultList[tmp_node] += (1-alpha)*former_resultList[i] / length;
                            }
                        }
                    }
                }
                resultList[s] += alpha;
                clock_t t_power_current = clock();
                // double current_time = t_power_current - t_power_start;
                // if(current_time > random_walk_this_round*rd_ratio) {
                //     cout << "K: " << k << " break: power too much" << endl;
                //     break;
                // }
            }
            clock_t t_power_end = clock();
            cout << "power time: " << (t_power_end - t_power_start) / (double) CLOCKS_PER_SEC << endl;

            batch_size--;
            delete[] res;
        }

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* PPF_k(int s, double epsilon, double& time, string graphd, int batch_size, double rd_ratio, int power_iterations) {
        clock_t t0 = clock();
        cout << "node: " << s << endl;
        double* resultList = new double[vert];
        double* former_resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        if(g.getOutSize(s) == 0){
            resultList[s] = 1;
            return resultList;
        }
        int K = 1000;
        double num_forests = log(vert)/epsilon/100*walk_forest_time_ratio;
        int num_forest_each_batch = num_forests / batch_size;

        double walk_num = vert/epsilon/100; // burn up random walks
        clock_t t_walk0_start = clock();
        for(double i = 0; i < walk_num; i++){
            int tempNode = s;
            resultList[tempNode] += alpha/walk_num;
            while(R.drand() > alpha){
                int length = g.getOutSize(tempNode);
                if(length == 0)
                    tempNode = tempNode;
                else{
                    int r = R.generateRandom() % length;
                    tempNode = g.getOutVert(tempNode, r);
                }
                resultList[tempNode] += alpha/walk_num;
            }
        }
        clock_t t_walk0_end = clock();
        cout << "burn up random walk time: " << (t_walk0_end - t_walk0_start) / (double) CLOCKS_PER_SEC << endl;

        double random_walk_this_round;

        while(batch_size>0) {
            clock_t t_build_start = clock();
            double* res = new double[vert]();
            for(int i=0; i<vert; i++) {
                if(resultList[i] != 0) {
                    if(g.getOutSize(i) == 0) {
                        res[i] += (1-alpha)/alpha*resultList[i];
                    } else {
                        for(int j=0; j<g.getOutSize(i); j++) {
                            int tmp_node = g.getOutVert(i,j);
                            res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                res[i] -= resultList[i]/alpha;
            }
            res[s] += 1.;
            clock_t t_build_end = clock();
            cout << "compute r time: " << (t_build_end - t_build_start) / (double) CLOCKS_PER_SEC << endl;
            clock_t t_alias_start = clock();
            double r_sum = 0.;
            double max_r = 0.;
            vector<pair<int, double> > aliasP;
            
            for(int i = 0; i < vert; i++){
                if(res[i] != 0) {
                    r_sum += abs(res[i]);
                    aliasP.push_back(pair<int, double>(i, abs(res[i])));
                    if(abs(res[i])>max_r) max_r = abs(res[i]);
                }
            }
            cout << "rsum: " << r_sum << endl;
            cout << "rmax: " << max_r << endl;
            Alias alias = Alias(aliasP);
            clock_t t_alias_end = clock();
            cout << "bulid alias time: " << (t_alias_end - t_alias_start) / (double) CLOCKS_PER_SEC << endl;

            clock_t t_walk_start = clock();
            // int forest_num = ceil((double)walk_num_each_batch/vert);
            // int forest_num = 10;
            cout << "sample rooted spanning forests: " << num_forest_each_batch << endl;
            for(int i=0; i<num_forest_each_batch; i++) {
                int* root = loop_erased_walk();
                if(graphd == "directed") {
                    for(int id=0; id<vert; id++) {
                        // if(g.getOutSize(i) == 0) continue; 
                        resultList[root[id]] += res[id]/(double)num_forest_each_batch; 
                    }
                } else if(graphd == "undirected") {
                    double* partition_acc = new double[vert]();
                    double* residual_over_partition = new double[vert]();
                    for(int id=0; id<vert; id++) {
                        partition_acc[root[id]] += g.getOutSize(id); 
                        residual_over_partition[root[id]] += res[id];
                    }
                    for(int id=0; id<vert; id++) {
                        // if(g.degree[id] == 0) continue;
                        resultList[id] += g.getOutSize(id)*residual_over_partition[root[id]]/partition_acc[root[id]]/(double)num_forest_each_batch;
                    }
                    delete[] partition_acc;
                    delete[] residual_over_partition;
                }
                delete[] root;
            }
            clock_t t_walk_end = clock();
            cout << "sampling spanning forests time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;

            random_walk_this_round = t_walk_end - t_walk_start;

            clock_t t_power_start = clock();
            for(int k=0; k<power_iterations; k++) {
                for(int i=0; i<vert; i++) {
                    former_resultList[i] = resultList[i];
                    resultList[i] = 0;
                }
                for(int i=0; i<vert; i++) {
                    if(former_resultList[i]!=0) {
                        int length = g.getOutSize(i);
                        if(length == 0) {
                            resultList[i] += (1-alpha)*former_resultList[i];
                        } else {
                            for(int j=0; j<length; j++) {
                                int tmp_node = g.getOutVert(i,j);
                                resultList[tmp_node] += (1-alpha)*former_resultList[i] / length;
                            }
                        }
                    }
                }
                resultList[s] += alpha;
                clock_t t_power_current = clock();
                // double current_time = t_power_current - t_power_start;
                // if(current_time > random_walk_this_round*rd_ratio) {
                //     cout << "K: " << k << " break: power too much" << endl;
                //     break;
                // }
            }
            clock_t t_power_end = clock();
            cout << "power time: " << (t_power_end - t_power_start) / (double) CLOCKS_PER_SEC << endl;

            batch_size--;
            delete[] res;
        }

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* pc_PPW_k(double epsilon, double& time, int batch_size, double rd_ratio, int power_iterations) {
        clock_t t0 = clock();
        cout << "pc: PPW"<< endl;
        double* resultList = new double[vert];
        double* former_resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        double walk_num = vert*log(vert)/epsilon/100;
        int K = 1000;
        int walk_num_each_batch = walk_num / batch_size;
        
        double random_walk_this_round;
        while(batch_size>0) {
            clock_t t_build_start = clock();
            double* res = new double[vert]();
            for(int i=0; i<vert; i++) {
                if(resultList[i] != 0) {
                    if(g.getOutSize(i) == 0) {
                        res[i] += (1-alpha)/alpha*resultList[i];
                    } else {
                        for(int j=0; j<g.getOutSize(i); j++) {
                            int tmp_node = g.getOutVert(i,j);
                            res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                res[i] -= resultList[i]/alpha;
            }
            for(int i=0; i<vert; i++) {
                res[i] += 1./vert;
            }
            clock_t t_build_end = clock();
            cout << "compute r time: " << (t_build_end - t_build_start) / (double) CLOCKS_PER_SEC << endl;
            clock_t t_alias_start = clock();
            double r_sum = 0.;
            double max_r = 0.;
            vector<pair<int, double> > aliasP;
            
            for(int i = 0; i < vert; i++){
                if(res[i] != 0) {
                    r_sum += abs(res[i]);
                    aliasP.push_back(pair<int, double>(i, abs(res[i])));
                    if(abs(res[i])>max_r) max_r = abs(res[i]);
                }
            }
            cout << "rsum: " << r_sum << endl;
            cout << "rmax: " << max_r << endl;
            Alias alias = Alias(aliasP);
            clock_t t_alias_end = clock();
            cout << "bulid alias time: " << (t_alias_end - t_alias_start) / (double) CLOCKS_PER_SEC << endl;

            clock_t t_walk_start = clock();
            for(int i = 0; i < walk_num_each_batch; i++) {
                int tempNode = alias.generateRandom(R);
                int sourceNode = tempNode;
                if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                while(R.drand() > alpha){
                    int length = g.getOutSize(tempNode);
                    if(length == 0){                    
                        tempNode = tempNode;                   
                    }
                    else{
                        int tempIndex = R.generateRandom() % length;
                        tempNode = g.getOutVert(tempNode, tempIndex);
                    }
                    if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                    if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                }
                // if(res[sourceNode]>0) resultList[tempNode] += r_sum / walk_num_each_batch;
                // if(res[sourceNode]<0) resultList[tempNode] -= r_sum / walk_num_each_batch;
            }
            clock_t t_walk_end = clock();
            cout << "random walk time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;
            cout << "random walk number: " << walk_num_each_batch << endl;

            random_walk_this_round = t_walk_end - t_walk_start;

            clock_t t_power_start = clock();
            for(int k=0; k<power_iterations; k++) {
                for(int i=0; i<vert; i++) {
                    former_resultList[i] = resultList[i];
                    resultList[i] = 0;
                }
                for(int i=0; i<vert; i++) {
                    if(former_resultList[i]>0) {
                        int length = g.getOutSize(i);
                        if(length == 0) {
                            resultList[i] += (1-alpha)*former_resultList[i];
                        } else {
                            for(int j=0; j<length; j++) {
                                int tmp_node = g.getOutVert(i,j);
                                resultList[tmp_node] += (1-alpha)*former_resultList[i] / length;
                            }
                        }
                    }
                }
                for(int i=0; i<vert; i++) {
                    resultList[i] += alpha/vert;
                }
                clock_t t_power_current = clock();
                // double current_time = t_power_current - t_power_start;
                // if(current_time > random_walk_this_round*rd_ratio) {
                //     cout << "K: " << k << " break: power too much" << endl;
                //     break;
                // }
            }
            clock_t t_power_end = clock();
            cout << "power time: " << (t_power_end - t_power_start) / (double) CLOCKS_PER_SEC << endl;

            batch_size--;
            delete[] res;
        }

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* pc_PPF_k(double epsilon, double& time, string graphd, int batch_size, double rd_ratio, int power_iterations) {
        clock_t t0 = clock();
        cout << "pc: PPF"<< endl;
        double* resultList = new double[vert];
        double* former_resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 1./vert;
        int K = 1000;
        double num_forests = log(vert)/epsilon/100*walk_forest_time_ratio;
        double num_forest_each_batch = num_forests / batch_size;

        double random_walk_this_round;

        while(batch_size>0) {
            clock_t t_build_start = clock();
            double* res = new double[vert]();
            for(int i=0; i<vert; i++) {
                if(resultList[i] != 0) {
                    if(g.getOutSize(i) == 0) {
                        res[i] += (1-alpha)/alpha*resultList[i];
                    } else {
                        for(int j=0; j<g.getOutSize(i); j++) {
                            int tmp_node = g.getOutVert(i,j);
                            res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                res[i] -= resultList[i]/alpha;
            }
            for(int i=0; i<vert; i++) {
                res[i] += 1./vert;
            }
            clock_t t_build_end = clock();
            cout << "compute r time: " << (t_build_end - t_build_start) / (double) CLOCKS_PER_SEC << endl;
            clock_t t_alias_start = clock();
            double r_sum = 0.;
            double max_r = 0.;
            vector<pair<int, double> > aliasP;
            
            for(int i = 0; i < vert; i++){
                if(res[i] != 0) {
                    r_sum += abs(res[i]);
                    aliasP.push_back(pair<int, double>(i, abs(res[i])));
                    if(abs(res[i])>max_r) max_r = abs(res[i]);
                }
            }
            cout << "rsum: " << r_sum << endl;
            cout << "rmax: " << max_r << endl;
            Alias alias = Alias(aliasP);
            clock_t t_alias_end = clock();
            cout << "bulid alias time: " << (t_alias_end - t_alias_start) / (double) CLOCKS_PER_SEC << endl;

            clock_t t_walk_start = clock();
            cout << "sample rooted spanning forests: " << num_forest_each_batch << endl;
            for(int i=0; i<num_forest_each_batch; i++) {
                int* root = loop_erased_walk();
                if(graphd == "directed") {
                    for(int id=0; id<vert; id++) { 
                        resultList[root[id]] += res[id]/(double)num_forest_each_batch; 
                    }
                } else if(graphd == "undirected") {
                    double* partition_acc = new double[vert]();
                    double* residual_over_partition = new double[vert]();
                    for(int id=0; id<vert; id++) {
                        partition_acc[root[id]] += g.getOutSize(id); 
                        residual_over_partition[root[id]] += res[id];
                    }
                    for(int id=0; id<vert; id++) {
                        resultList[id] += g.getOutSize(id)*residual_over_partition[root[id]]/partition_acc[root[id]]/(double)num_forest_each_batch;
                    }
                    delete[] partition_acc;
                    delete[] residual_over_partition;
                }
                delete[] root;
            }
            clock_t t_walk_end = clock();
            cout << "sampling spanning forests time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;

            random_walk_this_round = t_walk_end - t_walk_start;

            clock_t t_power_start = clock();
            for(int k=0; k<power_iterations; k++) {
                for(int i=0; i<vert; i++) {
                    former_resultList[i] = resultList[i];
                    resultList[i] = 0;
                }
                for(int i=0; i<vert; i++) {
                    if(former_resultList[i]>0) {
                        int length = g.getOutSize(i);
                        if(length == 0) {
                            resultList[i] += (1-alpha)*former_resultList[i];
                        } else {
                            for(int j=0; j<length; j++) {
                                int tmp_node = g.getOutVert(i,j);
                                resultList[tmp_node] += (1-alpha)*former_resultList[i] / length;
                            }
                        }
                    }
                }
                for(int i=0; i<vert; i++) {
                    resultList[i] += alpha/vert;
                }
                clock_t t_power_current = clock();
                // double current_time = t_power_current - t_power_start;
                // if(current_time > random_walk_this_round*rd_ratio) {
                //     cout << "K: " << k << " break: power too much" << endl;
                //     break;
                // }
            }
            clock_t t_power_end = clock();
            cout << "power time: " << (t_power_end - t_power_start) / (double) CLOCKS_PER_SEC << endl;

            batch_size--;
            delete[] res;
        }

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* PPW(int s, double epsilon, double& time, int batch_size, double rd_ratio) {
        clock_t t0 = clock();
        cout << "node: " << s << endl;
        double* resultList = new double[vert];
        double* former_resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        if(g.getOutSize(s) == 0){
            resultList[s] = 1;
            return resultList;
        }
        // double walk_num = vert*log(vert)/epsilon/100;
        double walk_num = alpha*log(epsilon*epsilon)/log(1-alpha)*vert*log(vert)/100;
        int K = 1000;
        int walk_num_each_batch = walk_num / batch_size;

        double random_walk_this_round;
        while(batch_size>0) {
            clock_t t_build_start = clock();
            double* res = new double[vert]();
            for(int i=0; i<vert; i++) {
                if(resultList[i] != 0) {
                    if(g.getOutSize(i) == 0) {
                        res[i] += (1-alpha)/alpha*resultList[i];
                    } else {
                        for(int j=0; j<g.getOutSize(i); j++) {
                            int tmp_node = g.getOutVert(i,j);
                            res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                res[i] -= resultList[i]/alpha;
            }
            res[s] += 1.;
            clock_t t_build_end = clock();
            cout << "compute r time: " << (t_build_end - t_build_start) / (double) CLOCKS_PER_SEC << endl;
            clock_t t_alias_start = clock();
            double r_sum = 0.;
            double max_r = 0.;
            vector<pair<int, double> > aliasP;
            
            for(int i = 0; i < vert; i++){
                if(res[i] != 0) {
                    r_sum += abs(res[i]);
                    aliasP.push_back(pair<int, double>(i, abs(res[i])));
                    if(abs(res[i])>max_r) max_r = abs(res[i]);
                }
            }
            cout << "rsum: " << r_sum << endl;
            cout << "rmax: " << max_r << endl;
            Alias alias = Alias(aliasP);
            clock_t t_alias_end = clock();
            cout << "bulid alias time: " << (t_alias_end - t_alias_start) / (double) CLOCKS_PER_SEC << endl;

            clock_t t_walk_start = clock();
            for(int i = 0; i < walk_num_each_batch; i++) {
                int tempNode = alias.generateRandom(R);
                int sourceNode = tempNode;
                if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                while(R.drand() > alpha){
                    int length = g.getOutSize(tempNode);
                    if(length == 0){                    
                        tempNode = tempNode;                   
                    }
                    else{
                        int tempIndex = R.generateRandom() % length;
                        tempNode = g.getOutVert(tempNode, tempIndex);
                    }
                    if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                    if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                }
                // if(res[sourceNode]>0) resultList[tempNode] += r_sum / walk_num_each_batch;
                // if(res[sourceNode]<0) resultList[tempNode] -= r_sum / walk_num_each_batch;
            }
            clock_t t_walk_end = clock();
            cout << "random walk time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;
            cout << "random walk number: " << walk_num_each_batch << endl;

            random_walk_this_round = t_walk_end - t_walk_start;

            clock_t t_power_start = clock();
            for(int k=0; k<K; k++) {
                for(int i=0; i<vert; i++) {
                    former_resultList[i] = resultList[i];
                    resultList[i] = 0;
                }
                for(int i=0; i<vert; i++) {
                    if(former_resultList[i]>0) {
                        int length = g.getOutSize(i);
                        if(length == 0) {
                            resultList[i] += (1-alpha)*former_resultList[i];
                        } else {
                            for(int j=0; j<length; j++) {
                                int tmp_node = g.getOutVert(i,j);
                                resultList[tmp_node] += (1-alpha)*former_resultList[i] / length;
                            }
                        }
                    }
                }
                resultList[s] += alpha;
                clock_t t_power_current = clock();
                double current_time = t_power_current - t_power_start;
                if(current_time > random_walk_this_round*rd_ratio) {
                    cout << "K: " << k << " break: power too much" << endl;
                    break;
                }
            }
            clock_t t_power_end = clock();
            cout << "power time: " << (t_power_end - t_power_start) / (double) CLOCKS_PER_SEC << endl;

            batch_size--;
            delete[] res;
        }

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* PPF(int s, double epsilon, double& time, string graphd, int batch_size, double rd_ratio) {
        clock_t t0 = clock();
        cout << "node: " << s << endl;
        double* resultList = new double[vert];
        double* former_resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        if(g.getOutSize(s) == 0){
            resultList[s] = 1;
            return resultList;
        }
        int K = 1000;
        // double num_forests = log(vert)/epsilon/100*walk_forest_time_ratio;
        double num_forests = alpha*log(vert)*log(epsilon*epsilon)/log(1-alpha)/100*walk_forest_time_ratio;
        int num_forest_each_batch = num_forests / batch_size;

        // double rmax = 1e-4;

        // clock_t t_push_start = clock();
        // double* pres = new double[vert]();
        // bool* isInQueue = new bool[vert];
        // for(int i=0; i<vert; i++) {
        //     isInQueue[i] = false;
        // }
        // queue<int> r_queue;
        // r_queue.push(s);
        // pres[s] = 1.0;
        // double r_sum = 1.0;
        // while(r_queue.size() > 0){
        //     int tempNode = r_queue.front();
        //     r_queue.pop();
        //     isInQueue[tempNode] = false;
        //     int tempOutSize = g.getOutSize(tempNode);
        //     double tempRes = pres[tempNode];
        //     resultList[tempNode] += alpha * tempRes;
        //     r_sum -= alpha * tempRes;
        //     // cout << r_sum << endl;
        //     pres[tempNode] = 0;
        //     if(tempOutSize == 0) {
        //         int updateNode = tempNode;
        //         pres[updateNode] += (1-alpha)*tempRes;
        //         if(!isInQueue[updateNode] && pres[updateNode]>rmax) {
        //             r_queue.push(updateNode);
        //             isInQueue[updateNode] = true;
        //         }
        //     } else {
        //         for(int i=0; i<tempOutSize; i++) {
        //             int updateNode = g.getOutVert(tempNode, i);
        //             pres[updateNode] += (1-alpha)*tempRes / tempOutSize;
        //             if(!isInQueue[updateNode] && pres[updateNode]>rmax) {
        //                 r_queue.push(updateNode);
        //                 isInQueue[updateNode] = true;
        //             }
        //         }
        //     }
        // }
        // clock_t t_push_end = clock();
        // cout << "push time: " << (t_push_end - t_push_start) / (double) CLOCKS_PER_SEC << endl;
        // cout << "current r_sum: " << r_sum << endl;

        double walk_num = vert/epsilon/100; // burn up random walks
        clock_t t_walk0_start = clock();
        for(double i = 0; i < walk_num; i++){
            int tempNode = s;
            resultList[tempNode] += alpha/walk_num;
            while(R.drand() > alpha){
                int length = g.getOutSize(tempNode);
                if(length == 0)
                    tempNode = tempNode;
                else{
                    int r = R.generateRandom() % length;
                    tempNode = g.getOutVert(tempNode, r);
                }
                resultList[tempNode] += alpha/walk_num;
            }
        }
        clock_t t_walk0_end = clock();
        cout << "burn up random walk time: " << (t_walk0_end - t_walk0_start) / (double) CLOCKS_PER_SEC << endl;

        double random_walk_this_round;

        while(batch_size>0) {
            clock_t t_build_start = clock();
            double* res = new double[vert]();
            for(int i=0; i<vert; i++) {
                if(resultList[i] != 0) {
                    if(g.getOutSize(i) == 0) {
                        res[i] += (1-alpha)/alpha*resultList[i];
                    } else {
                        for(int j=0; j<g.getOutSize(i); j++) {
                            int tmp_node = g.getOutVert(i,j);
                            res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                res[i] -= resultList[i]/alpha;
            }
            res[s] += 1.;
            clock_t t_build_end = clock();
            cout << "compute r time: " << (t_build_end - t_build_start) / (double) CLOCKS_PER_SEC << endl;
            clock_t t_alias_start = clock();
            double r_sum = 0.;
            double max_r = 0.;
            vector<pair<int, double> > aliasP;
            
            for(int i = 0; i < vert; i++){
                if(res[i] != 0) {
                    r_sum += abs(res[i]);
                    aliasP.push_back(pair<int, double>(i, abs(res[i])));
                    if(abs(res[i])>max_r) max_r = abs(res[i]);
                }
            }
            cout << "rsum: " << r_sum << endl;
            cout << "rmax: " << max_r << endl;
            Alias alias = Alias(aliasP);
            clock_t t_alias_end = clock();
            cout << "bulid alias time: " << (t_alias_end - t_alias_start) / (double) CLOCKS_PER_SEC << endl;

            clock_t t_walk_start = clock();
            // int forest_num = ceil((double)walk_num_each_batch/vert);
            // int forest_num = 10;
            cout << "sample rooted spanning forests: " << num_forest_each_batch << endl;
            for(int i=0; i<num_forest_each_batch; i++) {
                int* root = loop_erased_walk();
                if(graphd == "directed") {
                    for(int id=0; id<vert; id++) {
                        // if(g.getOutSize(i) == 0) continue; 
                        resultList[root[id]] += res[id]/(double)num_forest_each_batch; 
                    }
                } else if(graphd == "undirected") {
                    double* partition_acc = new double[vert]();
                    double* residual_over_partition = new double[vert]();
                    for(int id=0; id<vert; id++) {
                        partition_acc[root[id]] += g.getOutSize(id); 
                        residual_over_partition[root[id]] += res[id];
                    }
                    for(int id=0; id<vert; id++) {
                        // if(g.degree[id] == 0) continue;
                        resultList[id] += g.getOutSize(id)*residual_over_partition[root[id]]/partition_acc[root[id]]/(double)num_forest_each_batch;
                    }
                    delete[] partition_acc;
                    delete[] residual_over_partition;
                }
                delete[] root;
            }
            clock_t t_walk_end = clock();
            cout << "sampling spanning forests time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;

            random_walk_this_round = t_walk_end - t_walk_start;

            clock_t t_power_start = clock();
            for(int k=0; k<K; k++) {
                for(int i=0; i<vert; i++) {
                    former_resultList[i] = resultList[i];
                    resultList[i] = 0;
                }
                for(int i=0; i<vert; i++) {
                    if(former_resultList[i]!=0) {
                        int length = g.getOutSize(i);
                        if(length == 0) {
                            resultList[i] += (1-alpha)*former_resultList[i];
                        } else {
                            for(int j=0; j<length; j++) {
                                int tmp_node = g.getOutVert(i,j);
                                resultList[tmp_node] += (1-alpha)*former_resultList[i] / length;
                            }
                        }
                    }
                }
                resultList[s] += alpha;
                clock_t t_power_current = clock();
                double current_time = t_power_current - t_power_start;
                if(current_time > random_walk_this_round*rd_ratio) {
                    cout << "K: " << k << " break: power too much" << endl;
                    break;
                }
            }
            clock_t t_power_end = clock();
            cout << "power time: " << (t_power_end - t_power_start) / (double) CLOCKS_PER_SEC << endl;

            batch_size--;
            delete[] res;
        }

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* DPI(int s, double epsilon, double& time, int batch_size, double rd_ratio, double l1_error) {
        clock_t t0 = clock();
        cout << "node: " << s << endl;
        double* resultList = new double[vert];
        double* former_resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        if(g.getOutSize(s) == 0){
            resultList[s] = 1;
            return resultList;
        }
        double walk_num = vert*log(vert)/100;
        int K = 1000;
        int walk_num_each_batch = walk_num / batch_size;

        double random_walk_this_round;
        while(batch_size>0) {
            clock_t t_build_start = clock();
            double* res = new double[vert]();
            for(int i=0; i<vert; i++) {
                if(resultList[i] != 0) {
                    if(g.getOutSize(i) == 0) {
                        res[i] += (1-alpha)/alpha*resultList[i];
                    } else {
                        for(int j=0; j<g.getOutSize(i); j++) {
                            int tmp_node = g.getOutVert(i,j);
                            res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                res[i] -= resultList[i]/alpha;
            }
            res[s] += 1.;
            clock_t t_build_end = clock();
            cout << "compute r time: " << (t_build_end - t_build_start) / (double) CLOCKS_PER_SEC << endl;
            clock_t t_alias_start = clock();
            double r_sum = 0.;
            double max_r = 0.;
            vector<pair<int, double> > aliasP;
            
            for(int i = 0; i < vert; i++){
                if(res[i] != 0) {
                    r_sum += abs(res[i]);
                    aliasP.push_back(pair<int, double>(i, abs(res[i])));
                    if(abs(res[i])>max_r) max_r = abs(res[i]);
                }
            }
            cout << "rsum: " << r_sum << endl;
            cout << "rmax: " << max_r << endl;
            Alias alias = Alias(aliasP);
            clock_t t_alias_end = clock();
            cout << "bulid alias time: " << (t_alias_end - t_alias_start) / (double) CLOCKS_PER_SEC << endl;

            clock_t t_walk_start = clock();
            for(int i = 0; i < walk_num_each_batch; i++) {
                int tempNode = alias.generateRandom(R);
                int sourceNode = tempNode;
                if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                while(R.drand() > alpha){
                    int length = g.getOutSize(tempNode);
                    if(length == 0){                    
                        tempNode = tempNode;                   
                    }
                    else{
                        int tempIndex = R.generateRandom() % length;
                        tempNode = g.getOutVert(tempNode, tempIndex);
                    }
                    if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                    if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                }
                // if(res[sourceNode]>0) resultList[tempNode] += r_sum / walk_num_each_batch;
                // if(res[sourceNode]<0) resultList[tempNode] -= r_sum / walk_num_each_batch;
            }
            clock_t t_walk_end = clock();
            cout << "random walk time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;
            cout << "random walk number: " << walk_num_each_batch << endl;

            random_walk_this_round = t_walk_end - t_walk_start;

            double threshold0 = l1_error / g.m / 2;
            // double* res = new double[vert]();
            for(int i=0; i<vert; i++) {
                if(resultList[i] != 0) {
                    if(g.getOutSize(i) == 0) {
                        res[i] += (1-alpha)/alpha*resultList[i];
                    } else {
                        for(int j=0; j<g.getOutSize(i); j++) {
                            int tmp_node = g.getOutVert(i,j);
                            res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                res[i] -= resultList[i]/alpha;
            }
            res[s] += 1.;
            queue<int> active_vertices;
            bool* is_active = new bool[vert];
            for(int i=0; i<vert; i++) is_active[i] = false;
            int queue_size_threshold = vert / 4;
            for(int i=0; i<vert; i++) {
                if(abs(res[i])>g.getOutSize(i)*threshold0) {
                    active_vertices.push(i);
                    is_active[i] = true;
                }
            }
            cout << "queue size: " << active_vertices.size() << endl;
            cout << "n/4: " << queue_size_threshold << endl;
            while(active_vertices.size()<queue_size_threshold) {
                int current_node = active_vertices.front();
                active_vertices.pop();
                is_active[current_node] = false;
                if(abs(res[current_node])>g.getOutSize(current_node)*threshold0) {
                    resultList[current_node] += alpha*res[current_node];
                    for(int i=0; i<g.getOutSize(current_node); i++) {
                        int updateNode = g.getOutVert(current_node, i);
                        res[updateNode] += (1-alpha)*res[i]/g.getOutSize(current_node);
                        if(abs(res[updateNode])>g.getOutSize(updateNode)*threshold0 && is_active[updateNode]) {
                            active_vertices.push(updateNode);
                            is_active[updateNode] = true;
                        }
                    }
                }
            }

            clock_t t_power_start = clock();
            int num_epoch = 8;
            int epoch = 0;
            while(epoch<num_epoch) {
                double l1_error_this_epoch = pow(l1_error, (1.0 + epoch) / num_epoch);
                const double threshold = l1_error_this_epoch / g.m / 2;
                for(int k=0; k<K; k++) {
                    for(int i=0; i<vert; i++) res[i] = 0.;
                    for(int i=0; i<vert; i++) {
                        if(abs(resultList[i])>threshold*g.getOutSize(i)) {
                            if(g.getOutSize(i) == 0) {
                                res[i] += (1-alpha)/alpha*resultList[i];
                            } else {
                                for(int j=0; j<g.getOutSize(i); j++) {
                                    int tmp_node = g.getOutVert(i,j);
                                    res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                                }
                            }
                        }
                    }
                    for(int i=0; i<vert; i++) {
                        res[i] -= resultList[i]/alpha;
                    }
                    res[s] += 1.;
                    r_sum = 0.;
                    for(int i = 0; i < vert; i++){
                        if(res[i] != 0) {
                            r_sum += abs(res[i]);
                        }
                    }
                    if(r_sum<l1_error_this_epoch) {
                        cout << "number of iterations: " << k << endl;
                        cout << "l1_error is below: " << l1_error_this_epoch << endl;
                        break;
                    }
                    for(int i=0; i<vert; i++) {
                        resultList[i] += alpha*res[i];
                    }
                    clock_t t_power_current = clock();
                }
                epoch++;
            }
            clock_t t_power_end = clock();
            cout << "power time: " << (t_power_end - t_power_start) / (double) CLOCKS_PER_SEC << endl;

            batch_size--;
            delete[] res;
        }

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* DPI_naive(int s, double epsilon, double& time, int batch_size, double rd_ratio, double l1_error) {
        clock_t t0 = clock();
        cout << "node: " << s << endl;
        double* resultList = new double[vert];
        double* former_resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        if(g.getOutSize(s) == 0){
            resultList[s] = 1;
            return resultList;
        }
        double walk_num = vert*log(vert)/100;
        int K = 1000;
        int walk_num_each_batch = walk_num / batch_size;

        double random_walk_this_round;
        while(batch_size>0) {
            clock_t t_build_start = clock();
            double* res = new double[vert]();
            for(int i=0; i<vert; i++) {
                if(resultList[i] != 0) {
                    if(g.getOutSize(i) == 0) {
                        res[i] += (1-alpha)/alpha*resultList[i];
                    } else {
                        for(int j=0; j<g.getOutSize(i); j++) {
                            int tmp_node = g.getOutVert(i,j);
                            res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                res[i] -= resultList[i]/alpha;
            }
            res[s] += 1.;
            clock_t t_build_end = clock();
            cout << "compute r time: " << (t_build_end - t_build_start) / (double) CLOCKS_PER_SEC << endl;
            clock_t t_alias_start = clock();
            double r_sum = 0.;
            double max_r = 0.;
            vector<pair<int, double> > aliasP;
            
            for(int i = 0; i < vert; i++){
                if(res[i] != 0) {
                    r_sum += abs(res[i]);
                    aliasP.push_back(pair<int, double>(i, abs(res[i])));
                    if(abs(res[i])>max_r) max_r = abs(res[i]);
                }
            }
            cout << "rsum: " << r_sum << endl;
            cout << "rmax: " << max_r << endl;
            Alias alias = Alias(aliasP);
            clock_t t_alias_end = clock();
            cout << "bulid alias time: " << (t_alias_end - t_alias_start) / (double) CLOCKS_PER_SEC << endl;

            clock_t t_walk_start = clock();
            for(int i = 0; i < walk_num_each_batch; i++) {
                int tempNode = alias.generateRandom(R);
                int sourceNode = tempNode;
                if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                while(R.drand() > alpha){
                    int length = g.getOutSize(tempNode);
                    if(length == 0){                    
                        tempNode = tempNode;                   
                    }
                    else{
                        int tempIndex = R.generateRandom() % length;
                        tempNode = g.getOutVert(tempNode, tempIndex);
                    }
                    if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                    if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                }
                // if(res[sourceNode]>0) resultList[tempNode] += r_sum / walk_num_each_batch;
                // if(res[sourceNode]<0) resultList[tempNode] -= r_sum / walk_num_each_batch;
            }
            clock_t t_walk_end = clock();
            cout << "random walk time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;
            cout << "random walk number: " << walk_num_each_batch << endl;

            random_walk_this_round = t_walk_end - t_walk_start;

            clock_t t_power_start = clock();
            for(int k=0; k<K; k++) {
                for(int i=0; i<vert; i++) res[i] = 0.;
                for(int i=0; i<vert; i++) {
                    if(resultList[i] != 0) {
                        if(g.getOutSize(i) == 0) {
                            res[i] += (1-alpha)/alpha*resultList[i];
                        } else {
                            for(int j=0; j<g.getOutSize(i); j++) {
                                int tmp_node = g.getOutVert(i,j);
                                res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                            }
                        }
                    }
                }
                for(int i=0; i<vert; i++) {
                    res[i] -= resultList[i]/alpha;
                }
                res[s] += 1.;
                r_sum = 0.;
                for(int i = 0; i < vert; i++){
                    if(res[i] != 0) {
                        r_sum += abs(res[i]);
                    }
                }
                if(r_sum<l1_error) {
                    cout << "number of iterations: " << k << endl;
                    cout << "l1_error is below: " << l1_error << endl;
                    break;
                }
                for(int i=0; i<vert; i++) {
                    former_resultList[i] = resultList[i];
                    resultList[i] = 0;
                }
                for(int i=0; i<vert; i++) {
                    if(former_resultList[i]>0) {
                        int length = g.getOutSize(i);
                        if(length == 0) {
                            resultList[i] += (1-alpha)*former_resultList[i];
                        } else {
                            for(int j=0; j<length; j++) {
                                int tmp_node = g.getOutVert(i,j);
                                resultList[tmp_node] += (1-alpha)*former_resultList[i] / length;
                            }
                        }
                    }
                }
                resultList[s] += alpha;
                clock_t t_power_current = clock();
            }
            clock_t t_power_end = clock();
            cout << "power time: " << (t_power_end - t_power_start) / (double) CLOCKS_PER_SEC << endl;

            batch_size--;
            delete[] res;
        }

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    void compare_walk_and_forest_variance(string filename) {
        int zero_out_neighbor_count = 0;
        for(int i=0; i<vert; i++) {
            if(g.getOutSize(i)==0) zero_out_neighbor_count++;
        }
        cout << "0 outneighbor node count: " << zero_out_neighbor_count << endl;

        // compute co-root probability
        int num = 1000;
        double average_probability_sum = 0.;
        double forest_time = 0.;
        for(int i=0; i<num; i++) {
            clock_t t_start = clock();
            int* root = loop_erased_walk();
            clock_t t_end = clock();
            forest_time += t_end - t_start;
            int* component_size = new int[vert]();
            long probability_sum = 0;
            for(int id=0; id<vert; id++) {
                component_size[root[id]] += 1;
            }
            for(int id=0; id<vert; id++) {
                probability_sum += component_size[root[id]];
            }
            // cout << "probability sum: " << probability_sum << endl;
            average_probability_sum += probability_sum / (double)num;
            delete[] root;
            delete[] component_size;
        }
        
        cout << "average probability sum: " << average_probability_sum << endl;
        cout << "average probability sum / n: " << average_probability_sum / vert << endl;
        forest_time = forest_time / (double) CLOCKS_PER_SEC / (double) num;
        cout << "avarage forest time: " << forest_time << endl;
        walk_forest_time_ratio = average_walk_time / forest_time;
        cout << "t_walk / t_forest: " << average_walk_time / forest_time << endl;

        stringstream ssout;
        ssout << "ppr-answer/" << filename << "wf-ratio-" << alpha << ".txt";
        string outputFile = ssout.str();
        ofstream fout(outputFile);
        fout << walk_forest_time_ratio << endl;
        fout.close();
    }

    void compare_walk_and_forest_time() {
        double walk_num = 1000000.;
        clock_t t_start = clock();
        for(double i = 0; i < walk_num; i++){
            int tempNode = R.generateRandom() % vert;
            while(R.drand() > alpha){
                int length = g.getOutSize(tempNode);
                if(length == 0)
                    tempNode = tempNode;
                else{
                    int r = R.generateRandom() % length;
                    tempNode = g.getOutVert(tempNode, r);
                }
            }
        }
        clock_t t_end = clock();
        average_walk_time = (t_end - t_start) / (double) CLOCKS_PER_SEC / (double) walk_num * vert;
        cout << "average walk time: " << average_walk_time << endl;
    }

    void norm_of_pi() {
        double* gt = new double[vert]();
        stringstream ss;
        ss << "ppr-answer/" << target_filename << "/" << "PageRankCentrality-" << alpha << ".txt";
        string inputFile = ss.str();
        ifstream fin(inputFile);
        double real_sum = 0.;
        for(int i=0; i<vert; i++) {
            fin >> gt[i];
        }
        fin.close();
        double norm = 0;
        for(int i=0; i<vert; i++) {
            norm += gt[i]*gt[i];
        }
        cout << "dataset: " << target_filename << " norm_of_pi: " << norm << " n*norm_of_pi: " << norm*vert << endl;
    }

    void compare_l1_error(string graphd) {
        double* gt = new double[vert]();
        stringstream ss;
        ss << "ppr-answer/" << target_filename << "/" << "PageRankCentrality-" << alpha << ".txt";
        string inputFile = ss.str();
        ifstream fin(inputFile);
        double real_sum = 0.;
        for(int i=0; i<vert; i++) {
            fin >> gt[i];
            real_sum += gt[i];
        }
        fin.close();
        cout << "real sum: " << real_sum << endl;

        double walk_num = vert*log(vert);
        cout << "number of walks: " << walk_num << endl;
        double walk_time;
        double* resultList = new double[vert];

        resultList = pc_naive_monte_carlo(walk_num, walk_time);
        
        double mc_max_err_walk = 0.;
        double mc_l1_err_walk = 0.;

        for(int i = 0; i<vert; i++){
            double tmp_err;
            tmp_err = abs(gt[i] - resultList[i]);
            mc_l1_err_walk += tmp_err;
            if(mc_max_err_walk<tmp_err) mc_max_err_walk = tmp_err;
        }
        cout << "walk max error: " << mc_max_err_walk << endl;
        cout << "walk l1 error: " << mc_l1_err_walk << endl;

        double mc_max_err_forest = 0.;
        double mc_l1_err_forest = 0.;

        double forest_num = walk_num / (double)vert * walk_forest_time_ratio;
        cout << "number of forests: " << forest_num << endl;

        resultList = pc_naive_forest(forest_num, walk_time, graphd);

        for(int i = 0; i<vert; i++){
            double tmp_err;
            tmp_err = abs(gt[i] - resultList[i]);
            mc_l1_err_forest += tmp_err;
            if(mc_max_err_forest<tmp_err) mc_max_err_forest = tmp_err;
        }
        cout << "forest max error: " << mc_max_err_forest << endl;
        cout << "forest l1 error: " << mc_l1_err_forest << endl;
    }

    void single_source(int s, string method, string graphd, int power_iterations, int batch_size, double rd_ratio, double epsilon, double l1_error) {
        vector<int> realList = getRealTopK(s, vert);
        unordered_map<int, double> realMap = getRealTopKMap(s, vert);
        double* resultList = new double[vert];

        double mc_time;
        if(method == "naive") {
            resultList = naive_monte_carlo(s, epsilon, mc_time);
        } else if(method == "PW") {
            resultList = PPW(s, epsilon, mc_time, 1, rd_ratio);
        } else if(method == "PPW") {
            resultList = PPW(s, epsilon, mc_time, batch_size, rd_ratio);
        } else if(method == "PF") {
            resultList = PPF(s, epsilon, mc_time, graphd, 1, rd_ratio);
        } else if(method == "PPF") {
            resultList = PPF(s, epsilon, mc_time, graphd, batch_size, rd_ratio);
        } else if(method == "DPI") {
            resultList = DPI(s, epsilon, mc_time, batch_size, rd_ratio, l1_error);
        } else if(method == "DPI_naive") {
            resultList = DPI_naive(s, epsilon, mc_time, batch_size, rd_ratio, l1_error);
        }

        double mc_max_err = 0.;
        double mc_l1_err = 0.;
        double real_sum = 0.;

        for(int i = 0; i<vert; i++){
            double tmp_err;
            tmp_err = abs(realMap[i] - resultList[i]);
            mc_l1_err += tmp_err;
            if(mc_max_err<tmp_err) mc_max_err = tmp_err;
            real_sum += resultList[i];
        }

        cout << "max error: " << mc_max_err << endl;
        cout << "l1 error: " << mc_l1_err << endl;
        cout << "real sum: " << real_sum << endl;

        avg_mc_time +=  mc_time;
        avg_mc_max_err += mc_max_err;
        avg_mc_l1_err += mc_l1_err;
    }

    void pagerank_centrality(string filename, string algo, string method, string graphd, int power_iterations, int batch_size, double rd_ratio, double epsilon) {
        double* gt = new double[vert]();
        stringstream ss;
        ss << "ppr-answer/" << target_filename << "/" << "PageRankCentrality-" << alpha << ".txt";
        string inputFile = ss.str();
        ifstream fin(inputFile);
        for(int i=0; i<vert; i++) {
            fin >> gt[i];
        }
        fin.close();

        double* resultList = new double[vert];

        double mc_time;
        if(method == "MCW") {
            resultList = MCW(mc_time, epsilon);
        } else if(method == "MCF") {
            resultList = MCF(mc_time, epsilon, graphd);
        } else if(method == "FORA") {
            resultList = FORA(epsilon, mc_time);
        } else if(method == "PPW") {
            resultList = pc_PPW(epsilon, mc_time, batch_size, rd_ratio);
        } else if(method == "PPF") {
            resultList = pc_PPF(epsilon, mc_time, graphd, batch_size, rd_ratio);
        } else if(method == "PW") {
            resultList = pc_PPW(epsilon, mc_time, 1, rd_ratio);
        } else if(method == "PF") {
            resultList = pc_PPF(epsilon, mc_time, graphd, 1, rd_ratio);
        }

        double mc_max_err = 0.;
        double mc_l1_err = 0.;
        double real_sum = 0.;

        for(int i = 0; i<vert; i++){
            double tmp_err;
            tmp_err = abs(gt[i] - resultList[i]);
            mc_l1_err += tmp_err;
            if(mc_max_err<tmp_err) mc_max_err = tmp_err;
            real_sum += resultList[i];
        }
        cout << "max error: " << mc_max_err << endl;
        cout << "l1 error: " << mc_l1_err << endl;
        cout << "real sum: " << real_sum << endl;

        stringstream ssout;
        ssout << filename << "-" << method << "-" << algo << alpha << ".txt";
        string outputFile = ssout.str();
        ofstream fout(outputFile, ios::app);
        fout << mc_time << " " << mc_l1_err << endl;
        fout.close();
    }

    void pagerank_centrality_ground_truth(int K) {
        clock_t t_power_start = clock();
        double* former_resultList = new double[vert]();
        for(int i=0; i<vert; i++) {
            resultList[i] = 1./vert;
        }
        for(int k=0; k<K; k++) {
            for(int i=0; i<vert; i++) {
                former_resultList[i] = resultList[i];
                resultList[i] = 0;
            }
            for(int i=0; i<vert; i++) {
                if(former_resultList[i]>0) {
                    int length = g.getOutSize(i);
                    if(length == 0) {
                        resultList[i] += (1-alpha)*former_resultList[i];
                    } else {
                        for(int j=0; j<length; j++) {
                            int tmp_node = g.getOutVert(i,j);
                            resultList[tmp_node] += (1-alpha)*former_resultList[i] / length;
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                resultList[i] += alpha/vert;
            }
        }
        clock_t t_power_end = clock();
        cout << "PageRank Centrality ground-truth time: " << (t_power_end-t_power_start) / (double) CLOCKS_PER_SEC << endl;

        stringstream ss;
        ss << "ppr-answer/" << target_filename << "/" << "PageRankCentrality-" << alpha << ".txt";
        string outputFile = ss.str();
        ofstream fout(outputFile);
        fout.precision(17);
        for(int i=0; i<vert; i++) {
            fout << resultList[i] << endl;
        }
        fout.close();
    }

    void any_distribution_ground_truth(double* source_distribution, int t, int K) {
        clock_t t_power_start = clock();
        double* former_resultList = new double[vert]();
        for(int i=0; i<vert; i++) {
            resultList[i] = source_distribution[i];
        }
        for(int k=0; k<K; k++) {
            for(int i=0; i<vert; i++) {
                former_resultList[i] = resultList[i];
                resultList[i] = 0;
            }
            for(int i=0; i<vert; i++) {
                if(former_resultList[i]>0) {
                    int length = g.getOutSize(i);
                    if(length == 0) {
                        resultList[i] += (1-alpha)*former_resultList[i];
                    } else {
                        for(int j=0; j<length; j++) {
                            int tmp_node = g.getOutVert(i,j);
                            resultList[tmp_node] += (1-alpha)*former_resultList[i] / length;
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                resultList[i] += alpha*source_distribution[i];
            }
        }
        clock_t t_power_end = clock();
        cout << "any_distribution-ground-truth time: " << (t_power_end-t_power_start) / (double) CLOCKS_PER_SEC << endl;

        stringstream ss;
        ss << "ppr-answer/" << target_filename << "/" << t << "-" << alpha << ".source-txt";
        string outputFile = ss.str();
        ofstream fout(outputFile);
        fout.precision(17);
        for(int i=0; i<vert; i++) {
            fout << resultList[i] << endl;
        }
        fout.close();
    }

    void t_PowerMethod(vector<int> nodeList, int iterations){
        for(int i = 0; i < nodeList.size(); i++){
            int tempNode = nodeList[i];
            stringstream ss;
            ss << "ppr-answer/" << target_filename << "/" << tempNode << "-" << alpha << ".txt";
            string outputFile = ss.str();
            cout << "file: " << outputFile << endl;
            PowerMethodK(iterations, outputFile, tempNode, 500);
            cout << outputFile << "done!"  << endl;        
        }
    }

    void PowerMethodK(int iterations, string outputFile, int u, int k){
        unordered_map<int, double> map_residual;
        map_residual.clear();
        map_residual[u] = 1.0;

        int num_iter=0;
        double* map_ppr = new double[vert];
        for(int i = 0; i < vert; i++){
            map_ppr[i] = 0;
        }
        while( num_iter < iterations ){
            //cout << u << ": iter " << num_iter << endl;
            num_iter++;

            vector< pair<int,double> > pairs(map_residual.begin(), map_residual.end());
            map_residual.clear();
            for(auto &p: pairs){
                if(p.second > 0){
                    map_ppr[p.first] += alpha*p.second;
                    int out_deg = g.getOutSize(p.first);

                    double remain_residual = (1-alpha)*p.second;
                    if(out_deg==0){
                        map_residual[p.first] += remain_residual;
                    }
                    else{
                        double avg_push_residual = remain_residual / out_deg;
                        for(int i = 0; i < g.getOutSize(p.first); i++){
                            int next = g.getOutVert(p.first, i);
                            map_residual[next] += avg_push_residual;
                        }
                    }
                }
            }
        }
        ofstream fout(outputFile);
        vector<pair<int, double> > pprs;
        for(int j = 0; j < vert; j++){
            pprs.push_back(pair<int, double>(j, map_ppr[j]));
        }
        sort(pprs.begin(), pprs.end(), maxScoreCmp);
        for(int j = 0; j < vert; j++){
            if(pprs[j].second >= 0){
                fout << setprecision(16) << pprs[j].first << " " << setprecision(16) << pprs[j].second << "\n";
            }
        }
        fout.close();
        delete[] map_ppr;
    }

    //generate random query node
    void generateQueryNode(int nodeNum, ofstream& fout){
        for(int i = 0; i < nodeNum; i++){
            int tempNode = R.generateRandom() % vert;
            if(g.getOutSize(tempNode) == 0){
                i--;
                continue;   
            }
            fout << tempNode << endl;
        }
    }

    void generateSourceDistribution(ofstream& fout) {
        double* source_dis = new double[vert];
        double sum = 0.;
        for(int i=0; i<vert; i++) {
            source_dis[i] = R.drand_t();
            sum += source_dis[i];
        }
        for(int i=0; i<vert; i++) {
            source_dis[i] /= sum;
        }
        fout.precision(17);
        for(int i=0; i<vert; i++) {
            fout << source_dis[i] << endl;
        }
    }

    double* MCW(double& time, double epsilon) {
        cout << "naive monte carlo: PageRank Centrality" << endl;
        double* resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        double walk_num = vert*log(vert)/epsilon/epsilon;
        clock_t t0 = clock();
        for(double i = 0; i < walk_num; i++){
            int tempNode = R.generateRandom() % vert;
            while(R.drand() > alpha){
                int length = g.getOutSize(tempNode);
                if(length == 0)
                    tempNode = tempNode;
                else{
                    int r = R.generateRandom() % length;
                    tempNode = g.getOutVert(tempNode, r);
                }
            }
            resultList[tempNode] += 1;
        }
        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        for(int i = 0; i < vert; i++){
            resultList[i] /= (double) walk_num;
        }
        return resultList;
    }

    double* MCF(double& time, double epsilon, string graphd) {
        cout << "naive monte carlo forest: PageRank Centrality" << endl;
        double* resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        double* res = new double[vert]();
        for(int i=0; i<vert; i++) res[i] = 1./vert;

        clock_t t0 = clock();
        double forest_num = log(vert)/epsilon/epsilon*walk_forest_time_ratio;
        int forest_num_int = ceil(forest_num);
        for(int i=0; i<forest_num_int; i++) {
            int* root = loop_erased_walk();
            if(graphd == "directed") {
                for(int id=0; id<vert; id++) {
                    resultList[root[id]] += res[id]/(double)forest_num_int; 
                }
            } else if(graphd == "undirected") {
                double* partition_acc = new double[vert]();
                double* residual_over_partition = new double[vert]();
                for(int id=0; id<vert; id++) {
                    partition_acc[root[id]] += g.getOutSize(id); 
                    residual_over_partition[root[id]] += res[id];
                }
                for(int id=0; id<vert; id++) {
                    resultList[id] += g.getOutSize(id)*residual_over_partition[root[id]]/partition_acc[root[id]]/(double)forest_num_int;
                }
                delete[] partition_acc;
                delete[] residual_over_partition;
            }
            delete[] root;
        }
        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* pc_PPW(double epsilon, double& time, int batch_size, double rd_ratio) {
        clock_t t0 = clock();
        cout << "pc: PPW"<< endl;
        double* resultList = new double[vert];
        double* former_resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        double walk_num = alpha*log(epsilon*epsilon)/log(1-alpha)*vert*log(vert)/100;
        // double walk_num = vert*log(vert)/epsilon/100;
        int K = 1000;
        int walk_num_each_batch = walk_num / batch_size;
        
        double random_walk_this_round;
        while(batch_size>0) {
            clock_t t_build_start = clock();
            double* res = new double[vert]();
            for(int i=0; i<vert; i++) {
                if(resultList[i] != 0) {
                    if(g.getOutSize(i) == 0) {
                        res[i] += (1-alpha)/alpha*resultList[i];
                    } else {
                        for(int j=0; j<g.getOutSize(i); j++) {
                            int tmp_node = g.getOutVert(i,j);
                            res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                res[i] -= resultList[i]/alpha;
            }
            for(int i=0; i<vert; i++) {
                res[i] += 1./vert;
            }
            clock_t t_build_end = clock();
            cout << "compute r time: " << (t_build_end - t_build_start) / (double) CLOCKS_PER_SEC << endl;
            clock_t t_alias_start = clock();
            double r_sum = 0.;
            double max_r = 0.;
            vector<pair<int, double> > aliasP;
            
            for(int i = 0; i < vert; i++){
                if(res[i] != 0) {
                    r_sum += abs(res[i]);
                    aliasP.push_back(pair<int, double>(i, abs(res[i])));
                    if(abs(res[i])>max_r) max_r = abs(res[i]);
                }
            }
            cout << "rsum: " << r_sum << endl;
            cout << "rmax: " << max_r << endl;
            Alias alias = Alias(aliasP);
            clock_t t_alias_end = clock();
            cout << "bulid alias time: " << (t_alias_end - t_alias_start) / (double) CLOCKS_PER_SEC << endl;

            clock_t t_walk_start = clock();
            for(int i = 0; i < walk_num_each_batch; i++) {
                int tempNode = alias.generateRandom(R);
                int sourceNode = tempNode;
                if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                while(R.drand() > alpha){
                    int length = g.getOutSize(tempNode);
                    if(length == 0){                    
                        tempNode = tempNode;                   
                    }
                    else{
                        int tempIndex = R.generateRandom() % length;
                        tempNode = g.getOutVert(tempNode, tempIndex);
                    }
                    if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                    if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                }
                // if(res[sourceNode]>0) resultList[tempNode] += r_sum / walk_num_each_batch;
                // if(res[sourceNode]<0) resultList[tempNode] -= r_sum / walk_num_each_batch;
            }
            clock_t t_walk_end = clock();
            cout << "random walk time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;
            cout << "random walk number: " << walk_num_each_batch << endl;

            random_walk_this_round = t_walk_end - t_walk_start;

            clock_t t_power_start = clock();
            for(int k=0; k<K; k++) {
                for(int i=0; i<vert; i++) {
                    former_resultList[i] = resultList[i];
                    resultList[i] = 0;
                }
                for(int i=0; i<vert; i++) {
                    if(former_resultList[i]>0) {
                        int length = g.getOutSize(i);
                        if(length == 0) {
                            resultList[i] += (1-alpha)*former_resultList[i];
                        } else {
                            for(int j=0; j<length; j++) {
                                int tmp_node = g.getOutVert(i,j);
                                resultList[tmp_node] += (1-alpha)*former_resultList[i] / length;
                            }
                        }
                    }
                }
                for(int i=0; i<vert; i++) {
                    resultList[i] += alpha/vert;
                }
                clock_t t_power_current = clock();
                double current_time = t_power_current - t_power_start;
                if(current_time > random_walk_this_round*rd_ratio) {
                    cout << "K: " << k << " break: power too much" << endl;
                    break;
                }
            }
            clock_t t_power_end = clock();
            cout << "power time: " << (t_power_end - t_power_start) / (double) CLOCKS_PER_SEC << endl;

            batch_size--;
            delete[] res;
        }

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* pc_PPF(double epsilon, double& time, string graphd, int batch_size, double rd_ratio) {
        clock_t t0 = clock();
        cout << "pc: PPF"<< endl;
        double* resultList = new double[vert];
        double* former_resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        int K = 1000;
        double num_forests = alpha*log(epsilon*epsilon)/log(1-alpha)*log(vert)/100*walk_forest_time_ratio;
        // double num_forests = log(vert)/epsilon/100*walk_forest_time_ratio;
        double num_forest_each_batch = num_forests / batch_size;

        double random_walk_this_round;

        while(batch_size>0) {
            clock_t t_build_start = clock();
            double* res = new double[vert]();
            for(int i=0; i<vert; i++) {
                if(resultList[i] != 0) {
                    if(g.getOutSize(i) == 0) {
                        res[i] += (1-alpha)/alpha*resultList[i];
                    } else {
                        for(int j=0; j<g.getOutSize(i); j++) {
                            int tmp_node = g.getOutVert(i,j);
                            res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                res[i] -= resultList[i]/alpha;
            }
            for(int i=0; i<vert; i++) {
                res[i] += 1./vert;
            }
            clock_t t_build_end = clock();
            cout << "compute r time: " << (t_build_end - t_build_start) / (double) CLOCKS_PER_SEC << endl;
            clock_t t_alias_start = clock();
            double r_sum = 0.;
            double max_r = 0.;
            vector<pair<int, double> > aliasP;
            
            for(int i = 0; i < vert; i++){
                if(res[i] != 0) {
                    r_sum += abs(res[i]);
                    aliasP.push_back(pair<int, double>(i, abs(res[i])));
                    if(abs(res[i])>max_r) max_r = abs(res[i]);
                }
            }
            cout << "rsum: " << r_sum << endl;
            cout << "rmax: " << max_r << endl;
            Alias alias = Alias(aliasP);
            clock_t t_alias_end = clock();
            cout << "bulid alias time: " << (t_alias_end - t_alias_start) / (double) CLOCKS_PER_SEC << endl;

            clock_t t_walk_start = clock();
            cout << "sample rooted spanning forests: " << num_forest_each_batch << endl;
            for(int i=0; i<num_forest_each_batch; i++) {
                int* root = loop_erased_walk();
                if(graphd == "directed") {
                    for(int id=0; id<vert; id++) { 
                        resultList[root[id]] += res[id]/(double)num_forest_each_batch; 
                    }
                } else if(graphd == "undirected") {
                    double* partition_acc = new double[vert]();
                    double* residual_over_partition = new double[vert]();
                    for(int id=0; id<vert; id++) {
                        partition_acc[root[id]] += g.getOutSize(id); 
                        residual_over_partition[root[id]] += res[id];
                    }
                    for(int id=0; id<vert; id++) {
                        resultList[id] += g.getOutSize(id)*residual_over_partition[root[id]]/partition_acc[root[id]]/(double)num_forest_each_batch;
                    }
                    delete[] partition_acc;
                    delete[] residual_over_partition;
                }
                delete[] root;
            }
            clock_t t_walk_end = clock();
            cout << "sampling spanning forests time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;

            random_walk_this_round = t_walk_end - t_walk_start;

            clock_t t_power_start = clock();
            for(int k=0; k<K; k++) {
                for(int i=0; i<vert; i++) {
                    former_resultList[i] = resultList[i];
                    resultList[i] = 0;
                }
                for(int i=0; i<vert; i++) {
                    if(former_resultList[i]>0) {
                        int length = g.getOutSize(i);
                        if(length == 0) {
                            resultList[i] += (1-alpha)*former_resultList[i];
                        } else {
                            for(int j=0; j<length; j++) {
                                int tmp_node = g.getOutVert(i,j);
                                resultList[tmp_node] += (1-alpha)*former_resultList[i] / length;
                            }
                        }
                    }
                }
                for(int i=0; i<vert; i++) {
                    resultList[i] += alpha/vert;
                }
                clock_t t_power_current = clock();
                double current_time = t_power_current - t_power_start;
                if(current_time > random_walk_this_round*rd_ratio) {
                    cout << "K: " << k << " break: power too much" << endl;
                    break;
                }
            }
            clock_t t_power_end = clock();
            cout << "power time: " << (t_power_end - t_power_start) / (double) CLOCKS_PER_SEC << endl;

            batch_size--;
            delete[] res;
        }

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* pc_naive_monte_carlo(double walk_num, double& time) {
        cout << "naive monte carlo: PageRank Centrality" << endl;
        double* resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        clock_t t0 = clock();
        for(double i = 0; i < walk_num; i++){
            int tempNode = R.generateRandom() % vert;
            while(R.drand() > alpha){
                int length = g.getOutSize(tempNode);
                if(length == 0)
                    tempNode = tempNode;
                else{
                    int r = R.generateRandom() % length;
                    tempNode = g.getOutVert(tempNode, r);
                }
            }
            resultList[tempNode] += 1;
        }
        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        for(int i = 0; i < vert; i++){
            resultList[i] /= (double) walk_num;
        }
        return resultList;
    }

    double* pc_naive_forest(int forest_num, double& time, string graphd) {
        cout << "naive monte carlo forest: PageRank Centrality" << endl;
        double* resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        double* res = new double[vert]();
        for(int i=0; i<vert; i++) res[i] = 1./vert;

        clock_t t0 = clock();
        int forest_num_int = ceil(forest_num);
        for(int i=0; i<forest_num_int; i++) {
            int* root = loop_erased_walk();
            if(graphd == "directed") {
                for(int id=0; id<vert; id++) {
                    resultList[root[id]] += res[id]/(double)forest_num_int; 
                }
            } else if(graphd == "undirected") {
                double* partition_acc = new double[vert]();
                double* residual_over_partition = new double[vert]();
                for(int id=0; id<vert; id++) {
                    partition_acc[root[id]] += g.getOutSize(id); 
                    residual_over_partition[root[id]] += res[id];
                }
                for(int id=0; id<vert; id++) {
                    resultList[id] += g.getOutSize(id)*residual_over_partition[root[id]]/partition_acc[root[id]]/(double)forest_num_int;
                }
                delete[] partition_acc;
                delete[] residual_over_partition;
            }
            delete[] root;
        }
        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* FORA(double epsilon, double& time) {
        clock_t t0 = clock();
        cout << "FORA: pagerank centrality" << endl;
        double* resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;

        clock_t t_push_start = clock();
        double* res = new double[vert];
        bool* isInQueue = new bool[vert];
        queue<int> r_queue;
        double rmax = epsilon / vert;
        for(int i=0; i<vert; i++) {
            res[i] = 1./vert;
            r_queue.push(i);
            isInQueue[i] = true;
        }
        double r_sum = 1.0;
        while(r_queue.size() > 0){
            int tempNode = r_queue.front();
            r_queue.pop();
            isInQueue[tempNode] = false;
            int tempOutSize = g.getOutSize(tempNode);
            double tempRes = res[tempNode];
            resultList[tempNode] += alpha * tempRes;
            r_sum -= alpha * tempRes;
            // cout << r_sum << endl;
            res[tempNode] = 0;
            if(tempOutSize == 0) {
                int updateNode = tempNode;
                res[updateNode] += (1-alpha)*tempRes;
                if(!isInQueue[updateNode] && res[updateNode]>rmax) {
                    r_queue.push(updateNode);
                    isInQueue[updateNode] = true;
                }
            } else {
                for(int i=0; i<tempOutSize; i++) {
                    int updateNode = g.getOutVert(tempNode, i);
                    res[updateNode] += (1-alpha)*tempRes / tempOutSize;
                    if(!isInQueue[updateNode] && res[updateNode]>rmax) {
                        r_queue.push(updateNode);
                        isInQueue[updateNode] = true;
                    }
                }
            }
        }
        clock_t t_push_end = clock();
        cout << "push time: " << (t_push_end - t_push_start) / (double) CLOCKS_PER_SEC << endl;
        cout << "current r_sum: " << r_sum << endl;

        vector<pair<int, double> > aliasP;
        for(int i = 0; i < vert; i++){
            if(res[i] != 0) {
                aliasP.push_back(pair<int, double>(i, abs(res[i])));
            }
        }
        Alias alias = Alias(aliasP);

        double walk_num = vert/epsilon;

        clock_t t_walk_start = clock();
        for(double i = 0; i < walk_num; i++) {
            int tempNode = alias.generateRandom(R);
            int sourceNode = tempNode;
            if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num;
            if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num;
            while(R.drand() > alpha){
                int length = g.getOutSize(tempNode);
                if(length == 0){                    
                    tempNode = tempNode;                   
                }
                else{
                    int tempIndex = R.generateRandom() % length;
                    tempNode = g.getOutVert(tempNode, tempIndex);
                }
                if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num;
                if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num;
            }
        }
        clock_t t_walk_end = clock();
        cout << "random walk time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    void any_distribution(double* source_distribution, int t, string method, string graphd, int power_iterations, int batch_size, double rd_ratio, double epsilon) {
        stringstream ss;
        ss << "ppr-answer/" << target_filename << "/" << t << "-" << alpha << ".source-txt";
        string inputFile = ss.str();
        double* gt = new double[vert];
        ifstream fin(inputFile);
        for(int i=0; i<vert; i++) {
            fin >> gt[i];
        }
        fin.close();

        double* resultList = new double[vert];

        double mc_time;
        if(method == "MCW") {
            resultList = dis_MCW(source_distribution, mc_time, epsilon);
        } else if(method == "MCF") {
            resultList = dis_MCF(source_distribution, mc_time, epsilon, graphd);
        } else if(method == "FORA") {
            resultList = dis_FORA(source_distribution, epsilon, mc_time);
        } else if(method == "PPW") {
            resultList = dis_PPW(source_distribution, epsilon, mc_time, batch_size, rd_ratio);
        } else if(method == "PPF") {
            resultList = dis_PPF(source_distribution, epsilon, mc_time, graphd, batch_size, rd_ratio);
        } else if(method == "PW") {
            resultList = dis_PPW(source_distribution, epsilon, mc_time, 1, rd_ratio);
        } else if(method == "PF") {
            resultList = dis_PPF(source_distribution, epsilon, mc_time, graphd, 1, rd_ratio);
        }

        double mc_max_err = 0.;
        double mc_l1_err = 0.;
        double real_sum = 0.;

        for(int i = 0; i<vert; i++){
            double tmp_err;
            tmp_err = abs(gt[i] - resultList[i]);
            mc_l1_err += tmp_err;
            if(mc_max_err<tmp_err) mc_max_err = tmp_err;
            real_sum += resultList[i];
        }

        cout << "max error: " << mc_max_err << endl;
        cout << "l1 error: " << mc_l1_err << endl;
        cout << "real sum: " << real_sum << endl;

        avg_mc_time +=  mc_time;
        avg_mc_max_err += mc_max_err;
        avg_mc_l1_err += mc_l1_err;
    }

    double* dis_MCW(double* source_dis, double& time, double epsilon) {
        cout << "naive monte carlo: any distribution" << endl;
        double* resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        double walk_num = vert*log(vert)/epsilon/epsilon;
        vector<pair<int, double> > aliasP;
        for(int i = 0; i < vert; i++){
            aliasP.push_back(pair<int, double>(i, source_dis[i]));
        }
        Alias alias = Alias(aliasP);
        clock_t t0 = clock();
        for(double i = 0; i < walk_num; i++){
            int tempNode = alias.generateRandom(R);
            while(R.drand() > alpha){
                int length = g.getOutSize(tempNode);
                if(length == 0)
                    tempNode = tempNode;
                else{
                    int r = R.generateRandom() % length;
                    tempNode = g.getOutVert(tempNode, r);
                }
            }
            resultList[tempNode] += 1;
        }
        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        for(int i = 0; i < vert; i++){
            resultList[i] /= (double) walk_num;
        }
        return resultList;
    }

    double* dis_MCF(double* source_dis, double& time, double epsilon, string graphd) {
        cout << "naive monte carlo forest: any distribution" << endl;
        double* resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        double* res = new double[vert]();
        for(int i=0; i<vert; i++) res[i] = source_dis[i];

        clock_t t0 = clock();
        double forest_num = log(vert)/epsilon/epsilon*walk_forest_time_ratio;
        int forest_num_int = ceil(forest_num);
        for(int i=0; i<forest_num_int; i++) {
            int* root = loop_erased_walk();
            if(graphd == "directed") {
                for(int id=0; id<vert; id++) {
                    resultList[root[id]] += res[id]/(double)forest_num_int; 
                }
            } else if(graphd == "undirected") {
                double* partition_acc = new double[vert]();
                double* residual_over_partition = new double[vert]();
                for(int id=0; id<vert; id++) {
                    partition_acc[root[id]] += g.getOutSize(id); 
                    residual_over_partition[root[id]] += res[id];
                }
                for(int id=0; id<vert; id++) {
                    resultList[id] += g.getOutSize(id)*residual_over_partition[root[id]]/partition_acc[root[id]]/(double)forest_num_int;
                }
                delete[] partition_acc;
                delete[] residual_over_partition;
            }
            delete[] root;
        }
        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* dis_FORA(double* source_dis, double epsilon, double& time) {
        clock_t t0 = clock();
        cout << "FORA: any distribution" << endl;
        double* resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;

        clock_t t_push_start = clock();
        double* res = new double[vert];
        bool* isInQueue = new bool[vert];
        queue<int> r_queue;
        double rmax = epsilon / vert;
        for(int i=0; i<vert; i++) {
            res[i] = source_dis[i];
            if(res[i]>rmax) {
                r_queue.push(i);
                isInQueue[i] = true;
            }
        }
        double r_sum = 1.0;
        while(r_queue.size() > 0){
            int tempNode = r_queue.front();
            r_queue.pop();
            isInQueue[tempNode] = false;
            int tempOutSize = g.getOutSize(tempNode);
            double tempRes = res[tempNode];
            resultList[tempNode] += alpha * tempRes;
            r_sum -= alpha * tempRes;
            // cout << r_sum << endl;
            res[tempNode] = 0;
            if(tempOutSize == 0) {
                int updateNode = tempNode;
                res[updateNode] += (1-alpha)*tempRes;
                if(!isInQueue[updateNode] && res[updateNode]>rmax) {
                    r_queue.push(updateNode);
                    isInQueue[updateNode] = true;
                }
            } else {
                for(int i=0; i<tempOutSize; i++) {
                    int updateNode = g.getOutVert(tempNode, i);
                    res[updateNode] += (1-alpha)*tempRes / tempOutSize;
                    if(!isInQueue[updateNode] && res[updateNode]>rmax) {
                        r_queue.push(updateNode);
                        isInQueue[updateNode] = true;
                    }
                }
            }
        }
        clock_t t_push_end = clock();
        cout << "push time: " << (t_push_end - t_push_start) / (double) CLOCKS_PER_SEC << endl;
        cout << "current r_sum: " << r_sum << endl;

        vector<pair<int, double> > aliasP;
        for(int i = 0; i < vert; i++){
            if(res[i] != 0) {
                aliasP.push_back(pair<int, double>(i, abs(res[i])));
            }
        }
        Alias alias = Alias(aliasP);

        double walk_num = vert/epsilon;

        clock_t t_walk_start = clock();
        for(double i = 0; i < walk_num; i++) {
            int tempNode = alias.generateRandom(R);
            int sourceNode = tempNode;
            if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num;
            if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num;
            while(R.drand() > alpha){
                int length = g.getOutSize(tempNode);
                if(length == 0){                    
                    tempNode = tempNode;                   
                }
                else{
                    int tempIndex = R.generateRandom() % length;
                    tempNode = g.getOutVert(tempNode, tempIndex);
                }
                if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num;
                if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num;
            }
        }
        clock_t t_walk_end = clock();
        cout << "random walk time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* dis_PPW(double* source_dis, double epsilon, double& time, int batch_size, double rd_ratio) {
        clock_t t0 = clock();
        cout << "dis: PPW"<< endl;
        double* resultList = new double[vert];
        double* former_resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        double walk_num = vert*log(vert)/epsilon/100;
        int K = 1000;
        int walk_num_each_batch = walk_num / batch_size;
        
        double random_walk_this_round;
        while(batch_size>0) {
            clock_t t_build_start = clock();
            double* res = new double[vert]();
            for(int i=0; i<vert; i++) {
                if(resultList[i] != 0) {
                    if(g.getOutSize(i) == 0) {
                        res[i] += (1-alpha)/alpha*resultList[i];
                    } else {
                        for(int j=0; j<g.getOutSize(i); j++) {
                            int tmp_node = g.getOutVert(i,j);
                            res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                res[i] -= resultList[i]/alpha;
            }
            for(int i=0; i<vert; i++) {
                res[i] += source_dis[i];
            }
            clock_t t_build_end = clock();
            cout << "compute r time: " << (t_build_end - t_build_start) / (double) CLOCKS_PER_SEC << endl;
            clock_t t_alias_start = clock();
            double r_sum = 0.;
            double max_r = 0.;
            vector<pair<int, double> > aliasP;
            
            for(int i = 0; i < vert; i++){
                if(res[i] != 0) {
                    r_sum += abs(res[i]);
                    aliasP.push_back(pair<int, double>(i, abs(res[i])));
                    if(abs(res[i])>max_r) max_r = abs(res[i]);
                }
            }
            cout << "rsum: " << r_sum << endl;
            cout << "rmax: " << max_r << endl;
            Alias alias = Alias(aliasP);
            clock_t t_alias_end = clock();
            cout << "bulid alias time: " << (t_alias_end - t_alias_start) / (double) CLOCKS_PER_SEC << endl;

            clock_t t_walk_start = clock();
            for(int i = 0; i < walk_num_each_batch; i++) {
                int tempNode = alias.generateRandom(R);
                int sourceNode = tempNode;
                if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                while(R.drand() > alpha){
                    int length = g.getOutSize(tempNode);
                    if(length == 0){                    
                        tempNode = tempNode;                   
                    }
                    else{
                        int tempIndex = R.generateRandom() % length;
                        tempNode = g.getOutVert(tempNode, tempIndex);
                    }
                    if(res[sourceNode]>0) resultList[tempNode] += alpha*r_sum / walk_num_each_batch;
                    if(res[sourceNode]<0) resultList[tempNode] -= alpha*r_sum / walk_num_each_batch;
                }
                // if(res[sourceNode]>0) resultList[tempNode] += r_sum / walk_num_each_batch;
                // if(res[sourceNode]<0) resultList[tempNode] -= r_sum / walk_num_each_batch;
            }
            clock_t t_walk_end = clock();
            cout << "random walk time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;
            cout << "random walk number: " << walk_num_each_batch << endl;

            random_walk_this_round = t_walk_end - t_walk_start;

            clock_t t_power_start = clock();
            for(int k=0; k<K; k++) {
                for(int i=0; i<vert; i++) {
                    former_resultList[i] = resultList[i];
                    resultList[i] = 0;
                }
                for(int i=0; i<vert; i++) {
                    if(former_resultList[i]>0) {
                        int length = g.getOutSize(i);
                        if(length == 0) {
                            resultList[i] += (1-alpha)*former_resultList[i];
                        } else {
                            for(int j=0; j<length; j++) {
                                int tmp_node = g.getOutVert(i,j);
                                resultList[tmp_node] += (1-alpha)*former_resultList[i] / length;
                            }
                        }
                    }
                }
                for(int i=0; i<vert; i++) {
                    resultList[i] += alpha*source_dis[i];
                }
                clock_t t_power_current = clock();
                double current_time = t_power_current - t_power_start;
                if(current_time > random_walk_this_round*rd_ratio) {
                    cout << "K: " << k << " break: power too much" << endl;
                    break;
                }
            }
            clock_t t_power_end = clock();
            cout << "power time: " << (t_power_end - t_power_start) / (double) CLOCKS_PER_SEC << endl;

            batch_size--;
            delete[] res;
        }

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }

    double* dis_PPF(double* source_dis, double epsilon, double& time, string graphd, int batch_size, double rd_ratio) {
        clock_t t0 = clock();
        cout << "dis: PPF"<< endl;
        double* resultList = new double[vert];
        double* former_resultList = new double[vert];
        for(int i = 0; i < vert; i++)
            resultList[i] = 0;
        int K = 1000;
        double num_forests = log(vert)/epsilon/100*walk_forest_time_ratio;
        double num_forest_each_batch = num_forests / batch_size;

        double random_walk_this_round;

        while(batch_size>0) {
            clock_t t_build_start = clock();
            double* res = new double[vert]();
            for(int i=0; i<vert; i++) {
                if(resultList[i] != 0) {
                    if(g.getOutSize(i) == 0) {
                        res[i] += (1-alpha)/alpha*resultList[i];
                    } else {
                        for(int j=0; j<g.getOutSize(i); j++) {
                            int tmp_node = g.getOutVert(i,j);
                            res[tmp_node] += (1-alpha)/alpha*resultList[i] / g.getOutSize(i);
                        }
                    }
                }
            }
            for(int i=0; i<vert; i++) {
                res[i] -= resultList[i]/alpha;
            }
            for(int i=0; i<vert; i++) {
                res[i] += source_dis[i];
            }
            clock_t t_build_end = clock();
            cout << "compute r time: " << (t_build_end - t_build_start) / (double) CLOCKS_PER_SEC << endl;
            clock_t t_alias_start = clock();
            double r_sum = 0.;
            double max_r = 0.;
            vector<pair<int, double> > aliasP;
            
            for(int i = 0; i < vert; i++){
                if(res[i] != 0) {
                    r_sum += abs(res[i]);
                    aliasP.push_back(pair<int, double>(i, abs(res[i])));
                    if(abs(res[i])>max_r) max_r = abs(res[i]);
                }
            }
            cout << "rsum: " << r_sum << endl;
            cout << "rmax: " << max_r << endl;
            Alias alias = Alias(aliasP);
            clock_t t_alias_end = clock();
            cout << "bulid alias time: " << (t_alias_end - t_alias_start) / (double) CLOCKS_PER_SEC << endl;

            clock_t t_walk_start = clock();
            cout << "sample rooted spanning forests: " << num_forest_each_batch << endl;
            for(int i=0; i<num_forest_each_batch; i++) {
                int* root = loop_erased_walk();
                if(graphd == "directed") {
                    for(int id=0; id<vert; id++) { 
                        resultList[root[id]] += res[id]/(double)num_forest_each_batch; 
                    }
                } else if(graphd == "undirected") {
                    double* partition_acc = new double[vert]();
                    double* residual_over_partition = new double[vert]();
                    for(int id=0; id<vert; id++) {
                        partition_acc[root[id]] += g.getOutSize(id); 
                        residual_over_partition[root[id]] += res[id];
                    }
                    for(int id=0; id<vert; id++) {
                        resultList[id] += g.getOutSize(id)*residual_over_partition[root[id]]/partition_acc[root[id]]/(double)num_forest_each_batch;
                    }
                    delete[] partition_acc;
                    delete[] residual_over_partition;
                }
                delete[] root;
            }
            clock_t t_walk_end = clock();
            cout << "sampling spanning forests time: " << (t_walk_end - t_walk_start) / (double) CLOCKS_PER_SEC << endl;

            random_walk_this_round = t_walk_end - t_walk_start;

            clock_t t_power_start = clock();
            for(int k=0; k<K; k++) {
                for(int i=0; i<vert; i++) {
                    former_resultList[i] = resultList[i];
                    resultList[i] = 0;
                }
                for(int i=0; i<vert; i++) {
                    if(former_resultList[i]>0) {
                        int length = g.getOutSize(i);
                        if(length == 0) {
                            resultList[i] += (1-alpha)*former_resultList[i];
                        } else {
                            for(int j=0; j<length; j++) {
                                int tmp_node = g.getOutVert(i,j);
                                resultList[tmp_node] += (1-alpha)*former_resultList[i] / length;
                            }
                        }
                    }
                }
                for(int i=0; i<vert; i++) {
                    resultList[i] += alpha*source_dis[i];
                }
                clock_t t_power_current = clock();
                double current_time = t_power_current - t_power_start;
                if(current_time > random_walk_this_round*rd_ratio) {
                    cout << "K: " << k << " break: power too much" << endl;
                    break;
                }
            }
            clock_t t_power_end = clock();
            cout << "power time: " << (t_power_end - t_power_start) / (double) CLOCKS_PER_SEC << endl;

            batch_size--;
            delete[] res;
        }

        clock_t t1 = clock();
        time = (t1 - t0) / (double) CLOCKS_PER_SEC;
        cout << "MonteCarlo time: " << (t1 - t0) / (double) CLOCKS_PER_SEC << endl;
        return resultList;
    }
};

void ppr_t_PowerMethod(PPR* ppr, vector<int> nodeList, int iterations){
    return ppr->t_PowerMethod(nodeList, iterations);
}

void PPR::PowerMethodMulti(int iterations, int node_count, int num_thread){
    struct timeval t_start,t_end; 
    gettimeofday(&t_start, NULL); 
    long start = ((long)t_start.tv_sec)*1000+(long)t_start.tv_usec/1000; 
    string inputFile = data_path + "/" + target_filename + ".query";
    ifstream node_file(inputFile);
    vector<int> nodes;
    for(int i = 0; i < node_count; i++){
        int temp_node;
        node_file >> temp_node;
        if(g.getOutSize(temp_node) == 0){
            i--;
            cout << "illegal : " << temp_node << endl;
            continue;
        }
        nodes.push_back(temp_node);
    }
    node_file.close();
    if(node_count < num_thread){
        num_thread = node_count;
    }
    vector<thread> threads;
    for(int i = 0; i < num_thread-1; i++){
        vector<int> t_nodes;
        for(int j = 0; j < node_count / num_thread; j++){
            t_nodes.push_back(nodes[i * node_count / num_thread + j]);
        }
        threads.push_back(thread(ppr_t_PowerMethod, this, t_nodes, iterations));
    }
    vector<int> t_nodes;
    for(int j = 0; j < node_count / num_thread; j++){
        t_nodes.push_back(nodes[(num_thread-1) * node_count / num_thread + j]);
    }
    t_PowerMethod(t_nodes, iterations);
    for (int i = 0; i < num_thread - 1; i++){
        threads[i].join();
    }
    gettimeofday(&t_end, NULL); 
    long end = ((long)t_end.tv_sec)*1000+(long)t_end.tv_usec/1000; 
    int cost_time = end - start;

    cout << "cost: " << cost_time / (double) 1000 << endl;
}
#endif
