#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "ppr.h"
#include <unordered_set>    
#include <cstdlib>
#include <cstring>
#include <pybind11/pybind11.h>

namespace py=pybind11;

int check_inc(int i, int max) {
    if (i == max) {
        exit(1);
    }
    return i + 1;
}

bool maxCmp(const pair<int, double>& a, const pair<int, double>& b){
    return a.second > b.second;
}

vector<int> getRealTopK(int s, int k, string target_filename, int vert){
    stringstream ss;
    ss << "ppr-answer/" << target_filename << "/" << s << ".txt";
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

unordered_map<int, double> getRealTopKMap(int s, int k, string target_filename, int vert){
    unordered_map<int, double> answer_map;
    stringstream ss;
    ss << "ppr-answer/" << target_filename << "/" << s << ".txt";
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

//int main(int argc, char *argv[]){
//    int i = 1;
//    char *endptr;
//    string filename;
//    double alpha = 0.2;
//    int node_count = 20;
//    string algo = "TopPPR";
//    string method = "PW";
//    string graphd = "directed";
//    int power_iterations = 3;
//    double rd_ratio = 1.;
//    double samples = 0.004;
//    double epsilon = 0.5;
//    double l1_error = 1e-5;
//    int num_forests = 6;
//    int batch_size = 1;
//    while (i < argc) {
//        if (!strcmp(argv[i], "-d")) {
//            i = check_inc(i, argc);
//            filename = argv[i];
//        }
//        else if (!strcmp(argv[i], "-algo")) {
//            i = check_inc(i, argc);
//            algo = argv[i];
//        }
//        else if (!strcmp(argv[i], "-method")) {
//            i = check_inc(i, argc);
//            method = argv[i];
//        }
//        else if (!strcmp(argv[i], "-g")) {
//            i = check_inc(i, argc);
//            graphd = argv[i];
//        }
//        else if (!strcmp(argv[i], "-power_iterations")) {
//            i = check_inc(i, argc);
//            power_iterations = strtod(argv[i], &endptr);
//            cout << "power iterations:  " << power_iterations << endl;
//            if ((power_iterations < 0) && endptr) {
//                cerr << "Invalid power_iterations argument" << endl;
//                exit(1);
//            }
//        }
//        else if (!strcmp(argv[i], "-num_forests")) {
//            i = check_inc(i, argc);
//            num_forests = strtod(argv[i], &endptr);
//            cout << "num_forests:  " << num_forests << endl;
//            if ((num_forests < 0) && endptr) {
//                cerr << "Invalid num_forests argument" << endl;
//                exit(1);
//            }
//        }
//        else if (!strcmp(argv[i], "-samples")) {
//            i = check_inc(i, argc);
//            samples = strtod(argv[i], &endptr);
//            if ((samples < 0) && endptr) {
//                cerr << "Invalid samples argument" << endl;
//                exit(1);
//            }
//        }
//        else if (!strcmp(argv[i], "-rd_ratio")) {
//            i = check_inc(i, argc);
//            rd_ratio = strtod(argv[i], &endptr);
//            if ((rd_ratio < 0) && endptr) {
//                cerr << "Invalid rd_ratio argument" << endl;
//                exit(1);
//            }
//        }
//        else if (!strcmp(argv[i], "-epsilon")) {
//            i = check_inc(i, argc);
//            epsilon = strtod(argv[i], &endptr);
//            if ((rd_ratio < 0) && endptr) {
//                cerr << "Invalid epsilon argument" << endl;
//                exit(1);
//            }
//        }
//        else if (!strcmp(argv[i], "-l1_error")) {
//            i = check_inc(i, argc);
//            l1_error = strtod(argv[i], &endptr);
//            if ((l1_error < 0) && endptr) {
//                cerr << "Invalid l1_error argument" << endl;
//                exit(1);
//            }
//        }
//        else if (!strcmp(argv[i], "-batch_size")) {
//            i = check_inc(i, argc);
//            batch_size = strtod(argv[i], &endptr);
//            cout << "batch size: " << batch_size << endl;
//            if ((batch_size < 0) && endptr) {
//                cerr << "Invalid batch_size argument" << endl;
//                exit(1);
//            }
//        }
//        else if (!strcmp(argv[i], "-n")) {
//            i = check_inc(i, argc);
//            node_count = strtod(argv[i], &endptr);
//            if ((node_count < 0) && endptr) {
//                cerr << "Invalid node_count argument" << endl;
//                exit(1);
//            }
//        }
//        else if (!strcmp(argv[i], "-a")) {
//            i = check_inc(i, argc);
//            alpha = strtod(argv[i], &endptr);
//            if (((alpha < 0) || (alpha > 1)) && endptr) {
//                cerr << "Invalid alpha argument" << endl;
//                exit(1);
//            }
//        }
//        else {
//            cout << "invalid argument" << endl;
//            exit(1);
//        }
//        i++;
//    }
//
//    PPR ppr = PPR(filename, alpha);
//    if(algo == "GEN_QUERY"){
//        ofstream outFile("dataset/" + filename + ".query");
//        ppr.generateQueryNode(node_count, outFile);
//        outFile.close();
//    } else if(algo == "GEN_QUERY_DIS"){
//        for(int t=0; t<node_count; t++) {
//            stringstream queryname_ss;
//            queryname_ss << "dataset/" << filename << "/" << t << ".distribution_query";
//            string queryname = queryname_ss.str();
//            ofstream outFile(queryname);
//            ppr.generateSourceDistribution(outFile);
//            outFile.close();
//        }
//    }
//    else if(algo == "GEN_GROUND_TRUTH"){
//        string queryname = "dataset/" + filename + ".query";
//        if(!ppr.is_file_exist(queryname)){
//            cout << "please generate query file first" << endl;
//        }
//        else {
//            if(alpha<0.1) ppr.PowerMethodMulti(2000, node_count, 10);
//            if(alpha>0.1) ppr.PowerMethodMulti(100, node_count, 10);
//        }
//        // ppr.PowerMethodMulti(100, node_count, 10);/*  多线程PowerMethparameter: iteration loops, node size, thread num */
//    }
//    else if(algo == "GEN_GROUND_TRUTH_CENTRALITY"){
//        if(alpha<0.1) ppr.pagerank_centrality_ground_truth(2000);
//        if(alpha>0.1) ppr.pagerank_centrality_ground_truth(100);
//    }
//    else if(algo == "SS") {
//        cout << "method: " << method << endl;
//        string queryname = "dataset/" + filename + ".query";
//        if(!ppr.is_file_exist(queryname)){
//            cout << "please generate query file first" << endl;
//            return 0;
//        }
//        ifstream nodes_file("dataset/" + filename + ".query");
//        vector<int> test_nodes;
//        while(!nodes_file.eof()){
//            int temp_node;
//            nodes_file >> temp_node;
//            test_nodes.push_back(temp_node);
//        }
//        cout << "read done!" << endl;
//        int realCount = 0;
//        int totalQuery = test_nodes.size();
//        for(int t = 0; t < node_count; t++) {
//            int test_node = test_nodes[realCount++];
//            stringstream ss;
//            ss << "ppr-answer/" << filename << "/" << test_node << "-" << alpha << ".txt";
//            string infile = ss.str();
//            if(!ppr.is_file_exist(infile)){
//                cout << "node:" << test_node << " groundtruth file not found, please generate groundtruth first" << endl;
//                return 0;
//            }
//            if(ppr.g.getOutSize(test_node) == 0){
//                t--;
//                node_count--;
//                continue;
//            }
//            ppr.single_source(test_node, method, graphd, power_iterations, batch_size, rd_ratio, epsilon, l1_error);
//        }
//        cout << "avg time: " << ppr.avg_mc_time / (double) node_count << endl;
//        cout << "avg max error: " << ppr.avg_mc_max_err / (double) node_count << endl;
//        cout << "avg l1 error: " << ppr.avg_mc_l1_err / (double) node_count << endl;
//        stringstream ss;
//        ss << filename << "-" << method << "-" << algo << alpha << ".txt";
//        string outputFile = ss.str();
//        ofstream fout(outputFile, ios::app);
//        fout << ppr.avg_mc_time / (double) node_count << " " << ppr.avg_mc_l1_err / (double) node_count << endl;
//        fout.close();
//    } else if(algo == "COMPARE") {
//        ppr.compare_walk_and_forest_time();
//        ppr.compare_walk_and_forest_variance(filename);
//        ppr.compare_l1_error("directed");
//        if(graphd == "undirected") ppr.compare_l1_error("undirected");
//        // ppr.norm_of_pi();
//    } else if(algo == "PC") {
//        cout << "method: " << method << endl;
//        ppr.pagerank_centrality(filename, algo, method, graphd, power_iterations, batch_size, rd_ratio, epsilon);
//    } else if(algo == "DIS") {
//        for(int t = 0; t < node_count; t++) {
//            stringstream queryname_ss;
//            queryname_ss << "dataset/" << filename << "/" << t << ".distributoin_query";
//            string queryname = queryname_ss.str();
//            if(!ppr.is_file_exist(queryname)){
//                cout << "please generate query file first" << endl;
//                return 0;
//            }
//            ifstream nodes_file(queryname);
//            double* source_distribution = new double[ppr.vert];
//            for(int i=0; i<ppr.vert; i++) {
//                nodes_file >> source_distribution[i];
//            }
//            cout << "read done!" << endl;
//            stringstream ss;
//            ss << "ppr-answer/" << filename << "/" << t << "-" << alpha << ".source-txt";
//            string infile = ss.str();
//            if(!ppr.is_file_exist(infile)){
//                cout << "source distribution: " << t << " groundtruth file not found, please generate groundtruth first" << endl;
//                if(alpha<0.1) ppr.any_distribution_ground_truth(source_distribution, t, 2000);
//                if(alpha>0.1) ppr.any_distribution_ground_truth(source_distribution, t, 100);
//                return 0;
//            }
//            ppr.any_distribution(source_distribution, t, method, graphd, power_iterations, batch_size, rd_ratio, epsilon);
//        }
//        cout << "avg time: " << ppr.avg_mc_time / (double) node_count << endl;
//        cout << "avg max error: " << ppr.avg_mc_max_err / (double) node_count << endl;
//        cout << "avg l1 error: " << ppr.avg_mc_l1_err / (double) node_count << endl;
//        stringstream ss;
//        ss << filename << "-" << method << "-" << algo << alpha << ".txt";
//        string outputFile = ss.str();
//        ofstream fout(outputFile, ios::app);
//        fout << ppr.avg_mc_time / (double) node_count << " " << ppr.avg_mc_l1_err / (double) node_count << endl;
//        fout.close();
//    }
//    return 0;
//};

int algorithm(string algo,string filename,string datapath,double alpha,int node_count,string method,string graphd,int power_iterations,double rd_ratio,double samples,double epsilon,double l1_error,int num_forests,int batch_size)
    {
        char *endptr;
        PPR ppr = PPR(datapath, filename, alpha);
    if(algo == "GEN_QUERY"){
        ofstream outFile(datapath + "/" + filename + ".query");
        ppr.generateQueryNode(node_count, outFile);
        outFile.close(); 
    } else if(algo == "GEN_QUERY_DIS"){
        for(int t=0; t<node_count; t++) {
            stringstream queryname_ss;
            queryname_ss << datapath << "/" << filename << "/" << t << ".distribution_query";
            string queryname = queryname_ss.str();
            ofstream outFile(queryname);
            ppr.generateSourceDistribution(outFile);
            outFile.close();
        }
    }
    else if(algo == "GEN_GROUND_TRUTH"){
        string queryname = datapath + "/" + filename + ".query";
        if(!ppr.is_file_exist(queryname)){
            cout << "please generate query file first" << endl;
        }
        else {
            if(alpha<0.1) ppr.PowerMethodMulti(2000, node_count, 10);
            if(alpha>0.1) ppr.PowerMethodMulti(100, node_count, 10);
        }
        // ppr.PowerMethodMulti(100, node_count, 10);/*  多线程PowerMethparameter: iteration loops, node size, thread num */
    }
    else if(algo == "GEN_GROUND_TRUTH_CENTRALITY"){
        if(alpha<0.1) ppr.pagerank_centrality_ground_truth(2000);
        if(alpha>0.1) ppr.pagerank_centrality_ground_truth(100);
    }
    else if(algo == "SS") {
        cout << "method: " << method << endl;
        string queryname = datapath + "/" + filename + ".query";
        if(!ppr.is_file_exist(queryname)){
            cout << "please generate query file first" << endl;
            return 0;
        }
        ifstream nodes_file(datapath + "/" + filename + ".query");
        vector<int> test_nodes;
        while(!nodes_file.eof()){
            int temp_node;
            nodes_file >> temp_node;
            test_nodes.push_back(temp_node);
        }
        cout << "read done!" << endl;
        int realCount = 0;
        int totalQuery = test_nodes.size();
        for(int t = 0; t < node_count; t++) {
            int test_node = test_nodes[realCount++];
            stringstream ss;
            ss << "ppr-answer/" << filename << "/" << test_node << "-" << alpha << ".txt";
            string infile = ss.str();
            if(!ppr.is_file_exist(infile)){
                cout << "node:" << test_node << " groundtruth file not found, please generate groundtruth first" << endl;
                return 0;
            }
            if(ppr.g.getOutSize(test_node) == 0){
                t--;
                node_count--;
                continue;
            }
            ppr.single_source(test_node, method, graphd, power_iterations, batch_size, rd_ratio, epsilon, l1_error);
        }
        cout << "avg time: " << ppr.avg_mc_time / (double) node_count << endl;
        cout << "avg max error: " << ppr.avg_mc_max_err / (double) node_count << endl;
        cout << "avg l1 error: " << ppr.avg_mc_l1_err / (double) node_count << endl;
        stringstream ss;
        ss << filename << "-" << method << "-" << algo << alpha << ".txt";
        string outputFile = ss.str();
        ofstream fout(outputFile, ios::app);
        fout << ppr.avg_mc_time / (double) node_count << " " << ppr.avg_mc_l1_err / (double) node_count << endl;
        fout.close();
    } else if(algo == "COMPARE") {
        ppr.compare_walk_and_forest_time();
        ppr.compare_walk_and_forest_variance(filename);
        ppr.compare_l1_error("directed");
        if(graphd == "undirected") ppr.compare_l1_error("undirected");
        // ppr.norm_of_pi();
    } else if(algo == "PC") {
        cout << "method: " << method << endl;
        ppr.pagerank_centrality(filename, algo, method, graphd, power_iterations, batch_size, rd_ratio, epsilon);
    } else if(algo == "DIS") {
        for(int t = 0; t < node_count; t++) {
            stringstream queryname_ss;
            queryname_ss << datapath << "/" << filename << "/" << t << ".distributoin_query";
            string queryname = queryname_ss.str();
            if(!ppr.is_file_exist(queryname)){
                cout << "please generate query file first" << endl;
                return 0;
            }
            ifstream nodes_file(queryname);
            double* source_distribution = new double[ppr.vert];
            for(int i=0; i<ppr.vert; i++) {
                nodes_file >> source_distribution[i];
            }
            cout << "read done!" << endl;
            stringstream ss;
            ss << "ppr-answer/" << filename << "/" << t << "-" << alpha << ".source-txt";
            string infile = ss.str();
            if(!ppr.is_file_exist(infile)){
                cout << "source distribution: " << t << " groundtruth file not found, please generate groundtruth first" << endl;
                if(alpha<0.1) ppr.any_distribution_ground_truth(source_distribution, t, 2000);
                if(alpha>0.1) ppr.any_distribution_ground_truth(source_distribution, t, 100);
                return 0;
            }
            ppr.any_distribution(source_distribution, t, method, graphd, power_iterations, batch_size, rd_ratio, epsilon);
        }
        cout << "avg time: " << ppr.avg_mc_time / (double) node_count << endl;
        cout << "avg max error: " << ppr.avg_mc_max_err / (double) node_count << endl;
        cout << "avg l1 error: " << ppr.avg_mc_l1_err / (double) node_count << endl;
        stringstream ss;
        ss << filename << "-" << method << "-" << algo << alpha << ".txt";
        string outputFile = ss.str();
        ofstream fout(outputFile, ios::app);
        fout << ppr.avg_mc_time / (double) node_count << " " << ppr.avg_mc_l1_err / (double) node_count << endl;
        fout.close();
    }
    return 1;
    }

PYBIND11_MODULE(TopPPR, m)
{
    m.def("ppr_t_PowerMethod", &ppr_t_PowerMethod);
    py::class_<PPR>(m, "PPR")
        .def(py::init<std::string,std::string,double>())
        .def("is_file_exist", &PPR::is_file_exist)
        .def("getRealTopK", &PPR::getRealTopK)
        .def("getRealTopKMap", &PPR::getRealTopKMap)
        .def("getTopK", &PPR::getTopK)
        .def("naive_monte_carlo", &PPR::naive_monte_carlo)
        .def("PPW_k", &PPR::PPW_k)
        .def("PPF_k", &PPR::PPF_k)
        .def("pc_PPW_k", &PPR::pc_PPW_k)
        .def("pc_PPF_k", &PPR::pc_PPF_k)
        .def("PPW", &PPR::PPW)
        .def("PPF", &PPR::PPF)
        .def("DPI", &PPR::DPI)
        .def("DPI_naive", &PPR::DPI_naive)
        .def("compare_walk_and_forest_variance", &PPR::compare_walk_and_forest_variance)
        .def("compare_walk_and_forest_time", &PPR::compare_walk_and_forest_time)
        .def("norm_of_pi", &PPR::norm_of_pi)
        .def("compare_l1_error", &PPR::compare_l1_error)
        .def("single_source", &PPR::single_source)
        .def("pagerank_centrality", &PPR::pagerank_centrality)
        .def("pagerank_centrality_ground_truth", &PPR::pagerank_centrality_ground_truth)
        .def("any_distribution_ground_truth", &PPR::any_distribution_ground_truth)
        .def("t_PowerMethod", &PPR::t_PowerMethod)
        .def("PowerMethodK", &PPR::PowerMethodK)
        .def("generateQueryNodegenerateQueryNode", &PPR::generateQueryNode)
        .def("generateSourceDistribution", &PPR::generateSourceDistribution)
        .def("MCW", &PPR::MCW)
        .def("MCF", &PPR::MCF)
        .def("pc_PPW", &PPR::pc_PPW)
        .def("pc_PPF", &PPR::pc_PPF)
        .def("pc_naive_monte_carlo", &PPR::pc_naive_monte_carlo)
        .def("pc_naive_forest", &PPR::pc_naive_forest)
        .def("FORA", &PPR::FORA)
        .def("any_distribution", &PPR::any_distribution)
        .def("dis_MCW", &PPR::dis_MCW)
        .def("dis_MCF", &PPR::dis_MCF)
        .def("dis_FORA", &PPR::dis_FORA)
        .def("dis_PPW", &PPR::dis_PPW)
        .def("dis_PPF", &PPR::dis_PPF)
        .def("PowerMethodMulti", &PPR::PowerMethodMulti);
    m.def("getRealTopKMap", &getRealTopKMap);
    m.def("algorithm", &algorithm,py::arg("algo")="PW",py::arg("filename")="youtube_u",py::arg("datapath")="dataset/",py::arg("alpha")=0.2,py::arg("node_count")=20,py::arg("method")="PW",py::arg("graphd")="directed",py::arg("power_iterations")=3,py::arg("rd_ratio")=1,py::arg("samples")=0.004,py::arg("epsilon")=0.5,py::arg("l1_error")=1e-5,py::arg("num_forests")=6,py::arg("batch_size")=1);
}