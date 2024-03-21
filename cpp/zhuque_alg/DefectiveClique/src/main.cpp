#include <random>
#include "cliqueEnum.h"
#include "polyEnum.h"
#include <pybind11/pybind11.h>
namespace py=pybind11;
void random_graph(int n, double prob, int seed)
{
    assert(n > 0);
    assert(seed >= 0);
    int edges = 0;
    vector<vector<int>> adj;
    adj.resize(n);
    srand(seed);
    for (int i = 0; i < n; ++i) {
        for (int j = i+1;j < n; ++j) {
            if (rand() < RAND_MAX * prob) {
                adj[i].emplace_back(j);
                adj[j].emplace_back(i);
                edges += 2;
            }
        }
    }
    char outfile[1024];
    sprintf(outfile, "randomGraphs/random-%d-%f-%d.txt",n,prob,seed);
    FILE *out = fopen(outfile, "w");
    if (out == NULL) {
        printf("Failed to open %s\n", outfile);
        exit(1);
    }
    fprintf(out, "%d %d\n", n, edges/2);
    for (int i = 0; i < n; ++i) {
        for(auto v: adj[i]) {
            if (i < n) fprintf(out,"%d %d\n", i, v);
        }
    }
    fclose(out);

}

int main(int argc, char *argv[])
{
    Algorithm *mc = NULL;
    if (argc > 2 && strcmp(argv[1],"delay")==0) {
        mc = new PolyEnum();
    }
    else if (argc > 2 && strcmp(argv[1],"branch")==0) {
       mc = new CliqueEnum();
    }
    if (mc == NULL) return 1;
    int alg = 1;
    mc->read_graph(argv[2]);
    mc->setParameters(argc, argv);
    mc->run();
    delete mc;
    return 1;
}

char* concatenate(int num, const char* str) {
    std::string numStr = std::to_string(num);
    std::string result = str+numStr ;
    int length = strlen(result.c_str());
    char* mutableString = new char[length + 1];
    strcpy(mutableString, result.c_str());
    return mutableString;
}
int algorithm(string work,string data,int q,int k,int a,long r)
{
    int c_length = strlen(work.c_str());
    char* work_c = new char[c_length + 1];
    strcpy(work_c, work.c_str());
    Algorithm *mc = NULL;
    if(strstr(work_c,"delay"))
    {
        mc = new PolyEnum();
    }
    else if(strstr(work_c,"branch"))
    {
        mc = new CliqueEnum();
    }
    if (mc == NULL) return 1;
    int alg = 1;
    
    c_length = strlen(data.c_str());
    char* data_c = new char[c_length + 1];
    strcpy(data_c, data.c_str());
    char *argv[]={work_c,data_c,concatenate(q,"-q="),concatenate(k,"-k="),concatenate(a,"-a="),concatenate(r,"-r=")};
    int argc=5;
    mc->read_graph(argv[1]);
    mc->setParameters(argc, argv);
    mc->run();
    delete mc;
    return 1;
}
PYBIND11_MODULE(clique, m)
{
    m.def("algorithm", &algorithm,py::arg("work")="branch",py::arg("data")="datas/ca-GrQc.txt",py::arg("q")=10,py::arg("k")=1,py::arg("a")=6,py::arg("r")=LONG_MAX);
}