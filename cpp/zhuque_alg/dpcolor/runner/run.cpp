#include "../tools/getArgs.hpp"

#include "../graph/graph.hpp"

#include "../kclique/ccpath.hpp"
#include "../kclique/ccpathParallel.hpp"
#include "../kclique/pivoterMsk.hpp"
// #include "../tools/multiThreads.hpp"
#include "../kclique/ccParallel.hpp"
#include <pybind11/pybind11.h>
#include <cassert>
#include <string>
#include <iostream>
using std::string;

int run(int deb, string filePath, int N, double alpha, int k, int threads, string algo) {
    if (deb == 0) {
        Graph *g = new Graph();
        v_size n, maxK, tmp;
        double exCnt = 0.0;
        FILE *f = fopen((filePath + "s.txt").c_str(), "r");
        int err = fscanf(f, "%u", &n);
        if (err != 1) {
            printf("s.txt not exist\n");
            return 0;
        }
        if (~fscanf(f, "%u", &maxK)) {
            if (maxK >= k) {
                while (~fscanf(f, "%u-clique: %lf", &tmp, &exCnt)) {
                    if (tmp == k) break;
                }
            }
        }
        g->load(filePath + "edge.bin", filePath + "idx.bin", n);
        samplePlusExact *pt = new samplePlusExact(g);
        if (algo=="cc")
            pt->runCC(k, deb, exCnt, alpha, N);
        else if (algo=="ccpath")
            pt->runCCPath(k, deb, exCnt, alpha, N);
        else
            pt->run(k, deb, exCnt, alpha, N);
        delete g;
        delete pt;
    } else if (deb == 1) {
        int num = 9;
        string name[] = {
                "skitter/", "google/", "berkstan/", "stanford/",
                "amazon0601/", "gowalla/", "comlj/", "orkut/",
                "friender/"
        };
        bool ok[] = {
                false, false, false, false,
                false, false, false, false,
                true
        };
        double ks[500];
        for (int i = 0; i < num; i++) {
            if (!ok[i]) continue;
            Graph *g = new Graph();
            string filePath = "data/";
            filePath += name[i];
            v_size n, maxK, tmp;
            double exCnt = 0.0;
            FILE *f = fopen((filePath + "s.txt").c_str(), "r");
            int err = fscanf(f, "%u", &n);
            if (err != 1) {
                printf("s.txt not exist\n");
                return 0;
            }
            if (~fscanf(f, "%u", &maxK)) {
                while (~fscanf(f, "%u-clique: %lf", &tmp, &exCnt)) {
                    ks[tmp] = exCnt;
                }
            }
            g->load(filePath + "edge.bin", filePath + "idx.bin", n);
            samplePlusExact *pt = new samplePlusExact(g);
            pt->multiRunInit(15, 1.0);
            std::cout << filePath << std::endl;
            for (v_size k = 3; k <= 15; k++) {
                for (v_size tt = 5000000; tt <= 4000000 * 10; tt += 10000000)
                    pt->multiRun(k, 1, ks[k], tt);
            }
            delete g;
            delete pt;
        }
    } else if (deb == 3) {
        Graph *g = new Graph();
        v_size n, maxK, tmp;
        double exCnt = 0.0;
        FILE *f = fopen((filePath + "s.txt").c_str(), "r");
        int err = fscanf(f, "%u", &n);
        if (err != 1) {
            printf("s.txt not exist\n");
            return 0;
        }
        if (~fscanf(f, "%u", &maxK)) {
            if (maxK >= k) {
                while (~fscanf(f, "%u-clique: %lf", &tmp, &exCnt)) {
                    if (tmp == k) break;
                }
            }
        }
        fclose(f);
        g->load(filePath + "edge.bin", filePath + "idx.bin", n);
        if (algo=="cc") {
            ccParallel *pt = new ccParallel(g, threads);
            pt->run(k, deb, exCnt, alpha, N);
            delete pt;
        } else {
            samplePlusExactParallel *pt =
                    new samplePlusExactParallel(g, threads);
            pt->run(k, deb, exCnt, alpha, N);
            delete pt;
        }
        delete g;
    }
    return 0;
}

int makeCSR(string f1, string pEdgePath, string pIdxPath) {
    FILE *fpEdge = fopen(pEdgePath.c_str(), "wb");
    FILE *fpIdx = fopen(pIdxPath.c_str(), "wb");
    FILE *fp = fopen(f1.c_str(), "r");

    ui n, u, v;
    ui m;
    ui *pEdge = nullptr, *pIdx = nullptr, *pDeg = nullptr;
    ui i = 0, preU = 0, preV = 0;

    fscanf(fp, "%u %u", &n, &m);
    m *= 2;

    pDeg = new ui[n]();
    pEdge = new ui[m];
    pIdx = new ui[n + 1]();

    while(fscanf(fp, "%u %u", &u, &v) != EOF) {
        pDeg[u]++;
        pDeg[v]++;
    }

    pIdx[0] = 0;
    for(ui u = 1; u <= n; u++) {
        pIdx[u] = pIdx[u - 1] + pDeg[u - 1];
    }

    fseek(fp, 0, SEEK_SET);
    fscanf(fp, "%u%u", &u, &v);
    printf("%u %u\n", n, m);

    while(fscanf(fp, "%u %u", &u, &v) != EOF) {
        pEdge[pIdx[u]++] = v;
        pEdge[pIdx[v]++] = u;
    }

    pIdx[0] = 0;
    for(ui u = 1; u <= n; u++) {
        pIdx[u] = pIdx[u - 1] + pDeg[u - 1];
    }

    fclose(fp);

    fwrite(pIdx, 4, n + 1, fpIdx);
    fclose(fpIdx);

    fwrite(pEdge, 4, m, fpEdge);
    fclose(fpEdge);

    delete [] pIdx;
    delete [] pEdge;
    delete [] pDeg;
    return 0;
}

int changeToD(string pEdgePath, string pIdxPath, v_size vCnt_){
    Graph * g = new Graph();
    g->load(pEdgePath, pIdxPath, vCnt_);
    printf("load\n");
    g->changeToDegeneracy(pEdgePath, pIdxPath);
    delete g;
    return 0;
}

PYBIND11_MODULE(dpcolor, m){
m.def("run", &run,
pybind11::arg("deb")=0,pybind11::arg("filePath")="data/dblp/",pybind11::arg("N")=5000000,
pybind11::arg("alpha")=1.0,pybind11::arg("k")=10,pybind11::arg("threads")=10,pybind11::arg("algo")="");
m.def("makeCSR", &makeCSR,
pybind11::arg("f1")="./data/dblp/dblp.txt",pybind11::arg("pEdgePath")="./data/dblp/tmpedge.bin",
pybind11::arg("pIdxPath")="./data/dblp/tmpidx.bin");
m.def("changeToD", &changeToD,
pybind11::arg("pEdgePath")="./data/dblp/tmpedge.bin",pybind11::arg("pIdxPath")="./data/dblp/tmpidx.bin",
pybind11::arg("vCnt_")=425957);
}