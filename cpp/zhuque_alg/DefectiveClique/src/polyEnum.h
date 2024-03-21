#pragma once
#include "algorithms.h"
#include <cstring>

// #define NOTEXCL 0
// #define EXCL -1
// #define DISONE 1
// #define DISTWO 2

#define upair pair<int32, int32>
typedef unsigned long long ulong64;

class PolyEnum : public Algorithm
{
private:
    int32 k          = 2;
    int32 minsize    = 2;
    int32 alg        = 1;

    int32 MCliqueSize  = 0;
    ulong64 cliquenums = 0;
    ulong64 iterations = 0;

    vector<bool> inR;
    vector<CuckooHash>  cuhash;
    vector<vector<int>> fwdadj;
    int32 *recore, *redeg, **readj;
    vector<int> maps;
    vector<int> inRi;

    int32 maxcore = 0;

    clock_t global_time;
    long limited_results = LONG_MAX;
    long cur_out_size = 0;
public:
    PolyEnum(/* args */);
    ~PolyEnum();
    void run();
    void updates(int v, int nonbrs, vector<upair> &C, int32 s, int32 t, vector<upair> &Res);
    void updates1(int v, int nonbrs, vector<upair> &C, int32 s, int32 t, int32 tn, vector<upair> &Res);

    bool pivot(int32 rsize, vector<upair> &C, vector<upair> &X, int &pivotSize);
    inline bool is_nbr(int32 v, int32 u);
    
    void setParameters(int argc, char *argv[]) {
        for (int i = 1; i < argc; ++i) {
            char *p = argv[i];
            if (strstr(p, "-q=")) {
                minsize = atoi(p+3);
            }
            else if (strstr(p, "-k=")) {
                k = atoi(p+3);
            }
            else if (strstr(p, "-a=")) {
                alg = atoi(p+3);
            }
            else if (strstr(p, "-r=")) {
                limited_results = atoi(p+3);
                assert(limited_results>0);
            } 
        }
        printf("minsize=%d, k=%d, alg=%d\n", minsize, k, alg);
    }

    //Polynomial delay algorithm 
    void polyRelabelVertices(int32 *nodeset);
    void polyEnum();
    int  polyGetPI(vector<int> &R, int rsize, int nonbrs);
    bool checkAdd(vector<int> &R, int rsize, int nonbrs, int v, int &_nonbrs);
    void polyRecursion(vector<int> &R, int rsize, int nonbrs, int v);
    
    bool validExtendMax(vector<int> &R, int rsize, int nonbrs, vector<int> &_R, int &_rsize, int &_nonbrs,vector<upair> &CR);
    void recursiveInc(vector<int> &R, int rsize, int nonbrs, vector<upair> &C, vector<upair> &X, vector<upair> &CR);

    upair polyExtendMax(vector<int> &R, int rsize, int nonbrs);
    void polyInitRoot(int v, vector<int32> &R, vector<upair> &P, vector<bool> &visited);
    void polyGenerateAlSat(vector<int> &R, int rsize, int nonbrs, int v);
    void GenerateSubset(vector<int> &P,  int psize, int nnbrP, vector<upair> &S, int deep);

    void recursiveReduc(vector<int> &R, int rsize, int nonbrs, vector<upair> &C, vector<upair> &X, vector<upair> &CR, int deep);
    void reids();

};