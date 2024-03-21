#pragma once
#include "algorithms.h"
#include <cstring>

// #define NOTEXCL 0
// #define EXCL -1
// #define DISONE 1
// #define DISTWO 2

#define upair pair<int32, int32>
typedef unsigned long long ulong64;

class CliqueEnum : public Algorithm
{
private:
        int32 k          = 1;
        int32 minsize    = 2;
        int32 alg        = 5;

        int32 MCliqueSize  = 0;
        ulong64 cliquenums = 0;
        ulong64 iterations = 0;

        vector<bool> inR;
        vector<bool> executed;
        vector<CuckooHash>  cuhash;
        vector<vector<int>> fwdadj;
        int32 *recore, *redeg, **readj;

        vector<int> subdeg;
        int32 maximumGap = 0;
        vector<int> testC;

        clock_t global_time;
        long limited_results = LONG_MAX;
        long cur_out_size;

public:
    CliqueEnum(/* args */);
    ~CliqueEnum();
    void run();
    void recursion(vector<int32> &R, int32 rsize, int32 nonbrs, vector<upair> &C, vector<upair> &X, vector<int32> &P, int32 &psize);
    void updates(int v, int nonbrs, vector<upair> &C, int32 s, int32 t, vector<upair> &Res);

    void initIndexAndFwd(int32 *nodeset);
    void initRoot(int v, vector<int32> &R, vector<upair> &C, vector<upair> &X, vector<bool> &visited);
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
                printf("limited_results=%ld\n", limited_results);
            }
        }
        printf("minsize=%d, k=%d, alg=%d\n", minsize, k, alg);
    }

    //Basic branch algorithm
    void basicEnum();
    void basicEnum2d();
    void basicBranchN(vector<int> &R, int rsize, int nonbrs, vector<int> &C, int csize, vector<int> &X);
    void basicBranchW(vector<int> &R, int rsize, int nonbrs, vector<upair> &C, int csize, vector<upair> &X);
    void basicCORBranch(vector<int> &R, int rsize, int nnbrs, vector<upair> &C, int csize, vector<upair> &X);
    // void updateSet(int v, int nonbrs, vector<upair> &C, int32 s, int32 t, vector<upair> &Res);
    void initHash();

    //Pivot-based branch algorithm
    void pivotEnum();
    void pivotOrderingEnum();
    void combiNonbrs(vector<int> &R, int rsize, int nonbrs, vector<upair> &C, vector<upair> &X, vector<upair> &X1, int xsp);
    void combiNonbrs1(vector<int> &R, int rsize, int nonbrs, vector<upair> &C, vector<upair> &X, vector<upair> &X1, int xsp);
    void enum2d();
    void enum2dOrdering();
    void initRoot1(int v, vector<int32> &R, vector<upair> &C, vector<upair> &X, vector<bool> &visited);
    void TwoHopNbrs(int v, vector<int32> &R, vector<upair> &C, vector<upair> &X, vector<bool> &visited);

    void pivotEnumTwoPhases();
    void pivotBranch(vector<int> &R, int rsize, int nonbrs, vector<upair> &C, int csize, vector<upair> &X);
    void pivotBranch1(vector<int> &R, int rsize, int nonbrs, vector<upair> &C, int csize, vector<upair> &X);
    void pivotBranchUB(vector<int> &R, int rsize, int nonbrs, vector<upair> &C, int csize, vector<upair> &X, const int &mxsize);
    void twoHopNbrs(int v, vector<upair> &C,vector<upair> &X, vector<bool> &visited, vector<bool> &computed);

    // Branch enumeration
    void branchEnum();
    void branchClique(vector<int> &R, int rsize, vector<int> &C, int psize, vector<int> &X);

    void subgraph(vector<upair> &C);

    int testCnt1 = 0;
    int testCnt2 = 0;
    int testCnt3 = 0;
};