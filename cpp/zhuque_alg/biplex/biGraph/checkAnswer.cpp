#include "biGraph.hpp"
#include <string>
#include <iostream>
#include <vector>
#include <fstream> 
#include <algorithm>
#include <set>
using namespace std;

struct Node {
    vector<int> R[2];

    Node(int n0, int n1, int * r0, int * r1) {
        R[0].resize(n0);
        for(int i = 0; i < n0; i++) R[0][i] = r0[i];
        R[1].resize(n1);
        for(int i = 0; i < n1; i++) R[1][i] = r1[i];

        sort(R[0].begin(), R[0].end());
        sort(R[1].begin(), R[1].end());
    }

    bool operator == (const Node & t) const {
        if(R[0].size() != t.R[0].size()) return false;
        if(R[1].size() != t.R[1].size()) return false;
        for(int i = 0; i < R[0].size(); i++) {
            if(R[0][i] != t.R[0][i]) return false;
        }
        for(int i = 0; i < R[1].size(); i++) {
            if(R[1][i] != t.R[1][i]) return false;
        }

        return true;
    }

    bool operator < (const Node & t) const {
        if(R[0].size() < t.R[0].size()) {
            return false;
        }
        else if(R[0].size() == t.R[0].size()) {
            if(R[1].size() < t.R[1].size()) {
                return false;
            }
            else if(R[1].size() == t.R[1].size()) {
                for(int i = 0; i < R[0].size(); i++) {
                    if(R[0][i] < t.R[0][i]) return false;
                    else if (R[0][i] > t.R[0][i]) return true;
                }
                for(int i = 0; i < R[1].size(); i++) {
                    if(R[1][i] < t.R[1][i]) return false;
                    else if (R[1][i] > t.R[1][i]) return true;
                }

                return false;
            }

            return true;
        }

        return true;
    }
};

set<Node> st;

int main(int argc, char * argv[])
{
    // "../../../Biclique/data/movielens-10m_ti/movie.txt"

    biGraph * g = new biGraph(argv[2], std::stoi(argv[3]),
             "two");

    int R[2][100000];
    int n[2];

    auto check = [&]()->int {
        for(int i = 0; i < n[0]; i++) {
            int d = 0;
            for(int j = 0; j < n[1]; j++) {
                if(g->connect(R[0][i], R[1][j], 0)) d++;
            }

            if(d + 1 < n[1]) return -1;
        }

        for(int i = 0; i < n[1]; i++) {
            int d = 0;
            for(int j = 0; j < n[0]; j++) {
                if(g->connect(R[1][i], R[0][j], 1)) d++;
            }

            if(d + 1 < n[0]) return -2;
        }

        for(int t = 0, z = 1; t <= 1; t++, z = 0) {
            for(int i = 0; i < n[t]; i++) {
                int u = R[t][i];
                for(int j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                    int v = g->e[t][j];

                    bool inR = false;
                    for(int k = 0; k < n[z]; k++) if(R[z][k] == v) {
                        inR = true;
                        break;
                    }

                    if(inR) continue;

                    int d = 0;
                    for(int k = 0; k < n[t]; k++) {
                        if(g->connect(v, R[t][k], z)) {
                            d++;
                        }
                    }

                    if(d + 1 < n[t]) return -3;
                }
            }
        }

        return true;
    };

    string s;
    n[0] = n[1] = 0;
    int t = 0;

    fstream ccin;
    ccin.open(argv[1], ios::in);

    while(ccin >> s) {
// cout << s << endl;
// cout << s.substr(0, 3) << endl;
        if(s.substr(0, 3) == "R_l") {

            if(n[0] > 0 && check() < 0) {
                printf("error1, %d\n", check());
                return 0;
            }
            Node tmp(n[0], n[1], R[0], R[1]);
            if(st.find(tmp) != st.end()) {
                printf("mulitple\n");
                return 0;
            }
            st.insert(move(tmp));

            n[0] = n[1] = 0;
// cout << s.substr(4) << endl;

            R[0][n[0]++] = stoi(s.substr(4));
            t = 0;
        }
        else if(s.substr(0, 3) == "R_r"){
            R[1][n[1]++] = stoi(s.substr(4));
            t = 1;
        }
        else {

            R[t][n[t]++] = stoi(s);
        }
    }

    if(n[0] > 0 && !check()) {
        printf("error2\n");
        return 0;
    }

    printf("%u\n", st.size());

    delete g;

    return 0;
}

// g++ checkAnswer.cpp -std=c++11 -O3 -o checkAnswer -I ../tools/fastIO.hpp ../tools/listLinearHeap.hpp ../tools/hopstotchHash.hpp -march=native -mavx