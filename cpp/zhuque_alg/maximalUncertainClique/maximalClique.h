#include "algorithms.h"

class maximalClique : public Algorithm
{
private:
	FILE *outs;
	bool output;

	int maximum_clique_size;
	unsigned long long maximal_cliques;
	unsigned long long iterations, comcolors;
	int pruning;

	vector<bool> subset, checknodes;
	vector<int> counts, clset;
	vector<int> LB;

	vector<vector<pairs>> subgraph;
	vector<int> subgraphdeg, subgraphcolor;

public:
	maximalClique();
	~maximalClique();

	void setOutput(const char *str, double eta);
	void initSubgrah() {subgraph.resize(V), subgraphdeg.resize(V);}
	void setSubgraph(node *Can, list &X, vector<bool> &visited);
	void clearSetSubgraph(node *IN, list &Out);
	void printSubgraph();

	//pivot
	void hybrid_pivot_enumerate(int *RP, int *R, double q, int rsize, node *Can, int csize, list &X, int xsize);

	int generateCan(node *Can, int i, int csize, double qn, int u, node *CanN);
	int generateCan(list &Can, double qn, int u, list &CanN);

	bool CheckColors(node *Can, int csize, int curRsize, int curMaxsize);

	void run_pivot_enumerate();
    void run();
};