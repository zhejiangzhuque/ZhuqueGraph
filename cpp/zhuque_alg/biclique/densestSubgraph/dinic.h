#ifndef DINIC
#define DINIC

#include <vector>
#include <queue>
#include <cstring>
#include <cstdio>

#define INF 1e10

struct FlowGraph {
	int n, m;
	std::vector<int> head, vertex, nxt;
	std::vector<double> capcity;
	uint32_t cnt = 0;

	void init(int n, int m) {
		this->n = n;
		this->m = m;
		head.resize(n + 10, -1);
		vertex.resize(m * 2 + 10);
		nxt.resize(m * 2 + 10, -1);
		capcity.resize(m * 2 + 10);
		cnt = 0;
	}

	void addEdge(int u, int v, double w) {
		vertex[cnt] = v;
		capcity[cnt] = w;
		nxt[cnt] = head[u];
		head[u] = cnt++;
	}
};

struct Dinic
{
	std::vector<int> distance;//限制增广路不要重复走点
    std::vector<int> currentArc;
	int s, t;
	FlowGraph g;

	Dinic(int n, int m, int s, int t):s(s), t(t) {
		g.init(n, m);
		distance.resize(n + 1);
        currentArc.resize(n + 1);
	}

	// void print() {
	// 	printf("n %d m %d\n", g.n, g.m);
	// 	printf("s %d t %d\n", s, t);
	// 	for(int u = 1; u <= g.n; u++) {
	// 		for(int e = g.head[u]; ~e; e = g.nxt[e]) {
	// 			printf("%d %d %lld\n", u, g.vertex[e], g.capcity[e]);
	// 		}
	// 	}
	// }

	void addEdge(int u, int v, double w) {
		// if(g.cnt + 2 >= g.m * 2) {
		// 	printf("edges count error\n");fflush(stdout);
		// 	exit(-1);
		// }
		// if(u >= g.n || v >= g.n) {
		// 	printf("nodes count error, %d %d, g.n %d\n", u, v, g.n);fflush(stdout);
		// 	exit(-1);
		// }
		
// printf("addEdge %u %u, cnt %u, m*2 %d\n", u, v, g.cnt, g.m * 2);
// fflush(stdout);
		g.addEdge(u, v, w);
		g.addEdge(v, u, 0);//反向边
	}

	bool bfs() {
		std::fill(distance.begin(), distance.end(), -1);
		
        std::queue<int> que;
        que.push(s);
	    distance[s] = 0;

        while(!que.empty()) {
            int u = que.front(); que.pop();

            for(int e = g.head[u]; ~e; e = g.nxt[e]) {
                int v = g.vertex[e];
                if(std::abs(g.capcity[e]) < 1e-8 || distance[v] != -1) continue;

                distance[v] = distance[u] + 1;
                que.push(v);
            }
        }
        
        return distance[t] != -1;
	}

    double dfs(int u, double inFlow) {
        double outFlow = 0;
		for (int e = currentArc[u]; ~e && outFlow < inFlow; e = g.nxt[e]) {
            currentArc[u] = e;

            int v = g.vertex[e];
			if (g.capcity[e] == 0 || distance[v] != distance[u] + 1)
				continue;
            
            double minF = std::min(g.capcity[e], inFlow - outFlow);
			double pathFlow = v == t ? minF : dfs(v, minF);
            outFlow += pathFlow;
			
            g.capcity[e] -= pathFlow;
            g.capcity[e ^ 1] += pathFlow;
		}

		return outFlow;
	}

	double maxFlow() {
		double totalFlow = 0;

        while(bfs()) {
            memcpy(currentArc.data(), g.head.data(), sizeof(int) * (g.n + 1));
            totalFlow += dfs(s, INF);
        }
		
		return totalFlow;
	}

	uint32_t getEdges() {
		return g.cnt;
	}

	void copyCapacity(std::vector<double> & capcity) {
		// memcpy(capcity.data(), g.capcity.data(), sizeof(double)*g.cnt);
		for(uint32_t i = 0; i < g.cnt; i++) capcity[i] = g.capcity[i];
	}
	void updateEdgesFrom(uint32_t nowEdges, double mid, std::vector<double> & capcity) {
		// memcpy(g.capcity.data(), capcity.data(), sizeof(double)*capcity.size());
		for(uint32_t i = 0; i < capcity.size(); i++) g.capcity[i] = capcity[i];
		for(uint32_t i = nowEdges; i < g.cnt; i += 2) {
			g.capcity[i] = mid;
		}
		for(uint32_t i = nowEdges + 1; i < g.cnt; i += 2) {
			g.capcity[i] = 0.0;
		}
	}

	// void print(uint32_t threshold, std::vector<bool> & vis) {
	// 	printf("the densest sub graph:");
	// 	for(uint32_t u = 0; u < threshold; u++) {
	// 		if(distance[u] != -1) {
	// 			printf("%u ", u);
	// 			vis[u] = true;
	// 		}
	// 	}
	// 	printf("\n");
	// }

	void print(uint32_t threshold) {
		printf("the densest sub graph:");
		uint32_t sum = 0;
		for(uint32_t u = 0; u < threshold; u++) {
			if(distance[u] != -1) {
				printf("%u ", u);
				sum++;
			}
		}
		printf("\nthe densest subgraph size:%u\n", sum);
	}

	uint32_t SSize(uint32_t threshold) {
		uint32_t sum = 0;
		for(uint32_t u = 0; u < threshold; u++) {
			if(distance[u] != -1) sum++;
		}
		return sum;
	}

void getS(uint32_t threshold, std::vector<uint32_t> & nodes) {
	for(uint32_t u = 0; u < threshold; u++) if(distance[u] != -1) nodes.push_back(u);
}
	
void tmp() {
	uint32_t u = 7;
	for (int e = g.head[u]; ~e; e = g.nxt[e]) {
		if(e & 1) continue;
		int v = g.vertex[e];
		printf("%u-%d, %f\n", u, v, g.capcity[e]);
	}
}

double getCut(uint32_t threshold, std::vector<bool> & vis) {
	double sum[3] = {0.0, 0.0, 0.0};
	for(uint32_t u = 0; u < threshold; u++) {
		if(vis[u]) {
			for (int e = g.head[u]; ~e; e = g.nxt[e]) {
				if(e & 1) continue;
				int v = g.vertex[e];
				if(!vis[v]) {
					if(v < 8) sum[0] += g.capcity[e];
					if(8 <= v && v < 30) sum[1] += g.capcity[e];
					if(v == 33 || v == 32) sum[2] += g.capcity[e];
				}
			}
		}
	}

	for(uint32_t u = 0; u < threshold; u++) {
		if(vis[u]) {
			for (int e = g.head[u]; ~e; e = g.nxt[e]) {
				if(e & 1) continue;
				int v = g.vertex[e];
				if(!vis[v]) {
					if(8 <= v && v < 30) {
						printf("%u-%d %f\n", u, v, g.capcity[e]);
					}
				}
			}
		}
	}

	printf("%f %f %f %f\n", sum[0], sum[1], sum[2], sum[0]  + sum[1] + sum[2]);
}

void printEdges() {
	for(uint32_t i = 0; i < g.cnt; i++) {
		printf("i_%u, %u %.2f %u\n", i, g.vertex[i], g.capcity[i], g.nxt[i]);
	}

	
}
};

#endif