#ifndef BIGRAPH_HPP
#define BIGRAPH_HPP

#include <vector>
#include <algorithm>
#include <queue>
#include <iostream>

#include "fastIO.hpp"
#include "listLinearHeap.hpp"
#include "hopstotchHash.hpp"

#define unfilled -1
class CuckooHash
{
private:
	/* data */
	uint32_t capacity;
	uint32_t mask;
	uint32_t size;
	const uint32_t buff_size = sizeof(uint32_t);
	uint32_t *hashtable;

	void rehash(uint32_t **_table) {
		uint32_t oldcapacity = capacity;
		mask = mask == 0 ? 1 : ((mask << 1) | 1);
		capacity = (mask + 1) * buff_size;
		uint32_t *newhash = new uint32_t[capacity];
		memset((newhash), unfilled, sizeof(uint32_t) * capacity);
		for (uint32_t i = 0; i < oldcapacity; ++i){
			if ((*_table)[i] != unfilled) insert((*_table)[i], &newhash);
		}
		std::swap((*_table), newhash);
		delete[] newhash;
	}
    
	void insert(const uint32_t &_u, uint32_t **_table) {
		
		uint32_t hs = hash1(_u);
		for (uint32_t i = 0; i < buff_size; ++i) {
			if ((*_table)[hs * buff_size + i] == unfilled){
				(*_table)[hs * buff_size + i] = _u;
				return;
			}
		}
		hs = hash2(_u);
		for (uint32_t i = 0; i < buff_size; ++i) {
			if ((*_table)[hs * buff_size + i] == unfilled){
				(*_table)[hs * buff_size + i] = _u;
				return;
			}
		}

		bool use_hash1 = true;
		uint32_t u = _u;
		for (uint32_t i = 0; i < mask; ++i) {
			uint32_t replaced;
			if (use_hash1) hs = hash1(u);
			else hs = hash2(u);
			uint32_t j = 0;
			for (; j < buff_size; ++j) {
				if ((*_table)[hs * buff_size + j] == unfilled) break;
			}
			if (buff_size == j) {
				replaced = std::move((*_table)[hs * buff_size]);
				j = 1;
				for (; j < buff_size; j++) {
					(*_table)[hs * buff_size + j - 1] =
						std::move((*_table)[hs * buff_size + j]);
				}
				(*_table)[hs * buff_size + j - 1] = u;
			}
			else {
				replaced = std::move((*_table)[hs * buff_size + j]);
				(*_table)[hs * buff_size + j] = u;
			}
			use_hash1 = hs == hash2(replaced);
			u = std::move(replaced);
			if (u == unfilled) return;
		}
		rehash(_table);
		insert(u, _table);
	}

	uint32_t hash1(const uint32_t &x) { return x & mask;}
	uint32_t hash2(const uint32_t &x) { return ~x & mask;}

public:
	CuckooHash(/* args */) {
		capacity = 0;
		hashtable = NULL;
		mask = 0;
		size = 0;
	}
	~CuckooHash() {
		if (hashtable) delete[] hashtable;
	}

	void reserve(uint32_t _size) {
		if (capacity >= _size) return;
		mask = mask == 0 ? 1 : ((mask << 1) | 1);
		while (_size >= mask * buff_size) mask = (mask << 1) | 1;
		capacity = (mask + 1) * buff_size;
		if (hashtable) delete[] hashtable;
		hashtable = new uint32_t[capacity];
		memset(hashtable, unfilled, sizeof(uint32_t) * capacity);
	}

	void insert(const uint32_t &_u) {
		if (find(_u)) return;
		insert(_u, &hashtable);
		size++;
	}

	bool find(const uint32_t &_u) {
		uint32_t hs1 = hash1(_u);
		uint32_t hs2 = hash2(_u);
// if(!(buff_size == 4 && sizeof(uint32_t) == 4))
// printf("%u %u\n", buff_size, sizeof(uint32_t));fflush(stdout);
		// assert(buff_size == 4 && sizeof(uint32_t) == 4);
		__m128i cmp = _mm_set1_epi32(_u);
		__m128i b1 = _mm_load_si128((__m128i*)&hashtable[buff_size * hs1]);
		__m128i b2 = _mm_load_si128((__m128i*)&hashtable[buff_size * hs2]);
        __m128i flag = _mm_or_si128(_mm_cmpeq_epi32(cmp, b1), _mm_cmpeq_epi32(cmp, b2));

		return _mm_movemask_epi8(flag) != 0;
	}
	uint32_t getcapacity() {return capacity;}
	uint32_t getsize() {return size;}
	uint32_t getmask() {return mask;}
	uint32_t *gethashtable() {return hashtable;}
};


struct biGraph {
    uint32_t n1, n2, m, maxDu, maxDv;
    uint32_t core[2];
    uint32_t n[2];

    struct Edge{
        uint32_t u, v;
    };
    std::vector<Edge> edges;
    
    std::vector<uint32_t> pU, e1, pV, e2;
    std::vector<uint32_t> p[2];
    std::vector<uint32_t> e[2];

    std::vector<CuckooHash> cuhash[2]; 

    biGraph() {}

    biGraph(const std::string & filePath, int mode, const std::string & order) {
        if(mode == 1) read(filePath, order);
        else readWithOutUVM(filePath, order);

        // createHashTables();
        cuhash[0].resize(n[0]);
        cuhash[1].resize(n[1]);
        for(uint32_t t = 0; t <= 1; t++) {
            for(uint32_t i = 0; i < n[t]; i++) {
                uint32_t d = p[t][i + 1] - p[t][i];
                cuhash[t][i].reserve(d+1);
                for (uint32_t j = p[t][i]; j < p[t][i + 1];++j)
                    cuhash[t][i].insert(e[t][j]);
            }
        }

        // for(uint32_t i = 0; i < n1; i++) {
        //     uint32_t d = pU[i + 1] - pU[i];
        //     cuhash[0][i].reserve(d+1);
        //     for (uint32_t j = pU[i]; j < pU[i + 1];++j)
        //         cuhash[0][i].insert(e1[j]);
        // }
        // for(uint32_t i = 0; i < n2; i++) {
        //     uint32_t d = pV[i + 1] - pV[i];
        //     cuhash[1][i].reserve(d+1);
        //     for (uint32_t j = pV[i]; j < pV[i + 1];++j)
        //         cuhash[1][i].insert(e2[j]);
        // }
    }

    void read(const std::string & filePath, const std::string & order) {
        fastIO in(filePath, "r");

        n1 = in.getUInt();
        n2 = in.getUInt();
        m = in.getUInt();
   
        edges.resize(m);
        e1.resize(m);
        e2.resize(m);
        pU.resize(n1 + 5);
        pV.resize(n2 + 5);

        // uint32_t minL = n1*10, minR = n2*10;
        for(uint32_t i = 0; i < m; i++) {
            edges[i].u = in.getUInt();
            edges[i].v = in.getUInt();
        }
// printf("there\n");fflush(stdout);
        // changeToDegreeOrder();
// printf("there\n");fflush(stdout);
        // changeToCoreOrder();
        // rawOrder();
        if(order == "core") changeToCoreOrderVersion2();
        else if(order == "two") changeToTwoHopCoreOrder();
        else printf("order parameter error\n");

        n[0] = n1;
        n[1] = n2;
        p[0] = std::move(pU);
        p[1] = std::move(pV);
        e[0] = std::move(e1);
        e[1] = std::move(e2);
    }

    void readWithOutUVM(const std::string & filePath, const std::string & order) {
        fastIO in(filePath, "r");

        // n1 = in.getUInt();
        // n2 = in.getUInt();
        // m = in.getUInt();

        // printf("graph size:n1: %u n2:%u, m %u\n", n1, n2, m);fflush(stdout);
        m = 0;
        uint32_t minL = 1<<30, minR = 1<<30;
        uint32_t maxL = 0, maxR = 0;
        while(!in.empty()) {
            uint32_t u = in.getUInt();
            uint32_t v = in.getUInt();
// printf("%u %u ", u, v);fflush(stdout);
            edges.push_back(Edge{u, v});
            m++;
            minL = std::min(minL, u);
            minR = std::min(minR, v);
            maxL = std::max(maxL, u);
            maxR = std::max(maxR, v);
        }

        // printf("there\n");

        for(uint32_t i = 0; i < m; i++) {
            edges[i].u -= minL;
            edges[i].v -= minR;
        }

        n1 = maxL - minL + 1;
        n2 = maxR - minR + 1;
        // printf("graph size:n1: %u n2:%u, m %u\n", n1, n2, m);fflush(stdout);
// for(uint32_t i = 0; i < m; i++) {
//     printf("%u %u\n", edges[i].u, edges[i].v);
// }

        // edges.resize(m);
        e1.resize(m);
        e2.resize(m);
        pU.resize(n1 + 5);
        pV.resize(n2 + 5);

        // uint32_t minL = n1*10, minR = n2*10;
        // for(uint32_t i = 0; i < m; i++) {
        //     edges[i].u = in.getUInt();
        //     edges[i].v = in.getUInt();
        // }
// printf("there\n");fflush(stdout);
        // changeToDegreeOrder();
// printf("there\n");fflush(stdout);
        // changeToCoreOrder();
        // rawOrder();

        if(order == "core") changeToCoreOrderVersion2();
        else if(order == "two") changeToTwoHopCoreOrder();
        else printf("order parameter error\n");

        n[0] = n1;
        n[1] = n2;
        p[0] = std::move(pU);
        p[1] = std::move(pV);
        e[0] = std::move(e1);
        e[1] = std::move(e2);
    }

    void changeToTwoHopCoreOrder() {
        // printf("here two\n");fflush(stdout);
        std::vector<uint32_t> d1, d2;
        
        d1.resize(n1);
        d2.resize(n2);

        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        maxDu = 0;
        for(uint32_t i = 0; i < n1; i++) {
            maxDu = std::max(maxDu, d1[i]);
        }
        maxDv = 0;
        for(uint32_t i = 0; i < n2; i++) {
            maxDv = std::max(maxDv, d2[i]);
        }

        pU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        }
        for(uint32_t i = 0; i < m; i++) {
            e1[pU[edges[i].u]++] = edges[i].v; 
        }
        pU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        } 
        
        pV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }
        for(uint32_t i = 0; i < m; i++) {
            e2[pV[edges[i].v]++] = edges[i].u; 
        }
        pV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }

        uint32_t maxTwoHopDegreeU = 0, maxTwoHopDegreeV = 0;
        d1.clear();
        std::vector<uint32_t> stk;
        uint32_t n = std::max(n1, n2);
        std::vector<uint32_t> ids(n);
        std::vector<uint32_t> keys(n);
        std::vector<uint32_t> labelsL(n1);
        std::vector<uint32_t> labelsR(n2);
        uint32_t l = 0;
        
        for(uint32_t u = 0; u < n1; u++) ids[u] = u;
        for(uint32_t u = 0; u < n1; u++) {
            stk.clear();
            uint32_t tmp = pU[u + 1] - pU[u];

            for(uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                uint32_t v = e1[i];
                for(uint32_t j = pV[v]; j < pV[v + 1]; j++) {
                    uint32_t w = e2[j];
                    if(d1[w] == 0) {
                        d1[w] = 1;
                        stk.push_back(w);
                    }
                }
            }
            tmp += stk.size();

            maxTwoHopDegreeU = std::max(maxTwoHopDegreeU, tmp);
            keys[u] = tmp;

            for(auto w : stk) d1[w] = 0;
        }

        ListLinearHeap heap(n1, maxTwoHopDegreeU + 1);
        heap.init(n1, maxTwoHopDegreeU + 1, ids.data(), keys.data());
        core[0] = 0;
        core[1] = 0;

        for(uint32_t i = 0; i < n1; i++) {
            uint32_t u, degU;

            if(!heap.pop_min(u, degU)) printf("errorLheap\n");
            labelsL[u] = l++;
            core[0] = std::max(core[0], degU);

            stk.clear();
            for(uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                uint32_t v = e1[i];
                for(uint32_t j = pV[v]; j < pV[v + 1]; j++) {
                    uint32_t w = e2[j];
                    if(d1[w] == 0) {
                        d1[w] = 1;
                        stk.push_back(w);
                        heap.decrement(w, 1);
                    }
                }
            }

            for(auto w : stk) d1[w] = 0;
        }

        d2.clear();
        l = 0;
        
        for(uint32_t u = 0; u < n2; u++) ids[u] = u;
        for(uint32_t u = 0; u < n2; u++) {
            stk.clear();
            uint32_t tmp = pV[u + 1] - pV[u];

            for(uint32_t i = pV[u]; i < pV[u + 1]; i++) {
                uint32_t v = e2[i];
                for(uint32_t j = pU[v]; j < pU[v + 1]; j++) {
                    uint32_t w = e1[j];
                    if(d2[w] == 0) {
                        d2[w] = 1;
                        stk.push_back(w);
                    }
                }
            }
            tmp += stk.size();

            maxTwoHopDegreeV = std::max(maxTwoHopDegreeV, tmp);
            keys[u] = tmp;

            for(auto w : stk) d2[w] = 0;
        }

        // printf("core %u %u\n", maxTwoHopDegreeU, maxTwoHopDegreeV);

        ListLinearHeap rheap(n2, maxTwoHopDegreeV + 1);
        rheap.init(n2, maxTwoHopDegreeV + 1, ids.data(), keys.data());

        for(uint32_t i = 0; i < n2; i++) {
            uint32_t u, degU;

            if(!rheap.pop_min(u, degU)) printf("errorLheap\n");
            labelsR[u] = l++;
            core[1] = std::max(core[1], degU);

            stk.clear();
            for(uint32_t i = pV[u]; i < pV[u + 1]; i++) {
                uint32_t v = e2[i];
                for(uint32_t j = pU[v]; j < pU[v + 1]; j++) {
                    uint32_t w = e1[j];
                    if(d2[w] == 0) {
                        d2[w] = 1;
                        stk.push_back(w);
                        heap.decrement(w, 1);
                    }
                }
            }

            for(auto w : stk) d2[w] = 0;
        }

        for(uint32_t i = 0; i < m; i++) {
            edges[i].u = labelsL[edges[i].u];
            edges[i].v = labelsR[edges[i].v];
        }

        std::fill(d1.begin(), d1.begin() + n1, 0);
        std::fill(d2.begin(), d2.begin() + n2, 0);
        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);

        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        for(uint32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for(uint32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for(uint32_t i = 0; i < m; i++) {
            e1[ pU[edges[i].u]++ ] = edges[i].v;
        }
        for(uint32_t i = 0; i < m; i++) {
            e2[ pV[edges[i].v]++ ] = edges[i].u;
        }

        pU[0] = pV[0] = 0;
        for(uint32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for(uint32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for(uint32_t i = 0; i < n1; i++) {
            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
        }
        for(uint32_t i = 0; i < n2; i++) {
            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
        }
    }

    void coreReduction(int p, int q) {
        std::queue<uint32_t> qL, qR;
        std::vector<uint32_t> d1(n1 + 1), d2(n2 + 1);
        std::vector<uint32_t> labelsL(n1 + 1), labelsR(n2 + 1);
        std::vector<bool> visL(n1 + 1), visR(n2 + 1);
                
        for(uint32_t i = 0; i < n1; i++) {
            d1[i] = deg1(i);
            if(deg1(i) < q) {
                qL.push(i);
                visL[i] = true;
            }
        }
        for(uint32_t i = 0; i < n2; i++) {
            d2[i] = deg2(i);
            if(deg2(i) < p) {
                qR.push(i);
                visR[i] = true;
            }
        }

        while(!qL.empty() || !qR.empty()) {
            while(!qL.empty()) {
                uint32_t u = qL.front(); qL.pop();

                for(uint32_t i = 0; i < d1[u]; i++) {
                    uint32_t v = e1[pU[u] + i];
                    // if(d2[v] < q) continue;

                    for(uint32_t j = pV[v]; j < pV[v] + d2[v]; j++) {
                        if(e2[j] == u) {
                            --d2[v];
                            std::swap(e2[j], e2[pV[v] + d2[v]]);

                            if(d2[v] == p - 1 && !visR[v]) {
                                qR.push(v);
                                visR[v] = true;
                            }
                            break;
                        }
                    }
                }
            }

            while(!qR.empty()) {
                uint32_t v = qR.front(); qR.pop();

                for(uint32_t i = 0; i < d2[v]; i++) {
                    uint32_t u = e2[pV[v] + i];
                    // if(d1[u] < p) continue;

                    for(uint32_t j = pU[u]; j < pU[u] + d1[u]; j++) {
                        if(e1[j] == v) {
                            --d1[u];
                            std::swap(e1[j], e1[pU[u] + d1[u]]);

                            if(d1[u] == q - 1 && !visL[u]) {
                                qL.push(u);
                                visL[u] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }

        uint32_t pL = 1, pR = 1;
        for(uint32_t u = 0; u < n1; u++) {
            if(!visL[u]) labelsL[u] = pL++;
        }
        for(uint32_t v = 0; v < n2; v++) {
            if(!visR[v]) labelsR[v] = pR++;
        }
        
        uint32_t pm = 0;
        for(uint32_t u = 0; u < n1; u++) {
            if(visL[u]) continue;
            for(uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                uint32_t v = e1[i];
                if(!visR[v]) {
                    edges[pm].u = labelsL[u] - 1;
                    edges[pm].v = labelsR[v] - 1;
                    ++pm;
                } 
            }
        }
        m = pm;

        n1 = pL - 1;
        n2 = pR - 1;
        // printf("n1 %u, n2 %u, m %u\n", n1, n2, m);
        
        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);

        std::fill(d1.begin(), d1.begin() + n1 + 1, 0);
        std::fill(d2.begin(), d2.begin() + n2 + 1, 0);
        changeToDegreeOrder();

        // std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
        // std::fill(pV.begin(), pV.begin() + n2 + 1, 0);

        // std::fill(d1.begin(), d1.begin() + n1 + 1, 0);
        // std::fill(d2.begin(), d2.begin() + n2 + 1, 0);

        // for(uint32_t i = 0; i < m; i++) {
        //     ++d1[edges[i].u];
        //     ++d2[edges[i].v];
        // }

        // maxDu = 0;
        // for(uint32_t i = 0; i < n1; i++) {
        //     maxDu = std::max(maxDu, d1[i]);
        // }

        // maxDv = 0;
        // for(uint32_t i = 0; i < n2; i++) {
        //     maxDv = std::max(maxDv, d2[i]);
        // }

        // for(uint32_t i = 0; i < n1; i++) {
        //     pU[i + 1] = pU[i] + d1[i];
        // }
        // for(uint32_t i = 0; i < n2; i++) {
        //     pV[i + 1] = pV[i] + d2[i];
        // }

        // for(uint32_t i = 0; i < m; i++) {
        //     e1[ pU[edges[i].u]++ ] = edges[i].v;
        // }
        // for(uint32_t i = 0; i < m; i++) {
        //     e2[ pV[edges[i].v]++ ] = edges[i].u;
        // }

        // pU[0] = pV[0] = 0;
        // for(uint32_t i = 0; i < n1; i++) {
        //     pU[i + 1] = pU[i] + d1[i];
        // }
        // for(uint32_t i = 0; i < n2; i++) {
        //     pV[i + 1] = pV[i] + d2[i];
        // }

        // for(uint32_t i = 0; i < n1; i++) {
        //     std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
        // }
        // for(uint32_t i = 0; i < n2; i++) {
        //     std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
        // }

        // printf("core reduction\n");
    }

    void changeToCoreOrderVersion2() {
// printf("here core\n");fflush(stdout);
        std::vector<uint32_t> d1, d2;
        
        d1.resize(n1);
        d2.resize(n2);

        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        maxDu = 0;
        for(uint32_t i = 0; i < n1; i++) {
            maxDu = std::max(maxDu, d1[i]);
        }
        maxDv = 0;
        for(uint32_t i = 0; i < n2; i++) {
            maxDv = std::max(maxDv, d2[i]);
        }

        pU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        }
        for(uint32_t i = 0; i < m; i++) {
            e1[pU[edges[i].u]++] = edges[i].v; 
        }
        pU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        } 
        
        pV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }
        for(uint32_t i = 0; i < m; i++) {
            e2[pV[edges[i].v]++] = edges[i].u; 
        }
        pV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }
// printf("there\n");fflush(stdout);
        ListLinearHeap lheap(n1, maxDu + 1), rheap(n2, maxDv + 1);
        uint32_t n = std::max(n1, n2);
        std::vector<uint32_t> ids(n);
        std::vector<uint32_t> keys(n);
        std::vector<uint32_t> labelsL(n1);
        std::vector<uint32_t> labelsR(n2);
        uint32_t l1 = 0, l2 = 0;

        for(uint32_t i = 0; i < n1; i++) {
            ids[i] = i;
            keys[i] = d1[i] + 1;
        }
        lheap.init(n1, maxDu + 1, ids.data(), keys.data());
        for(uint32_t i = 0; i < n2; i++) {
            ids[i] = i;
            keys[i] = d2[i] + 1;
// printf("%u %u\n", ids[i], keys[i]);
        }
        rheap.init(n2, maxDv + 1, ids.data(), keys.data());
        
// printf("there\n");fflush(stdout);
        // uint32_t minN = std::min(n1, n2);
        core[0] = 0;
        core[1] = 0;
        for(uint32_t i = 0; i < n1 + n2; i++) {
            uint32_t u, degU = n2 + 11;
            uint32_t v, degV = n1 + 11;
// printf("    %u\n", i);fflush(stdout);
            if(!lheap.empty() && !lheap.pop_min(u, degU)) printf("errorLheap\n");
// printf("%u %u\n", u, degU);fflush(stdout);
            if(!rheap.empty() && !rheap.pop_min(v, degV)) printf("errorRheap\n");
// printf("%u %u\n", v, degV);fflush(stdout);
            if(degU <= degV) {
                if(degV != n1 + 11)
                    rheap.insert(v, degV);
// printf("rm L %u\n", u);      
                for(uint32_t j = pU[u]; j < pU[u + 1]; j++) {
                    ui d= rheap.decrement(e1[j]);
// printf("afer %u %u\n", e1[j], d);
                }
                labelsL[u] = l1++;
                core[0] = std::max(core[0], degU);
            }
            else {
                if(degU != n2 + 11)
                    lheap.insert(u, degU);
// printf("rm R %u\n", v);   
                for(uint32_t j = pV[v]; j < pV[v + 1]; j++) {
                    lheap.decrement(e2[j]);
                }
                labelsR[v] = l2++;
                core[1] = std::max(core[1], degV);
            }
        }
// printf("there\n");fflush(stdout);
        printf("core %u %u\n", core[0], core[1]);
        

        // if(n1 < n2) {
        //     for(uint32_t j = n1; j < n2; j++) {
        //         uint32_t v, degV;
        //         if(!rheap.pop_min(v, degV)) printf("errorRheap\n");
        //         labelsR[v] = j;
        //     }
        // }
        // else if(n1 > n2) {
        //     for(uint32_t j = n2; j < n1; j++) {
        //         uint32_t u, degU;
        //         if(!lheap.pop_min(u, degU)) printf("errorRheap\n");
        //         labelsL[u] = j;
        //     }
        // }

    

        for(uint32_t i = 0; i < m; i++) {
            edges[i].u = labelsL[edges[i].u];
            edges[i].v = labelsR[edges[i].v];
        }

        std::fill(d1.begin(), d1.begin() + n1, 0);
        std::fill(d2.begin(), d2.begin() + n2, 0);
        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);

        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        for(uint32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for(uint32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for(uint32_t i = 0; i < m; i++) {
            e1[ pU[edges[i].u]++ ] = edges[i].v;
        }
        for(uint32_t i = 0; i < m; i++) {
            e2[ pV[edges[i].v]++ ] = edges[i].u;
        }

        pU[0] = pV[0] = 0;
        for(uint32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for(uint32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for(uint32_t i = 0; i < n1; i++) {
            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
        }
        for(uint32_t i = 0; i < n2; i++) {
            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
        }

    }

    void changeToCoreOrder() {
        std::vector<uint32_t> d1, d2;
        
        d1.resize(n1);
        d2.resize(n2);

        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        maxDu = 0;
        for(uint32_t i = 0; i < n1; i++) {
            maxDu = std::max(maxDu, d1[i]);
        }
        maxDv = 0;
        for(uint32_t i = 0; i < n2; i++) {
            maxDv = std::max(maxDv, d2[i]);
        }

        pU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        }
        for(uint32_t i = 0; i < m; i++) {
            e1[pU[edges[i].u]++] = edges[i].v; 
        }
        pU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        } 
        
        pV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }
        for(uint32_t i = 0; i < m; i++) {
            e2[pV[edges[i].v]++] = edges[i].u; 
        }
        pV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }

        ListLinearHeap lheap(n1, maxDu + 1), rheap(n2, maxDv + 1);
        uint32_t n = std::max(n1, n2);
        std::vector<uint32_t> ids(n);
        std::vector<uint32_t> keys(n);
        std::vector<uint32_t> labelsL(n1);
        std::vector<uint32_t> labelsR(n2);

        for(uint32_t i = 0; i < n1; i++) {
            ids[i] = i;
            keys[i] = d1[i] + 1;
        }
        lheap.init(n1, maxDu + 1, ids.data(), keys.data());
        for(uint32_t i = 0; i < n2; i++) {
            ids[i] = i;
            keys[i] = d2[i] + 1;
        }
        rheap.init(n2, maxDv + 1, ids.data(), keys.data());
        
        uint32_t minN = std::min(n1, n2);
        for(uint32_t i = 0; i < minN; i++) {
            uint32_t u, degU;
            uint32_t v, degV;

            if(!lheap.pop_min(u, degU)) printf("errorLheap\n");
            if(!rheap.pop_min(v, degV)) printf("errorRheap\n");
    // printf("i %u, %u %u, %u %u\n", i, u, degU, v, degV);
    // fflush(stdout);
            labelsL[u] = i;
            labelsR[v] = i;
            for(uint32_t j = pU[u]; j < pU[u + 1]; j++) {
                rheap.decrement(e1[j]);
            }
            for(uint32_t j = pV[v]; j < pV[v + 1]; j++) {
                lheap.decrement(e2[j]);
            }
        }
        if(n1 < n2) {
            for(uint32_t j = n1; j < n2; j++) {
                uint32_t v, degV;
                if(!rheap.pop_min(v, degV)) printf("errorRheap\n");
                labelsR[v] = j;
            }
        }
        else if(n1 > n2) {
            for(uint32_t j = n2; j < n1; j++) {
                uint32_t u, degU;
                if(!lheap.pop_min(u, degU)) printf("errorRheap\n");
                labelsL[u] = j;
            }
        }

        for(uint32_t i = 0; i < m; i++) {
            edges[i].u = labelsL[edges[i].u];
            edges[i].v = labelsR[edges[i].v];
        }

        std::fill(d1.begin(), d1.begin() + n1, 0);
        std::fill(d2.begin(), d2.begin() + n2, 0);
        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);

        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        for(uint32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for(uint32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for(uint32_t i = 0; i < m; i++) {
            e1[ pU[edges[i].u]++ ] = edges[i].v;
        }
        for(uint32_t i = 0; i < m; i++) {
            e2[ pV[edges[i].v]++ ] = edges[i].u;
        }

        pU[0] = pV[0] = 0;
        for(uint32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for(uint32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for(uint32_t i = 0; i < n1; i++) {
            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
        }
        for(uint32_t i = 0; i < n2; i++) {
            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
        }

        // print();
        // fflush(stdout);
    }

    void rawOrder() {
        std::vector<uint32_t> d1, d2;

        d1.resize(n1);
        d2.resize(n2);

        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        maxDu = 0;
        for(uint32_t i = 0; i < n1; i++) {
            maxDu = std::max(maxDu, d1[i]);
        }
        maxDv = 0;
        for(uint32_t i = 0; i < n2; i++) {
            maxDv = std::max(maxDv, d2[i]);
        }

        pU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        }
        for(uint32_t i = 0; i < m; i++) {
            e1[pU[edges[i].u]++] = edges[i].v; 
        }
        pU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        } 
        
        pV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }   
        for(uint32_t i = 0; i < m; i++) {
            e2[pV[edges[i].v]++] = edges[i].u; 
        }
        pV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }

        for(uint32_t i = 0; i < n1; i++) {
            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
        }
        for(uint32_t i = 0; i < n2; i++) {
            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
        }
    }

    void changeToDegreeOrder() {
        std::vector<uint32_t> d1, d2;
        std::vector<uint32_t> label1, label2;

        d1.resize(std::max(n1, n2) + 1);
        d2.resize(std::max(n1, n2) + 1);
        label1.resize(n1);
        label2.resize(n2);
        
        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }
        
        maxDu = 0;
        for(uint32_t i = 0; i < n1; i++) {
            maxDu = std::max(maxDu, d1[i]);
        }

        maxDv = 0;
        for(uint32_t i = 0; i < n2; i++) {
            maxDv = std::max(maxDv, d2[i]);
        }
        
        for(uint32_t i = 0; i < n1; i++) {
// printf("%u-%u-%u\n", i, d1[i] + 1, pU[d1[i] + 1]);
            pV[d1[i] + 1]++;
// printf("%u-%u-%u\n", i, d1[i] + 1, pU[d1[i] + 1]);
        }
// for(uint32_t i = 0; i <= maxDu; i++) {
//     printf("%d-%u\n", i, pU[i]);
// }
// -f data\exam2.txt -p 4 -q 4 -pm
        for(uint32_t i = 0; i < maxDu; i++) {
            pV[i + 1] += pV[i];
        }
        
// printf("labels:");
        for(uint32_t i = 0; i < n1; i++) {
            label1[i] = pV[d1[i]]++;
// printf("%u-%u ", i,label1[i]);
        }
// printf("\n");
// printf("there\n");fflush(stdout);
        for(uint32_t i = 0; i < n2; i++) {
            pU[d2[i] + 1]++;
        }
        for(uint32_t i = 0; i < maxDv; i++) {
            pU[i + 1] += pU[i];
        }
        
        for(uint32_t i = 0; i < n2; i++) {
            label2[i] = pU[d2[i]]++;
        }

        for(uint32_t i = 0; i < m; i++) {
            edges[i].u = label1[edges[i].u];
            edges[i].v = label2[edges[i].v];
        }

        std::fill(d1.begin(), d1.begin() + std::max(n1, n2) + 1, 0);
        std::fill(d2.begin(), d2.begin() + std::max(n1, n2) + 1, 0);
        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);
        // memset(buffer, 0, sizeof(uint32_t) * bufferSize);
        for(uint32_t i = 0; i < m; i++) {
// printf("i%u:%u-%u ", i, edges[i].u, edges[i].v);fflush(stdout);
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }
// printf("\n");
// printf("there2\n");fflush(stdout);
        for(uint32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for(uint32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for(uint32_t i = 0; i < m; i++) {
            e1[ pU[edges[i].u]++ ] = edges[i].v;
        }
        for(uint32_t i = 0; i < m; i++) {
            e2[ pV[edges[i].v]++ ] = edges[i].u;
        }
// printf("there3\n");fflush(stdout);
        pU[0] = pV[0] = 0;
        for(uint32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for(uint32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for(uint32_t i = 0; i < n1; i++) {
            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
        }
        for(uint32_t i = 0; i < n2; i++) {
            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
        }
// printf("thereed\n");
// print();
        // delete [] buffer;
    }

    void swapUV() {
        std::swap(n[0], n[1]);
        std::swap(maxDu, maxDv);
        
        p[0].swap(p[1]);
        e[0].swap(e[1]);
    }

    uint32_t deg1(uint32_t u) {
        return p[0][u + 1] - p[0][u];
    }
    uint32_t deg2(uint32_t v) {
        return p[1][v + 1] - p[1][v];
    }

    // bool connectUV(uint32_t u, uint32_t v) {
    //     return std::binary_search(
    //         e[0].begin() + p[0][u], e[0].begin() + p[0][u + 1],
    //         v
    //     );
    // }
    bool connectUV(uint32_t u, uint32_t v) {
        return cuhash[0][u].find(v);
    }
    // bool connectVU(uint32_t v, uint32_t u) {
    //     return std::binary_search(
    //         e[1].begin() + p[1][v], e[1].begin() + p[1][v + 1],
    //         u
    //     );
    // }
    bool connectVU(uint32_t v, uint32_t u) {
        return cuhash[1][v].find(u);
    }
    // bool connect(uint32_t u, uint32_t v, uint32_t t) {
    //     return std::binary_search(
    //         e[t].begin() + p[t][u], e[t].begin() + p[t][u + 1], v
    //     );
    // }

    bool connect(uint32_t u, uint32_t v, uint32_t t) {
        return cuhash[t][u].find(v);
    }
    

    void print() {
        printf("U:\n");
        for(uint32_t u = 0; u < n[0]; u++) {
            printf("%u:", u);
            for(uint32_t i = p[0][u]; i < p[0][u + 1]; i++) {
                printf("%u ", e[0][i]);
            }
            printf("\n");
        }
        printf("V:\n");
        for(uint32_t v = 0; v < n[1]; v++) {
            printf("%u:", v);
            for(uint32_t i = p[1][v]; i < p[1][v + 1]; i++) {
                printf("%u ", e[1][i]);
            }
            printf("\n");
        }
    }

    uint32_t deg(uint32_t u, uint32_t t) {
        return p[t][u + 1] - p[t][u];
    }


//SIMD
    std::vector<hopstotchHash> hashTables[2];

    void createHashTables() {
        for(uint32_t t = 0; t <= 1; t++) {
            hashTables[t].resize(n[t]);
        
            for(uint32_t u = 0; u < n[t]; u++) {
                // if(deg(u, t) > H) {//H is the bucket size of hopstotch hash
                //     hashTables[t][u].build(e[t].data() + p[t][u], deg(u, t));
                // }
                hashTables[t][u].build(e[t].data() + p[t][u], deg(u, t));
            }
        }
    }

    // void createHashTableR() {
    //     hashTablesR.resize(n2);
        
    //     for(uint32_t v = 0; v < n2; v++) {
    //         if(deg2(v) > H) {
    //             hashTablesR[v].build(e2.data() + pV[v], deg2(v));
    //         }
    //     }
    // }

    // bool connectUVFast(uint32_t u, uint32_t v) {
    //     // return connectUV(u, v);
    //     if(hashTablesL[u].n > 0) return hashTablesL[u].contain(v);
    //     else {
    //         return connectUV(u, v);

    //         // if(deg1(u) < 8) return connectUV(u, v);
    //         // return ccSIMD(e1.data(), p[0][u], v) || std::binary_search(
    //         //         e1.begin() + p[0][u] + 8, e1.begin() + p[0][v + 1], v);
    //     }
    // }

    // bool connectVUFast(uint32_t v, uint32_t u) {
    //     // return connectVU(v, u);
    //     if(hashTablesR[v].n > 0) return hashTablesR[v].contain(u);
    //     else {
    //         return connectVU(v, u);

    //         // if(deg2(v) < 8) return connectVU(v, u);
    //         // return ccSIMD(e2.data(), pV[v], u) || std::binary_search(
    //         //         e2.begin() + pV[v] + 8, e2.begin() + pV[v + 1], u);
    //     }
    // }
    // bool connectRaw(uint32_t u, uint32_t v, uint32_t t) {
    //     return std::binary_search(
    //         e[t].begin() + p[t][u], e[t].begin() + p[t][u + 1], v
    //     );
    // }
    // bool connect(uint32_t u, uint32_t v, uint32_t t) {
    //     if(hashTables[t][u].n > 0) return hashTables[t][u].contain(v);
    //     else return connectRaw(u, v, t);
    // }

//     bool ccSIMD(uint32_t * startAddress, uint32_t seekSize, uint32_t v) {
//         __m256i eightV = _mm256_set1_epi32(v);
        
//         auto address = reinterpret_cast<const __m256i *>(startAddress + seekSize);
//         __m256i eightNextV = _mm256_loadu_si256(address);

//         __m256i cmpRes = _mm256_cmpeq_epi32(eightV, eightNextV);
//         auto msk = _mm256_movemask_epi8(cmpRes);
// #ifdef _WIN32
//         if(__popcnt(msk) > 0) return true;
// #else
//         if(_popcnt32(msk) > 0) return true;
// #endif
//         return false;
//     }
};

// print();
//         ListLinearHeap lheap(n1, maxDu + 1), rheap(n2, maxDv + 1);
//         uint32_t n = std::max(n1, n2);
//         std::vector<uint32_t> ids(n);
//         std::vector<uint32_t> keys(n);
//         std::vector<uint32_t> labelsL(n1);
//         std::vector<uint32_t> labelsR(n2);
//         std::vector<bool> removedL(n1), removedR(n2);
//         std::vector<uint32_t> d1(n1), d2(n2);
//         uint32_t pL = 1, pR = 1;

//         uint32_t id1 = 0, id2 = 0;
        
//         bool ok = false;
// printf("%u %u\n", n1, n2);
//         for(uint32_t i = 0; i < n1; i++) {
//             ids[i] = i;
//             d1[i] = keys[i] = p[0][i + 1] - p[0][i];
//         }
// // for(int i = n1 - 10; i < n1; i++) {
// //     printf("%u\n", keys[i]);
// // }
// // for(int i = n1/2 - 10; i < n1/2; i++) {
// //     printf("%u\n", keys[i]);
// // }
//         lheap.init(n1, maxDu + 1, ids.data(), keys.data());
//         for(uint32_t i = 0; i < n2; i++) {
//             ids[i] = i;
//             d2[i] = keys[i] = pV[i + 1] - pV[i];
//         }
//         rheap.init(n2, maxDv + 1, ids.data(), keys.data());

//         uint32_t minN = std::min(n1, n2);
// uint32_t debug = 1410904;
//         uint32_t i = 0;
//         for(; i < minN; i++) {
//             uint32_t u, degU;
//             uint32_t v, degV;

//             if(!lheap.pop_min(u, degU)) printf("errorLheap\n");
//             if(!rheap.pop_min(v, degV)) printf("errorRheap\n");

//             removedL[u] = removedR[v] = true;

//             if(ok) {
//                 labelsL[u] = pL++;
//                 labelsR[v] = pR++;
//                 continue;
//             }

//             if(degU < p && degV < q) {
//                 for(uint32_t j = p[0][u]; j < p[0][u + 1]; j++) {
//                     if(!removedR[e1[j]]) {
//                         rheap.decrement(e1[j]);
//                         d2[e1[j]]--;
//                     }
//                 }
//                 for(uint32_t j = pV[v]; j < pV[v + 1]; j++) {
//                     if(!removedL[e2[j]]) {
//                         lheap.decrement(e2[j]);
//                         d1[e2[j]]--;
//                     }
//                 }
//             }
//             else if(degU < p && degV >= q) {
//                 for(uint32_t j = p[0][u]; j < p[0][u + 1]; j++) {
//                     if(e1[j] == v) {
//                         degV--;
//                         d2[v]--;
//                     }
//                     if(!removedR[e1[j]]) {
//                         rheap.decrement(e1[j]);
//                         d2[e1[j]]--;
//                     }
//                 }
                
//                 if(degV >= q) {
//                     removedR[v] = false;
//                     rheap.insert(v, degV);
//                 }
//             }
//             else if(degU >= p && degV < q) {
//                 for(uint32_t j = pV[v]; j < pV[v + 1]; j++) {
//                     if(u == e2[j]) {
//                         degU--;
//                         d1[u]--;
//                     }
//                     if(!removedL[e2[j]]) {
//                         lheap.decrement(e2[j]);
//                         d1[e2[j]]--;
//                         // lheap.decrement(e2[j]);
//                         // d1[e2[j]]--;
//                     }
//                 }

//                 if(degU >= p) {
//                     removedL[u] = false;
//                     lheap.insert(u, degU);
//                 }
//             }
//             else {
// // printf("ok:%u %u, %u %u\n", u, degU, v, degV);
// // fflush(stdout);
//                 ok = true;
//                 labelsL[u] = pL++;
//                 labelsR[v] = pR++;
//             }
//         }

//         if(!ok) {
//             n1 = n2 = m = 0;
//             return;
//         }

// printf("there %d %d\n", pL, pR);

//         if(n1 < n2) {
//             for(uint32_t j = n1; j < n2; j++) {
//                 uint32_t v, degV;
//                 if(!rheap.pop_min(v, degV)) printf("errorRheap\n");
//                 labelsR[v] = pR++;
//             }
//         }
//         else if(n1 > n2) {
//             for(uint32_t j = n2; j < n1; j++) {
//                 uint32_t u, degU;
//                 if(!lheap.pop_min(u, degU)) printf("errorRheap\n");
//                 labelsL[u] = pL++;
//             }
//         }


#endif