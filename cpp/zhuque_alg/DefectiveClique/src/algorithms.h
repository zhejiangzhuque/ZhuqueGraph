#pragma once
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <queue>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cassert>
#include <cmath>
#include <utility>
#include <immintrin.h>
#include <climits>

using namespace std;

typedef int int32;
typedef char char8;
typedef bool boolean;
typedef vector<int32> vectori;
typedef vector<boolean> vectorb;
typedef unsigned int uint32;

class Algorithm
{
protected:
	int32 n, m, md;
    int32 *deg, *datas, **adj;

    vectori core, topcore, colors;
    
    int32 mincliquesize;
    int32 algorithm;

    /* data */
public:
    Algorithm(/* args */);
    virtual ~Algorithm();
    virtual void run() {}
    virtual void setParameters(int argc, char **argv) {}

    void read_graph(const char *str);

    // void read_graph(const char *str);
    void scalability(bool randomv, float scal);
    void testprintGraph();

    // void setparemeters(double eta, int alg, int k) {
    //     this->eta = eta; this->algorithm = alg;
	// 	mincliquesize = k;
    //     printf("eta=%lf, alg=%d, mincliquesize=%d\n", eta, alg, k);
    // }

    int core_decompsition(int *nodeset, int nodesize);
	int coloring(int *nodeset, int nodesize);
};

class Bitset
{
private:
    /* data */
public:

    uint32 size, capacity;
    uint32 *buff;
    uint32 _type_size = sizeof(uint32) + 1;
    uint32 _mask = (1 << _type_size) - 1;

    Bitset(/* args */) {
        size = capacity = 0;
        buff = NULL;
    }
    ~Bitset() {
       if (buff) delete[] buff; 
    }

    void allocacte(uint32 _cap) {
        _cap = (_cap >> _type_size) + 1;
        if (capacity >= _cap) return;
        if (buff) delete[] buff;
        capacity = _cap;
        buff = new uint32[_cap];
    }
    void init(uint32 _size) {
        _size = (_size >> _type_size) + 1;
        size = _size;
        memset(buff, uint32(0), sizeof(uint32) * size);
    }
    void insert(uint32 id) {
        //assert((id >> _type_size) < capacity);
        buff[id >> _type_size] ^= 1 << (id & _mask);
    }
    bool find(uint32 id) {
        return buff[id >> _type_size] >> (id & _mask) & 1;
    }
    bool empty() { 
        return size == 0;
    }
};


#define unfilled -1
class CuckooHash
{
private:
	/* data */
	int32 capacity;
	int32 mask;
	int32 size;
	int32 buff_size = sizeof(int32);
	int32 *hashtable;

	void rehash(int32 **_table) {
		int32 oldcapacity = capacity;
		mask = mask == 0 ? 1 : ((mask << 1) | 1);
		capacity = (mask + 1) * buff_size;
		int32 *newhash = new int32[capacity];
		memset((newhash), unfilled, sizeof(int32) * capacity);
		for (int32 i = 0; i < oldcapacity; ++i){
			if ((*_table)[i] != unfilled) insert((*_table)[i], &newhash);
		}
		swap((*_table), newhash);
		delete[] newhash;
	}
	void insert(const int32 &_u, int32 **_table) {
		
		int32 hs = hash1(_u);
		for (int32 i = 0; i < buff_size; ++i) {
			if ((*_table)[hs * buff_size + i] == unfilled){
				(*_table)[hs * buff_size + i] = _u;
				return;
			}
		}
		hs = hash2(_u);
		for (int32 i = 0; i < buff_size; ++i) {
			if ((*_table)[hs * buff_size + i] == unfilled){
				(*_table)[hs * buff_size + i] = _u;
				return;
			}
		}

		bool use_hash1 = true;
		int32 u = _u;
		for (int32 i = 0; i < mask; ++i) {
			int32 replaced;
			if (use_hash1) hs = hash1(u);
			else hs = hash2(u);
			int32 j = 0;
			for (; j < buff_size; ++j) {
				if ((*_table)[hs * buff_size + j] == unfilled) break;
			}
			if (buff_size == j) {
				replaced = move((*_table)[hs * buff_size]);
				j = 1;
				for (; j < buff_size; j++) {
					(*_table)[hs * buff_size + j - 1] =
						move((*_table)[hs * buff_size + j]);
				}
				(*_table)[hs * buff_size + j - 1] = u;
			}
			else {
				replaced = move((*_table)[hs * buff_size + j]);
				(*_table)[hs * buff_size + j] = u;
			}
			use_hash1 = hs == hash2(replaced);
			u = move(replaced);
			if (u == unfilled) return;
		}
		rehash(_table);
		insert(u, _table);
	}

	int32 hash1(const int32 &x) { return x & mask;}
	int32 hash2(const int32 &x) { return ~x & mask;}

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

	void reserve(int32 _size) {
		if (capacity >= _size) return;
		mask = mask == 0 ? 1 : ((mask << 1) | 1);
		while (_size >= mask * buff_size) mask = (mask << 1) | 1;
		capacity = (mask + 1) * buff_size;
		if (hashtable) delete[] hashtable;
		hashtable = new int32[capacity];
		memset(hashtable, unfilled, sizeof(int32) * capacity);
	}

	void insert(const int32 &_u) {
		if (find(_u)) return;
		insert(_u, &hashtable);
		size++;
	}

	bool find(const int32 &_u) {
		int32 hs1 = hash1(_u);
		int32 hs2 = hash2(_u);

		assert(buff_size == 4 && sizeof (int32) == 4);
		__m128i cmp = _mm_set1_epi32(_u);
		__m128i b1 = _mm_load_si128((__m128i*)&hashtable[buff_size * hs1]);
		__m128i b2 = _mm_load_si128((__m128i*)&hashtable[buff_size * hs2]);
		__m128i flag = _mm_or_si128(_mm_cmpeq_epi32(cmp, b1), _mm_cmpeq_epi32(cmp, b2));

		return _mm_movemask_epi8(flag) != 0;
	}
	int32 getcapacity() {return capacity;}
	int32 getsize() {return size;}
	int32 getmask() {return mask;}
	int32 *gethashtable() {return hashtable;}
};

typedef struct _Node
{
    int size, exsize; 
    vector<int> nodes;
    vector<int> counts;

	_Node() {
		size 	= 0;
		exsize  = 0;
	}
	_Node(int n) {
		nodes.resize(n);
		counts.resize(n);
		size 	= n;
		exsize  = 0;
	}
	void resize(int n) {
		nodes.resize(n);
		counts.resize(n);
		size 	= n;
		exsize  = 0;
	}
	void reserve(int n) {
		nodes.reserve(n);
		counts.reserve(n);
		size 	= 0;
		exsize  = 0;
	}
	void cinsert(int x, int c) {
		nodes.emplace_back(x);
		counts.emplace_back(c);
		size++;
	}
	void einsert(int x, int c) {
		nodes.emplace_back(x);
		counts.emplace_back(c);
		size++;
		exsize++;
	}
	void cinsert(int id, int x, int c) {
		nodes[id]  = x;
		counts[id] = c;
		size++;
	}
	void einsert(int id, int x, int c) {
		nodes[id]  = x;
		counts[id] = c;
		size++;
		exsize++;
	}
	void clear() {
		size	= 0;
		exsize  = 0;
		nodes.clear();
		counts.clear();
	}
} Node;