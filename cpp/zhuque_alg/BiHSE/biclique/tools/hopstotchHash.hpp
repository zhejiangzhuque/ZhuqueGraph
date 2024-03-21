#ifndef HOPSTOTCH
#define HOPSTOTCH

#ifndef _WIN32
#include <x86intrin.h>
#else
#include <intrin.h>
#endif

#include <stdint.h>
#include <immintrin.h>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <utility>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

#define EMPTY 0xffffffff

constexpr uint32_t L2_CACHE_LINE_SIZE = 64;
constexpr uint32_t H = L2_CACHE_LINE_SIZE / sizeof(uint32_t);

#ifdef __AVX512__
#define setFunction _mm512_set1_epi32
#define cmpFunction _mm512_cmpeq_epi32_mask
#define loadFunction _mm512_loadu_si512
#endif

#ifdef __AVX512__
const __m256i eightEMPTY = _mm256_set1_epi32(EMPTY);
const __m512i sixteenEMPTY = _mm512_set1_epi32(EMPTY);
#endif

class hopstotchHash {
public:
    uint32_t * v;
    uint32_t n; //邻接点数量
    uint32_t roundN, preZero;
    uint32_t hashKey = 1e9 + 7;
    
public:
#ifndef __SSE__
#define __SSE__
#endif

#ifdef 	__SSE__
#ifndef __AVX512__
    bool containSIMD(uint32_t u) {
        uint32_t hashV = hash(u);
        uint32_t p = hashV;

        __m256i eightU = _mm256_set1_epi32(u);
        
        auto address = reinterpret_cast<const __m256i *>(v + p);
        // while(true) {
        __m256i eightNextV = _mm256_loadu_si256(address);
        __m256i cmpRes = _mm256_cmpeq_epi32(eightU, eightNextV);
        auto msk = _mm256_movemask_epi8(cmpRes);
#ifdef _WIN32
        if(__popcnt(msk) > 0) return true;
#else
        if(_popcnt32(msk) > 0) return true;
#endif
        eightNextV = _mm256_loadu_si256(address + 1);
        cmpRes = _mm256_cmpeq_epi32(eightU, eightNextV);
        msk = _mm256_movemask_epi8(cmpRes);
#ifdef _WIN32
        if(__popcnt(msk) > 0) return true;
#else
        if(_popcnt32(msk) > 0) return true;
#endif 
    
        return false;
    }
#else
    bool containSIMD(uint32_t u) {
        uint32_t hashV = hash(u);
        uint32_t p = hashV;

        __m256i eightU = _mm256_set1_epi32(u);
        
        // while(true) {
        __m256i eightNextV = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(v + p));
        __mmask8 cmpRes = _mm256_cmpeq_epi32_mask(eightU, eightNextV);
        if(cmpRes) return true;    
        
        eightNextV = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(v + p + 8));
        cmpRes = _mm256_cmpeq_epi32_mask(eightU, eightNextV);
        if(cmpRes) return true;    

        return false;
    }
#endif
#endif

    ~hopstotchHash() {
        delete [] v;
    }
    hopstotchHash():n(0) {}

    void build(uint32_t * adjList, uint32_t n_) {
        n = n_;
        // preZero = __builtin_clz(2 * n - 1);
        // if(n == 0) return;

        if(n > 0) {
            preZero = 0;
            while((1u<<preZero) <= 2*n-1) preZero++;
            preZero = 32 - preZero;

            roundN = 1u << (32 - preZero);
            if(roundN < H) roundN = H;
            v = new uint32_t[roundN + H];
        }
        else {
            preZero = 32;
            roundN = 0;
            v = new uint32_t[H];
        }
        
// printf("zero:%u roundN%u, n%u\n", preZero, roundN, n);fflush(stdout);
        memset(v, EMPTY, sizeof(uint32_t) * (roundN + H));
        for(uint32_t i = 0; i < n; i++) {
            uint32_t u = adjList[i];
            if(!insert(u)) {
                bool f = reBuild();
                if(!f) {
                    printf("error build hash table\n");
                    fflush(stdout);
                    exit(-1);
                }
                else insert(u);
            }
        }
        // memcpy(v + roundN, v, sizeof(uint32_t)*(H-1));
// printf("zero:%u roundN%u\n", preZero, roundN);
        for(uint32_t i = 0; i < H - 1; i++) {
            v[i+roundN] = v[i];
        }

        for(uint32_t i = 0; i < n; i++) {
            uint32_t u = adjList[i];
            if(!contain(u)) {
                printf("error build hash table 2\n");
                fflush(stdout);
                exit(-1);
            }
        }
        // for(uint32_t i = 0; i < H - 1; i++) {
        //     if(v[i] != v[i + roundN]) {
        //         printf("error st and ed of hashTable\n");fflush(stdout);
        //     }
        // }
    }

    bool insert(uint32_t u) {
        uint32_t hashV = hash(u);
        uint32_t p = hashV;
        uint32_t i = 0;
        for(; i < roundN; i++) {
            if(v[p] == EMPTY) break;
            p = (p + 1) % roundN;
            // if(v[p] == u) return true;
        }

        while((p - hashV + roundN) % roundN >= H) {
            uint32_t t = (p - H + 1 + roundN) % roundN;
            bool f = false;

            for(uint32_t i = t; i != p; i = (i + 1) % roundN) {
                if((p - hash(v[i]) + roundN) % roundN < H) {
                    v[p] = v[i];
                    v[i] = EMPTY;
                    p = i;
                    f = true;
                    break;
                }
            }

            if(!f) return false;
        }

        v[p] = u;
        return true;
    }

    bool containNormal(uint32_t u) {
        uint32_t hashV = hash(u);

        for(uint32_t i = hashV; i < hashV + H; i++) {
            // uint32_t j = i % roundN;
            if(v[i] == EMPTY) return false;
            if(v[i] == u) return true;
        }

        return false;
    }

    bool contain(uint32_t u) {

#ifdef __AVX512__
        return containSIMD512(u);
#else
        return containSIMD(u);
        // return containNormal(u);
#endif
    }
#ifdef __AVX512__
    bool containSIMD512(uint32_t u) {
        uint32_t hashV = hash(u);
        uint32_t p = hashV;

        __m512i sixteenU = _mm512_set1_epi32(u);
        
        // while(true) {
        __m512i sixteenNextV = _mm512_loadu_si512(reinterpret_cast<const __m512i *>(v + p));
        __mmask16 cmpRes = _mm512_cmpeq_epi32_mask(sixteenU, sixteenNextV);
        if(cmpRes) return true;    
        // if(_mm512_cmpeq_epi32_mask(sixteenNextV, sixteenEMPTY)) return false;

        //     p += parallelism;
        // }
        
        return false;
    }
#endif
    bool reBuild() {
        uint32_t * tmpV = v;
        v = new uint32_t[roundN];
        memset(v, 0xff, sizeof(uint32_t) * roundN);

        uint32_t rebuildTimes = 0;
        while(rebuildTimes < 100) {
            bool f = true;
            for(uint32_t i = rebuildTimes; i < roundN + rebuildTimes; i++) {
                uint32_t j = i % roundN;
                if(tmpV[j] != EMPTY) {
                    if(!insert(tmpV[j])) {
                        f = false;
                        break;
                    }
                }
            }

            if(f) break;
            else {
                rebuildTimes++;
                memset(v, 0xff, sizeof(uint32_t) * roundN);
            }
        }
        
        delete [] tmpV;

        if(rebuildTimes == 100) return false;

        return true;
    }
//h(k) = (A∗k mod 2^w) >> (w − r), w = 32, A常数，r=log(roundN)
    uint32_t hash(uint32_t t) {
        uint32_t tmp = hashKey * t;
        uint32_t ans = tmp >> preZero;
        // assert(ans < roundN);
        return ans;
    }
};

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

#endif