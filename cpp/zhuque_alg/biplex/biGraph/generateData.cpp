#include "biGraph.hpp"

int main(int argc, char * argv[])
{
    // "../../../Biclique/data/movielens-10m_ti/movie.txt"
    biGraph * g = new biGraph(argv[1], std::stoi(argv[2]),
             "two");

    printf("%u %u %u\n", g->n[0] + g->n[1], g->n[1], g->m);

    for(uint32_t i = 0; i < g->n[0]; i++) {
        printf("%u", i);
        for(uint32_t j = g->p[0][i]; j < g->p[0][i + 1]; j++) {
            uint32_t v = g->e[0][j];

            printf(" %u", v + g->n[0]);
        }

        printf("\n");
    }

    for(uint32_t i = 0; i < g->n[1]; i++) {
        printf("%u", g->n[0] + i);
        for(uint32_t j = g->p[1][i]; j < g->p[1][i + 1]; j++) {
            uint32_t v = g->e[1][j];

            printf(" %u", v);
        }

        printf("\n");
    }

    delete g;

    return 0;
}

// g++ generateData.cpp -std=c++11 -O3 -o generateData -I ../tools/fastIO.hpp ../tools/listLinearHeap.hpp ../tools/hopstotchHash.hpp -march=native -mavx