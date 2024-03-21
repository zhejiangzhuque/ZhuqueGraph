#include "listLinearHeap.hpp"

int main()
{
    printf("shas\n");
    int n = 10, m = 10;
    ListLinearHeap p(n, m);
    unsigned int x[100], y[100];

    for(unsigned int i = 0; i < n; i++) {
        x[i] = i;
        y[i] = i;
    }
    p.init(n, m, x, y);
// printf("there\n");
    while(!p.empty()) {
        unsigned a, b;
        p.pop_min(a, b);

        if(a != 5) {
            printf("d %u\n",  p.decrement(5));
        }
        printf("x%u, %u\n", a, b);
    }

    return 0;
}