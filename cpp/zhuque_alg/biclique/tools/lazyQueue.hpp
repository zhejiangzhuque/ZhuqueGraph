#ifndef LAZYQUEUE
#define LZAYQUEUE

#include <queue>
#include <vector>

template<class T>
struct lazyQueue {
    std::priority_queue<T, std::vector<T>, std::greater<T>> qa, qb;

    void push(T v) {
        qa.push(v);
    }

    void remove(T v) {
        qb.push(v);
    }

    T top() {
        flush();
        return qa.top();
    }

    void pop() {
        flush();
        qa.pop();
    }

    void flush() {
        while(!qa.empty() && !qb.empty() && qa.top() == qb.top()) {
            qa.pop();
            qb.pop();
        }
    }

    uint32_t size() {
        flush();
        return qa.size() - qb.size();
    }
};

#endif
