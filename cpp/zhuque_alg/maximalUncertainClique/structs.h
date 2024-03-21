
#ifndef _INCLUDE_STRUCT_H
#define _INCLUDE_STRUCT_H

#include <cassert>
#include <cstdio>

typedef struct _node {
	int u;
	int pre, next;
	double p;
} node;
typedef struct _list {
	int size, capacity;
	int head, tail;
	node *Nodes;
} list;

static void node_push_back(node *nodes, const int pid, const int &u, const double &p);

static void list_init_size(list &listnodes, int n);
static void list_init_size(list &listnodes);
static void list_clear(list &listnodes);
static void list_free(list &listnodes);
static void list_push_back(list &listnodes, const int &u, const double &p);
static void list_insert(list &listnodes, int insertpos, int u, double p);

void node_push_back(node *nodes, int pid, const int &u, const double &p)
{
	nodes[pid].u = u;
	nodes[pid].p = p;
	nodes[pid].next = -1;
	nodes[pid].pre = pid - 1;
	nodes[pid-1].next = pid;
}

void list_init_size(list &listnodes, int n)
{
	listnodes.capacity = n;
	listnodes.head = -1;
	listnodes.tail = -1;
	listnodes.size = 0;
	listnodes.Nodes = new node[n];
}
void list_init_size(list &listnodes)
{
	listnodes.capacity = 0;
	listnodes.head = -1;
	listnodes.tail = -1;
	listnodes.size = 0;
	listnodes.Nodes = NULL;
}
void list_clear(list &listnodes)
{
	listnodes.head = -1;
	listnodes.tail = -1;
	listnodes.size = 0;
}
void list_free(list &listnodes)
{
	listnodes.capacity = 0;
	listnodes.head = -1;
	listnodes.tail = -1;
	listnodes.size = 0;
	if (listnodes.Nodes !=NULL) 
		delete[] listnodes.Nodes;
}
void list_push_back(list &listnodes, const int &u, const double &p)
{
	assert(listnodes.size < listnodes.capacity);
	listnodes.Nodes[listnodes.size].u = u;
	listnodes.Nodes[listnodes.size].p = p;
	listnodes.Nodes[listnodes.size].pre = listnodes.tail;
	listnodes.Nodes[listnodes.size].next = -1;
	if (listnodes.head == -1) listnodes.head = listnodes.size;
	if (listnodes.tail != -1) listnodes.Nodes[listnodes.tail].next = listnodes.size;
	listnodes.tail = listnodes.size++;
}

void list_insert(list &listnodes, int insertpos, int u, double p)
{
	assert(listnodes.size < listnodes.capacity);
    if (insertpos < 0 || listnodes.head == -1) 
        list_push_back(listnodes, u, p);
	else if (listnodes.head == insertpos) {
		listnodes.Nodes[listnodes.size].u = u;
		listnodes.Nodes[listnodes.size].p = p;
		listnodes.Nodes[listnodes.size].pre = -1;
		listnodes.Nodes[listnodes.size].next = listnodes.head;
		listnodes.Nodes[listnodes.head].pre = listnodes.size;
		listnodes.head = listnodes.size++;
	}
	else {
		listnodes.Nodes[listnodes.size].u = u;
		listnodes.Nodes[listnodes.size].p = p;
		listnodes.Nodes[listnodes.size].pre = listnodes.Nodes[insertpos].pre;
		listnodes.Nodes[listnodes.size].next = insertpos;
        listnodes.Nodes[listnodes.Nodes[insertpos].pre].next = listnodes.size;
		listnodes.Nodes[insertpos].pre = listnodes.size++;
	}
}

typedef struct _epro {
	double p;
	int u;
} epro;

static bool decreasefun(const epro &a, const epro &b)
{
    if (a.p > b.p + 1e-16) return true;
    else if (a.p + 1e-16 < b.p) return false;
    else return a.u > b.u;
}

// static bool decreasefun(double &a, double &b)
// {
//     if (a > b + 1e-16) return true;
//     else return false;
// }
#endif
