#ifndef TREE_H
#define TREE_H

#include <string>
using namespace std;

struct Vertex {
    bool isLeaf;
    string name;
    Vertex() : isLeaf(false), name("-") {}
    Vertex(string name) : isLeaf(true), name(name){};
};

struct Edge {
    int v1, v2;
    Edge(int v1 = -1, int v2 = -1) : v1(v1), v2(v2) {}
    int has(int id) {
        if (v1 == id)
            return v2;
        if (v2 == id)
            return v1;
        return -1;
    }
};

#endif