#ifndef TREE_H
#define TREE_H

#include "mutationmodels.h"
#include <vector>
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
        if (v1 == id) return v2;
        if (v2 == id) return v1;
        return -1;
    }
};

class Tree {
    private:
        uint64_t orig;
        int species, seqlen;
        double p_mutate;
        mutation_model smm;
        vector<vector<int> > sequences;
        vector<Vertex> vertices;
        vector<Edge> edges;
        
        void recursiveNewick(string& newick, int id, bool* visited);

    public:
        Tree(uint64_t seed, const int species, const int seqlen, double p_mutate, mutation_model smm) : 
            orig(seed), species(species), seqlen(seqlen), p_mutate(p_mutate), smm(smm),
            vertices(2 * species - 2), edges(2 * species - 3) {
                sequences = vector<vector<int>>(2 * species - 2, vector<int>(seqlen));
            }
        void generateTopology();
        void dfsSequenceGen();
        void writeToNexus();
        string toNewick();
};

#endif
