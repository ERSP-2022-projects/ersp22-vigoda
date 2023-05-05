#ifndef TREE_H
#define TREE_H

#include "mutationmodels.h"
#include <vector>
using namespace std;

struct Vertex {
    bool isLeaf;
    string name;
    int edges[3] = {-1,-1,-1};
    Vertex() : isLeaf(false), name("-") {}
    Vertex(string name) : isLeaf(true), name(name) {};
};

struct Edge {
    int v1, v2;
    double branchLength;
    Edge(int v1 = -1, int v2 = -1, double branchLength = 0.1) : v1(v1), v2(v2), branchLength(branchLength) {}
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
        mutation_model smm;
        vector<vector<int> > sequences;
        vector<Vertex> vertices;
        vector<Edge> edges;
        
        void makeVertices();
        void splitEdge(int e, int v);
        void recursiveNewick(string& newick, int id, bool* visited);

    public:
        Tree(uint64_t seed, const int species, const int seqlen, double p_mutate = 0.1, mutation_model smm = jc69) : 
            orig(seed), species(species), seqlen(seqlen), smm(smm),
            vertices(2 * species - 2) {
                sequences = vector<vector<int>>(2 * species - 2, vector<int>(seqlen));
                edges = vector<Edge>(2 * species - 3, Edge(p_mutate));
                makeVertices();
            }
        void generateRandomTopology();
        void generateTopology(vector<int> edgearr);
        void copyTopology(Tree& tree);
        void setBranchLengths(vector<double> lengths);
        void randomizeBranchLengths(double min, double max);
        void dfsSequenceGen();
        void mixSequences(vector<Tree> trees);
        void writeToNexus();
        string toNewick();
        void printVE();
};

#endif
