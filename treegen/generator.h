#ifndef GENERATOR_H
#define GENERATOR_H

#include "rng.h"
#include "consts.h"
#include <iostream>
using namespace std;

struct Vertex {
    bool isLeaf;
    string name;
    Vertex() : isLeaf(false), name("-") {}
    Vertex(string name) : isLeaf(true), name(name) {};
};

struct Edge {
    int v1, v2;
    Edge(int v1 = -1, int v2 = -1) : v1(v1), v2(v2) {}
    int has(int id);
};

class Generator {

    private: 
        uint64_t seed, orig;
        Vertex vertices[2*SPECIES-2];
        Edge edges[2*SPECIES-3];
        int order[2*SPECIES-2];
        int sequences[2*SPECIES-2][SEQUENCE_LENGTH];
        void random_sequence(int id);
        void mutate_sequence(int ancestor, int descendant);
        int roll_nucleotide(int orig);
        char convert_nucleotide(int id);

    public:
        Generator(uint64_t seed = time(0));
        void generate_tree();

};

#endif