#include "tree.h"
#include <fstream>
#include <filesystem>
#include <stack>
using namespace std;

void Tree::generateTopology() {
    uint64_t seed = init(orig + 1);
    vertices[0] = Vertex(SPECIES_NAMES[0]);
    vertices[1] = Vertex(SPECIES_NAMES[1]);
    edges[0] = Edge(0, 1);
    for (int i = 2; i < species; i++) {
        vertices[2 * i - 2] = Vertex(SPECIES_NAMES[i]);
        vertices[2 * i - 1] = Vertex();
        edges[2 * i - 3] = Edge(2 * i - 2, 2 * i - 1);
        int e = nextInt(&seed, 2 * i - 4);
        edges[2 * i - 2] = Edge(2 * i - 1, edges[e].v2);
        edges[e].v2 = 2 * i - 1;
    }
}

void Tree::logTreeTxt() {
    ofstream file;
    file.open("results/tree_" + to_string(orig) + ".txt");
    for (int i = 0; i < 2 * species - 3; i++) {
        Edge e = edges[i];
        if ((vertices[e.v1].isLeaf))
            file << vertices[e.v1].name;
        else
            file << e.v1;
        file << ", ";
        if ((vertices[e.v2].isLeaf))
            file << vertices[e.v2].name;
        else
            file << e.v2;
        file << endl;
    }
    file.close();
}

void Tree::dfsSequenceGen() {
    uint64_t seed = init(orig + 2);
    stack<int> s;
    bool visited[2 * species - 2];
    for (int i = 0; i < 2 * species - 2; i++)
        visited[i] = false;
    s.push(0); // starting node doesn't matter
    for (int i = 0; i < seqlen; i++)
        sequences[0][i] = nextInt(&seed, 4); // randomize 1st node
    while (!s.empty()) {
        int id = s.top();
        visited[id] = true;
        s.pop();
        for (int i = 0; i < 2 * species - 3; i++) {
            int v = edges[i].has(id);
            if (v != -1 && !visited[v]) {
                s.push(v);
                // mutate the sequence
                // note id = ancestor, v = descendant
                for (int i = 0; i < seqlen; i++) {
                    sequences[v][i] = sequences[id][i];
                    if (nextFloat(&seed) < p_mutate) {
                        sequences[v][i] += roll_nucleotide(smm, &seed);
                        sequences[v][i] %= 4;
                    }
                }
            }
        }
    }
}

void Tree::writeToNexus() {
    // filesystem::create_directory("results/" + to_string(orig));
    ofstream file;
    file.open("results/" + to_string(orig) + "_data.nex");
    file << "begin data;" << endl;
    file << "dimensions ntax=" << species << " nchar=" << seqlen << ";" << endl;
    file << "format datatype=dna interleave=no gap=-;" << endl;
    file << "matrix" << endl;
    for (int v = 0; v < 2 * species - 2; v++) {
        if (!vertices[v].isLeaf)
            continue;
        file << vertices[v].name << "\t";
        for (int i = 0; i < seqlen; i++)
            file << convert_nucleotide(sequences[v][i]);
        file << endl;
    }
    file << ";" << endl;
    file << "end;" << endl;
    file.close();
}
