#include "generator.h"
#include "consts.h"
#include <stack>
#include <fstream>
using namespace std;

int Edge::has(int id) {
    if (v1 == id) return v2;
    if (v2 == id) return v1;
    return -1;
}

Generator::Generator(uint64_t seed) {
    this->seed = init(seed);
    this->orig = seed;
}

void Generator::generate_tree() {
    // generate topology
    vertices[0] = Vertex(SPECIES_NAMES[0]);
    vertices[1] = Vertex(SPECIES_NAMES[1]);
    edges[0] = Edge(0, 1);
    for (int i = 2; i < SPECIES; i++) {
        vertices[2*i-2] = Vertex(SPECIES_NAMES[i]);
        vertices[2*i-1] = Vertex();
        edges[2*i-3] = Edge(2*i-2, 2*i-1);
        int e = nextInt(&seed, 2*i-4);
        edges[2*i-2] = Edge(2*i-1, edges[e].v2);
        edges[e].v2 = 2*i-1;
    }

    // log tree info
    ofstream file;
    file.open("results/tree_" + to_string(orig) + ".txt");
    // file << "root: " << vertices[start].name << endl;
    // file << "edges: " << endl;
    for (int i = 0; i < 2*SPECIES-3; i++) {
        Edge e = edges[i];
        if ((vertices[e.v1].isLeaf)) file << vertices[e.v1].name;
        else file << e.v1;
        file << ", ";
        if ((vertices[e.v2].isLeaf)) file << vertices[e.v2].name;
        else file << e.v2;
        file << endl;
    }
    

    // dfs sequence generation
    stack<int> s;
    bool visited[2*SPECIES-2] = {false};
    int start = nextInt(&seed, SPECIES);
    if (start > 1) start = 2*start-2;
    s.push(start);
    random_sequence(start);
    while(!s.empty()) {
        int id = s.top();
        visited[id] = true;
        s.pop();
        for (int i = 0; i < 2*SPECIES-3; i++) {
            int v = edges[i].has(id);
            if (v != -1 && !visited[v]) {
                s.push(v);
                mutate_sequence(id, v);
                // write this edge
                // Edge e = edges[i];
                // if ((vertices[e.v1].isLeaf)) file << vertices[e.v1].name;
                // else file << e.v1;
                // file << ", ";
                // if ((vertices[e.v2].isLeaf)) file << vertices[e.v2].name;
                // else file << e.v2;
                // file << endl;
            }
        }
    }
    
    // write the data
    file.close();
    file.open("results/data_" + to_string(orig) + ".nex");
    file << "begin data;" << endl;
    file << "dimensions ntax=" << SPECIES << " nchar=" << SEQUENCE_LENGTH << ";" << endl;
    file << "format datatype=dna interleave=no gap=-;" << endl;
    file << "matrix" << endl;
    for (int v = 0; v < 2*SPECIES-2; v++) {
        if (vertices[v].isLeaf) {
            file << vertices[v].name << "\t";
            for (int i = 0; i < SEQUENCE_LENGTH; i++)
                file << convert_nucleotide(sequences[v][i]);
            file << endl;
        }
    }
    file << ";" << endl;
    file << "end;" << endl;
    file.close();
}

void Generator::random_sequence(int id) {
    for (int i = 0; i < SEQUENCE_LENGTH; i++)
        sequences[id][i] = nextInt(&seed, 4);
}

void Generator::mutate_sequence(int ancestor, int descendant) {
    for (int i = 0; i < SEQUENCE_LENGTH; i++) {
        sequences[descendant][i] = 
            roll_nucleotide(sequences[ancestor][i]);
    }
}

int Generator::roll_nucleotide(int orig) {
    if (nextFloat(&seed) < P_MUTATE) return nextInt(&seed, 4);
    else return orig;
}

char Generator::convert_nucleotide(int id) {
    if (id == 0) return 'A';
    if (id == 1) return 'C';
    if (id == 2) return 'G';
    if (id == 3) return 'T';
    return '-';
}