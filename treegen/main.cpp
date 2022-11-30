#include <iostream>
#include <fstream>
#include <stack>
#include "rng.h"
using namespace std;

static const string SPECIES_NAMES[10] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"};

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

char convert_nucleotide(int id);

int main(int argc, char** argv) {
    uint64_t seed =     (argc >= 2) ? stoll(argv[1]) : time(0);
    int species =       (argc == 5) ? stoi(argv[2]) : 10; // default # species = 10
    int seq_length =    (argc == 5) ? stoi(argv[3]) : 1000; // default sequence length = 1000
    double p_mutate =   (argc == 5) ? stod(argv[4]) : 0.2; // default site mutation probability = 0.2

    seed %= 10000; // reduce seed mod 10000 for simplification
    string orig = to_string(seed);
    seed = init(seed);
    Vertex vertices[2*species-2];
    Edge edges[2*species-3];
    int order[2*species-2];
    int sequences[2*species-2][seq_length];

    // generate topology
    vertices[0] = Vertex(SPECIES_NAMES[0]);
    vertices[1] = Vertex(SPECIES_NAMES[1]);
    edges[0] = Edge(0, 1);
    for (int i = 2; i < species; i++) {
        vertices[2*i-2] = Vertex(SPECIES_NAMES[i]);
        vertices[2*i-1] = Vertex();
        edges[2*i-3] = Edge(2*i-2, 2*i-1);
        int e = nextInt(&seed, 2*i-4);
        edges[2*i-2] = Edge(2*i-1, edges[e].v2);
        edges[e].v2 = 2*i-1;
    }

    // log tree info
    ofstream file;
    file.open("results/tree_" + orig + ".txt");
    // file << "root: " << vertices[start].name << endl;
    // file << "edges: " << endl;
    for (int i = 0; i < 2*species-3; i++) {
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
    bool visited[2*species-2];
    for (int i = 0; i < 2*species-2; i++) visited[i] = false;
    int start = nextInt(&seed, species);
    if (start > 1) start = 2*start-2;
    s.push(start);
    for (int i = 0; i < seq_length; i++) sequences[start][i] = nextInt(&seed, 4); // randomize start node
    while(!s.empty()) {
        int id = s.top();
        visited[id] = true;
        s.pop();
        for (int i = 0; i < 2*species-3; i++) {
            int v = edges[i].has(id);
            if (v != -1 && !visited[v]) {
                s.push(v);
                // mutate the sequence
                for (int i = 0; i < seq_length; i++) {
                    sequences[v][i] = // v = descendant
                        (nextFloat(&seed) < p_mutate) ? nextInt(&seed, 4) : sequences[id][i];
                }
            }
        }
    }
    
    // write the data
    file.close();
    file.open("results/data_" + orig + ".nex");
    file << "begin data;" << endl;
    file << "dimensions ntax=" << species << " nchar=" << seq_length << ";" << endl;
    file << "format datatype=dna interleave=no gap=-;" << endl;
    file << "matrix" << endl;
    for (int v = 0; v < 2*species-2; v++) {
        if (vertices[v].isLeaf) {
            file << vertices[v].name << "\t";
            for (int i = 0; i < seq_length; i++)
                file << convert_nucleotide(sequences[v][i]);
            file << endl;
        }
    }
    file << ";" << endl;
    file << "end;" << endl;
    file.close();
    
    return 0;
}

int Edge::has(int id) {
    if (v1 == id) return v2;
    if (v2 == id) return v1;
    return -1;
}

char convert_nucleotide(int id) {
    if (id == 0) return 'A';
    if (id == 1) return 'C';
    if (id == 2) return 'G';
    if (id == 3) return 'T';
    return '-';
}