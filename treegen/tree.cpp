#include "tree.h"
#include <stack>
using namespace std;

void Tree::generateTopology() {
    uint64_t seed = init(orig + 1);
    vertices[0] = Vertex("1");
    vertices[1] = Vertex("2");
    edges[0] = Edge(0, 1);
    for (int i = 2; i < species; i++) {
        vertices[2 * i - 2] = Vertex(to_string(i+1));
        vertices[2 * i - 1] = Vertex();
        edges[2 * i - 3] = Edge(2 * i - 2, 2 * i - 1);
        int e = nextInt(&seed, 2 * i - 4);
        edges[2 * i - 2] = Edge(2 * i - 1, edges[e].v2);
        edges[e].v2 = 2 * i - 1;
    }
}

void Tree::dfsSequenceGen(double p_mutate, mutation_model smm) {
    uint64_t seed = init(orig + 2);
    stack<int> s;
    bool visited[2 * species - 2];
    for (int i = 0; i < 2 * species - 2; i++)
        visited[i] = false;
    s.push(0); // starting node doesn't matter
    for (int i = 0; i < seqlen; i++) sequences[0][i] = nextInt(&seed, 4); // randomize 1st node
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

void Tree::dropData(double pct_missing) {
    uint64_t seed = init(orig + 3);
    for (int node = 0; node < 2 * species - 2; node++) {
        if (!vertices[node].isLeaf) continue;
        for (int c = 0; c < seqlen; c++) {
            if (nextFloat(&seed) < pct_missing) sequences[node][c] = -1;
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
        if (!vertices[v].isLeaf) continue;
        file << vertices[v].name << "\t";
        for (int i = 0; i < seqlen; i++) file << convert_nucleotide(sequences[v][i]);
        file << endl;
    }
    file << ";" << endl;
    file << "end;" << endl;
    file.close();
}

string Tree::toNewick() {
    string newick = "";
    bool visited[2 * species - 2];
    for (int i = 0; i < 2 * species - 2; i++) visited[i] = false;
    recursiveNewick(newick, 3, visited);
    return newick;
}

void Tree::recursiveNewick(string& newick, int id, bool* visited) {
        visited[id] = true;
        if (vertices[id].isLeaf) {
            newick += vertices[id].name;
        }
        else {
            newick += "(";
            bool first = false;
            for (int i = 0; i < 2 * species - 3; i++) {
                int v = edges[i].has(id);
                if (v != -1 && !visited[v]) {
                    if (first) newick += ",";
                    else first = true;
                    recursiveNewick(newick, v, visited);
                }
            }
            newick += ")";
        }
}
