#include "tree.h"
#include <fstream>
#include <stack>
using namespace std;

void Tree::makeVertices() {
    for (int i = 0; i < species; i++) {
        vertices[i] = Vertex(to_string(i+1));
    }
}

void Tree::generateRandomTopology() {
    uint64_t seed = init(orig + 1);
    edges[0] = Edge(0, 1);
    for (int i = 2; i < species; i++) {
        edges[2 * i - 3] = Edge(i, species+i-2);
        int e = nextInt(&seed, 2 * i - 4);
        edges[2 * i - 2] = Edge(species+i-2, edges[e].v2);
        edges[e].v2 = species+i-2;
    }
}

void Tree::generateTopology(vector<int> edgearr) {
    if (edgearr.size() != 4 * species - 6) {
        cerr << "# of edges does not match" << endl;
        return;
    }
    for (int i = 0; i < 2 * species - 3; i++) {
        edges[i] = Edge(edgearr[2*i]-1,edgearr[2*i+1]-1);
    }
}

void Tree::copyTopology(Tree& tree) {
    if (tree.species != species) {
        cerr << "# of leaves/edges does not match" << endl;
        return;
    }
    this->edges = tree.edges;
    // for (int i = 0; i < 2 * species - 3; i++) {
    //     edges[i] = Edge(tree.edges[i].v1, tree.edges[i].v2);
    // }
}

void Tree::setBranchLengths(vector<double> lengths) {
    if (lengths.size() != 2 * species - 3) {
        cerr << "# of branch lengths does not match # of edges" << endl;
        return;
    }
    for (int i = 0; i < 2 * species - 3; i++) {
        edges[i].branchLength = lengths[i];
    }
}

void Tree::randomizeBranchLengths(double min, double max) {
    uint64_t seed = init(orig + 3);
    for (int i = 0; i < 2 * species - 3; i++) {
        edges[i].branchLength = nextFloat(&seed) * (max - min) + min;
    }
}

void Tree::dfsSequenceGen() {
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
                    if (nextFloat(&seed) < edges[i].branchLength) {
                        sequences[v][i] += roll_nucleotide(smm, &seed);
                        sequences[v][i] %= 4;
                    }
                }
            }
        }
    }
}

void Tree::mixSequences(vector<Tree> trees) {
    uint64_t seed = init(orig + 3);
    for (int c = 0; c < seqlen; c++) {
        Tree use = trees[nextInt(&seed, trees.size())];
        for (int i = 0; i < 2 * species - 2; i++) {
            sequences[i][c] = use.sequences[i][c];
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
    for (int v = 0; v < species; v++) {
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
    recursiveNewick(newick, species, visited);
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

void Tree::printEdges() {
    for (Edge e : edges) {
        cout << "(" << e.v1 << ", " << e.v2 << "), " << e.branchLength << endl;
    }
}