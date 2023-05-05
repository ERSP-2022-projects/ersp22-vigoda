#include "tree.h"
#include <fstream>
#include <stack>
using namespace std;

void Tree::makeVertices() {
    for (int i = 0; i < species; i++) {
        vertices[i] = Vertex(to_string(i+1));
    }
}

void Tree::splitEdge(int e, int v) {
    int in_node = species + v - 2;
    int ext_edge = 2 * v - 3;
    int in_edge = 2 * v - 2;
    int moving_node = edges[e].v2;
    edges[ext_edge] = Edge(v, in_node);
    edges[in_edge] = Edge(in_node, moving_node);
    edges[e].v2 = in_node;
    vertices[v].edges[0] = ext_edge;
    vertices[in_node].edges[0] = e;
    vertices[in_node].edges[1] = in_edge;
    vertices[in_node].edges[2] = ext_edge;
    for (int i = 0; i < 3; i++)
        if (vertices[moving_node].edges[i] == e) vertices[moving_node].edges[i] = in_edge;
}

void Tree::generateRandomTopology() {
    uint64_t seed = init(orig + 1);
    edges[0] = Edge(0, 1);
    vertices[0].edges[0] = 0;
    vertices[1].edges[0] = 0;
    for (int i = 2; i < species; i++) {
        int e = nextInt(&seed, 2 * i - 4); // randomly select an edge to split
        splitEdge(e, i);
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
    this->vertices = tree.vertices;
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
    for (int i = 0; i < 2 * species - 2; i++) visited[i] = false;
    s.push(0); // starting node doesn't matter
    for (int i = 0; i < seqlen; i++) sequences[0][i] = nextInt(&seed, 4); // randomize 1st node
    while (!s.empty()) {
        int id = s.top();
        visited[id] = true;
        s.pop();
        for (int e = 0; e < (vertices[id].isLeaf ? 1 : 3); e++) {
            int edge = vertices[id].edges[e];
            int v = edges[edge].has(id);
            if (v != -1 && !visited[v]) {
                s.push(v);
                // mutate the sequence
                // note id = ancestor, v = descendant
                for (int i = 0; i < seqlen; i++) {
                    sequences[v][i] = sequences[id][i];
                    if (nextFloat(&seed) < edges[edge].branchLength) {
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
        for (int i = 0; i < 2 * species - 2; i++) {
            sequences[i][c] = trees[nextInt(&seed, trees.size())].sequences[i][c];
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
            for (int i = 0; i < (vertices[id].isLeaf ? 1 : 3); i++) {
            int v = edges[vertices[id].edges[i]].has(id);
                if (v != -1 && !visited[v]) {
                    if (first) newick += ",";
                    else first = true;
                    recursiveNewick(newick, v, visited);
                }
            }
            newick += ")";
        }
}

void Tree::printVE() {
    cout << "vertices" << endl;
    for (int i = 0; i < 2 * species - 2; i++) {
        Vertex v = vertices[i];
        cout << i << ": (" << v.edges[0] << ", " << v.edges[1] << ", " << v.edges[2] << ", " << ")" << endl;
    }
    cout << "edges" << endl;
    for (int i = 0; i < 2 * species - 3; i++) {
        Edge e = edges[i];
        cout << i << ": (" << e.v1 << ", " << e.v2 << "), " << e.branchLength << endl;
    }
}