#include "tree.h"
using namespace std;

int main() {
    int species = 10;
    int seqlen = 1000;
    uint64_t seed = 57279;
    // uint64_t seed = (time(0) % 100000);
    Tree tree1(seed, species, seqlen, 0.1, jc69);
    Tree tree2(seed+1, species, seqlen, 0.1, jc69);
    Tree tree3(seed, species, seqlen, 0.1, jc69);
    tree1.generateRandomTopology();
    tree2.generateRandomTopology();
    tree1.dfsSequenceGen();
    tree2.dfsSequenceGen();
    tree3.mixSequences(tree1, tree2);
    cout << tree1.toNewick() << endl;
    cout << tree2.toNewick() << endl;
    tree3.writeToNexus();
    // cout << tree3.toNewick() << endl;

    return 0;
}