#include "tree.h"
using namespace std;

int main() {
    int species = 5;
    int seqlen = 1000;
    // uint64_t seed = 57279;
    uint64_t seed = (time(0) % 100000);
    Tree tree1(seed, species, seqlen);
    Tree tree2(seed, species, seqlen);
    Tree tree3(seed, species, seqlen);
    tree1.generateRandomTopology();
    tree1.randomizeBranchLengths(0.01,0.5);
    // tree1.printEdges();
    tree2.copyTopology(tree1);
    tree2.randomizeBranchLengths(0.01,0.5);
    tree1.dfsSequenceGen();
    tree2.dfsSequenceGen();
    tree3.mixSequences({tree1, tree2});
    // cout << tree1.toNewick() << endl;
    // cout << tree2.toNewick() << endl;
    tree3.writeToNexus();
    // cout << tree3.toNewick() << endl;

    return 0;
}