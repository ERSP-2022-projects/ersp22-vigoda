#include "tree.h"
using namespace std;

int main() {
    int species = 10;
    int seqlen = 20000;
    // uint64_t seed = 45346;
    uint64_t seed = (time(0) % 100000);
    Tree tree1(seed, species, seqlen);
    Tree tree2(seed, species, seqlen);
    Tree tree3(seed, species, seqlen);
    tree1.generateRandomTopology();
    // cout << "1" << endl;
    tree1.randomizeBranchLengths(0.005,0.5);
    // cout << "2" << endl;
    // tree1.printEdges();
    tree2.copyTopology(tree1);
    // cout << "3" << endl;
    tree2.randomizeBranchLengths(0.005,0.5);
    // cout << "4" << endl;
    tree1.dfsSequenceGen();
    // cout << "5" << endl;
    // cout << "" << endl;
    tree2.dfsSequenceGen();
    // cout << "6" << endl;
    tree3.mixSequences({tree1, tree2});
    // cout << "7" << endl;
    tree3.writeToNexus();
    cout << "exported results for " << to_string(seed) << endl;
    cout << tree1.toNewick() << endl;
    // cout << tree2.toNewick() << endl;
    // cout << tree3.toNewick() << endl;

    return 0;
}