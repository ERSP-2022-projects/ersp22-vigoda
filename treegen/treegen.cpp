#include "tree.h"
#include "NewickTree.h"
using namespace std;

int main(int argc, char **argv) {
    // initializations + cmd line reading
    uint64_t seed = (time(0) % 100000); // default seed = time % 1M
    int species = 10;                   // default # species = 10
    int seqlen = 1000;                  // default sequence length = 1000
    double p_mutate = 0.2;              // default mutation probability = 0.2
    mutation_model smm = jc69;          // default site mutation model = jc69
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        int eq = arg.find('=');
        if (eq == -1) {
            cerr << "argument \'" << arg << "\" is not a valid parameter assignment" << endl;
            continue;
        }
        string param = arg.substr(0, eq);
        string after = arg.substr(eq + 1);
        if (param == "seed")
            seed = stoll(after);
        else if (param == "species" || param == "numspecies")
            species = stoi(after);
        else if (param == "seqlen")
            seqlen = stoi(after);
        else if (param == "p_mutate")
            p_mutate = stod(after);
        else if (param == "smm" || param == "mutation_model")
            smm = stomm(after);
        else
            cerr << param << " is not a valid parameter" << endl;
    }

    Tree tree(seed, species, seqlen, p_mutate, smm);
    tree.generateTopology();
    tree.dfsSequenceGen();
    tree.writeToNexus();

    NewickTree newickFormatted(tree.getVertices(), tree.getEdges(), 3);
    newickFormatted.printNewick(true);
    newickFormatted.exportNewick("results/" + to_string(seed) +"_tree.txt", true);
    cout << "exported results to " << to_string(seed) << endl;

}
