#include "tree.h"
using namespace std;

int main(int argc, char **argv)
{

    // initializations + cmd line reading
    uint64_t seed = (time(0) % 1000000); // default seed = time % 1M
    int species = 10;                    // default # species = 10
    int seqlen = 1000;                   // default sequence length = 1000
    double p_mutate = 0.2;               // default mutation probability = 0.2
    mutation_model smm = jc69;           // default site mutation model = jc69
    double bl_lower_bound = 0.01;
    double bl_upper_bound = 0.5;

    for (int i = 1; i < argc; i++)
    {
        string arg = argv[i];
        int eq = arg.find('=');
        if (eq == -1)
        {
            cerr << "argument \'" << arg << "\" is not a valid parameter assignment" << endl;
            continue;
        }
        string param = arg.substr(0, eq);
        string after = arg.substr(eq + 1);
        if (param == "seed")
            seed = stoll(after);
        else if (param == "ntax" || param == "species" || param == "numspecies" || param == "leaves" || param == "n")
            species = stoi(after);
        else if (param == "nchar" || param == "seqlen" || param == "length" || param == "N")
            seqlen = stoi(after);
        else if (param == "p" || param == "p_mutate")
            p_mutate = stod(after);
        else if (param == "smm" || param == "mutation_model")
            smm = stomm(after);
        else
            cerr << param << " is not a valid parameter" << endl;
    }

    Tree tree1(seed, species, seqlen, p_mutate, smm);
    Tree tree2(seed, species, seqlen, p_mutate, smm);
    Tree tree3(seed, species, seqlen, p_mutate, smm);
    tree1.generateRandomTopology();
    tree1.randomizeBranchLengths(bl_lower_bound, bl_upper_bound);
    // tree1.printEdges();
    tree2.copyTopology(tree1);
    tree2.randomizeBranchLengths(bl_lower_bound, bl_upper_bound);
    tree1.dfsSequenceGen();
    tree2.dfsSequenceGen();
    tree3.mixSequences({tree1, tree2});
    // cout << tree1.toNewick() << endl;
    // cout << tree2.toNewick() << endl;
    tree3.writeToNexus();
    cout << tree3.toNewick() << endl;
    cout << "exported results to " << to_string(seed) << endl;
    return 0;
}