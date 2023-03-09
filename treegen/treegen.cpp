#include "tree.h"
#include "NewickTree.h"
#include <iostream>
#include <random>  // for random_device and uniform_int_distribution
#include <cstdint> // for uint64_t
#include <ctime>   // for time
using namespace std;

int main(int argc, char **argv)
{
    // initializations + cmd line reading
    std::random_device rd;
    std::uniform_int_distribution<uint64_t> dist(0, 99999999999999);
    uint64_t seed = dist(rd);  // default seed = time % 1M
    int species = 10;          // default # species = 10
    int seqlen = 1000;         // default sequence length = 1000
    double p_mutate = 0.2;     // default mutation probability = 0.2
    mutation_model smm = jc69; // default site mutation model = jc69
    double b_length = 0.1;
    string filepath = "";
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
        else if (param == "bl" || param == "branch_length")
        {
            b_length = stod(after);
        }
        else if (param == "f" || param == "filepath")
        {
            filepath = after;
        }
        else
            cerr << param << " is not a valid parameter" << endl;
    }

    NewickTree tree1(species, seed);
    tree1.setBranchLengths(b_length);
    tree1.printNewick(true);
    tree1.generateSequences(seqlen);
    tree1.exportNexus("results/" + filepath + to_string(seed) + ".nex");
    tree1.exportNewick("results/" + filepath + to_string(seed) + "_tree.txt", true);
    cout << "exported results to " << to_string(seed) << endl;
}