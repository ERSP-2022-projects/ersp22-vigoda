#include "tree.h"
using namespace std;

int main(int argc, char **argv) {
    // initializations + cmd line reading
    uint64_t seed = (time(0) % 100000); // default seed = time % 1M
    int species = 10;                   // default # species = 10
    int seqlen = 1000;                  // default sequence length = 1000
    double p_mutate = 0.2;              // default mutation probability = 0.2
    mutation_model smm = jc69;          // default site mutation model = jc69
    bool missing_data = false;          // default no missing data
    double missing_pct = 0;              // ^
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
        else if (param == "ntax" || param == "species" || param == "numspecies" || param == "leaves" || param == "n")
            species = stoi(after);
        else if (param == "nchar" || param == "seqlen" || param == "length" || param == "N")
            seqlen = stoi(after);
        else if (param == "p" || param == "p_mutate")
            p_mutate = stod(after);
        else if (param == "smm" || param == "mutation_model")
            smm = stomm(after);
        else if (param == "drop" || param == "missing" || param == "pctmissing") {
            missing_data = true;
            missing_pct = stod(after);
        }
        else
            cerr << param << " is not a valid parameter" << endl;
    }

    Tree tree(seed, species, seqlen);
    tree.generateTopology();
    tree.dfsSequenceGen(p_mutate, smm);
    if (missing_data) tree.dropData(missing_pct);
    tree.writeToNexus();

    ofstream outfile;
    outfile.open("sumt.txt", ofstream::app);
    outfile << endl;
    outfile << "seed: " << seed << endl;
    outfile << "pct missing: " << missing_pct << endl;
    outfile << "original tree: " << tree.toNewick() << endl;
    outfile << "trials: " << endl;
    outfile.close();

    cout << "execute /Users/katytsao/Documents/GitHub/ersp22-vigoda/treegen/results/" << seed << "_data.nex" << endl;
    cout << "mcmc ngen=5000 printfreq=1000 diagnfreq=100 samplefreq=20" << endl;
    cout << "./resultreader.out " << seed << " 50" << endl;

    return 0;
}
