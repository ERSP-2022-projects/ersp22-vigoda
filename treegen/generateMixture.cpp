#include "NewickTree.h"
#include <iostream>
#include <filesystem>
#include <random>
using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " <tree1 file> <tree2 file>" << endl;
        return 1;
    }

    random_device rd;
    srand(rd());
    int orig = rand() % static_cast<int>(1e6);
    int numLeafs = 10;
    int sequenceLength = 1000;
    // NewickTree tree1 = NewickTree(numLeafs);
    // NewickTree tree2 = NewickTree(numLeafs);
    NewickTree tree1;
    NewickTree tree2;
    tree1.importNewick(argv[1]);
    tree2.importNewick(argv[2]);
    tree1.generateSequences(sequenceLength);
    tree2.generateSequences(sequenceLength);
    filesystem::create_directory("results/" + to_string(orig));
    filesystem::create_directory("results/" + to_string(orig) + "/mixture");
    cout << "storing results in results/" << to_string(orig) << endl;
    tree1.exportNewick("results/" + to_string(orig) + "_tree1.txt");
    tree2.exportNewick("results/" + to_string(orig) + "_tree2.txt");
    tree1.exportNexus("results/" + to_string(orig) + "_tree1.nex", true);
    tree2.exportNexus("results/" + to_string(orig) + "_tree2.nex", true);
    map<string, string> mixture = NewickTree::mixtureModel(vector<NewickTree *>{&tree1, &tree2}, vector<double>{0.5, 0.5});
    NewickTree::toNexus(
        "results/" + to_string(orig) + "/mixture/mixture.nex", mixture, vector<NewickTree *>{&tree1, &tree2}, true);
}
