#include "NewickTree.h"
#include <iostream>
#include <filesystem>
#include <random>
using namespace std;

int main(int argc, char *argv[])
{

    random_device rd;
    srand(rd());
    int orig = rand() % static_cast<int>(1e6);

    string foldername = to_string(orig);
    if (argc <= 1)
    {
        foldername = foldername;
    }
    else
    {
        foldername = argv[1];
    }
    int numLeafs = 10;
    int sequenceLength = 1000;
    NewickTree tree1 = NewickTree(numLeafs);
    NewickTree tree2 = NewickTree(numLeafs);
    // NewickTree tree1;
    // NewickTree tree2;
    // tree1.importNewick("tree1.txt");
    // tree2.importNewick("tree2.txt");
    tree1.generateSequences(sequenceLength);
    tree2.generateSequences(sequenceLength);
    filesystem::create_directory("results/" + foldername);
    filesystem::create_directory("results/" + foldername + "/mixture");
    cout << "storing results in results/" << foldername << endl;
    tree1.exportNewick("results/" + foldername + "/tree1.txt");
    tree2.exportNewick("results/" + foldername + "/tree2.txt");
    tree1.exportNexus("results/" + foldername + "/tree1.nex", true);
    tree2.exportNexus("results/" + foldername + "/tree2.nex", true);
    map<string, string> mixture = NewickTree::mixtureModel(vector<NewickTree *>{&tree1, &tree2}, vector<double>{0.5, 0.5});
    NewickTree::toNexus(
        "results/" + foldername + "/mixture/mixture.nex", mixture, vector<NewickTree *>{&tree1, &tree2}, true);
}
