#include <iostream>
#include <fstream>
#include <filesystem>
#include <random>
#include "NewickTree.h"
using namespace std;

int main(int argc, char *argv[])
{

    int numLeafs = 10;
    int sequenceLength = 4000;

    if (argc > 1)
    {
        numLeafs = stoi(argv[1]);
    }
    if (argc > 2)
    {
        sequenceLength = stoi(argv[2]);
    }
    string parentFolder = "results/" + to_string(numLeafs) + "leaves";
    if (!filesystem::exists(parentFolder))
    {
        filesystem::create_directory(parentFolder);
        filesystem::copy("automation/alltopos.py", parentFolder + "/");
    }

    random_device rd;
    srand(rd());
    int orig = rand() % static_cast<int>(1e6);

    string foldername = to_string(orig);

    while (filesystem::exists(parentFolder + "/" + foldername))
    {
        orig = rand() % static_cast<int>(1e6);
        foldername = to_string(orig);
    }

    NewickTree tree1(numLeafs);
    NewickTree tree2(tree1);

        filesystem::create_directory(parentFolder + "/" + foldername);
    cout << "storing results in " + parentFolder + "/" << foldername << endl;
}