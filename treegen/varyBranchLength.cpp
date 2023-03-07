#include <iostream>
#include <fstream>
#include <filesystem>
#include <random>
#include "NewickTree.h"
using namespace std;

struct config
{
    string name;
    double internalLength;
    double terminalLength;
    config() {}
    config(string name, double internalLength, double terminalLength)
    {
        this->name = name;
        this->internalLength = internalLength;
        this->terminalLength = terminalLength;
    }
};

int main(int argc, char *argv[])
{

    int numLeafs = 10;
    int sequenceLength = 4000;

    if (argc <= 1)
    {
        numLeafs = numLeafs;
    }
    else
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

    NewickTree generatingTree(numLeafs);
    filesystem::create_directory(parentFolder + "/" + foldername);
    cout << "storing results in " + parentFolder + "/" << foldername << endl;
    double sl = 0.1;
    double ll = 1.0;
    vector<config> configs = {
        config("allShort", sl, sl),
        config("allLong", ll, ll),
        config("internalShort", sl, ll),
        config("internalLong", ll, sl)};
    filesystem::copy("automation/allconfig.sh", parentFolder + "/" + foldername + "/");
    filesystem::copy("automation/aggregateTopo.py", parentFolder + "/" + foldername + "/");

    for (config setting : configs)
    {
        string folderpath = parentFolder + "/" + foldername + "/" + setting.name;
        filesystem::create_directory(folderpath);
        generatingTree.setInternalLengths(setting.internalLength);
        generatingTree.setTerminalLengths(setting.terminalLength);
        generatingTree.generateSequences(sequenceLength);
        generatingTree.exportNexus(folderpath + "/" + setting.name + ".nex", true);
        generatingTree.exportNewick(folderpath + "/generatingTree.txt");
        filesystem::copy("automation/parseTree.out", folderpath + "/");
        filesystem::copy("automation/treeCorrectness.R", folderpath + "/");
        ofstream fout;
        for (int i = 1; i <= 5; i++)
        {
            fout.open(folderpath + "/batch" + to_string(i) + ".txt");
            fout << "set autoclose=yes nowarn=yes" << endl;
            fout << "execute " + setting.name << ".nex" << endl;
            fout << "mcmc ngen=1000000 stoprule=yes stopval=0.01 samplefreq=100 diagnfreq=100 printfreq=100 file=output" << i << ".nex" << endl;
            fout << "sumt" << endl;
            fout << "quit" << endl;
            fout.close();
        }

        fout.open(folderpath + "/makefile");
        fout << "clean: " << endl;
        fout << "\trm *.nex.*" << endl
             << endl;
        fout << "strip:" << endl;
        fout << "\tmkdir runinfo\n";
        fout << "\tmv *.run1.* *.run2.* *.con.tre *.tstat *.vstat *.parts *.ckp *.ckp~ runinfo 2>/dev/null";
        fout.close();
        fout.open(folderpath + "/runbatch.sh");
        fout << "#!/ bin / bash\n"
             << "mb<batch1.txt> log1.txt\n"
             << "mb<batch2.txt> log2.txt\n"
             << "mb<batch3.txt> log3.txt\n"
             << "mb<batch4.txt> log4.txt\n"
             << "mb<batch5.txt> log5.txt\n"
             << "make strip\n"
             << "./parseTree.out\n"
             << "Rscript treeCorrectness.R";
        fout.close();
    }
}