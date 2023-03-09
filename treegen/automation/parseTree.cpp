#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <exception>
using namespace std;
string highProbTree(string filename)
{
    ifstream fin;
    fin.open(filename);
    string line;
    while (getline(fin, line))
    {
        if (line.find("tree tree_1") != string::npos)
        {
            int treeStart = line.find("(");
            return line.substr(treeStart);
        }
    }
    return "";
}

bool isNumber(string word)
{
    for (char a : word)
    {
        if (!isdigit(a))
        {
            return false;
        }
    }
    return true;
}

int iterationsToConvergence(string filename)
{
    ifstream fin;
    fin.open(filename);
    string line;
    int iters = -1;
    while (getline(fin, line))
    {
        if (line.size() == 0)
        {
            continue;
        }
        int end = line.find("\t");
        if (end == string::npos)
        {
            continue;
        }
        if (!isNumber(line.substr(0, end)))
        {
            continue;
        }
        iters = stoi(line.substr(0, end));
    }
    return iters;
}

int parseSeed(string filename)
{
    ifstream fin;
    int seed = -1;
    fin.open(filename);
    string line;
    while (getline(fin, line))
    {
        if (line.find("Seed = ") != string::npos)
        {

            seed = stoi(line.substr(line.find("Seed = ") + 7, string::npos));
            break;
        }
    }
    fin.close();
    if (seed == -1)
    {
        throw invalid_argument(filename + " does not contain 'Seed = '");
    }
    return seed;
}

int parseSwapSeed(string filename)
{
    ifstream fin;
    int swapSeed = -1;
    fin.open(filename);
    string line;
    while (getline(fin, line))
    {
        if (line.find("Swapseed = ") != string::npos)
        {
            swapSeed = stoi(line.substr(line.find("Swapseed = ") + 11, string::npos));
            break;
        }
    }
    fin.close();
    if (swapSeed == -1)
    {
        throw invalid_argument(filename + " does not contain 'Swapseed = '");
    }
    return swapSeed;
}
int main()
{
    ofstream fout;
    for (int i = 1; i <= 5; i++)
    {
        fout.open("output" + to_string(i) + "tr.txt");
        fout << highProbTree("output" + to_string(i) + ".nex.trprobs");
        fout.close();
    }
    fout.open("iter_summary.txt");
    for (int i = 1; i <= 5; i++)
    {
        // cout << "tree " << i << ": " << iterationsToConvergence("output" + to_string(i) + ".nex.mcmc");
        fout << iterationsToConvergence("output" + to_string(i) + ".nex.mcmc");
        fout << endl;
    }
    fout.close();
    fout.open("seed_summary.txt");
    for (int i = 1; i <= 5; i++)
    {
        fout << parseSeed("log" + to_string(i) + ".txt") << " " << parseSwapSeed("log" + to_string(i) + ".txt");
        fout << endl;
    }
    fout.close();
}