#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main(int argc, char **argv) {
    if (argc != 3) {
        cout << "use: ./resultreader [seed] [ntax]" << endl;
        return 0;
    }
    string seed = argv[1];
    int ntax = atoi(argv[2]);
    ofstream outfile;
    outfile.open("sumt.txt", ofstream::app);

    ifstream infile;
    infile.open("results/" + seed + "_data.nex.mcmc");
    string line;
    for (int i = 0; i < 6; i++) getline(infile, line); // text lines
    outfile << "  1st asdsf: " << flush;
    while (getline(infile, line)) {
        int gen = stoi(line.substr(0, line.find('\t')));
        if (stod(line.substr(line.find_last_of('\t')+1)) < 0.01) {
            outfile << gen << flush;
            break;
        }
    }
    outfile << endl;
    infile.close();

    infile.open("results/" + seed + "_data.nex.trprobs");
    for (int i = 0; i < 9; i++) getline(infile, line); // text lines
    for (int i = 0; i < ntax; i++) getline(infile, line); // leaves lines
    // outfile << "  trees:" << endl;
    while (getline(infile,line) && line[0] != 'e') {
        outfile << " " << line << endl;
    }
    // outfile << "max p: " << line.substr(line.find('=') + 2, 5) << endl; // probability
    // int spaceindex = line.find_last_of(' ')+1;
    // outfile << "tree: " << line.substr(spaceindex, line.length()-spaceindex-1) << endl; // tree
    infile.close();

    // outfile << endl;
    outfile.close();
    return 0;
}
