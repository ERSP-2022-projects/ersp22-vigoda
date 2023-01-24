#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <unordered_map>
#include <time.h>
#include <fstream>
using namespace std;

const int TAXA = 10; // set number of taxa
const int LENGTH = 1000; // set length of dna sequence

// Only use for setting the initial char, choose U.A.R
char getRandomChar() {
    double num = (double)rand() / RAND_MAX;
    if (num < 0.25) { return 'A'; }
    else if (num < 0.5) { return 'G'; }
    else if (num < 0.75) { return 'T'; }
    else { return 'C'; }
}

// Set char based on Jukes Cantor model. Input substitution probability
char getRandomCharJC(double subProb, char c) {
    vector<char> chars = {'A','G','T','C'};
    for (int i=0; i<chars.size(); ++i) {
        if (chars.at(i) == c) {
            chars.erase(chars.begin() + i);
        }
    }
    double num = (double)rand() / RAND_MAX;
    if (num < subProb) { return chars.at(0); }
    else if (num < 2 * subProb) { return chars.at(1); }
    else if (num < 3 * subProb) { return chars.at(2); }
    else { return c; }
}

// Better tree topology generator
vector<vector<int>> createTopology(int n) {
    vector<vector<int>> adjList(2 * n - 1); // total # of nodes will be 2n-1
    vector<int> leaves;
    leaves.push_back(0);
    int count = 0;
    while (leaves.size() < n) {
        // randomly pick a leaf to split
        int index = rand() % (int)leaves.size(); 
        int parent = leaves.at(index);
        leaves.erase(leaves.begin() + index); // O(n) x_x
        // add left child
        leaves.push_back(++count);
        adjList.at(parent).push_back(count);
        adjList.at(count).push_back(parent);
        // add right child
        leaves.push_back(++count);
        adjList.at(parent).push_back(count);
        adjList.at(count).push_back(parent);
    }
    return adjList;
}

// Input: Tree topology, starting node, sequence length
// Output: DNA sequences for all leaf nodes
// Uses one distribution
vector<string> generateData(vector<vector<int>>& adjList, int length) {
    int start = 0;
    vector<string> res(adjList.size(), "");
    for (int i=0; i<length; ++i) {
        queue<int> queue;
        vector<bool> visited(adjList.size(), false);
        queue.push(start);
        char c = getRandomChar();
        res.at(start).push_back(c);
        while (!queue.empty()) {
            int curr = queue.front();
            queue.pop();
            visited.at(curr) = true;
            for (int j=0; j<adjList.at(curr).size(); ++j) {
                int neighbor = adjList.at(curr).at(j);
                if (!visited.at(neighbor)) {
                    c = getRandomCharJC(0.1, res.at(curr).back()); // modify substituion rates here, rate = [0, 0.25]
                    res.at(neighbor).push_back(c);
                    queue.push(neighbor);
                }
            }
        }
    }
    return res;
}

// Input: Tree topology, starting node, sequence length
// Output: DNA sequences for all leaf nodes
// Uses two distributions
vector<string> generateData2(vector<vector<int>>& adjList, int length) {
    vector<double> probs = {0.1, 0.12}; // modify substitution rates here
    int start = 0;
    vector<string> res(adjList.size(), "");
    for (int i=0; i<length; ++i) {
        queue<int> queue;
        vector<bool> visited(adjList.size(), false);
        queue.push(start);
        char c = getRandomChar();
        res.at(start).push_back(c);
        while (!queue.empty()) {
            int curr = queue.front();
            queue.pop();
            visited.at(curr) = true;
            for (int j=0; j<adjList.at(curr).size(); ++j) {
                int neighbor = adjList.at(curr).at(j);
                if (!visited.at(neighbor)) {
                    int index = rand() % 2;
                    c = getRandomCharJC(probs[index], res.at(curr).back());
                    res.at(neighbor).push_back(c);
                    queue.push(neighbor);
                }
            }
        }
    }
    return res;
}

// Extracts <leaf node, DNA> from tree topology and DNA sequences
unordered_map<int, string> getLeaves(vector<vector<int>>& adjList, vector<string>& data) {
    unordered_map<int, string> res;
    for (int i=0; i<adjList.size(); ++i) {
        if (adjList.at(i).size() == 1) {
            res.insert({i, data.at(i)});
        }
    }
    return res;
}

// For testing purposes, gets the average number of substitutions from starting node to leaves
void getDifferenceRate(vector<string>& res, int start) {
    string original = res.at(start);
    vector<int> diffs(TAXA, 0);
    vector<string> leaves;
    for (int i = 0; i<TAXA; ++i) {
        leaves.push_back(res.at(res.size() - 1 - i));
    }
    double sum = 0.0;
    for (int i=0; i<TAXA; ++i) {
        for (int j=0; j<LENGTH; ++j) {
            if (leaves.at(i).at(j) != original.at(j)) {
                ++diffs.at(i);
            }
        }
        sum += diffs.at(i);
    }
    double avgDiffs = sum / TAXA;
    cout << "Average Diffs: " << avgDiffs << endl;
}

// Creates and writes to file
void createFile(unordered_map<int, string>& leaves) {
    fstream myFile;
    myFile.open("test-tree.nex", ios::out);
    if (myFile.is_open()) {
        myFile << "#NEXUS" << endl << endl;
        myFile << "begin data;" << endl;
        myFile << "\tdimensions ntax=" << TAXA << " nchar=" << LENGTH << ";" << endl;
        myFile << "\tformat datatype=dna interleave=no gap=-;" << endl;
        myFile << "\tmatrix" << endl;
        for (auto it : leaves) {
            myFile << "Node" << it.first << "\t" << it.second << endl;
        }
        myFile << "\t;" << endl;
        myFile << "end;" << endl;
        myFile.close();
        return;
    }
    else {
        cout << "error creating file" << endl;
    }
}

// Prints adjacency list for testing
void printTopology(vector<vector<int>>& adjList) {
    for (int i=0; i<adjList.size(); ++i) {
        cout << i << ": ";
        for (int j=0; j<adjList.at(i).size(); ++j) {
            cout << adjList.at(i).at(j) << " ";
        }
        cout << endl;
    }
}

int main() {
    srand(time(NULL));
    vector<vector<int>> adjList = createTopology(TAXA);
    vector<string> res = generateData(adjList, LENGTH);
    unordered_map<int, string> leaves = getLeaves(adjList, res);
    createFile(leaves);
    cout << "TOPOLOGY:" << endl;
    printTopology(adjList);
    // cout << "ALL DNA:" << endl;
    // for (int i=0; i<res.size(); ++i) {
    //     cout << i << ": " << res.at(i) << endl;
    // }
    return 0;
}