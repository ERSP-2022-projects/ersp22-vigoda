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
    vector<char> chars= {'A','G','T','C'};
    int index = rand() % 4;
    return chars[index];
}

// Set char based on Jukes Cantor model. Input substitution probability
char getRandomCharJC(double subProb, char c) {
    vector<char> chars = {'A','G','T','C'};
    for (int i=0; i<chars.size(); ++i) {
        if (chars[i] == c) {
            chars.erase(chars.begin() + i);
        }
    }
    double num = (double)rand() / RAND_MAX;
    if (num < subProb) { return chars[0]; }
    else if (num < 2 * subProb) { return chars[1]; }
    else if (num < 3 * subProb) { return chars[2]; }
    else { return c; }
}

// TODO: Add K2P substituion model

// Tree topology generator
vector<vector<int>> createTopology(int n) {
    vector<vector<int>> adjList(2 * n - 1); // total # of nodes will be 2n-1
    vector<int> leaves;
    leaves.push_back(0);
    int count = 0;
    while (leaves.size() < n) {
        // randomly pick a leaf to split
        int index = rand() % (int)leaves.size(); 
        int parent = leaves[index];
        leaves.erase(leaves.begin() + index); // O(n) x_x
        // add left child
        leaves.push_back(++count);
        adjList[parent].push_back(count);
        adjList[count].push_back(parent);
        // add right child
        leaves.push_back(++count);
        adjList[parent].push_back(count);
        adjList[count].push_back(parent);
    }
    return adjList;
}

// Input: Tree topology via an adjacency list
// Output: Maps each leaf node to a taxa label
unordered_map<int,int> getLeafMap(vector<vector<int>>& adjList) {
    unordered_map<int,int> res;
    int label = 0;
    for (int i=0; i<adjList.size(); ++i) {
        if (adjList[i].size() == 1) {
            res.insert({i, label++});
        }
    }
    return res;
}

// Input: Vector of k tree topologies 
// Output: Vector of leaf node dna sequences
vector<string> generateData(vector<vector<vector<int>>>& treeTops, vector<unordered_map<int,int>>& leafMaps) {
    vector<string> res(TAXA, "");
    for (int i=0; i<LENGTH; ++i) {
        int treeIndex = rand() % (int)treeTops.size(); // randomly picks a tree topology
        queue<int> queue;
        vector<bool> visited(treeTops[0].size(), false);
        queue.push(0);
        visited[0] = true;
        char c = getRandomChar();
        while (!queue.empty()) {
            int curr = queue.front();
            queue.pop();
            for (int j=0; j<treeTops[treeIndex][curr].size(); ++j) {
                int neigh = treeTops[treeIndex][curr][j];
                if (!visited[neigh]) {
                    visited[neigh] = true;
                    c = getRandomCharJC(0.1, c);
                    queue.push(neigh);
                    if (treeTops[treeIndex][neigh].size() == 1) { // if leaf node, append char to DNA sequence
                        int label = leafMaps[treeIndex][neigh];
                        res[label].push_back(c);
                    }
                }
            }
        }
    }
    return res;
}

// Creates and writes to .nex file
void createFile(vector<string>& leafData) {
    fstream myFile;
    myFile.open("tree.nex", ios::out);
    if (myFile.is_open()) {
        myFile << "#NEXUS" << endl << endl;
        myFile << "begin data;" << endl;
        myFile << "\tdimensions ntax=" << TAXA << " nchar=" << LENGTH << ";" << endl;
        myFile << "\tformat datatype=dna interleave=no gap=-;" << endl;
        myFile << "\tmatrix" << endl;
        for (int i=0; i<leafData.size(); ++i) {
            myFile << "Taxa" << i+1 << "\t" << leafData[i] << endl;
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

// DFS for Newick format
void newickDfs(vector<vector<int>>& adjList, unordered_map<int,int>& leafMap, string& res, int cur, int par) {
    if (leafMap.find(cur) == leafMap.end()) {
        vector<int> children;
        for (auto i : adjList[cur]) {
            if (i != par) {
                children.push_back(i);
            }
        }
        res += "(";
        newickDfs(adjList, leafMap, res, children[0], cur);
        res += ",";
        newickDfs(adjList, leafMap, res, children[1], cur);
        res += ")";
    }
    else {
        res += to_string(leafMap[cur]+1);
    }   
}

// Creates Newick format
string createNewick(vector<vector<int>>& adjList, unordered_map<int,int>& leafMap) {
    string res = "";
    newickDfs(adjList, leafMap, res, 0, -1);
    return res;
}

// For testing, prints adjacency list
void printTopology(vector<vector<int>>& adjList) {
    cout << "TOPOLOGY:" << endl;
    for (int i=0; i<adjList.size(); ++i) {
        cout << i << ": ";
        for (int j=0; j<adjList[i].size(); ++j) {
            cout << adjList[i][j] << " ";
        }
        cout << endl;
    }
}

// For testing, prints leaf map
void printLeafMap(unordered_map<int,int>& leafMap) {
    for (auto i : leafMap) {
        cout << "(" << i.first << ", " << i.second << ")  ";
    }
    cout << endl;
}

// Main function, generates leaf data using a k-mixture of tree topologies
void generate(int k) {
    vector<vector<vector<int>>> treeTops(k); // stores tree topologies
    vector<unordered_map<int,int>> leafMaps(k); // stores leaf map for each tree topology
    vector<string> newickTrees(k); // stores newick format for each tree topology
    for (int i=0; i<k; ++i) {
        treeTops[i] = createTopology(TAXA);
        leafMaps[i] = getLeafMap(treeTops[i]);
        newickTrees[i] = createNewick(treeTops[i], leafMaps[i]);
    }
    vector<string> leafData = generateData(treeTops, leafMaps);
    createFile(leafData);
    for (int i=0; i<k; ++i) {
        cout << "Newick Format " << i+1 << ": " << newickTrees[i] << endl; 
    }
    // TESTING:
    // for (int i=0; i<k; ++i) {
    //     cout << endl;
    //     printTopology(treeTops[i]);
    //     cout << endl;
    //     printLeafMap(leafMaps[i]);
    // }
}

int main() {
    srand(time(NULL));
    generate(2); // MR BAYES FAILED :P
    return 0;
}