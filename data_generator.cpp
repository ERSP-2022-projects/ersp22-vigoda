#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <time.h>
using namespace std;

int TAXA = 10; // set number of taxa
int LENGTH = 1000; // set length of dna sequence

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

// Input: Number of taxa
// Output: Adjacency list representing tree topology
vector<vector<int>> createTopology(int n) {
    vector<vector<int>> adjList(2 * n - 1); // based on my tree setup, total # of nodes is 2n-1. 
    int count = 0;
    queue<int> queue;
    queue.push(count++);
    while (queue.size() < n) {
        int curr = queue.front();
        queue.pop();
        // insert left neighbor
        adjList.at(curr).push_back(count);
        adjList.at(count).push_back(curr);
        queue.push(count++);
        // insert right neighbor
        adjList.at(curr).push_back(count);
        adjList.at(count).push_back(curr);
        queue.push(count++);
    }
    return adjList;
}

// Input: Tree topology, starting node, sequence length
// Output: DNA sequences for all leaf nodes
vector<string> generateData(vector<vector<int>>& adjList, int start, int length) {
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
                    c = getRandomCharJC(0.01, res.at(curr).back()); // modify substituion rates here, rate = [0, 0.25]
                    res.at(neighbor).push_back(c);
                    queue.push(neighbor);
                }
            }
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

int main() {
    srand(time(NULL));
    vector<vector<int>> adjList = createTopology(TAXA);
    int start = rand() % (adjList.size() - TAXA);
    vector<string> res = generateData(adjList, start, LENGTH);
    for (int i=res.size() - TAXA; i<res.size(); ++i) {
        cout << "Species " << i - TAXA + 1 << ": " << res.at(i) << endl;
        cout << endl;
    }
    getDifferenceRate(res, start);
    // DEBUG lmao sobbing and shitting rn
    // for (int i=0; i<adjList.size(); ++i) {
    //     cout << i << ": ";
    //     for (int j=0; j<adjList.at(i).size(); ++j) {
    //         cout << adjList.at(i).at(j) << " ";
    //     }
    //     cout << endl;
    // }
    return 0;
}