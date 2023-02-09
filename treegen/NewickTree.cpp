#include "NewickTree.h"
#include <queue>
#include <utility>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stack>
#include <map>
#include <set>
#include <random>
#include <stdexcept>
using namespace std;

NewickTree::NewickTree()
{
    root = nullptr;
    sequenceLength = 0;
    transitionMatrix = map<char, map<char, double>>();
}
NewickTree::NewickTree(string rootName) : NewickTree()
{
    root = new NewickTree::TreeNode(rootName, nullptr);
}
NewickTree::NewickTree(vector<Vertex> vertices, vector<Edge> edges, int rootIndex) : NewickTree()
{
    vector<string> names(vertices.size());
    for (int i = 0; i < names.size(); i++)
    {
        names[i] = vertices[i].name;
    }
    vector<vector<int>> adj(vertices.size());
    for (Edge e : edges)
    {
        adj[e.v1].push_back(e.v2);
        adj[e.v2].push_back(e.v1);
    }
    // should check that adj and names are same size and rootIndex < size later
    adjToTree(adj, names, rootIndex);
}
NewickTree::NewickTree(vector<vector<int>> adj, vector<string> names, int rootIndex) : NewickTree()
{
    // should check that adj and names are same size and rootIndex < size later
    sequenceLength = 0;
    transitionMatrix = map<char, map<char, double>>();
    adjToTree(adj, names, rootIndex);
}
NewickTree::NewickTree(int numLeafs) : NewickTree()
{
    vector<vector<int>> adjList(2 * numLeafs - 1); // total # of nodes will be 2n-1
    vector<int> leaves;
    leaves.push_back(0);
    int count = 0;
    while (leaves.size() < numLeafs)
    {
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

    vector<string> names(adjList.size(), "-");
    for (int i = 0; i < leaves.size(); i++)
    {
        names[leaves[i]] = to_string(i + 1);
    }
    adjToTree(adjList, names, 0);
}
void NewickTree::adjToTree(vector<vector<int>> adj, vector<string> names, int rootIndex)
{
    vector<bool> visited(adj.size());
    root = new NewickTree::TreeNode(names[rootIndex], nullptr);
    queue<pair<NewickTree::TreeNode *, int>> next;
    visited[rootIndex] = true;
    next.push(make_pair(root, rootIndex));
    while (next.size() > 0)
    {
        NewickTree::TreeNode *currNode = next.front().first;
        int currIndex = next.front().second;
        next.pop();
        for (int neighborIndex : adj[currIndex])
        {
            if (!visited[neighborIndex])
            {
                // add as child of curr, mark as visited, add to queue
                NewickTree::TreeNode *neighborNode = currNode->addChild(names[neighborIndex]);
                visited[neighborIndex] = true;
                next.push(make_pair(neighborNode, neighborIndex));
            }
        }
    }
}
void NewickTree::deleteDFS(NewickTree::TreeNode *start)
{
    if (start == nullptr)
    {
        return;
    }
    for (NewickTree::TreeNode *child : start->children)
    {
        deleteDFS(child);
    }
    delete start;
}
NewickTree::~NewickTree()
{
    deleteDFS(root);
}
// add leafOnly option
void NewickTree::newickDFS(NewickTree::TreeNode *start, vector<string> &symbols, bool leafOnly)
{
    if (!start->isLeaf())
    {
        symbols.push_back("(");
        for (int i = 0; i < start->childrenCount(); i++)
        {
            NewickTree::TreeNode *child = start->children[i];
            if (i > 0)
            {
                symbols.push_back(",");
            }
            newickDFS(child, symbols, leafOnly);
        }
        symbols.push_back(")");
    }
    if (!leafOnly || start->isLeaf())
        symbols.push_back(start->name);
}
int NewickTree::getSequenceLength()
{
    return sequenceLength;
}

void NewickTree::printNewick(bool leafOnly)
{
    vector<string> symbols;
    newickDFS(root, symbols, leafOnly);
    for (string symbol : symbols)
    {
        cout << symbol;
    }
    cout << ";";
    cout << endl;
}

void NewickTree::exportNewick(string filename, bool leafOnly)
{
    ofstream treefile;
    treefile.open(filename);
    vector<string> symbols;
    newickDFS(root, symbols, leafOnly);
    for (string symbol : symbols)
    {
        treefile << symbol;
    }
    treefile << ";";
    treefile.close();
}
NewickTree::TreeNode *NewickTree::importDFS(TreeNode *parent, string &newickString, int &pos)
{
    // pos should be position before current node
    // cout << "running at " << pos << endl;
    int nextEnd = newickString.find_first_of("(,)", pos + 1);
    // cout << "found nextEnd at " << nextEnd << " which is " << newickString.at(nextEnd) << endl;
    string currName;
    if (nextEnd == pos + 1)
    {
        currName = "-";
    }
    else
    {
        currName = newickString.substr(pos + 1, nextEnd - pos - 1);
        reverse(currName.begin(), currName.end());
    }
    TreeNode *currNode = new TreeNode(currName, parent);
    // if (parent != nullptr)
    //     cout << "creating " << currName << " under " << parent->name << endl;
    pos = nextEnd;
    if (newickString.at(nextEnd) == ')' || parent == nullptr && pos < newickString.size())
    {
        while (newickString.at(pos) != '(')
        {
            currNode->addChild(importDFS(currNode, newickString, pos));
        }
        pos += 1;
    }
    return currNode;
}
void NewickTree::importNewick(string filename)
{
    deleteDFS(root);
    ifstream treefile;
    treefile.open(filename);
    string newickString;
    getline(treefile, newickString);
    if (*newickString.rbegin() == ';')
    {
        newickString.pop_back();
    }
    reverse(newickString.begin(), newickString.end());
    int pos = -1;
    cout << newickString << endl;
    root = importDFS(nullptr, newickString, pos);
    return;
}
int NewickTree::getLeafCount()
{
    queue<TreeNode *> next;
    int leafCount = 0;
    if (root == nullptr)
    {
        return 0;
    }
    next.push(root);
    while (next.size() > 0)
    {
        TreeNode *currNode = next.front();
        next.pop();
        if (currNode->isLeaf())
        {
            leafCount++;
        }
        for (TreeNode *child : currNode->children)
        {
            next.push(child);
        }
    }
    return leafCount;
}
void NewickTree::setTransition(map<char, map<char, double>> transitionMatrix)
{

    if (transitionMatrix.size() == 0)
    {
        double p_mutate = 0.3;
        transitionMatrix = map<char, map<char, double>>{
            {'A', {{'A', 1.0 - p_mutate}, {'T', p_mutate / 3.0}, {'G', p_mutate / 3.0}, {'C', p_mutate / 3.0}}},
            {'T', {{'A', p_mutate / 3.0}, {'T', 1.0 - p_mutate}, {'G', p_mutate / 3.0}, {'C', p_mutate / 3.0}}},
            {'G', {{'A', p_mutate / 3.0}, {'T', p_mutate / 3.0}, {'G', 1.0 - p_mutate}, {'C', p_mutate / 3.0}}},
            {'C', {{'A', p_mutate / 3.0}, {'T', p_mutate / 3.0}, {'G', p_mutate / 3.0}, {'C', 1.0 - p_mutate}}}};
    }
    this->transitionMatrix = transitionMatrix;
}
void NewickTree::generateSequences(int sequenceLength)
{
    if (transitionMatrix.size() != symbols.size())
    {
        setTransition();
        cerr << "dimensions of transitionMatrix are incorrect, setting transitionMatrix to default (mutation probability of 0.3 equally split)" << endl;
    }
    for (pair<char, map<char, double>> mapping : transitionMatrix)
    {
        if (mapping.second.size() != symbols.size())
        {
            setTransition();
            cerr << "dimensions of transitionMatrix are incorrect, setting transitionMatrix to default (mutation probability of 0.3 equally split)" << endl;
            break;
        }
    }
    this->sequenceLength = sequenceLength;
    random_device rd;
    mt19937 mt(rd());
    default_random_engine gen;
    uniform_real_distribution<double> distribution(0.0,
                                                   1.0);
    string rootSeq = "";
    for (int i = 0; i < sequenceLength; i++)
    {
        char newChar = symbols[mt() % symbols.size()];
        rootSeq.push_back(newChar);
    }
    root->sequence = rootSeq;
    queue<TreeNode *> next;
    next.push(root);
    while (next.size() > 0)
    {
        TreeNode *currNode = next.front();
        next.pop();
        for (TreeNode *child : currNode->children)
        {
            string childSeq = "";
            for (int i = 0; i < currNode->sequence.size(); i++)
            {
                char nucleotide = currNode->sequence.at(i);
                double randChoice = distribution(gen);
                map<char, double> probabilities = transitionMatrix[nucleotide];
                double cummSum = 0.0;
                for (pair<char, double> candidate : probabilities)
                {
                    cummSum += candidate.second;
                    if (cummSum >= randChoice)
                    {
                        // cout << nucleotide << " " << candidate.first << " " << cummSum << " " << candidate.second << endl;
                        childSeq.push_back(candidate.first);
                        break;
                    }
                }
                // deals with floating point errors that cause transition probabilities to not sum up to 1 (should also manually check that they sum to 1 in each row)
                if (childSeq.size() < i + 1)
                {
                    childSeq.push_back(probabilities.rbegin()->first);
                }
            }
            child->sequence = childSeq;
            next.push(child);
        }
        // if (currNode->isLeaf())
        // {
        //     cout << currNode->name << " " << currNode->sequence << endl;
        // }
    }
}
void NewickTree::toNexus(string filename, map<string, string> sequences, vector<NewickTree *> trees, bool numericSort)
{
    ofstream file;
    file.open(filename);
    file << "begin data;" << endl;
    file << "dimensions ntax=" << sequences.size() << " nchar=" << sequences.begin()->second.size() << ";" << endl;
    file << "format datatype=dna interleave=no gap=-;" << endl;
    file << "matrix" << endl;
    // allows you to process
    vector<pair<string, string>> ordered_sequences;
    for (pair<string, string> sequence : sequences)
    {
        ordered_sequences.push_back(sequence);
    }
    if (numericSort)
    {
        sort(ordered_sequences.begin(), ordered_sequences.end(), ([](const pair<string, string> &a, const pair<string, string> &b)
                                                                  { return stoi(a.first) < stoi(b.first); }));
    }
    for (pair<string, string> sequence : ordered_sequences)
    {
        file << sequence.first << "\t" << sequence.second << endl;
    }

    file << ";" << endl;
    file << "end;" << endl;
    if (trees.size() > 0)
    {
        file << "begin trees;" << endl;
        for (int i = 0; i < trees.size(); i++)
        {
            NewickTree *treePointer = trees[i];
            file << "tree tree" << i + 1 << " = ";
            vector<string> symbols;
            treePointer->newickDFS(treePointer->root, symbols, true);
            for (string symbol : symbols)
            {
                file << symbol;
            }
            file << ";" << endl;
        }
        file << "end;" << endl;
    }
    file.close();
}
map<string, string> NewickTree::getSequences()
{
    queue<TreeNode *> next;
    map<string, string> sequences;
    if (root == nullptr)
    {
        return sequences;
    }
    next.push(root);
    while (next.size() > 0)
    {
        TreeNode *currNode = next.front();
        next.pop();
        if (currNode->isLeaf())
        {
            sequences[currNode->name] = currNode->sequence;
        }
        for (TreeNode *child : currNode->children)
        {
            next.push(child);
        }
    }
    return sequences;
}
void NewickTree::exportNexus(string filename, bool numericSort)
{
    toNexus(filename, getSequences(), vector<NewickTree *>{this}, numericSort);
}
map<string, string> NewickTree::mixtureModel(vector<NewickTree *> trees, vector<double> weights = vector<double>())
{
    // validate that trees is not empty
    if (trees.size() == 0)
    {
        throw invalid_argument("NewickTree::mixtureModel 'trees' must be non empty vector of type NewickTree");
    }
    // validate that weights is the same length as trees or is 0
    if (weights.size() != 0 && weights.size() != trees.size())
    {
        throw invalid_argument("NewickTree::mixtureModel 'weights' vector must be either empty or of the same dimensions as 'trees'");
    }
    // validate that all trees have the same leaves
    int leafCount = trees[0]->getLeafCount();
    for (NewickTree *tree : trees)
    {
        if (tree->getLeafCount() != leafCount)
        {
            throw invalid_argument("NewickTree::mixtureModel all NewickTree's in 'trees' vector must have the same number of leaf nodes");
        }
    }
    vector<map<string, string>> sequenceList(trees.size());
    for (int i = 0; i < trees.size(); i++)
    {
        sequenceList[i] = trees[i]->getSequences();
    }
    set<string> keySet;
    for (pair<string, string> keyPair : sequenceList[0])
    {
        keySet.insert(keyPair.first);
    }
    for (int i = 0; i < trees.size(); i++)
    {
        for (pair<string, string> keyPair : sequenceList[i])
        {
            if (keySet.find(keyPair.first) == keySet.end())
            {
                throw invalid_argument("NewickTree::mixtureModel all NewickTree's in 'trees' vector must have identically labeled leaf nodes");
            }
        }
    }

    // validate that all trees have equal sequence lengths
    int sequenceLength = trees[0]->getSequenceLength();
    for (NewickTree *tree : trees)
    {
        if (tree->getSequenceLength() != sequenceLength)
        {
            throw invalid_argument("NewickTree::mixtureModel all NewickTree's in 'trees' vector must have the same sequence length");
        }
    }

    if (weights.size() == 0)
    {
        weights = vector<double>(trees.size());
        for (int i = 0; i < weights.size(); i++)
        {
            weights[i] = 1.0 / weights.size();
        }
    }
    double epsilon = 1e-7;
    double weightSum = 0.0;
    for (double weight : weights)
    {
        weightSum += weight;
    }
    if (abs(1.0 - weightSum) > epsilon)
    {
        throw invalid_argument("NewickTree::mixtureModel 'weights' must sum up to 1.0");
    }
    default_random_engine gen;
    uniform_real_distribution<double> distribution(0.0, 1.0);

    map<string, string> mixedSequence;
    for (string key : keySet)
    {
        mixedSequence[key] = "";
    }
    for (int i = 0; i < sequenceLength; i++)
    {
        // pick random tree for this index
        double randChoice = distribution(gen);
        int randIndex = -1;
        double cummSum = 0.0;
        for (int j = 0; j < weights.size(); j++)
        {
            cummSum += weights[j];
            if (cummSum >= randChoice)
            {
                randIndex = j;
                break;
            }
        }
        if (randIndex == -1)
        {
            cerr << "Did not properly select random index for position " << i << " in mixture model. Defaulting to " << weights.size() - 1 << "." << endl;
            randIndex = weights.size() - 1;
        }
        for (string key : keySet)
        {
            mixedSequence[key].push_back(sequenceList[randIndex][key].at(i));
        }
    }
    return mixedSequence;
}