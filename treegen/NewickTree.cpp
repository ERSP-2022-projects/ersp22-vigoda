#include "NewickTree.h"
#include <queue>
#include <utility>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stack>
#include <map>
using namespace std;

NewickTree::NewickTree()
{
    root = nullptr;
}
NewickTree::NewickTree(string rootName)
{
    root = new NewickTree::TreeNode(rootName, nullptr);
}
NewickTree::NewickTree(vector<Vertex> vertices, vector<Edge> edges, int rootIndex)
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
                NewickTree::TreeNode *neighborNode = new NewickTree::TreeNode(names[neighborIndex], currNode);
                currNode->children.push_back(neighborNode);
                visited[neighborIndex] = true;
                next.push(make_pair(neighborNode, neighborIndex));
            }
        }
    }
}
NewickTree::NewickTree(vector<vector<int>> adj, vector<string> names, int rootIndex)
{
    // should check that adj and names are same size and rootIndex < size later
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
                NewickTree::TreeNode *neighborNode = new NewickTree::TreeNode(names[neighborIndex], currNode);
                currNode->children.push_back(neighborNode);
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
void NewickTree::printNewick(bool leafOnly)
{
    vector<string> symbols;
    newickDFS(root, symbols, leafOnly);
    for (string symbol : symbols)
    {
        cout << symbol;
    }
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
    treefile.close();
}
NewickTree::TreeNode *NewickTree::importDFS(TreeNode *parent, string &newickString, int &pos)
{
    // pos should be position before current node
    cout << "running at " << pos << endl;
    int nextEnd = newickString.find_first_of("(,)", pos + 1);
    cout << "found nextEnd at " << nextEnd << " which is " << newickString.at(nextEnd) << endl;
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
    if (parent != nullptr)
        cout << "creating " << currName << " under " << parent->name << endl;
    pos = nextEnd;
    if (newickString.at(nextEnd) == ')')
    {
        while (newickString.at(pos) != '(')
        {
            currNode->children.push_back(importDFS(currNode, newickString, pos));
        }
        pos += 1;
    }
    cout << "returning at pos " << pos << endl;
    return currNode;
}
void NewickTree::importNewick(string filename)
{
    cout << "got hereee!" << endl;
    deleteDFS(root);
    cout << "got here!!" << endl;
    ifstream treefile;
    treefile.open(filename);
    string newickString;
    getline(treefile, newickString);
    reverse(newickString.begin(), newickString.end());
    int pos = -1;
    cout << "got here!" << endl;
    cout << newickString << endl;
    root = importDFS(nullptr, newickString, pos);
    return;
}