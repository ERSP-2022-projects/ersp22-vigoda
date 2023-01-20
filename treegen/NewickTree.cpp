#include "NewickTree.h"
#include <queue>
#include <utility>
#include <iostream>
using namespace std;

NewickTree::NewickTree() {}
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
// add leafOnly option
void NewickTree::newickDFS(NewickTree::TreeNode *start, vector<string> &symbols)
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
            newickDFS(child, symbols);
        }
        symbols.push_back(")");
    }
    symbols.push_back(start->name);
}
void NewickTree::printNewick(bool leafOnly)
{
    vector<string> symbols;
    newickDFS(root, symbols);
    for (string symbol : symbols)
    {
        cout << symbol;
    }
    cout << endl;
}

void NewickTree::exportNewick(string filename, bool leafOnly)
{
    return;
}
