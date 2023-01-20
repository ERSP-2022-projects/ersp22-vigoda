#ifndef NEWICK_TREE
#define NEWICK_TREE

#include "tree.h"
#include <string>
#include <vector>
using namespace std;

class NewickTree
{
    struct TreeNode
    {
        string name;
        TreeNode *parent;
        vector<TreeNode *> children;
        TreeNode() {}
        TreeNode(string n, TreeNode *p) : name(n), parent(p) {}
        void addChild(string n)
        {
            children.push_back(new TreeNode(n, this));
        }
        int childrenCount()
        {
            return children.size();
        }
        bool isLeaf()
        {
            return (children.size() == 0);
        }
    };

private:
    TreeNode *root;
    void newickDFS(TreeNode *start, vector<string> &symbols);

public:
    NewickTree();
    NewickTree(string rootName);
    NewickTree(vector<Vertex> vertices, vector<Edge> edges, int rootIndex);
    NewickTree(vector<vector<int>> adj, vector<string> names, int rootIndex = 0);
    //~NewickTree();
    void printNewick(bool leafOnly = true);
    void exportNewick(string filename, bool leafOnly = true);
};

#endif