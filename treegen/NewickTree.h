#ifndef NEWICKTREE_H
#define NEWICKTREE_H

#include "tree.h"
#include <string>
#include <vector>
#include <utility>
#include <map>
using namespace std;

class NewickTree
{
    struct TreeNode
    {
        string name;
        string sequence;
        TreeNode *parent;
        vector<TreeNode *> children;
        TreeNode() : name(""), sequence(""), parent(nullptr), children(vector<TreeNode *>(0)) {}
        TreeNode(string n, TreeNode *p) : name(n), sequence(""), parent(p), children(vector<TreeNode *>(0)) {}
        TreeNode *addChild(string n)
        {
            TreeNode *newNode = new TreeNode(n, this);
            children.push_back(newNode);
            return newNode;
        }
        TreeNode *addChild(TreeNode *child)
        {
            children.push_back(child);
            return child;
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
    int sequenceLength;
    TreeNode *root;
    map<char, map<char, double>> transitionMatrix;
    inline static vector<char> symbols = {'A', 'T', 'G', 'C'};
    void adjToTree(vector<vector<int>> adj, vector<string> names, int rootIndex);
    void newickDFS(TreeNode *start, vector<string> &symbols, bool leafOnly);
    void deleteDFS(TreeNode *start);
    void toNexus(string filename, map<string, string> sequences);
    TreeNode *importDFS(TreeNode *parent, string &newickString, int &pos);

public:
    NewickTree();
    NewickTree(string rootName);
    NewickTree(vector<Vertex> vertices, vector<Edge> edges, int rootIndex);
    NewickTree(vector<vector<int>> adj, vector<string> names, int rootIndex = 0);
    ~NewickTree();
    int getSequenceLength();
    int getLeafCount();
    void generateSequences(int sequenceLength);
    void printNewick(bool leafOnly = true);
    void exportNewick(string filename, bool leafOnly = true);
    void importNewick(string filename);
    void exportNexus(string filename);
    void setTransition(map<char, map<char, double>> transitionMatrix = map<char, map<char, double>>());
    map<string, string> getSequences();
    static map<string, string> mixtureModel(vector<NewickTree> trees, vector<double> weights);
};

#endif