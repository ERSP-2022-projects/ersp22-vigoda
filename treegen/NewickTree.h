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
    static int createdCount;
    int sequenceLength;
    TreeNode *root;
    map<char, map<char, double>> transitionMatrix;
    inline static vector<char> symbols = {'A', 'T', 'G', 'C'};
    void adjToTree(vector<vector<int>> adj, vector<string> names, int rootIndex);
    void newickDFS(TreeNode *start, vector<string> &symbols, bool leafOnly);
    void deleteDFS(TreeNode *start);
    TreeNode *importDFS(TreeNode *parent, string &newickString, int &pos);

public:
    /**
     * Constructs an empty tree with a null 'root'
     */
    NewickTree();
    /**
     * Constructs a tree that only has a 'root' named 'rootName'
     *
     * @param rootName The name of the root of the tree
     */
    NewickTree(string rootName);
    /**
     * Constructs a random tree with the given number of leaf Nodes. Leaves (tips) will be labeled [1, 2, ..., numLeafs].
     *
     * @param numLeafs The number of leaf nodes in this random tree
     */
    NewickTree(int numLeafs);
    /**
     * Constructs a tree from a Vertex-Edge representation (see tree.h for details).
     *
     * @param vertices A vector of all the Vertices in a graph with their name
     * @param edges A vector of all the Edges in the graph
     * @param rootIndex The index of the root Vertex in 'vertices'.
     */
    NewickTree(vector<Vertex> vertices, vector<Edge> edges, int rootIndex);
    /**
     * Constructs a tree from an adjacency list representation. Note that naming of leaf and internal nodes is not done automatically and users
     * should take care to name their nodes appropriately.
     *
     * @param adj Adjacency list representation of desired tree.
     * @param names A vector of the name of each vertex. Index-aligned with 'adj'.
     * @param rootIndex The index of the root vertex in 'adj' and 'names'.
     */
    NewickTree(vector<vector<int>> adj, vector<string> names, int rootIndex = 0);
    /**
     * Deletes the NewickTree by doing a post-order traversal of each TreeNode starting from the root and deleting each TreeNode pointer.
     */
    ~NewickTree();
    /**
     * @return The length of character sequences in this tree. (0 if no sequences have been generated).
     */
    int getSequenceLength();
    /**
     * @return The number of leaf nodes in this tree.
     */
    int getLeafCount();
    /**
     * Generates a genetic sequence of length 'sequenceLength' on each TreeNode in the NewickTree. This process starts by generating a random
     * random sequence for 'root' and then changing each mutation site according to the probabilities set in 'transitionMatrix'. By default
     * the probability of any site mutating is 0.3 divided evenly among all of the other nucleotides. Other transition matrices can be set using
     * setTransition.
     *
     * @param sequenceLength The length of the genetic sequence to be generated (must be greater than 0).
     */
    void generateSequences(int sequenceLength);
    /**
     * Prints current tree topology in Newick format to the console.
     *
     * @param leafOnly Only print the leaf nodes in the tree (True by default).
     */
    void printNewick(bool leafOnly = true);
    /**
     * Stores the current tree topology in Newick format to a file.
     *
     * @param filename The file to store the tree topology in.
     * @param leafOnly Only print the leaf nodes in the tree (True by default).
     */
    void exportNewick(string filename, bool leafOnly = true);
    /**
     * Imports a tree topology from Newick format
     *
     * @param filename The file to import from.
     */
    void importNewick(string filename);
    /**
     * Exports the genetic sequences and current tree topology in Nexus format.
     *
     * @param filename The file to store the Nexus file in.
     */
    void exportNexus(string filename, bool numericSort = false);
    /**
     * Sets the transition matrix for generating genetic sequences (0.3 mutation rate Jukes-Cantor by default).
     *
     * @param transitionMatrix Map of nucleotides to their transition probability to any other nucleotide.
     * (e.g. transitionMatrix['A']['G']) is probability of Adenine mutating to Guanine between neighboring nodes.
     */
    void setTransition(map<char, map<char, double>> transitionMatrix = map<char, map<char, double>>());
    /**
     * Gets genetic sequences associated with each leaf node.
     *
     * @return Map of the name of each leaf node to its associated genetic sequence.
     */
    map<string, string> getSequences();
    /**
     * Converts a map of genetic sequences to Nexus format and stores it in a file.
     *
     * @param filename The file to store the Nexus file in.
     * @param sequences A map of the name of each leaf node to its associated genetic sequence.
     * @param trees (Optional) A list of the tree topologies to be stored in the Nexus format in NewickTree format.
     * @param numericSort (Optional) Set to true if the name of all leaf nodes can be converted to an integer and you would like them ordered
     * in ascending numeric order in the Nexus file.
     */

    static void toNexus(string filename, map<string, string> sequences, vector<NewickTree *> trees = vector<NewickTree *>(), bool numericSort = false);
    /**
     * Constructs a mixture model out of multiple tree topologies.
     *
     * @param trees A vector of NewickTree * representing each of the tree topologies in the mixture.
     * @param weights The index-aligned weight associated with each tree in 'trees'. Must sum up to 1. (Equally divided by default).
     */
    static map<string, string> mixtureModel(vector<NewickTree *> trees, vector<double> weights);
};

#endif