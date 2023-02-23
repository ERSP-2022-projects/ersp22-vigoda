from Bio import Phylo
from io import StringIO
import re
from ete3 import Tree


def validate(trprobsfilePath, txtFilePath):
    mcmcTree = displayAndReturnTree(trprobsfilePath,display= False)
   # trueTree =  

def generateNewickTree(nleafs):
    return treeToNewick(generate_random_tree(nleafs))

def generate_random_tree(nleafs,model="yule"):
    tree = Tree()
    tree.populate(nleafs,model)
    return tree

def treeToNewick(tree):
    return tree.write(format=1) #format 1 is newick format

def calculate_rf_distance(tree1, tree2):
    t1 = Tree(tree1)
    t2 = Tree(tree2)

    # Calculate the Robinson-Foulds distance
    rf_distance = t1.robinson_foulds(t2, unrooted_trees=True)[0]

    return rf_distance

    
def displayAndReturnTree(filepath, display = True):

    with open(filepath, "r") as handle:
        text = handle.read()

    # extract the translate block and parse it into a dictionary
    translate_block = re.search(r"translate(.*?);", text, re.DOTALL).group(1)
    label_dict = {}
    for match in re.finditer(r"(\d+) ([^,\n]+)", translate_block):
        label_dict[match.group(1)] = match.group(2)

    # extract the tree string and replace labels with taxon names
    tree_str = re.search(r"tree tree_1 .*?\((.*?)\);", text, re.DOTALL).group(1)
    for label, name in label_dict.items():
        tree_str = re.sub(r"\b%s\b" % label, name, tree_str)

    # create a Phylo tree object and draw it
    tree = Phylo.read(StringIO(tree_str + ";"), "newick")
    tree_str = "(" + tree_str + ");"
    if display:
        Phylo.draw(tree)
    
    return tree_str


def displayTree(tree_str):
    """Takes in a tree in newick string format and displays the tree 
    """
    tree = Phylo.read(StringIO(tree_str), "newick")
    Phylo.draw(tree) 

# def calculate_distance(tree1, tree2):
#     t1 = Phylo.read(StringIO(tree1), "newick")
#     t2 = Phylo.read(StringIO(tree2), "newick")
    
#     if t2 is None:
#         print("Error: No second tree provided.")
#         return None
    
#     return t1.distance(t2)




def main():
    pass
if __name__ == "__main__":
    main()