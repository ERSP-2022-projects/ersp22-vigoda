import random
from ete3 import Tree


class PhylogeneticTree:
    def __init__(self, nTaxa, nChars):
        self.nTaxa = nTaxa
        self.nChars = nChars
        self.newickTree = self._generate_random_tree(nTaxa)

    def _generate_random_tree(self, nTaxa):
        # Generate a random tree with the specified number of taxa
        tree = Tree()
        leaf_names = [str(i + 1) for i in range(nTaxa)]
        tree.populate(len(leaf_names), names_library=leaf_names)

        # Convert the tree to Newick format and return it as a string
        return tree.write(format=1)


# Example usage
nTaxa = 4
nChars = 10
phylo_tree = PhylogeneticTree(nTaxa, nChars)
print(phylo_tree.newickTree)
