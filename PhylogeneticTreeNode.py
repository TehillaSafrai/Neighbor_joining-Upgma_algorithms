class PhylogeneticTreeNode:
    """Represents a node in the phylogenetic tree.
        max_length = the maximum length of the node in the phylogenetic tree from the present
        name = the name of the node in the phylogenetic tree
        left = the node's left child'
        right = the node's right child'
        size = the node number of children, grandchildren and so on, in the phylogenetic tree
        left_length = the length of the left child of the node in the phylogenetic tree
        right_length = the length of the right child of the node in the phylogenetic tree
    """

    def __init__(self, name: str, left=None, right=None):
        self.max_length = 0
        self.name = name
        self.left = left
        self.right = right
        self.size = 1
        self.left_length = 0.0
        self.right_length = 0.0

