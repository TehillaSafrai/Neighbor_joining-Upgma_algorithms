import argparse
import os
from typing import List
import numpy as np
from Bio import SeqIO  # pip install biopython
from PhylogeneticTreeNode import PhylogeneticTreeNode
from CalculateDistanceMatrix import calculate_dist_matrix
from Bio import Phylo
from io import StringIO
import matplotlib.pyplot as plt


def parse_fasta_file(file_path: str):
    """
    Parses a FASTA file (plain or gzipped) and returns a mapping of sequence identifiers to nucleotide sequences.
    Parameters:
        file_path (str): The path to the FASTA file.

    Returns:
        dict: A dictionary with sequence IDs as keys and DNA sequences as values.
    """
    sequences = {}

    with open(file_path, 'r') as file_handle:
        for record in SeqIO.parse(file_handle, "fasta"):
            sequences[record.id] = str(record.seq)

    return sequences


def create_new_node(labels, min_i, min_j, node_list: List[PhylogeneticTreeNode],
                    counter) -> PhylogeneticTreeNode:
    """
    Node creation utility that:
    - Generates new internal tree node
    - Updates node lists and labels
    - Maintains tree structure integrity
    - Returns new PhylogeneticTreeNode
    """
    node_name = f'U{counter}' if counter is not None else 'U'
    new_node = PhylogeneticTreeNode(node_name)

    # Find or create child nodes
    new_node.left = next((node for node in node_list if node.name == labels[min_i]),
                         PhylogeneticTreeNode(labels[min_i]))
    new_node.right = next((node for node in node_list if node.name == labels[min_j]),
                          PhylogeneticTreeNode(labels[min_j]))
    new_node.size = new_node.left.size + new_node.right.size

    # Update lists (remove old nodes, add new node)
    indices = [labels.index(new_node.left.name), labels.index(new_node.right.name)]
    max_idx, min_idx = max(indices), min(indices)

    # Remove nodes and labels in correct order (larger index first)
    labels.pop(max_idx)
    node_list.pop(max_idx)
    labels.pop(min_idx)
    node_list.pop(min_idx)

    # Add new node
    labels.append(node_name)
    node_list.append(new_node)
    return new_node


def update_dist_matrix_and_labels(n, min_i, min_j, dist_matrix,
                                  distance_updating_func, new_node):
    """
    Matrix update function that:
    - Recalculates distances after node joins
    - Applies specified distance function
    - Maintains matrix consistency
    - Returns updated distance matrix
    """
    # step 3: calculate the new dists from other seqs.
    new_dist_mat = np.zeros((n - 1, n - 1))
    # change changes dist - all that collaids with the the node
    for k in range(n):
        if k != min_i and k != min_j:
            new_dist_k = distance_updating_func(dist_matrix, min_i, min_j, k, new_node)
            new_index = k if k < min(min_j, min_i) else k - 1 if k < max(min_i, min_j) else k - 2
            # min_index = min(min_i, min_j)
            new_dist_mat[new_index, n - 2] = new_dist_k
            new_dist_mat[n - 2, new_index] = new_dist_k

    # copy unchanged dis
    for i in range(n):
        for j in range(n):
            if i != min_i and i != min_j and j != min_i and j != min_j:
                new_i = i if i < min(min_i, min_j) else i - 1 if i < max(min_i, min_j) else i - 2
                new_j = j if j < min(min_i, min_j) else j - 1 if j < max(min_i, min_j) else j - 2
                new_dist_mat[new_i, new_j] = dist_matrix[i, j]

    return new_dist_mat


def distance_updating_func_upgma(dist_matrix, min_i, min_j, k, new_node):
    """
    UPGMA implementation:
    - Uses ultrametric tree assumption
    - Iteratively joins closest clusters
    - Calculates weighted average distances
    - Returns complete phylogenetic tree
    """
    a_size = new_node.left.size
    b_size = new_node.right.size
    return (dist_matrix[min_i, k] * a_size + dist_matrix[min_j, k] * b_size) / (a_size + b_size)


def distance_updating_func_joining_neighbors(dist_matrix, min_i, min_j, k, new_node):
    """
    Neighbor Joining implementation:
    - Handles non-ultrametric trees
    - Uses Q-matrix for joins
    - Calculates branch lengths dynamically
    - Returns optimized phylogenetic tree
    """
    return 0.5 * (dist_matrix[min_i, k] + dist_matrix[min_j, k] - dist_matrix[min_i, min_j])


def upgma(dist_matrix, labels):
    """Implement UPGMA algorithm for phylogenetic tree construction."""
    n = len(labels)
    node_list = [PhylogeneticTreeNode(label) for label in labels]
    newNodeCounter = 0

    while n > 1:
        # Find minimum distance
        dist_matrix = np.array(dist_matrix, dtype=float)
        np.fill_diagonal(dist_matrix, float('inf'))
        min_i, min_j = np.unravel_index(np.argmin(dist_matrix), dist_matrix.shape)

        # Create new node
        branch_length = dist_matrix[min_i, min_j] / 2
        new_node = create_new_node(labels, min_i, min_j, node_list, newNodeCounter)
        new_node.max_length = branch_length
        if new_node.left.name.startswith('U'):
            new_node.left_length = branch_length - new_node.left.max_length
        else:
            new_node.left_length = branch_length
        if new_node.right.name.startswith('U'):
            new_node.right_length = branch_length - new_node.right.max_length
        else:
            new_node.right_length = branch_length
        newNodeCounter += 1
        # Update distance matrix
        dist_matrix = update_dist_matrix_and_labels(n, min_i, min_j, dist_matrix,
                                                    distance_updating_func_upgma, new_node)
        n -= 1

    return node_list[0]  # Return root node


def change_format_to_newick(node: PhylogeneticTreeNode) -> str:
    """
    Convert a phylogenetic tree to Newick format.
    Args:
        node: Root node of the phylogenetic tree
    Returns:
        String representation of the tree in Newick format
    """
    if node is None:
        return ""

    # Handle leaf nodes
    if node.left is None and node.right is None:
        return f"{node.name}"  # Leaf nodes must have lengths if parent assigned them

    # Process child nodes
    parts = []

    # Add left subtree if it exists
    if node.left:
        left_str = change_format_to_newick(node.left)
        # Ensure every branch has a length, default to 0 if none
        length = node.left_length if node.left_length is not None else 0
        left_str += f":{length:.6f}"
        parts.append(left_str)

    # Add right subtree if it exists
    if node.right:
        right_str = change_format_to_newick(node.right)
        # Ensure every branch has a length, default to 0 if none
        length = node.right_length if node.right_length is not None else 0
        right_str += f":{length:.6f}"
        parts.append(right_str)

    # Join the parts and omit internal node names
    newick = f"({','.join(parts)})"

    return newick


def neighbor_joining(dist_matrix, labels, nodes_list):
    """Implement NJ algorithm for phylogenetic tree construction."""
    n = len(labels)
    new_nodes_counter = 0

    while n != 2:

        # calculate row sums
        R = np.sum(dist_matrix, axis=1)

        # calculate Q matrix
        Q = np.zeros_like(dist_matrix, dtype=float)
        for i in range(n):
            for j in range(n):
                if i != j:
                    Q[i, j] = (n - 2) * dist_matrix[i, j] - R[i] - R[j]

        # find minimum in Q
        # np.fill_diagonal(Q, np.inf)
        np.fill_diagonal(Q, float('inf'))
        min_i, min_j = np.unravel_index(np.argmin(Q), Q.shape)

        # Get all indices where Q equals the minimum value

        # calculate branch length, new node name Z
        diz = 0.5 * dist_matrix[min_i, min_j] + 1 / (2 * (n - 2)) * (R[min_i] - R[min_j])
        djz = 0.5 * dist_matrix[min_i, min_j] + 1 / (2 * (n - 2)) * (R[min_j] - R[min_i])

        # create new node for the two nodes
        new_node = create_new_node(labels, min_i, min_j, nodes_list, new_nodes_counter)
        new_nodes_counter += 1
        new_node.left_length = max(0, diz)
        new_node.right_length = max(0, djz)

        dist_matrix = update_dist_matrix_and_labels(n, labels, min_i, min_j, dist_matrix,
                                                    distance_updating_func_joining_neighbors, new_node)
        n -= 1

    node = PhylogeneticTreeNode(f'U{new_nodes_counter}')
    node.left = nodes_list[0]
    node.right = nodes_list[1]
    node.left_length = dist_matrix[0, 1] / 2
    node.right_length = dist_matrix[0, 1] / 2
    return node


def gen_newick_tree(fasta_path: str, algo: str):
    """
    Reads a FASTA file, computes a distance matrix, and generates a Newick-format UPGMA tree.

    Parameters:
        fasta_path (str): Path to the FASTA file.
        algo (str): Which algorithm will be used for tree construction. Either "nj" or "upgma".

    Prints:
        str: The Newick-format tree string. With the original leaf names, and without internal nodes names.
        for example: ((C:0.018,F:0.018):0.135,(D:0.125,(A:0.085,(B:0.072,E:0.072):0.0127):0.041):0.028);
    """
    sequences = parse_fasta_file(fasta_path)

    dist_matrix, labels = calculate_dist_matrix(sequences)

    tree = [PhylogeneticTreeNode(label) for label in labels]
    if algo.lower() == "upgma":
        tree = upgma(dist_matrix, labels)
    elif algo.lower() == "nj":
        tree = neighbor_joining(dist_matrix, labels, tree)
    else:
        raise ValueError("Algorithm must be either 'nj' or 'upgma'")

    newick_tree = change_format_to_newick(tree)
    base_name = os.path.splitext(os.path.basename(fasta_path))[0]
    output_file = f"{base_name}_{algo.lower()}.newick"
    export_newick_to_file(newick_tree, output_file)
    tree = Phylo.read(StringIO(newick_tree + ';'), "newick")

    # Create a Matplotlib figure for the tree
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)

    plt.title(f'Phylogenetic Tree {base_name} {algo.lower()}', fontsize=16)
    # plt.show()


def export_newick_to_file(newick_str: str, output_file: str):
    """
    Exports a Newick format tree string to a file.

    Parameters:
        newick_str (str): The Newick format string of the tree.
        output_file (str): Path to the output file.

    Returns:
        None
    """
    try:
        with open(output_file, 'w') as file:
            file.write(newick_str + ";\n")
        print(f"Newick tree successfully written to {output_file}")
    except IOError as e:
        print(f"Error writing to file {output_file}: {e}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_path", type=str,
                        help="Path to FASTA file containing multiple DNA sequences.",
                        required=True)
    parser.add_argument("--algo", type=str, help="Either nj or upgma.",
                        required=True)

    args = parser.parse_args()
    gen_newick_tree(args.fasta_path, args.algo)
