import pickle
from SubgroupRes import SubgroupRes
from Candidate import Candidate

NODE_DISTANCE = 3  # define here the desired node distance threshold

###################################################################################################################
"""
The main function to be used. The function removed amiRNAs that are designed to far target genes, and retains those 
whose targets are not too far.
:param dic_amiRNAs: dictionary of candidate AmiRNAs: key = AmiRNA, value = its attributes (should be different for sgRNAs).
:param subgroups: a list of subgroups of genes (representing internal nodes) sorted in a post order.
:param res: A path to the results of amiRNAs per internal node: key = subgroup of genes (internal node), 
        value = AmiRNAs
:return:A dictionary where key = a subgroup of genes, representing an internal node that was empty before. value: a dictionary
        # of AmiRNAs and their attributes
"""


def remove_candidates_with_distant_genes(node, subgroupres, max_gap_distance=3):
    """
    I will not call this function if the number of genes is smaller than the distance
    This function filters out genes that target genes that are too distanct in the tree.
    :param node: The internal node in the gene tree
    :param subgroupres: a SubgroupRes object containing all sgRNAs that were designed for the node in the tree
    :return:
    """
    gene_pair2distance_dict = create_gene2distance_dict(
        node)  # Calculate the distance between all genes in the internal node.
    for candidate in subgroupres.candidates_list:
        filtered_candidates_list = []
        genes_targeted = set(candidate.genes_score_dict.keys())

        list_of_clusters = []
        get_gene_clusters(genes_targeted, node,
                          list_of_clusters)  # Divide the genes targeted by the sgRNA into monophyletic groups.

        is_near = decide_if_near(list_of_clusters, gene_pair2distance_dict, max_gap_distance)

        if is_near:
            filtered_candidates_list.append(candidate)
    subgroupres.candidates_list = filtered_candidates_list  # Update the list of candidates
    return


def decide_if_near(list_of_clusters, gene_pair2distance_dict, max_gap_distance=3):
    for i, cluster_i in enumerate(list_of_clusters):
        for cluster_j in list_of_clusters[i + 1:]:
            near = False
            for gene_i in cluster_i:
                for gene_j in cluster_j:
                    # If at least two genes from each cluster have a distance smaller or equal to the threshold,
                    # The two clusters will be considered as close enough
                    if gene_pair2distance_dict[tuple(sorted([gene_i, gene_j]))] <= max_gap_distance:
                        near = True
                        break
                if near:
                    break
            # If at least one pair of clusters is too spaced, the genes targeted by the sgRNA are not considered near.
            if not near:
                return False
    return True


def find_distance_from_lca(gene_i, gene_j, node, gene2node_dict):
    """
    This function finds the distance of gene i from the least common ancestor of genes i and j
    :param gene2node_dict: A dictionary
    :param node:
    :param gene_i:
    :param gene_j:
    :return:
    """
    node_i = find_leaf(node, gene_i, gene2node_dict)  # Find the leaf node containing the gene
    dist_i = 0
    while gene_j not in set(node_i.node_leaves):  # From each node, go to the least common-ancestor of genes i & j
        node_i = node_i.parent
        dist_i += 1
    return dist_i  # Return the distance


def find_leaf(node, gene, gene2node_dict):
    """
    This function finds the leaf node of a gene within a subtree.
    :param node: A tree node.
    :param gene: The name of the gene
    :param gene2node_dict: A dictionary of gene -> node
    :return: A CladeNew object of the node with the single gene. If it's not in gene2dist_from_root_dict, it is updated
    """
    if gene in gene2node_dict:  # If the gene is already in the dictionary, return the value.
        return gene2node_dict[gene]
    gene2node_dict[gene] = node
    while gene2node_dict[gene].node_leaves != [gene]:  # Go down in the tree until the gene-node is reached.
        if gene in gene2node_dict[gene].clades[0].node_leaves:
            gene2node_dict[gene] = gene2node_dict[gene].clades[0]
        else:
            gene2node_dict[gene] = gene2node_dict[gene].clades[1]
    return gene2node_dict[gene]


def create_gene2distance_dict(node):
    """
    This function creates a dictionary that stores the distances between each pair of genes.
    :param node: A tree node
    :return: A dictionary (genei,genej) -> the distance between genei and gene j
    """
    gene_pair2distance_dict = {}
    gene2node_dict = {}  # Create a dictionary that stores the node
    for i, gene_i in enumerate(node.node_leaves):
        for gene_j in node.node_leaves[i + 1:]:
            if gene_i == gene_j or tuple(sorted([gene_i, gene_j])) in gene_pair2distance_dict:
                continue  # Continue if the genes are identical or if the distance was already calculated.

            dist_i = find_distance_from_lca(gene_i, gene_j, node,
                                            gene2node_dict)  # Find the distance of each gene from the least common ancestor
            dist_j = find_distance_from_lca(gene_j, gene_i, node, gene2node_dict)
            gene_pair2distance_dict[tuple(sorted([gene_i, gene_j]))] = dist_i + dist_j - 1
    return gene_pair2distance_dict


def get_gene_clusters(genes_targeted: set, node, list_of_clusters):
    """
    This recursive function divides the genes targeted by the sgRNA into monophyletic groups
    :param genes_targeted: The genes targeted by the sgRNA
    :param node: A tree node
    :param list_of_clusters: A list of monophyletic clusters
    :return: None
    """
    node_leaves_set = set(node.node_leaves)
    # If the genes targeted are the same as genes in the node, append the monophyletic group to the list
    if genes_targeted == node_leaves_set:
        list_of_clusters.append(node_leaves_set)
        return
    # If not, it means that the clusters are deeper in the tree.
    intersection_left = genes_targeted.intersection(set(node_leaves_set.clades[0].node_leaves))
    intersection_right = genes_targeted.intersection(set(node_leaves_set.clades[1].node_leaves))
    # If at least one gene appears in both the targeted genes set and the left/right clade, recursively call the function,
    # with the left/right child as the new node, and the intersection as the new set of genes targeted
    if intersection_left:
        get_gene_clusters(intersection_left, node.clades[0], list_of_clusters)
    if intersection_right:
        get_gene_clusters(intersection_right, node.clades[1], list_of_clusters)
    return
